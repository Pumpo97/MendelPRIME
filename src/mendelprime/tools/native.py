"""Native Exomiser and LIRICAL execution and complete output parsing."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import logging
import math
import os
import re
import tempfile
import time
import uuid
from collections import defaultdict
from collections.abc import Mapping
from pathlib import Path

import yaml

from mendelprime.domain import MendelPrimeError, Patient, Ranking, Source, validate_gene, variant_key
from mendelprime.runtime.process import run_process
from mendelprime.runtime.deployment import check_java_distribution
from mendelprime.runtime.storage import open_regular, private_directory, stat_fingerprint
from mendelprime.tools.vcf import (
    inspect_vcf,
    materialize_samples,
    read_observations,
    read_pedigree,
    resolve_family,
)

LOGGER = logging.getLogger(__name__)


def kinship_observation(executable, vcf_path, seed, work_root):
    """Run the preserved QC-only KING protocol; never infer roles from column order."""
    family, case_key = dict(seed["family"]), seed["profile"]["case_key"]
    if family["mode"] != "trio":
        return None
    roles = dict(family["roles"])
    record = {
        "status": "unavailable",
        "first_degree_both": False,
        "informative_sites_sufficient": False,
        "roles": roles,
        "interpretation": "Statistical kinship is not laboratory chain-of-custody confirmation; per-site quality is assessed separately.",
    }
    if executable and vcf_path:
        try:
            path = Path(vcf_path)
            with open_regular(path) as (stream, before):
                vcf_digest = hashlib.file_digest(stream, "sha256").hexdigest()
                vcf_before = stat_fingerprint(before)
            with open_regular(executable) as (stream, tool_before):
                tool_digest = hashlib.file_digest(stream, "sha256").hexdigest()
            if set(inspect_vcf(path)["samples"]) != set(roles.values()):
                raise ValueError("sample mismatch")
            original = {r["variant_key"]: r for r in seed["observations"]["observations"]}
            observed = read_observations(path, family, tuple(original), case_key=case_key)
            if any(
                r["samples"] != original[r["variant_key"]]["samples"] for r in observed["observations"]
            ) or len(observed["observations"]) != len(original):
                raise ValueError("frozen observations changed")
            with tempfile.TemporaryDirectory(prefix="mendelprime-kinship-", dir=work_root) as temporary:
                root = Path(temporary)
                parents, sex = root / "parents.tsv", root / "sex.tsv"
                parents.write_text(
                    f"{roles['proband']}\t{roles['father']}\t{roles['mother']}\n{roles['father']}\t0\t0\n{roles['mother']}\t0\t0\n"
                )
                sex.write_text(
                    f"#IID\tSEX\n{roles['proband']}\t{ {'MALE': 1, 'FEMALE': 2}.get(family['sex'], 0) }\n{roles['father']}\t1\n{roles['mother']}\t2\n"
                )
                prefix = root / "qc"
                command = [
                    str(executable),
                    "--vcf",
                    str(path),
                    "--update-parents",
                    str(parents),
                    "--update-sex",
                    str(sex),
                    "--snps-only",
                    "just-acgt",
                    "--max-alleles",
                    "2",
                    "--autosome",
                    "--make-king-table",
                    "counts",
                    "--threads",
                    "2",
                    "--var-filter",
                    "--vcf-min-gq",
                    "20",
                    "--vcf-min-dp",
                    "10",
                    "--out",
                    str(prefix),
                ]
                receipt = run_process(
                    command=command, root=root, timeout_seconds=180, artifacts=(prefix.with_suffix(".kin0"),)
                )
                if receipt["error_code"]:
                    raise MendelPrimeError("kinship_" + receipt["error_code"], category="resource")
                if stat_fingerprint(path.stat()) != vcf_before or stat_fingerprint(
                    executable.stat()
                ) != stat_fingerprint(tool_before):
                    raise MendelPrimeError("kinship_input_changed", category="resource")
                with prefix.with_suffix(".kin0").open() as stream:
                    pairs = {
                        frozenset((r["#IID1"], r["IID2"])): {
                            "nsnp": int(r["NSNP"]),
                            "hethet": int(r["HETHET"]),
                            "ibs0": int(r["IBS0"]),
                            "kinship": float(r["KINSHIP"]),
                        }
                        for r in csv.DictReader(stream, delimiter="\t")
                    }
                rows = {
                    name: pairs[frozenset((roles[a], roles[b]))]
                    for name, a, b in (
                        ("proband_father", "proband", "father"),
                        ("proband_mother", "proband", "mother"),
                        ("father_mother", "father", "mother"),
                    )
                }
                children, parents_row = (
                    (rows["proband_father"], rows["proband_mother"]),
                    rows["father_mother"],
                )
                first = all(0.177 <= r["kinship"] < 0.354 for r in children)
                sufficient = all(r["nsnp"] >= 10000 for r in children)
                ibs = parents_row["nsnp"] > 0 and all(
                    r["nsnp"] > 0 and r["ibs0"] / r["nsnp"] < parents_row["ibs0"] / parents_row["nsnp"]
                    for r in children
                )
                record.update(
                    vcf_sha256=vcf_digest,
                    executable_sha256=tool_digest,
                    vcf_stat=list(vcf_before),
                    status="kinship_supports_recorded_trio"
                    if first and sufficient and ibs
                    else "needs_review",
                    first_degree_both=first,
                    informative_sites_sufficient=sufficient,
                    child_ibs0_below_parent_pair=ibs,
                    relationships=rows,
                    filters={"PASS": True, "DP_when_present": 10, "GQ_when_present": 20},
                    method=next(iter((root / "stdout.log").read_text().splitlines()),
                                next(iter((root / "stderr.log").read_text().splitlines()), "sha256:" + tool_digest)),
                )
                LOGGER.info(
                    "[%s/%s] post kinship status=%s elapsed=%.2fs",
                    *case_key,
                    record["status"],
                    receipt["seconds"],
                )
        except Exception as error:
            record["error_code"] = (
                error.code if isinstance(error, MendelPrimeError) else "kinship_evidence_unavailable"
            )
    text = json.dumps(record, sort_keys=True, separators=(",", ":"))
    text += "\n" + "\n".join(
        f"{k}={str(record[k]).lower()}"
        for k in ("status", "first_degree_both", "informative_sites_sufficient")
    )
    return Source(
        provider="plink2-king",
        kind="patient_observation",
        locator="method:family/king-qc",
        version="preserved-20260925-kinship-rules-v1",
        content=text,
        case_key=case_key,
    )


_HGVS = re.compile(
    r"(?P<transcript>(?:ENS[TNR]|N[MRP]_|X[MRP]_)[^:]*):"
    r"(?:(?P<hgvs_c>c\.[^:]+))?(?::(?P<hgvs_p>p\.[^:]+))?"
)
_NUMBER = r"(?:[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?|[-+]?(?:nan|inf(?:inity)?))"
_LIRICAL_VARIANT = re.compile(
    rf"(?P<chromosome>[^:\s]+):(?P<position>\d+)(?P<reference>[^>\s]+)>"
    rf"(?P<alternate>\S+)\s+(?P<annotation>.+?)\s+"
    rf"pathogenicity:(?P<pathogenicity>{_NUMBER})\s+\[(?P<genotype>[^]]+)\]",
    re.IGNORECASE,
)
_EXOMISER_GENE_HEADERS = frozenset({"#RANK", "GENE_SYMBOL"})
_EXOMISER_VARIANT_HEADERS = frozenset({"GENE_SYMBOL", "CONTIG", "START", "REF", "ALT"})
_LIRICAL_HEADERS = frozenset({"rank", "diseaseCurie", "entrezGeneId"})
_GENE_MAP_HEADERS = frozenset({"entrez_id", "symbol"})
_FALLBACK_RANK = 10**9
ARTIFACTS = {
    "exomiser": ("exomiser.genes.tsv", "exomiser.variants.tsv"),
    "lirical": ("lirical.tsv", "lirical.json"),
}

_FREQUENCY_SOURCES = (
    "UK10K",
    "GNOMAD_E_AFR",
    "GNOMAD_E_AMR",
    "GNOMAD_E_EAS",
    "GNOMAD_E_NFE",
    "GNOMAD_E_SAS",
    "GNOMAD_G_AFR",
    "GNOMAD_G_AMR",
    "GNOMAD_G_EAS",
    "GNOMAD_G_NFE",
    "GNOMAD_G_SAS",
)
_REMOVED_EFFECTS = (
    "FIVE_PRIME_UTR_EXON_VARIANT",
    "FIVE_PRIME_UTR_INTRON_VARIANT",
    "THREE_PRIME_UTR_EXON_VARIANT",
    "THREE_PRIME_UTR_INTRON_VARIANT",
    "NON_CODING_TRANSCRIPT_EXON_VARIANT",
    "NON_CODING_TRANSCRIPT_INTRON_VARIANT",
    "CODING_TRANSCRIPT_INTRON_VARIANT",
    "UPSTREAM_GENE_VARIANT",
    "DOWNSTREAM_GENE_VARIANT",
    "INTERGENIC_VARIANT",
    "REGULATORY_REGION_VARIANT",
)
_INHERITANCE_MODES = {
    "AUTOSOMAL_DOMINANT": 0.1,
    "AUTOSOMAL_RECESSIVE_HOM_ALT": 0.1,
    "AUTOSOMAL_RECESSIVE_COMP_HET": 2,
    "X_DOMINANT": 0.1,
    "X_RECESSIVE_HOM_ALT": 0.1,
    "X_RECESSIVE_COMP_HET": 2,
    "MITOCHONDRIAL": 0.2,
}
_SEX_CODES = {"UNKNOWN": "0", "MALE": "1", "FEMALE": "2"}


def _error(field, code, category="response"):
    return MendelPrimeError(code, category=category, field=field)


def _check(condition, field, code, category="response"):
    if not condition:
        raise _error(field, code, category)


def _issue(code, artifact, record_number=None, gene=None):
    return {
        "code": code,
        "artifact": artifact,
        "record_number": record_number,
        "gene": gene,
    }


def _artifact(path):
    try:
        with open_regular(path) as (stream, before):
            raw = stream.read()
            after = os.fstat(stream.fileno())
    except OSError:
        raise _error(f"artifact.{path.name}", "unreadable_regular_file") from None
    _check(stat_fingerprint(before) == stat_fingerprint(after), path.name, "source_changed_during_read")
    try:
        content = raw.decode("utf-8")
    except UnicodeDecodeError:
        raise _error(f"artifact.{path.name}", "utf8_required") from None
    _check(bool(content), f"artifact.{path.name}", "nonempty_text_required")
    return {
        "name": path.name,
        "content": content,
        "content_sha256": hashlib.sha256(raw).hexdigest(),
    }


def _tsv_rows(content, artifact, required_headers, *, skip_prefix=""):
    semantic = content.removeprefix("\ufeff")
    if skip_prefix:
        semantic = "".join(
            line for line in semantic.splitlines(keepends=True) if not line.startswith(skip_prefix)
        )
    reader = csv.reader(io.StringIO(semantic, newline=""), delimiter="\t")
    headers = next((values for values in reader if values and any(values)), None)
    _check(headers is not None, f"{artifact}.headers", "header_required")
    _check(headers and all(headers), f"{artifact}.headers", "nonempty_headers_required")
    _check(len(headers) == len(set(headers)), f"{artifact}.headers", "duplicate_header")
    _check(required_headers <= set(headers), f"{artifact}.headers", "required_headers_missing")
    rows, malformed, record_number = [], [], 0
    for values in reader:
        if not values or all(value == "" for value in values):
            continue
        record_number += 1
        if len(values) > len(headers):
            malformed.append(record_number)
            continue
        padded = [*values, *(None for _ in range(len(headers) - len(values)))]
        rows.append((record_number, dict(zip(headers, padded, strict=True))))
    return rows, tuple(malformed)


def _semantic_json(value):
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def _source(provider, artifact, number, version, case_key, gene, content, key=None, fragment=None):
    locator = f"method:{provider}/{artifact}#record:{number}"
    if fragment is not None:
        locator += f"/variant:{fragment}"
    return Source(
        provider=provider,
        kind="tool_observation",
        locator=locator,
        version=version,
        content=content,
        genes=(gene,),
        variant_keys=() if key is None else (key,),
        case_key=case_key,
    )


def _optional(value):
    if value is None:
        return None
    text = str(value)
    return text if text else None


def _int(value):
    text = str(value or "")
    if text in {"", "."}:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def _float(value):
    text = str(value or "")
    if text in {"", "."}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _probability(value):
    text = str(value or "").strip()
    if text.endswith("%"):
        parsed = _float(text[:-1])
        return None if parsed is None else parsed / 100
    return _float(text)


def _bool(value):
    text = str(value or "").lower()
    if text in {"1", "true", "yes"}:
        return True
    if text in {"0", "false", "no"}:
        return False
    return None


def _split(value):
    return tuple(dict.fromkeys(item.strip() for item in str(value or "").split(",") if item.strip()))


def _rank(value):
    parsed = _int(value)
    return parsed if parsed not in {None, 0} else _FALLBACK_RANK


def _gene(value):
    gene = str(value or "").strip().upper()
    try:
        validate_gene(gene)
    except MendelPrimeError:
        return None
    return gene


def _entrez(value):
    return str(value or "").strip().removeprefix("NCBIGene:")


def _canonical_chromosome(chromosome):
    canonical = chromosome.removeprefix("chr")
    return "MT" if canonical == "M" else canonical


def _variant_identity(assembly, chromosome, position, reference, alternate):
    if position is None:
        return None
    try:
        return variant_key(assembly, chromosome, position, reference.upper(), alternate.upper())
    except MendelPrimeError:
        return None


def _validate_parse_request(case_key, assembly, version, gene_limit):
    _check(
        type(case_key) is tuple and len(case_key) == 2,
        "native.case_key",
        "case_key_required",
        "contract",
    )
    _check(assembly in {"GRCh37", "GRCh38"}, "native.assembly", "assembly_required", "contract")
    _check(
        isinstance(version, str) and bool(version.strip()) and version == version.strip(),
        "native.version",
        "canonical_version_required",
        "contract",
    )
    _check(
        type(gene_limit) is int and gene_limit > 0,
        "native.gene_limit",
        "positive_integer_required",
        "contract",
    )


def _native_gene(gene, records, variants=(), unparsed=()):
    return {
        "gene": gene,
        "records": tuple(records),
        "variants": tuple(variants),
        "unparsed_variants": tuple(unparsed),
    }


def _tool_result(case_key, tool_id, version, assembly, ranking, genes, artifacts, issues):
    return {
        "case_key": case_key,
        "tool_id": tool_id,
        "version": version,
        "assembly": assembly,
        "ranking": Ranking(genes=ranking),
        "genes": tuple(genes),
        "artifacts": tuple(artifacts),
        "issues": tuple(issues),
    }


def _variant_record(gene, key, chromosome, position, reference, alternate, match, genotype, details, source):
    return {
        "gene": gene,
        "key": key,
        "chromosome": chromosome,
        "position": position,
        "reference": reference,
        "alternate": alternate,
        "transcript": match.group("transcript") if match else None,
        "hgvs_c": match.group("hgvs_c") if match else None,
        "hgvs_p": match.group("hgvs_p") if match else None,
        "genotype_raw": genotype,
        "details": details,
        "source": source,
    }


def _exomiser_gene_record(row, gene, number, case_key, version, artifact):
    return {
        "source": _source("exomiser", artifact, number, version, case_key, gene, _semantic_json(row)),
        "rank": _int(row.get("#RANK")),
        "moi": _optional(row.get("MOI")),
        "p_value": _float(row.get("P-VALUE")),
        "combined_score": _float(row.get("EXOMISER_GENE_COMBINED_SCORE")),
        "phenotype_score": _float(row.get("EXOMISER_GENE_PHENO_SCORE")),
        "variant_score": _float(row.get("EXOMISER_GENE_VARIANT_SCORE")),
        "human_phenotype_score": _float(row.get("HUMAN_PHENO_SCORE")),
        "mouse_phenotype_score": _float(row.get("MOUSE_PHENO_SCORE")),
        "fish_phenotype_score": _float(row.get("FISH_PHENO_SCORE")),
        "walker_score": _float(row.get("WALKER_SCORE")),
        "phive_score": _float(row.get("PHIVE_ALL_SPECIES_SCORE")),
        "omim_score": _float(row.get("OMIM_SCORE")),
        "human_phenotype_evidence": row.get("HUMAN_PHENO_EVIDENCE") or "",
        "mouse_phenotype_evidence": row.get("MOUSE_PHENO_EVIDENCE") or "",
        "fish_phenotype_evidence": row.get("FISH_PHENO_EVIDENCE") or "",
        "human_ppi_evidence": row.get("HUMAN_PPI_EVIDENCE") or "",
        "mouse_ppi_evidence": row.get("MOUSE_PPI_EVIDENCE") or "",
        "fish_ppi_evidence": row.get("FISH_PPI_EVIDENCE") or "",
    }


def _exomiser_variant(row, gene, number, case_key, assembly, version, artifact):
    chromosome = row.get("CONTIG") or ""
    position = _int(row.get("START"))
    reference, alternate = row.get("REF") or "", row.get("ALT") or ""
    key = _variant_identity(assembly, _canonical_chromosome(chromosome), position, reference, alternate)
    match = _HGVS.search(row.get("HGVS") or "")
    source = _source("exomiser", artifact, number, version, case_key, gene, _semantic_json(row), key)
    details = {
        "moi": _optional(row.get("MOI")),
        "consequence": _optional(row.get("FUNCTIONAL_CLASS")),
        "quality": _float(row.get("QUAL")),
        "filter": _optional(row.get("FILTER")),
        "contributing": _bool(row.get("CONTRIBUTING_VARIANT")),
        "whitelist": _bool(row.get("WHITELIST_VARIANT")),
        "variant_score": _float(row.get("EXOMISER_VARIANT_SCORE")),
        "acmg_classification": _optional(row.get("EXOMISER_ACMG_CLASSIFICATION")),
        "acmg_evidence": _split(row.get("EXOMISER_ACMG_EVIDENCE")),
        "acmg_disease_id": _optional(row.get("EXOMISER_ACMG_DISEASE_ID")),
        "acmg_disease_name": _optional(row.get("EXOMISER_ACMG_DISEASE_NAME")),
        "clinvar_id": _optional(row.get("CLINVAR_VARIATION_ID")),
        "clinvar_classification": _optional(row.get("CLINVAR_PRIMARY_INTERPRETATION")),
        "clinvar_stars": _int(row.get("CLINVAR_STAR_RATING")),
        "loeuf": _float(row.get("GENE_CONSTRAINT_LOEUF")),
        "loeuf_lower": _float(row.get("GENE_CONSTRAINT_LOEUF_LOWER")),
        "loeuf_upper": _float(row.get("GENE_CONSTRAINT_LOEUF_UPPER")),
        "max_frequency_source": _optional(row.get("MAX_FREQ_SOURCE")),
        "max_frequency_percent": _float(row.get("MAX_FREQ")),
        "all_frequencies": _optional(row.get("ALL_FREQ")),
        "max_pathogenicity_source": _optional(row.get("MAX_PATH_SOURCE")),
        "max_pathogenicity_score": _float(row.get("MAX_PATH")),
        "all_pathogenicity_scores": _optional(row.get("ALL_PATH")),
    }
    variant = _variant_record(
        gene,
        key,
        chromosome,
        position,
        reference,
        alternate,
        match,
        _optional(row.get("GENOTYPE")),
        details,
        source,
    )
    return variant, None if key is not None else _issue(
        "unsupported_variant_identity", artifact, number, gene
    )


def variant_is_candidate(variant):
    """Return native computation eligibility without discarding the raw observation."""
    if variant["source"].provider != "lirical":
        return True
    score = variant.get("details", {}).get("pathogenicity")
    if type(score) is int:
        return score > 0
    return type(score) is float and math.isfinite(score) and score > 0


def parse_exomiser(root, case_key, assembly, version, gene_limit=100):
    """Parse all Exomiser gene and variant rows without projection limits."""
    _validate_parse_request(case_key, assembly, version, gene_limit)
    gene_artifact = _artifact(root / "exomiser.genes.tsv")
    variant_artifact = _artifact(root / "exomiser.variants.tsv")
    gene_rows, malformed_genes = _tsv_rows(
        gene_artifact["content"], gene_artifact["name"], _EXOMISER_GENE_HEADERS
    )
    variant_rows, malformed_variants = _tsv_rows(
        variant_artifact["content"], variant_artifact["name"], _EXOMISER_VARIANT_HEADERS
    )
    issues = [
        *(_issue("malformed_row", gene_artifact["name"], number) for number in malformed_genes),
        *(_issue("malformed_row", variant_artifact["name"], number) for number in malformed_variants),
    ]
    grouped = defaultdict(list)
    for number, row in gene_rows:
        gene = _gene(row.get("GENE_SYMBOL"))
        if gene is None:
            issues.append(_issue("invalid_gene", gene_artifact["name"], number))
            continue
        grouped[gene].append((number, row))
        if _int(row.get("#RANK")) is None:
            issues.append(_issue("invalid_rank", gene_artifact["name"], number, gene))
    _check(bool(grouped), "exomiser.ranking", "empty_ranking")
    ordered = tuple(
        sorted(
            grouped,
            key=lambda gene: (min(_rank(row.get("#RANK")) for _, row in grouped[gene]), gene),
        )
    )
    variants = defaultdict(list)
    for number, row in variant_rows:
        gene = _gene(row.get("GENE_SYMBOL"))
        if gene is None:
            issues.append(_issue("invalid_gene", variant_artifact["name"], number))
            continue
        if gene not in grouped:
            issues.append(_issue("variant_gene_unranked", variant_artifact["name"], number, gene))
            continue
        variant, issue = _exomiser_variant(
            row, gene, number, case_key, assembly, version, variant_artifact["name"]
        )
        variants[gene].append(variant)
        if issue is not None:
            issues.append(issue)
    genes = []
    for gene in ordered:
        records = tuple(
            _exomiser_gene_record(row, gene, number, case_key, version, gene_artifact["name"])
            for number, row in sorted(
                grouped[gene], key=lambda item: (_rank(item[1].get("#RANK")), item[1].get("MOI") or "")
            )
        )
        genes.append(_native_gene(gene, records, variants.get(gene, ())))
    return _tool_result(
        case_key,
        "exomiser",
        version,
        assembly,
        ordered[:gene_limit],
        genes,
        (gene_artifact, variant_artifact),
        issues,
    )


def _entrez_map(path, wanted):
    artifact = _artifact(path)
    rows, _ = _tsv_rows(artifact["content"], artifact["name"], _GENE_MAP_HEADERS)
    output = {}
    for _, row in rows:
        entrez = _entrez(row.get("entrez_id"))
        if entrez not in wanted:
            continue
        gene = _gene(row.get("symbol"))
        if not entrez or gene is None:
            continue
        if entrez in output and output[entrez]["gene"] != gene:
            raise _error("lirical.gene_map", "conflicting_entrez_mapping")
        output[entrez] = {
            "gene": gene,
            "source": Source(
                provider="hgnc",
                kind="annotation",
                locator=f"resource:hgnc/{artifact['name']}#entrez:{entrez}",
                version=f"sha256:{artifact['content_sha256']}",
                content=_semantic_json(row),
                genes=(gene,),
            ),
        }
    return output


def _lirical_json_records(content, artifact):
    try:
        payload = json.loads(content.removeprefix("\ufeff"))
    except json.JSONDecodeError:
        raise _error(f"{artifact}.json", "invalid_json") from None
    _check(isinstance(payload, dict), f"{artifact}.json", "object_required")
    raw_results = payload.get("analysisResults", ())
    _check(isinstance(raw_results, (list, tuple)), f"{artifact}.analysisResults", "array_required")
    output, issues = {}, []
    for number, item in enumerate(raw_results, 1):
        if not isinstance(item, dict):
            issues.append(_issue("invalid_json_record", artifact, number))
            continue
        disease_id = str(item.get("diseaseId") or "").strip()
        if not disease_id:
            issues.append(_issue("missing_disease_id", artifact, number))
            continue
        genotype = item.get("genotypeLR") if isinstance(item.get("genotypeLR"), dict) else {}
        gene_data = genotype.get("geneId") if isinstance(genotype.get("geneId"), dict) else {}
        record = {
            "item": item,
            "number": number,
            "gene_id": _entrez(gene_data.get("id")),
            "gene_symbol": str(gene_data.get("symbol") or "").strip().upper(),
            "match_type": str(genotype.get("matchType") or ""),
            "explanation": str(genotype.get("explanation") or ""),
        }
        previous = output.get(disease_id)
        if previous is not None:
            same_gene = (
                previous["gene_id"] == record["gene_id"]
                if previous["gene_id"] and record["gene_id"]
                else previous["gene_symbol"] == record["gene_symbol"]
            )
            _check(same_gene, "lirical.json.disease", "conflicting_disease_gene")
            issues.append(_issue("duplicate_json_disease", artifact, number))
        output[disease_id] = record
    return output, tuple(issues)


def _lirical_gene(row, json_record, entrez_map):
    tsv_id = _entrez(row.get("entrezGeneId"))
    json_id = json_record["gene_id"] if json_record is not None else ""
    _check(
        not (tsv_id and json_id and tsv_id != json_id),
        "lirical.gene_id",
        "tsv_json_gene_id_conflict",
    )
    mapping = entrez_map.get(tsv_id or json_id)
    mapped = mapping["gene"] if mapping is not None else ""
    json_symbol = json_record["gene_symbol"] if json_record is not None else ""
    issue = None
    if json_symbol:
        try:
            validate_gene(json_symbol, field="lirical.json.gene_symbol")
        except MendelPrimeError:
            json_symbol, issue = "", "invalid_json_gene"
    if mapped and json_symbol and mapped != json_symbol:
        issue = "gene_symbol_alias"
    return json_symbol or mapped or None, issue


def _lirical_disease_record(
    row,
    json_record,
    gene,
    number,
    case_key,
    version,
    tsv_artifact,
    json_artifact,
    mapping_source,
):
    json_source = None
    if json_record is not None:
        json_source = _source(
            "lirical",
            json_artifact,
            json_record["number"],
            version,
            case_key,
            gene,
            _semantic_json(json_record["item"]),
        )
    return {
        "source": _source("lirical", tsv_artifact, number, version, case_key, gene, _semantic_json(row)),
        "json_source": json_source,
        "gene_mapping_source": mapping_source,
        "rank": _int(row.get("rank")),
        "disease_id": _optional(row.get("diseaseCurie")),
        "disease_name": _optional(row.get("diseaseName")),
        "pretest_probability": _optional(row.get("pretestprob")),
        "posttest_probability": _probability(row.get("posttestprob")),
        "composite_lr": _float(row.get("compositeLR")),
        "genotype_match_type": _optional(json_record["match_type"]) if json_record else None,
        "genotype_explanation": _optional(json_record["explanation"]) if json_record else None,
    }


def _lirical_variants(text, disease_id, gene, number, case_key, assembly, version, artifact):
    parsed, rejected, issues = [], [], []
    for fragment_number, fragment in enumerate(re.split(r";\s*", text.strip()), 1):
        if not fragment:
            continue
        match = _LIRICAL_VARIANT.fullmatch(fragment.strip())
        if match is None:
            rejected.append(
                _source(
                    "lirical",
                    artifact,
                    number,
                    version,
                    case_key,
                    gene,
                    fragment,
                    fragment=fragment_number,
                )
            )
            issues.append(_issue("unparsed_variant_fragment", artifact, number, gene))
            continue
        chromosome = match.group("chromosome")
        position = int(match.group("position"))
        reference, alternate = match.group("reference"), match.group("alternate")
        key = _variant_identity(assembly, _canonical_chromosome(chromosome), position, reference, alternate)
        hgvs = _HGVS.search(match.group("annotation"))
        source = _source(
            "lirical",
            artifact,
            number,
            version,
            case_key,
            gene,
            fragment,
            key,
            fragment_number,
        )
        parsed.append(
            _variant_record(
                gene,
                key,
                chromosome,
                position,
                reference,
                alternate,
                hgvs,
                match.group("genotype"),
                {
                    "pathogenicity": _float(match.group("pathogenicity")),
                    "disease_id": disease_id,
                    "raw_annotation": fragment,
                },
                source,
            )
        )
        if key is None:
            issues.append(_issue("unsupported_variant_identity", artifact, number, gene))
    return parsed, rejected, issues


def parse_lirical(root, case_key, assembly, version, gene_map, gene_limit=100):
    """Parse all LIRICAL TSV/JSON records and native variant fragments."""
    _validate_parse_request(case_key, assembly, version, gene_limit)
    tsv_artifact = _artifact(root / "lirical.tsv")
    json_artifact = _artifact(root / "lirical.json")
    rows, malformed = _tsv_rows(
        tsv_artifact["content"], tsv_artifact["name"], _LIRICAL_HEADERS, skip_prefix="!"
    )
    issues = [_issue("malformed_row", tsv_artifact["name"], number) for number in malformed]
    json_records, json_issues = _lirical_json_records(json_artifact["content"], json_artifact["name"])
    issues.extend(json_issues)
    wanted = set()
    for _, row in rows:
        record = json_records.get(row.get("diseaseCurie") or "")
        wanted.add(_entrez(row.get("entrezGeneId")) or (record["gene_id"] if record else ""))
    entrez_map = _entrez_map(gene_map, wanted - {""})
    grouped, disease_genes = defaultdict(list), {}
    for number, row in rows:
        disease_id = _optional(row.get("diseaseCurie"))
        if disease_id is None:
            issues.append(_issue("missing_disease_id", tsv_artifact["name"], number))
            continue
        json_record = json_records.get(disease_id)
        gene, gene_issue = _lirical_gene(row, json_record, entrez_map)
        if gene_issue is not None and json_record is not None:
            issues.append(_issue(gene_issue, json_artifact["name"], json_record["number"]))
        if gene is None:
            issues.append(_issue("unmapped_gene", tsv_artifact["name"], number))
            continue
        previous = disease_genes.setdefault(disease_id, gene)
        _check(previous == gene, "lirical.disease_gene", "disease_gene_conflict")
        grouped[gene].append((number, row, json_record))
        if _int(row.get("rank")) is None:
            issues.append(_issue("invalid_rank", tsv_artifact["name"], number, gene))
    _check(bool(grouped), "lirical.ranking", "empty_ranking")
    ordered = tuple(
        sorted(
            grouped,
            key=lambda gene: (
                min(_rank(row.get("rank")) for _, row, _ in grouped[gene]),
                gene,
            ),
        )
    )
    genes = []
    for gene in ordered:
        records = sorted(
            grouped[gene],
            key=lambda item: (_rank(item[1].get("rank")), item[1].get("diseaseCurie") or ""),
        )
        disease_records, variants, unparsed = [], [], []
        for number, row, json_record in records:
            identifier = _entrez(row.get("entrezGeneId")) or (json_record["gene_id"] if json_record else "")
            mapping = entrez_map.get(identifier)
            mapping_source = (
                mapping["source"]
                if mapping and (json_record is None or _gene(json_record["gene_symbol"]) is None)
                else None
            )
            disease_records.append(
                _lirical_disease_record(
                    row,
                    json_record,
                    gene,
                    number,
                    case_key,
                    version,
                    tsv_artifact["name"],
                    json_artifact["name"],
                    mapping_source,
                )
            )
            parsed, rejected, variant_issues = _lirical_variants(
                row.get("variants") or "",
                row.get("diseaseCurie") or "",
                gene,
                number,
                case_key,
                assembly,
                version,
                tsv_artifact["name"],
            )
            variants.extend(parsed)
            unparsed.extend(rejected)
            issues.extend(variant_issues)
        genes.append(_native_gene(gene, disease_records, variants, unparsed))
    return _tool_result(
        case_key,
        "lirical",
        version,
        assembly,
        ordered[:gene_limit],
        genes,
        (tsv_artifact, json_artifact),
        issues,
    )


def _fresh_private_directory(directory):
    _check(
        isinstance(directory, Path) and directory.is_absolute(),
        "directory",
        "absolute_path_required",
        "contract",
    )
    try:
        directory.mkdir(mode=0o700, parents=True, exist_ok=False)
        directory.chmod(0o700)
    except FileExistsError:
        raise _error("directory", "already_exists", "execution") from None
    except OSError:
        raise _error("directory", "private_directory_creation_failed", "execution") from None


def _selected_samples(family):
    roles = family["roles"]
    return tuple(sample for sample in (roles["proband"], roles.get("mother"), roles.get("father")) if sample)


def _pedigree_lines(family):
    pedigree = family.get("pedigree")
    supplied = {item["sample_id"]: item for item in pedigree["members"]} if pedigree else {}
    selected = _selected_samples(family)
    proband_member = supplied.get(family["roles"]["proband"])
    selected_member = next((supplied[sample] for sample in selected if sample in supplied), None)
    family_member = proband_member or selected_member
    family_id = family_member["family_id"] if family_member else "family"
    lines = []
    for sample in selected:
        member = supplied.get(sample)
        is_proband = sample == family["roles"]["proband"]
        father = family["roles"].get("father") if is_proband else None
        mother = family["roles"].get("mother") if is_proband else None
        sex = family["sex"] if is_proband else member["sex"] if member else "UNKNOWN"
        affected = member["affected"] if member else 2 if is_proband else 0
        lines.append(
            "\t".join(
                (
                    member["family_id"] if member else family_id,
                    sample,
                    father or "0",
                    mother or "0",
                    _SEX_CODES[sex],
                    str(affected),
                )
            )
            + "\n"
        )
    return tuple(lines)


def prepare_input(patient, family, directory):
    """Create private effective VCF/PED files without modifying source inputs."""
    _check(
        isinstance(patient, Patient) and patient.genomic is not None,
        "patient",
        "genomic_patient_required",
        "contract",
    )
    _fresh_private_directory(directory)
    vcf, ped = directory / "effective.vcf", directory / "effective.ped"
    materialize_samples(patient.genomic["vcf_path"], vcf, _selected_samples(family))
    try:
        with ped.open("x", encoding="utf-8", newline="\n") as handle:
            handle.writelines(_pedigree_lines(family))
        ped.chmod(0o600)
        vcf.chmod(0o600)
    except FileExistsError:
        raise _error("pedigree", "already_exists", "execution") from None
    except OSError:
        raise _error("pedigree", "materialization_failed", "execution") from None
    return {"vcf": vcf, "ped": ped}


def exomiser_analysis(patient, family, prepared):
    """Return the Exomiser 14 analysis configuration as structured YAML data."""
    return {
        "analysis": {
            "genomeAssembly": "hg19" if family["assembly"] == "GRCh37" else "hg38",
            "vcf": str(prepared["vcf"].resolve()),
            "ped": str(prepared["ped"].resolve()),
            "proband": family["roles"]["proband"],
            "hpoIds": [term["hpo_id"] for term in patient.phenotypes if term["present"]],
            "inheritanceModes": dict(_INHERITANCE_MODES),
            "analysisMode": "PASS_ONLY",
            "frequencySources": list(_FREQUENCY_SOURCES),
            "pathogenicitySources": ["REVEL", "MVP"],
            "steps": [
                {"failedVariantFilter": {}},
                {"variantEffectFilter": {"remove": list(_REMOVED_EFFECTS)}},
                {"frequencyFilter": {"maxFrequency": 2}},
                {"pathogenicityFilter": {"keepNonPathogenic": True}},
                {"inheritanceFilter": {}},
                {"omimPrioritiser": {}},
                {"hiPhivePrioritiser": {}},
            ],
        },
        "outputOptions": {
            "outputContributingVariantsOnly": False,
            "numGenes": 100,
            "outputDirectory": ".",
            "outputFileName": "exomiser",
            "outputFormats": ["TSV_GENE", "TSV_VARIANT"],
        },
    }


def _assembly_resources(config, assembly):
    _check(isinstance(config, Mapping), "genotype", "mapping_required", "resource")
    configured = config.get("assemblies", ())
    selected = next((item for item in configured if item.get("assembly") == assembly), None)
    _check(selected is not None, "genotype.assemblies", "assembly_resources_not_configured", "resource")
    return selected


def _available(path, field, *, directory=False, executable=False):
    _check(
        isinstance(path, Path) and path.is_absolute(),
        field,
        "absolute_path_required",
        "resource",
    )
    try:
        present = path.is_dir() if directory else path.is_file() and path.stat().st_size > 0
        mode = os.R_OK | os.X_OK if directory or executable else os.R_OK
        accessible = os.access(path, mode)
    except OSError:
        present = accessible = False
    code = "readable_directory_required" if directory else "readable_nonempty_file_required"
    _check(present and accessible, field, code, "resource")
    if executable:
        _check(os.access(path, os.X_OK), field, "executable_required", "resource")


def _validate_resources(config, family, tool_id):
    selected = _assembly_resources(config, family["assembly"])
    _available(config.get("java"), "genotype.java", executable=True)
    if tool_id == "exomiser":
        files = ((config, "exomiser_jar"), (config, "exomiser_properties"))
        directories = ((selected, "exomiser_genome_data"), (selected, "exomiser_phenotype_data"))
    else:
        files = (
            (config, "lirical_jar"),
            (config, "gene_map"),
            (selected, "lirical_variants_db"),
            (selected, "lirical_clinvar_db"),
        )
        directories = ((config, "lirical_data"), (selected, "lirical_exomiser_data"))
    for values, name in files:
        _available(values.get(name), f"genotype.{name}")
    for values, name in directories:
        _available(values.get(name), f"genotype.{name}", directory=True)
    check_java_distribution(config[f"{tool_id}_jar"], config[f"{tool_id}_version"])
    return selected


def exomiser_command(analysis, config, assembly=None):
    if assembly is None:
        _check(len(config["assemblies"]) == 1, "genotype.assemblies", "assembly_required", "resource")
        assembly = config["assemblies"][0]["assembly"]
    selected = _assembly_resources(config, assembly)
    build = "hg19" if assembly == "GRCh37" else "hg38"
    return [
        str(config["java"]),
        f"-Xmx{config['max_heap'].lower()}",
        "-jar",
        str(config["exomiser_jar"]),
        "--analysis",
        str(analysis),
        f"--spring.config.location={config['exomiser_properties']}",
        # Spring command-line properties override stale paths in the supplied template.
        f"--exomiser.data-directory={selected['exomiser_genome_data'].parent}",
        f"--exomiser.{build}.data-directory={selected['exomiser_genome_data']}",
        f"--exomiser.{build}.data-version={selected['exomiser_release']}",
        f"--exomiser.phenotype.data-directory={selected['exomiser_phenotype_data']}",
        f"--exomiser.phenotype.data-version={selected['phenotype_release']}",
    ]


def lirical_command(patient, family, prepared, root, config):
    selected = _assembly_resources(config, family["assembly"])
    suffix = "19" if family["assembly"] == "GRCh37" else "38"
    assembly = "hg19" if family["assembly"] == "GRCh37" else "hg38"
    present = ",".join(term["hpo_id"] for term in patient.phenotypes if term["present"])
    command = [
        str(config["java"]),
        f"-Xmx{config['max_heap'].lower()}",
        "-jar",
        str(config["lirical_jar"]),
        "prioritize",
        "-d",
        str(config["lirical_data"]),
        "--assembly",
        assembly,
        "--vcf",
        str(prepared["vcf"]),
        "--sample-id",
        family["roles"]["proband"],
        "-p",
        present,
        "--sex",
        family["sex"],
        "--transcript-db",
        config["lirical_transcript_db"],
        f"-ed{suffix}",
        str(selected["lirical_exomiser_data"]),
        f"-e{suffix}",
        str(selected["lirical_variants_db"]),
        f"-c{suffix}",
        str(selected["lirical_clinvar_db"]),
        "-o",
        str(root),
        "-x",
        "lirical",
        "-f",
        "tsv",
        "-f",
        "json",
        "--display-all-variants",
    ]
    absent = ",".join(term["hpo_id"] for term in patient.phenotypes if not term["present"])
    if absent:
        command.extend(("-n", absent))
    return command


def _command(tool_id, patient, family, prepared, root, config):
    if tool_id == "exomiser":
        analysis = root / "analysis.yml"
        with analysis.open("x", encoding="utf-8", newline="\n") as handle:
            yaml.safe_dump(
                exomiser_analysis(patient, family, prepared),
                handle,
                allow_unicode=False,
                explicit_start=True,
                sort_keys=False,
            )
        analysis.chmod(0o600)
        return exomiser_command(analysis, config, family["assembly"])
    return lirical_command(patient, family, prepared, root, config)


def _fingerprint(path, field):
    try:
        with open_regular(path) as (_, status):
            return stat_fingerprint(status)
    except OSError:
        raise _error(field, "source_unavailable_during_analysis", "input") from None


def _input_fingerprints(patient):
    genomic = patient.genomic
    inputs = [("vcf", genomic["vcf_path"])]
    if genomic.get("pedigree_path") is not None:
        inputs.append(("pedigree", genomic["pedigree_path"]))
    return tuple((field, path, _fingerprint(path, field)) for field, path in inputs)


def _require_inputs_unchanged(fingerprints):
    for field, path, expected in fingerprints:
        try:
            current = _fingerprint(path, field)
        except MendelPrimeError:
            current = None
        if current != expected:
            raise _error(field, "source_changed_during_analysis", "input")


def _new_workspace(work_root, case_key):
    _check(
        isinstance(work_root, Path) and work_root.is_absolute(),
        "work_root",
        "absolute_path_required",
        "contract",
    )
    try:
        private_directory(work_root)
    except (OSError, TypeError):
        raise _error("work_root", "private_directory_required", "execution") from None
    digest = hashlib.sha256(f"{case_key[0]}\0{case_key[1]}".encode()).hexdigest()[:24]
    workspace = work_root / f"{digest}-{uuid.uuid4().hex}"
    try:
        os.mkdir(workspace, 0o700)
        private_directory(workspace, create=False)
    except FileExistsError:
        raise _error("work_root", "workspace_already_exists", "execution") from None
    except (OSError, TypeError):
        raise _error("work_root", "private_workspace_creation_failed", "execution") from None
    return workspace


def _tool_sources(tool):
    sources = []
    for gene in tool["genes"]:
        for record in gene["records"]:
            sources.extend(
                record[name]
                for name in ("source", "json_source", "gene_mapping_source")
                if record.get(name) is not None
            )
        sources.extend(variant["source"] for variant in gene["variants"])
        sources.extend(gene.get("unparsed_variants", ()))
    return tuple(sources)


def family_sources_for(tools, observations):
    """Bind exact family observations to each native gene/allele scope."""
    by_key = {item["variant_key"]: item for item in observations["observations"]}
    scoped = {
        (gene["gene"], variant["key"]): by_key[variant["key"]]
        for tool in tools
        for gene in tool["genes"]
        for variant in gene["variants"]
        if variant["key"] in by_key
    }
    grouped = {}
    for (gene, _), observation in scoped.items():
        grouped.setdefault((gene, observation["vcf_variant_key"]), []).append(observation)
    sources = []
    for (gene, raw_key), family_observations in grouped.items():
        first = family_observations[0]
        content = {
            name: first[name]
            for name in (
                "alternate_index",
                "samples",
                "quality",
                "filter_text",
                "raw_record",
                "vcf_variant_key",
            )
        }
        content["vcf_variant_key"] = dict(
            zip(
                ("assembly", "chromosome", "position", "reference", "alternate"),
                first["vcf_variant_key"],
                strict=True,
            )
        )
        content["target_mappings"] = [
            {
                "variant_key": ":".join(map(str, item["variant_key"])),
                "method": item["match_method"],
            }
            for item in family_observations
        ]
        keys = tuple(dict.fromkeys((raw_key, *(item["variant_key"] for item in family_observations))))
        sources.append(
            Source(
                provider="vcf",
                kind="patient_observation",
                locator=f"input:vcf#allele/{':'.join(map(str, raw_key))}/{gene}",
                version="family-observation-v1",
                genes=(gene,),
                case_key=observations["case_key"],
                variant_keys=keys,
                content=json.dumps(content, ensure_ascii=False, sort_keys=True),
            )
        )
    return tuple(sources)


def _execute_tool(patient, family, prepared, workspace, config, tool_id, fingerprints):
    version = config[f"{tool_id}_version"]
    root = workspace / tool_id
    artifacts = tuple(root / name for name in ARTIFACTS[tool_id])
    LOGGER.info("Native genotype tool start tool=%s case_key=%s/%s", tool_id, *patient.key)
    started = time.monotonic()
    try:
        try:
            _validate_resources(config, family, tool_id)
        except MendelPrimeError as error:
            outcome = {"error_code": error.code, "seconds": 0.0, "returncode": None}
        else:
            _fresh_private_directory(root)
            outcome = run_process(
                command=_command(tool_id, patient, family, prepared, root, config),
                root=root,
                timeout_seconds=config["timeout_seconds"],
                artifacts=artifacts,
            )
        outcome = outcome or {"error_code": "", "seconds": 0.0, "returncode": None}
        code, seconds, returncode = (outcome[name] for name in ("error_code", "seconds", "returncode"))
    except Exception:
        code, seconds, returncode = "native_execution_failed", round(time.monotonic() - started, 6), None
    execution = {
        "tool_id": tool_id,
        "version": version,
        "status": "failed" if code else "ok",
        "error_code": code,
        "seconds": seconds,
        "returncode": returncode,
        "artifacts": artifacts,
    }
    LOGGER.info(
        "Native genotype tool end tool=%s status=%s duration_seconds=%.6f case_key=%s/%s",
        tool_id,
        execution["status"],
        execution["seconds"],
        *patient.key,
    )
    _require_inputs_unchanged(fingerprints)
    if execution["status"] != "ok":
        return (
            execution,
            None,
            {
                "tool_id": tool_id,
                "code": execution["error_code"] or "native_execution_failed",
            },
        )
    try:
        arguments = root, patient.key, family["assembly"], version
        tool = (
            parse_exomiser(*arguments)
            if tool_id == "exomiser"
            else parse_lirical(*arguments, config["gene_map"])
        )
        valid = (
            isinstance(tool, dict)
            and tool.get("case_key") == patient.key
            and tool.get("tool_id") == tool_id
            and tool.get("version") == version
            and tool.get("assembly") == family["assembly"]
        )
        _check(valid, "native.result", "source_identity_mismatch", "response")
        return execution, tool, None
    except Exception:
        return execution, None, {"tool_id": tool_id, "code": "native_parse_failed"}


def acquire(patient, genotype_config, work_root):
    """Run Exomiser then LIRICAL with source-local failures and return plain records."""
    _check(
        isinstance(patient, Patient) and patient.genomic is not None,
        "patient",
        "genomic_patient_required",
        "contract",
    )
    _check(
        any(term["present"] for term in patient.phenotypes), "phenotypes", "positive_hpo_required", "input"
    )
    fingerprints = _input_fingerprints(patient)
    header = inspect_vcf(patient.genomic["vcf_path"])
    pedigree = (
        read_pedigree(patient.genomic["pedigree_path"])
        if patient.genomic.get("pedigree_path") is not None
        else None
    )
    family = resolve_family(patient, header, pedigree)
    _assembly_resources(genotype_config, family["assembly"])
    workspace = _new_workspace(work_root, patient.key)
    prepared = prepare_input(patient, family, workspace / "input")
    _require_inputs_unchanged(fingerprints)
    tools, executions, failures = [], [], []
    for tool_id in ARTIFACTS:
        execution, tool, failure = _execute_tool(
            patient, family, prepared, workspace, genotype_config, tool_id, fingerprints
        )
        executions.append(execution)
        if tool is not None:
            tools.append(tool)
        if failure is not None:
            failures.append(failure)
    targets = tuple(
        dict.fromkeys(
            variant["key"]
            for tool in tools
            for gene in tool["genes"]
            for variant in gene["variants"]
            if variant["key"] is not None
        )
    )
    observations = read_observations(patient.genomic["vcf_path"], family, targets, case_key=patient.key)
    _require_inputs_unchanged(fingerprints)
    family_sources = family_sources_for(tuple(tools), observations)
    all_sources = tuple(
        {
            source.source_id: source
            for source in (
                *(source for tool in tools for source in _tool_sources(tool)),
                *family_sources,
            )
        }.values()
    )
    evidence = {
        "case_key": patient.key,
        "tools": tuple(tools),
        "observations": observations,
        "family_sources": family_sources,
        "sources": all_sources,
    }
    degraded = (
        failures
        or any(item["status"] != "ok" for item in executions)
        or any(tool["issues"] for tool in tools)
    )
    return {
        "family": family,
        "evidence": evidence,
        "executions": tuple(executions),
        "failures": tuple(failures),
        "status": "degraded" if degraded else "ok",
    }


__all__ = [
    "acquire",
    "exomiser_analysis",
    "exomiser_command",
    "lirical_command",
    "parse_exomiser",
    "parse_lirical",
    "prepare_input",
    "variant_is_candidate",
]
