"""Shared source-bound VCF preparation resources for ranking and optional reports."""

import configparser
import csv
import gzip
import hashlib
import io
import json
import os
import re
import sqlite3
import stat
import zlib
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
from datetime import date
from pathlib import Path, PurePosixPath, PureWindowsPath
from threading import RLock

from mendelprime.domain import MendelPrimeError, Source, validate_gene, variant_key, variant_text
from mendelprime.runtime.storage import private_directory, read_regular_bounded, stat_fingerprint
from mendelprime.tools.pvs1_driver import _consequence, preflight_source, resource_gene_aliases, run_pvs1

from . import literature
from .annotation import annotate, merged_refseq_annotation, population, predictions, public_resource
from .reference import ReferenceGenome

_MAX_MANIFEST = 1024 * 1024
_MAX_GUIDELINE = 32 * 1024 * 1024
_MAX_COMPRESSED = 16 * 1024 * 1024
_MAX_TEXT = 128 * 1024 * 1024
_SHA256 = re.compile(r"[0-9a-f]{64}", re.ASCII)
_GENE = re.compile(r"[A-Z0-9][A-Z0-9._-]*", re.ASCII)
_CODON = re.compile(r"(?:p\.\(?)?[A-Z][a-z]{2}(\d+)[A-Z][a-z]{2}")
_DOI = "10.1126/scitranslmed.adz4172"
_DOI_URL = f"https://doi.org/{_DOI}"
_EXPECTED_PROTOCOL = "bb12f60eb217a996bb5cfccb18a3241cacaaf993e2e16d6dafc55ff90d4cebc8"
_EXPECTED_MANIFEST = "04da36613e0d5ad47c034885326421eef71a228bda61781009206661384cd152"
# The public revision changes provenance paths, not the protocol or table data.
_PUBLIC_MANIFEST = "1d8954053c619fedd076ae399305fe815bac8e405b5fb3eb13272ce085840527"
_MANIFEST_SCHEMA = "mgagent_variant_knowledge_manifest_v1"
_KNOWLEDGEBASE_SCHEMA = "ai_cura_knowledgebase_v1"
_GENERIC_S8 = "Any"
_GENERIC_S10 = {
    "Genes not listed above with loss-of-function being the molecular mechanism",
    "Genes not listed above with associated gene being an enzyme/subunit of enzyme/transcriptional factor/transporter",
}
_KNOWLEDGEBASE_KEYS = {
    "schema_version",
    "source",
    "gene_disease_consistency",
    "pp1_recommendations",
    "pp4_knowledgebase",
    "ps4_thresholds",
    "ps2_pm6_strength_thresholds",
}
_CONTEXT_HEADERS = {
    "mane": (
        "NCBI_GeneID",
        "Ensembl_Gene",
        "HGNC_ID",
        "symbol",
        "name",
        "RefSeq_nuc",
        "RefSeq_prot",
        "Ensembl_nuc",
        "Ensembl_prot",
        "MANE_status",
        "GRCh38_chr",
        "chr_start",
        "chr_end",
        "chr_strand",
    ),
    "gnomad-constraint": (
        "gene",
        "transcript",
        "obs_mis",
        "exp_mis",
        "oe_mis",
        "mu_mis",
        "possible_mis",
        "obs_mis_pphen",
        "exp_mis_pphen",
        "oe_mis_pphen",
        "possible_mis_pphen",
        "obs_syn",
        "exp_syn",
        "oe_syn",
        "mu_syn",
        "possible_syn",
        "obs_lof",
        "mu_lof",
        "possible_lof",
        "exp_lof",
        "pLI",
        "pNull",
        "pRec",
        "oe_lof",
        "oe_syn_lower",
        "oe_syn_upper",
        "oe_mis_lower",
        "oe_mis_upper",
        "oe_lof_lower",
        "oe_lof_upper",
        "constraint_flag",
        "syn_z",
        "mis_z",
        "lof_z",
        "oe_lof_upper_rank",
        "oe_lof_upper_bin",
        "oe_lof_upper_bin_6",
        "n_sites",
        "classic_caf",
        "max_af",
        "no_lofs",
        "obs_het_lof",
        "obs_hom_lof",
        "defined",
        "p",
        "exp_hom_lof",
        "classic_caf_afr",
        "classic_caf_amr",
        "classic_caf_asj",
        "classic_caf_eas",
        "classic_caf_fin",
        "classic_caf_nfe",
        "classic_caf_oth",
        "classic_caf_sas",
        "p_afr",
        "p_amr",
        "p_asj",
        "p_eas",
        "p_fin",
        "p_nfe",
        "p_oth",
        "p_sas",
        "transcript_type",
        "gene_id",
        "transcript_level",
        "cds_length",
        "num_coding_exons",
        "gene_type",
        "gene_length",
        "exac_pLI",
        "exac_obs_lof",
        "exac_exp_lof",
        "exac_oe_lof",
        "brain_expression",
        "chromosome",
        "start_position",
        "end_position",
    ),
}
_CONTEXT_FIELD = {"mane": "symbol", "gnomad-constraint": "gene"}
_CONTEXT_URL = {
    "mane": "https://www.ncbi.nlm.nih.gov/refseq/MANE/",
    "gnomad-constraint": "https://gnomad.broadinstitute.org/",
}
_CONTEXT_SCOPE = {
    "mane": "Transcript identity context. GRCh38 coordinates are not a GRCh37 projection.",
    "gnomad-constraint": "Gene/transcript-level constraint, not an exact-allele pathogenicity assertion.",
}
_PROTEIN_COLUMNS = (
    "gene",
    "transcript",
    "coding_hgvs",
    "protein_hgvs",
    "codon",
    "variation_id",
    "clinical_significance",
    "review_status",
    "phenotype_ids",
    "phenotype_list",
)


_REPORT_FIELDS = {
    "mane_summary", "mane_release", "ai_cura_protocol", "ai_cura_workbook",
    "clingen_guideline_manifest", "clingen_dosage_table", "genereviews_index",
    "post_evidence_catalog", "kinship_executable",
}


def resource_settings(settings, *, report):
    """Exclude inactive resources from validation and ranking cache identities."""
    return {
        name: value for name, value in settings.items()
        if (report or name not in _REPORT_FIELDS)
        and (report or settings.get("ranking_vep_enabled", False) or not name.startswith("vep_"))
    }


def open_resources(settings_dict, clinical_callback, work_root, *, purpose="report", workers=4):
    """Validate resources once and return the nine reader actions plus ``close``."""
    if purpose not in {"ranking", "report"}:
        raise MendelPrimeError("ranking_or_report_required", field="assessment_resources.purpose")
    settings = settings_dict.get("explanation_resources", settings_dict)
    if not isinstance(settings, dict):
        raise MendelPrimeError("settings_object_required", field="explanation_resources")
    reporting = purpose == "report"
    settings = resource_settings(settings, report=reporting)
    if not callable(clinical_callback):
        raise MendelPrimeError("variant_knowledge_required", field="clinical")
    if not isinstance(work_root, Path) or not work_root.is_absolute():
        raise MendelPrimeError("absolute_path_required", field="work_root")
    try:
        private_directory(work_root, create=False)
    except (OSError, TypeError):
        raise MendelPrimeError("private_directory_required", field="work_root") from None
    ranking_vep = settings.get("ranking_vep_enabled", False)
    if type(ranking_vep) is not bool:
        raise MendelPrimeError("boolean_required", field="explanation_resources.ranking_vep_enabled")
    _check_resources(settings, vep_enabled=reporting or ranking_vep, report=reporting)
    if purpose == "report":
        for name in (
            "clingen_dosage_table",
            "kinship_executable",
            "genereviews_index",
            "post_evidence_catalog",
        ):
            if settings.get(name) is not None:
                _require_file(settings[name], f"explanation_resources.{name}")
        if settings.get("kinship_executable") and not os.access(settings["kinship_executable"], os.X_OK):
            raise MendelPrimeError("executable_required", field="explanation_resources.kinship_executable")
    protocol = _load_protocol(settings, report=reporting)
    clingen = _load_clingen_guidance(settings) if reporting else ()
    reference = ReferenceGenome(settings["reference_fasta"], assembly=settings["assembly"])
    try:
        context_state = {
            "resources": {
                **({"mane": (settings["mane_summary"], settings["mane_release"])} if reporting else {}),
                "gnomad-constraint": (
                    settings["gnomad_constraint"],
                    settings["constraint_release"],
                ),
            },
            "tables": {},
            "lock": RLock(),
        }

        reference_reader = {"resolve": reference.resolve, "equivalent": reference.equivalent}

        def clinical(gene, keys):
            return clinical_callback(gene, keys, reference=reference_reader)

        def guidance(genes, *, ranking=False):
            selected = set(genes)
            sources = [
                source for source in protocol["ranking_rows" if ranking else "rows"]
                if source.source_id in protocol["generic_source_ids"] or selected.intersection(source.genes)
            ]
            if reporting and not ranking:
                sources = [protocol["protocol_source"], *sources, *clingen]
            return tuple({source.source_id: source for source in sources}.values())

        actions = {
            "resolve": reference.resolve,
            "annotate": lambda targets: annotate(settings, work_root, targets),
            "population": lambda gene, key: population(settings, reference, gene, key),
            "predictions": lambda gene, key: predictions(settings, gene, key),
            "pvs1": lambda targets: run_pvs1(settings, work_root, targets),
            "gene_context": lambda gene: gene_context(context_state, gene),
            "protein_context": lambda gene, **kwargs: protein_context(settings, gene, **kwargs),
            "clinical": clinical,
            "guidance": guidance,
        }
        ranking_actions = {
            **actions,
            "annotate": lambda targets: annotate(settings, work_root, targets, enabled=ranking_vep),
            "gene_context": lambda gene: gene_context(context_state, gene, ranking=True),
            "guidance": lambda genes: guidance(genes, ranking=True),
            "resource_identity": lambda: _resource_identity(
                resource_settings(settings, report=False), vep_enabled=ranking_vep
            ),
        }
        if purpose == "report":
            post_state = {
                "settings": settings,
                "root": work_root,
                "reference": reference,
                "pubmed": literature.new_pubmed(work_root.parent / "pubmed"),
                "workers": workers,
                "tables": {},
                "lock": RLock(),
            }
            actions["annotate"] = lambda targets: _report_annotate(post_state, targets)
            actions["pvs1"] = lambda targets: _report_pvs1(post_state, targets)
            actions["population"] = lambda gene, key: population(settings, reference, gene, key, details=True)
            actions["complete_report_evidence"] = lambda seed, review, selected: _complete_report_evidence(
                post_state, seed, review, selected
            )
        return {
            "actions": actions if purpose == "report" else ranking_actions,
            "ranking_actions": ranking_actions,
            "reference": reference,
        }
    except BaseException:
        reference.close()
        raise


def _report_transcripts(state):
    """Share the exact-version transcript registry across report annotation and PVS1."""
    settings = state["settings"]
    with state["lock"]:
        if "transcripts" not in state["tables"]:
            config = configparser.ConfigParser()
            config.read(settings["autopvs1_root"] / "config.ini")
            state["tables"]["gene_aliases"] = resource_gene_aliases(settings["autopvs1_root"])
            aliases = state["tables"]["gene_aliases"][0]
            path = Path(config["HG19"]["transcript"])
            if not path.is_absolute():
                path = settings["autopvs1_root"] / path
            transcripts = {}
            for line in path.read_text().splitlines():
                fields = line.split("\t")
                if len(fields) >= 13 and fields[1].startswith("NM_"):
                    address = (
                        aliases.get(fields[12], fields[12]),
                        fields[1].split(".")[0],
                        fields[2].removeprefix("chr"),
                    )
                    transcripts.setdefault(address, set()).add(fields[1])
            state["tables"]["transcripts"] = {k: tuple(sorted(v)) for k, v in transcripts.items()}
        return state["tables"]["transcripts"]


def _report_annotate(state, targets):
    _report_transcripts(state)
    rows = annotate(state["settings"], state["root"], targets, gene_aliases=state["tables"]["gene_aliases"][0])

    def resolve(row):
        target = row["target"]
        if row.get("error_code") not in {"vep_exact_transcript_not_found", "vep_gene_not_found"} or not str(
            target["transcript_id"]
        ).startswith("NM_"):
            return row
        try:
            annotation = merged_refseq_annotation(
                state["settings"],
                state["root"],
                {"gene": target["gene"], "key": target["variant_key"], "transcript": target["transcript_id"]},
                _report_transcripts(state),
                gene_aliases=state["tables"].get("gene_aliases", ({}, "")),
            )
        except Exception as error:
            return {
                **row,
                "error_code": error.code
                if isinstance(error, MendelPrimeError)
                else "post_refseq_query_failed",
            }
        return {**row, "status": "annotated", "error_code": "", "annotation": annotation}

    with ThreadPoolExecutor(max_workers=min(state["workers"], len(rows)) or 1) as pool:
        return tuple(pool.map(resolve, rows))


def _report_pvs1(state, targets):
    """Keep the ranking adapter unchanged; resolve report-only RefSeq representations."""
    settings, root = state["settings"], state["root"]
    if not any("aliases" in target for target in targets):
        return run_pvs1(settings, root, targets)
    transcripts = _report_transcripts(state)

    def resolve(target):
        representations = (target, *target.get("aliases", ()))
        preferred = target.get("preferred_transcripts", ())
        exact = {
            (row.get("transcript"), row.get("hgvs_c"), row.get("hgvs_p")): row
            for row in representations
            if row.get("transcript")
            in transcripts.get(
                (target["gene"], str(row.get("transcript", "")).split(".")[0], target["key"][1]), ()
            )
            and re.fullmatch(r"c\.\d+[ACGT]>[ACGT]", row.get("hgvs_c") or "")
            and re.fullmatch(r"p\.\(?[A-Z](?:[a-z]{2})?\d+(?:Ter|\*)\)?", row.get("hgvs_p") or "")
            and len(target["key"][3]) == len(target["key"][4]) == 1
        }
        if len(exact) == 1 and not preferred:
            row = next(iter(exact.values()))
            return (
                {
                    **row,
                    "details": {**row.get("details", {}), "consequence": "stop_gained"},
                    "verify_native_snv": True,
                    "retain_context": True,
                },
                None,
                "",
            )
        possible = any(
            _consequence(str(v.get("details", {}).get("consequence") or ""))
            or re.search(r"fs|Ter|\*", str(v.get("hgvs_p") or ""))
            or re.search(r"[+-][12][ACGT]", str(v.get("hgvs_c") or ""))
            for v in representations
        )
        if not possible:
            return target, None, ""
        if not preferred and target.get("mane_candidates"):
            candidates = target["mane_candidates"]
            if len(candidates) != 1:
                return None, None, "post_refseq_missing_or_ambiguous"
            preferred = tuple(candidates)
        try:
            annotation = merged_refseq_annotation(
                settings,
                root,
                target,
                transcripts,
                gene_aliases=state["tables"].get("gene_aliases", ({}, "")),
                **({"preferred": preferred} if preferred else {}),
            )
            consequence = "&".join(annotation["consequence_terms"])
            prepared = {
                **target,
                "retain_context": True,
                "transcript": annotation["transcript_id"],
                "hgvs_c": annotation["hgvsc"].split(":", 1)[-1],
                "hgvs_p": annotation["hgvsp"].split(":", 1)[-1],
                "source": annotation["source"],
                "details": {**target.get("details", {}), "consequence": consequence},
                "verify_native_snv": consequence == "stop_gained"
                and len(target["key"][3]) == len(target["key"][4]) == 1,
            }
            return prepared, annotation, ""
        except Exception as error:
            return (
                None,
                None,
                error.code if isinstance(error, MendelPrimeError) else "post_refseq_query_failed",
            )

    with ThreadPoolExecutor(max_workers=min(state["workers"], len(targets)) or 1) as pool:
        resolved = tuple(pool.map(resolve, targets))
    assessed = iter(run_pvs1(settings, root, tuple(row[0] for row in resolved if row[0] is not None)))
    output = []
    for original, (prepared, annotation, error) in zip(targets, resolved, strict=True):
        if prepared is None:
            output.append({"target": original, "status": "unknown", "error_code": error, "annotation": None})
            continue
        outcome = next(assessed)
        value = outcome.get("annotation")
        if value and value.get("native_hgvs_verification") and annotation is None:
            proof = value["native_hgvs_verification"]
            source = Source(
                provider="local-hgvs",
                kind="annotation",
                locator=f"method:exact-refseq-hgvs/{variant_text(original['key'])}",
                version="pyhgvs+local-codon-v1;" + str(value.get("algorithm_assets_sha256", "")),
                content=json.dumps(
                    {"gene": original["gene"], "native_source_id": prepared["source"].source_id, **proof},
                    sort_keys=True,
                ),
                genes=(original["gene"],),
                variant_keys=(original["key"],),
                case_key=prepared["source"].case_key,
            )
            annotation = {
                "transcript_id": prepared["transcript"],
                "hgvsc": prepared["transcript"] + ":" + prepared["hgvs_c"],
                "hgvsp": prepared["hgvs_p"],
                "consequence_terms": ("stop_gained",),
                "source": source,
                "protein_id": None,
                "exon": None,
                "intron": None,
                "canonical": False,
            }
        if value is not None and annotation is not None:
            value = {
                **value,
                "resolved_annotation": annotation,
                "observed_transcript": original.get("transcript"),
            }
        output.append({**outcome, "target": original, "annotation": value, "resolved_annotation": annotation})
    return tuple(output)


def _genereviews_documents(state, gene):
    """Read the official gene/chapter index and original mechanism sections."""
    from bs4 import BeautifulSoup

    path = state["settings"].get("genereviews_index")
    if not path:
        return ()
    with state["lock"]:
        if "genereviews" not in state["tables"]:
            raw = read_regular_bounded(path, _MAX_MANIFEST)[0].decode("utf-8")
            rows = csv.DictReader(io.StringIO(raw.lstrip("#")), delimiter="\t")
            index = {}
            for row in rows:
                if re.fullmatch(r"NBK\d+", row["NBK_id"]):
                    index.setdefault(row["genesymbol"], []).append((row["NBK_id"], row["GR_shortname"]))
            state["tables"]["genereviews"] = index
    documents = []
    for accession, shortname in state["tables"]["genereviews"].get(gene, ()):
        url = f"https://www.ncbi.nlm.nih.gov/sites/books/{accession}/?report=reader"
        raw = public_resource(url, state["root"] / "public")
        soup = BeautifulSoup(raw, "html.parser")
        section = soup.find(id=shortname + ".Molecular_Pathogenesis")
        inheritance = soup.find(id=shortname + ".Mode_of_Inheritance")
        table = soup.find(id="__" + shortname + ".molgen.TB_lrgtbl__")
        if section is None or inheritance is None or table is None:
            raise ValueError("GeneReviews mechanism sections unavailable")
        text = " ".join(section.get_text(" ", strip=True).split())
        mode_text = " ".join(inheritance.get_text(" ", strip=True).split())
        if not re.search(r"\b" + re.escape(gene) + r"\b", text):
            raise ValueError("GeneReviews gene mismatch")
        mechanism = next(
            (
                " ".join(p.get_text(" ", strip=True).split())
                for p in section.find_all("p")
                if p.get_text(" ", strip=True).startswith("Mechanism of disease causation.")
            ),
            "",
        )
        relevant = re.findall(r"transcript (NM_\d+) is the most highly expressed", text)
        chapter_modes = [
            mode
            for phrase, mode in (
                ("autosomal dominant", "AD"),
                ("autosomal recessive", "AR"),
                ("X-linked", "XL"),
            )
            if phrase.casefold() in mode_text.casefold()
        ]
        record = {
            "record_type": "genereviews_mechanism_v1",
            "gene": gene,
            "disease_label": soup.find("h1").get_text(" ", strip=True),
            "disease_ids": re.findall(r"\b\d{6}\b", table.get_text(" ", strip=True)),
            "inheritance": chapter_modes[0] if len(chapter_modes) == 1 else "unknown",
            "inheritance_original": mode_text,
            "evidence_summary": text,
            "mechanism_original": mechanism,
            "relevant_transcripts": relevant,
            "revision_original": next(
                (
                    " ".join(p.get_text(" ", strip=True).split())
                    for p in soup.find_all("p")
                    if "Last Revision:" in p.get_text(" ", strip=True)
                ),
                "",
            ),
        }
        statement = Source(
            provider="genereviews",
            kind="gene_disease",
            locator=url,
            version="original-mechanism-v1;sha256=" + hashlib.sha256(raw).hexdigest(),
            content=json.dumps(record, ensure_ascii=False, sort_keys=True),
            genes=(gene,),
        )
        documents.append(
            _document(
                statement,
                "GENEREVIEWS:" + accession,
                record["disease_label"],
                statement.content,
                "gene_disease_context",
            )
        )
    return tuple(documents)


def _family_series(record, original, article, identity):
    """Validate curated table rows against their original study, without grading them."""
    fields = {
        "variant_label",
        "variant_table_quote",
        "phenotype_header_quote",
        "phenotype_columns",
        "family_rows",
        "recruitment_quote",
        "independence_quote",
        "phenotype_limitations_quote",
        "methods_quote",
        "parent_analysis_quote",
        "de_novo_quote",
        "reported_de_novo_count",
    }
    absent_de_novo = record.get("reported_de_novo_count") is None
    if set(record) != fields or any(
        (
            record[k] is not None
            if absent_de_novo and k in {"parent_analysis_quote", "de_novo_quote"}
            else not record[k] or record[k] not in original
        )
        for k in fields
        if k.endswith("_quote")
    ):
        raise ValueError("family ledger original passage mismatch")
    if record["variant_label"] not in identity.content or not record["variant_table_quote"].startswith(
        record["variant_label"] + " "
    ):
        raise ValueError("family ledger variant identity mismatch")
    match = re.search(r"(\d+(?:,\s*\d+)*)$", record["variant_table_quote"])
    identifiers = tuple(match[1].replace(" ", "").split(",")) if match else ()
    columns = record["phenotype_columns"]
    if (
        not identifiers
        or len(set(identifiers)) != len(identifiers)
        or not columns
        or len(set(columns)) != len(columns)
    ):
        raise ValueError("family ledger table columns or identifiers invalid")
    header, cursor = re.sub(r"[^a-z]", "", record["phenotype_header_quote"].lower()), 0
    for column in columns:
        label = re.sub(r"[^a-z]", "", column.lower())
        position = header.find(label, cursor)
        if not label or position < 0:
            raise ValueError("family phenotype columns do not match original header order")
        cursor = position + len(label)
    families = []
    for quote in record["family_rows"]:
        if quote not in original.splitlines():
            raise ValueError("family ledger row is not an original table row")
        values = quote.split()
        number, affected = values[:2]
        sporadic = values[2:3] == ["(S)"]
        cells = values[3:] if sporadic else values[2:]
        if (
            number not in identifiers
            or not affected.isdigit()
            or int(affected) < 1
            or len(cells) != len(columns)
        ):
            raise ValueError("family ledger row shape mismatch")
        for cell in cells:
            if cell == "Unknown":
                continue
            parts = cell.split("/")
            if (
                len(parts) != 2
                or not all(p.isdigit() for p in parts)
                or not 0 <= int(parts[0]) <= int(parts[1]) <= int(affected)
            ):
                raise ValueError("family phenotype ratio invalid")
        families.append(
            {
                "family_id": number,
                "affected_members": int(affected),
                "sporadic_marker": sporadic,
                "phenotypes": dict(zip(columns, cells, strict=True)),
                "original_table_row": quote,
            }
        )
    if tuple(f["family_id"] for f in families) != identifiers:
        raise ValueError("family ledger does not match the variant table")
    count = record["reported_de_novo_count"]
    words = "zero one two three four five six seven eight nine ten".split()
    unrelated = record["independence_quote"].lower()
    reported = (record["de_novo_quote"] or "").lower()

    def number_words(n):
        return (str(n), words[n]) if n < len(words) else (str(n),)

    if not any(f"{n} unrelated families" in unrelated for n in number_words(len(families))):
        raise ValueError("unrelated family count is not supported by original text")
    if count is not None and (
        type(count) is not int
        or not 0 <= count <= len(families)
        or not any(
            f"{n} of {total} unrelated families" in reported
            for n in number_words(count)
            for total in number_words(len(families))
        )
        or "de novo" not in record["parent_analysis_quote"].lower()
    ):
        raise ValueError("reported de novo subset is not supported by original text")
    return {
        "schema_version": "publication_family_series_v1",
        "gene": article.genes[0],
        "variant_key": variant_text(article.variant_keys[0]),
        "article_source_id": article.source_id,
        "identity_source_id": identity.source_id,
        "publication": article.locator,
        "variant_table_original": record["variant_table_quote"],
        "phenotype_header_original": record["phenotype_header_quote"],
        "families": families,
        "within_publication_family_count": len(families),
        "recruitment_original": record["recruitment_quote"],
        "methods_original": record["methods_quote"],
        "phenotype_limitations_original": record["phenotype_limitations_quote"],
        "reported_de_novo": {
            "count": count,
            "family_ids": None,
            "family_id_mapping": "not_reported",
            "parent_variant_analysis": "reported",
            "parental_relationships": "not_documented",
            "original": [record["parent_analysis_quote"], record["de_novo_quote"]],
        }
        if count is not None
        else None,
        "cross_publication_deduplication": "not_performed",
        "assigned_ps2_pm6_points": None,
    }


def _curated_documents(state, gene, targets):
    """Load reviewed original excerpts, never precomputed criteria or grades."""
    path = state["settings"].get("post_evidence_catalog")
    if not path:
        return ()
    with state["lock"]:
        if "post_catalog" not in state["tables"]:
            catalog = json.loads(read_regular_bounded(path, _MAX_GUIDELINE)[0])
            if (
                set(catalog) != {"schema_version", "records"}
                or catalog["schema_version"] != "verified_original_evidence_v1"
            ):
                raise ValueError("unsupported original evidence catalog")
            state["tables"]["post_catalog"] = catalog["records"]
    keys = {_key(t["primary"]["key"], "report.variant_key") for t in targets}
    documents = []
    for record in state["tables"]["post_catalog"]:
        if record["gene"] != gene:
            continue
        key = _key(tuple(record["variant_key"]), "catalog.variant_key")
        if key not in keys:
            continue
        if set(record) not in (
            {"gene", "variant_key", "identity", "article"},
            {"gene", "variant_key", "identity", "article", "family_series"},
        ):
            raise ValueError("unexpected original evidence fields")
        resolved = state["reference"].resolve(key)
        if (
            resolved["status"] != "sequence_verified"
            or resolved["normalized_hgvs"] != record["identity"]["genomic_hgvs"]
        ):
            raise ValueError("original evidence genomic identity mismatch")
        originals = []
        for kind, item in (("annotation", record["identity"]), ("literature", record["article"])):
            fields = {
                "provider",
                "url",
                "title",
                "version",
                "original_file",
                "original_sha256",
                "text_file",
                "text_sha256",
                "quotes",
            }
            if set(item) != (fields | {"genomic_hgvs"} if kind == "annotation" else fields):
                raise ValueError("unexpected original evidence metadata")
            original_file = (path.parent / item["original_file"]).resolve()
            if not original_file.is_relative_to(path.parent.resolve()):
                raise ValueError("original evidence path outside resource")
            original_bytes = read_regular_bounded(original_file, _MAX_TEXT)[0]
            if hashlib.sha256(original_bytes).hexdigest() != item["original_sha256"]:
                raise ValueError("original evidence document changed")
            file = (path.parent / item["text_file"]).resolve()
            if not file.is_relative_to(path.parent.resolve()):
                raise ValueError("original evidence path outside resource")
            raw = read_regular_bounded(file, _MAX_TEXT)[0]
            if hashlib.sha256(raw).hexdigest() != item["text_sha256"]:
                raise ValueError("original evidence text changed")
            original = raw.decode("utf-8")
            if not item["quotes"] or any(not q or q not in original for q in item["quotes"]):
                raise ValueError("original evidence quote mismatch")
            content = "\n\n".join(item["quotes"])
            if kind == "annotation" and (
                record["identity"]["genomic_hgvs"] not in content or gene not in content
            ):
                raise ValueError("original allele mapping quote missing")
            statement = Source(
                provider=item["provider"],
                kind=kind,
                locator=item["url"],
                version=item["version"] + ";text-sha256=" + item["text_sha256"],
                content=content,
                genes=(gene,),
                variant_keys=(key,),
            )
            documents.append(
                _document(
                    statement,
                    "VERIFIED:" + statement.source_id,
                    item["title"],
                    content,
                    "verified_original_evidence",
                )
            )
            originals.append(statement)
        if "family_series" in record:
            fact = _family_series(record["family_series"], original, originals[1], originals[0])
            content = json.dumps(fact, ensure_ascii=False, sort_keys=True)
            ledger = Source(
                provider="publication-family-ledger",
                kind="annotation",
                locator=originals[1].locator.rstrip("/") + "/#family-series",
                version="publication-family-series-v1;sha256=" + hashlib.sha256(content.encode()).hexdigest(),
                content=content,
                genes=(gene,),
                variant_keys=(key,),
            )
            from mendelprime.tools.acmg import validate_family_series_source

            validate_family_series_source(ledger, {s.source_id: s for s in originals})
            documents.append(
                _document(
                    ledger,
                    "FAMILIES:" + ledger.source_id,
                    "Original publication family ledger",
                    content,
                    "publication_family_series",
                )
            )
    return tuple(documents)


def _curation_documents(state, gene, contexts):
    """Read explicit original curation statements, not previous PVS1 verdicts."""
    from bs4 import BeautifulSoup

    documents, issues = [], []
    for doc in contexts:
        source = doc["source"]
        if source.provider != "clingen":
            continue
        try:
            record = json.loads(source.content)
        except ValueError:
            continue
        row = record.get("row", {})
        if record.get("record_type") != "clingen_gene_validity" or row.get("GENE SYMBOL") != gene:
            continue
        url = row.get("ONLINE REPORT", "")
        if not url.startswith("https://search.clinicalgenome.org/kb/gene-validity/"):
            continue
        try:
            with state["lock"]:
                raw = public_resource(url, state["root"] / "public")
            soup = BeautifulSoup(raw, "html.parser")
            label = soup.find(string=lambda text: bool(text) and text.strip() == "Evidence Summary:")
            cell = label.find_parent("td").find_next_sibling("td") if label else None
            text = " ".join(cell.get_text(" ", strip=True).split()) if cell else ""
            if not text:
                raise ValueError("missing evidence summary")
            record = {
                "gene": gene,
                "disease_id": row.get("DISEASE ID (MONDO)", ""),
                "disease_label": row.get("DISEASE LABEL", ""),
                "inheritance": row.get("MOI", ""),
                "classification": row.get("CLASSIFICATION", ""),
                "evidence_summary": text,
            }
            statement = Source(
                provider="clingen",
                kind="gene_disease",
                locator=url,
                version="evidence-summary;sha256=" + hashlib.sha256(raw).hexdigest(),
                content=json.dumps(record, ensure_ascii=False, sort_keys=True),
                genes=(gene,),
            )
            documents.append(
                _document(
                    statement,
                    "CLINGEN:MECHANISM:" + statement.source_id,
                    "ClinGen original evidence summary",
                    statement.content,
                    "gene_disease_context",
                )
            )
        except MendelPrimeError as error:
            issues.append(f"clingen:{gene}:{error.code}")
        except (AttributeError, ValueError):
            issues.append(f"clingen:{gene}:curation_summary_invalid")
        except Exception:
            issues.append(f"clingen:{gene}:curation_summary_failed")
    path = state["settings"].get("clingen_dosage_table")
    if path:
        try:
            with state["lock"]:
                if "dosage" not in state["tables"]:
                    lines = path.read_text().splitlines()
                    header = next(i for i, line in enumerate(lines) if line.startswith("#Gene Symbol\t"))
                    state["tables"]["dosage"] = {
                        r["Gene Symbol"]: r
                        for r in csv.DictReader(
                            [lines[header].lstrip("#"), *lines[header + 1 :]], delimiter="\t"
                        )
                    }
                row = state["tables"]["dosage"].get(gene)
        except (OSError, UnicodeError, ValueError, StopIteration):
            row = None
            issues.append(f"clingen:{gene}:dosage_table_unavailable")
        if row:
            statement = Source(
                provider="clingen",
                kind="gene_disease",
                locator=f"https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen/clingen_gene.cgi?sym={gene}",
                version="GRCh37 dosage;" + str(row.get("Date Last Evaluated", "")),
                genes=(gene,),
                content=json.dumps(
                    {"gene": gene, "record_type": "clingen_dosage", "row": row}, sort_keys=True
                ),
            )
            documents.append(
                _document(
                    statement,
                    "CLINGEN:DOSAGE:" + gene,
                    "ClinGen dosage curation",
                    statement.content,
                    "gene_disease_context",
                )
            )
    return tuple(documents), tuple(issues)


def _complete_report_evidence(state, seed, review, selected):
    from mendelprime.tools.native import kinship_observation

    knowledge, evidence = dict(seed["variant_knowledge"]), dict(seed["evidence"])
    registry, issues = {s.source_id: s for s in seed["sources"]}, []
    qc = kinship_observation(
        state["settings"].get("kinship_executable"), review.get("vcf_path"), seed, state["root"]
    )
    if qc is not None:
        registry[qc.source_id] = qc
        if "status=kinship_supports_recorded_trio" not in qc.content:
            detail = json.loads(qc.content.split("\n", 1)[0])
            issues.append("kinship:" + detail.get("error_code", "relationship_verification_unavailable"))

    def prepare(gene):
        current = knowledge.get(gene, {})
        curated, failures = _curation_documents(state, gene, current.get("gene_context", ()))
        targets = tuple(item for item in selected if item["primary"]["gene"] == gene)
        try:
            curated = (*curated, *_genereviews_documents(state, gene))
            originals = _curated_documents(state, gene, targets)
        except Exception as error:
            raise MendelPrimeError(
                "post_original_evidence_invalid", category="resource", field=gene
            ) from error
        return {
            "gene": gene,
            "current": current,
            "curated": curated,
            "targets": targets,
            "originals": originals,
            "failures": failures,
        }

    def complete(job):
        gene, current, curated = job["gene"], job["current"], job["curated"]
        expanded = literature.complete_report_documents(
            state["pubmed"],
            gene,
            evidence.get(gene),
            job["targets"],
            current.get("citation_seeds", ()),
            seed["profile"]["case_key"],
            fetched=fetched,
            contexts=(*current.get("gene_context", ()), *curated),
        )
        documents = {d["source"].source_id: d for d in (*expanded["documents"], *job["originals"])}
        expanded = {**expanded, "documents": tuple(documents.values())}
        return (
            gene,
            {**current, "gene_context": (*current.get("gene_context", ()), *curated)},
            expanded,
            job["failures"],
        )

    with ThreadPoolExecutor(max_workers=min(state["workers"], len(seed["fixed_genes"])) or 1) as pool:
        jobs = tuple(pool.map(prepare, seed["fixed_genes"]))
        pmids = []
        for job in jobs:
            current, gene = job["current"], job["gene"]
            existing = {
                literature._variant_article_key(d).removeprefix("PMID:")
                for d in evidence.get(gene, {}).get("documents", ())
                if d["source"].provider == "pubmed"
            }
            plan = literature.report_citations(
                gene,
                job["targets"],
                current.get("citation_seeds", ()),
                (*current.get("gene_context", ()), *job["curated"]),
            )
            pmids.extend(pmid for pmid in plan if pmid not in existing)
        fetched = literature.fetch_many(state["pubmed"], tuple(dict.fromkeys(pmids))) if pmids else ({}, {})
        for gene, current, expanded, failures in pool.map(complete, jobs):
            knowledge[gene], evidence[gene] = current, expanded
            issues.extend(failures)
            _register_report_sources(registry, current, expanded)
    return {
        **seed,
        "variant_knowledge": knowledge,
        "evidence": evidence,
        "sources": tuple(registry.values()),
    }, tuple(issues)


def _register_report_sources(registry, current, expanded):
    for document in (
        *current.get("gene_context", ()),
        *expanded.get("documents", ()),
        *expanded.get("omim_context", ()),
    ):
        registry[document["source"].source_id] = document["source"]
    for source in expanded.get("citation_resolution", ()):
        registry[source.source_id] = source


def close_resources(state):
    """Idempotently close an ``open_resources`` state."""
    reference = state.get("reference") if isinstance(state, dict) else None
    if reference is not None:
        reference.close()


def _check_resources(settings, *, vep_enabled=True, report=True):
    if report:
        for name in ("mane_summary", "mane_release", "ai_cura_protocol", "ai_cura_workbook", "clingen_guideline_manifest"):
            if not settings.get(name):
                raise MendelPrimeError("required_for_report", field=f"explanation_resources.{name}")
    if vep_enabled:
        for name in ("vep_command", "vep_cache", "vep_fasta", "vep_cache_version", "vep_tool_version"):
            if not settings.get(name):
                raise MendelPrimeError("required_for_vep", field=f"explanation_resources.{name}")
    files = (
        "reference_fasta",
        "vep_fasta",
        "clinvar_protein_disease_db",
        "gnomad_sites",
        "gnomad_coverage",
        "gnomad_constraint",
        "spliceai_snv",
        "spliceai_indel",
        "revel",
        "mane_summary",
        "ai_cura_protocol",
        "ai_cura_workbook",
        "ai_cura_knowledgebase",
        "ai_cura_manifest",
        "clingen_guideline_manifest",
    )
    for name in files:
        if name in _REPORT_FIELDS and not report:
            continue
        if name == "vep_fasta" and not vep_enabled:
            continue
        _require_file(settings.get(name), f"explanation_resources.{name}")
    for name in ("vep_cache", "autopvs1_root") if vep_enabled else ("autopvs1_root",):
        _require_directory(settings.get(name), f"explanation_resources.{name}")
    if vep_enabled:
        _require_directory(
            settings["vep_cache"]
            / "homo_sapiens"
            / f"{settings['vep_cache_version']}_{settings['assembly']}",
            "explanation_resources.vep_cache_release",
        )
    _require_fasta_indexes(settings["reference_fasta"], "explanation_resources.reference_fasta")
    if vep_enabled:
        _require_fasta_indexes(settings["vep_fasta"], "explanation_resources.vep_fasta")
    for name in ("gnomad_sites", "gnomad_coverage", "spliceai_snv", "spliceai_indel", "revel"):
        _require_index(settings[name], f"explanation_resources.{name}", (".tbi", ".csi"))
    if vep_enabled:
        executable = Path(settings["vep_command"][0])
        if executable.is_absolute():
            _require_file(executable, "explanation_resources.vep_command[0]")
            if not os.access(executable, os.X_OK):
                raise MendelPrimeError("executable_required", field="explanation_resources.vep_command[0]")
    error = preflight_source(settings["autopvs1_root"])
    if error:
        raise MendelPrimeError(error, field="explanation_resources.autopvs1_root")


def _resource_identity(settings, *, vep_enabled):
    """Track configured versions and file identities without reading shared databases."""
    paths, values = set(), {}
    for name, value in settings.items():
        if name.startswith("vep_") and not vep_enabled:
            continue
        if isinstance(value, Path):
            paths.add(value)
        else:
            values[name] = value
    root = settings.get("autopvs1_root")
    if isinstance(root, Path):
        config_path = root / "config.ini"
        paths.add(config_path)
        paths.update(root.rglob("*.py"))
        paths.update(
            root / name for name in ("maxentpy/data/score3_matrix.txt", "maxentpy/data/score5_matrix.txt")
        )
        config = configparser.ConfigParser()
        try:
            with config_path.open(encoding="utf-8") as handle:
                config.read_file(handle)
            for section in config:
                for value in config[section].values():
                    path = Path(value)
                    paths.add(path if path.is_absolute() else root / path)
        except (OSError, UnicodeError, configparser.Error):
            values["autopvs1_config"] = "unavailable"
    fingerprints = {}
    for path in sorted(paths):
        for target in (path, *(Path(f"{path}{suffix}") for suffix in (".fai", ".gzi", ".tbi", ".csi"))):
            try:
                fingerprints[str(target)] = stat_fingerprint(target.stat())
            except OSError:
                fingerprints[str(target)] = None
    return {"settings": values, "files": fingerprints, "vep_enabled": vep_enabled}


def _require_file(path, field):
    try:
        valid = (
            isinstance(path, Path) and path.is_file() and path.stat().st_size > 0 and os.access(path, os.R_OK)
        )
    except OSError:
        valid = False
    if not valid:
        raise MendelPrimeError("readable_regular_file_required", field=field)


def _require_directory(path, field):
    try:
        valid = isinstance(path, Path) and path.is_dir() and os.access(path, os.R_OK | os.X_OK)
    except OSError:
        valid = False
    if not valid:
        raise MendelPrimeError("readable_directory_required", field=field)


def _require_index(path, field, suffixes):
    try:
        valid = any(
            Path(f"{path}{suffix}").is_file() and Path(f"{path}{suffix}").stat().st_size > 0
            for suffix in suffixes
        )
    except OSError:
        valid = False
    if not valid:
        raise MendelPrimeError("index_required", field=field)


def _require_fasta_indexes(path, field):
    _require_index(path, field, (".fai",))
    if path.name.endswith((".gz", ".bgz")):
        _require_index(path, field, (".gzi",))


def _read_limited(path, limit, field):
    try:
        with path.open("rb") as stream:
            value = stream.read(limit + 1)
    except OSError:
        raise MendelPrimeError("resource_unavailable", field=field) from None
    if not value or len(value) > limit:
        raise MendelPrimeError("resource_size_invalid", field=field)
    return value


def _load_clingen_guidance(settings):
    field = "explanation_resources.clingen_guideline_manifest"
    try:
        manifest = json.loads(
            _read_limited(settings["clingen_guideline_manifest"], _MAX_MANIFEST, field).decode()
        )
    except (UnicodeError, json.JSONDecodeError, RecursionError):
        raise MendelPrimeError("schema_mismatch", field=field) from None
    if not isinstance(manifest, dict) or not isinstance(manifest.get("resources"), list):
        raise MendelPrimeError("schema_mismatch", field=field)
    root, sources = settings["clingen_guideline_manifest"].parent.resolve(), []
    for index, value in enumerate(manifest["resources"]):
        item_field = f"{field}.resources[{index}]"
        required = {"name", "url", "sha256", "text_file", "text_sha256"}
        if (
            not isinstance(value, dict)
            or not required.issubset(value)
            or any(not isinstance(value.get(name), str) or not value[name] for name in required)
        ):
            raise MendelPrimeError("schema_mismatch", field=item_field)
        if not _SHA256.fullmatch(value["sha256"]) or not _SHA256.fullmatch(value["text_sha256"]):
            raise MendelPrimeError("sha256_required", field=item_field)
        text_path = (root / value["text_file"]).resolve()
        if not text_path.is_relative_to(root):
            raise MendelPrimeError("resource_escape_forbidden", field=f"{item_field}.text_file")
        content = _read_limited(text_path, _MAX_GUIDELINE, f"{item_field}.text_file")
        if hashlib.sha256(content).hexdigest() != value["text_sha256"]:
            raise MendelPrimeError("hash_mismatch", field=f"{item_field}.text_sha256")
        try:
            text = content.decode()
        except UnicodeError:
            raise MendelPrimeError("utf8_text_required", field=f"{item_field}.text_file") from None
        sources.append(
            Source(
                provider="clingen_svi",
                kind="guideline",
                locator=value["url"],
                version=f"{value['name']};sha256={value['sha256']}",
                content=text,
            )
        )
    if not sources:
        raise MendelPrimeError("nonempty_resources_required", field=field)
    return tuple(sources)


def gene_context(state, gene, *, ranking=False):
    validate_gene(gene)
    documents, outcomes = [], []
    with state["lock"]:
        for provider, (path, release) in state["resources"].items():
            if ranking and provider == "mane":
                continue
            try:
                table = state["tables"].get(provider)
                if table is None:
                    table = _load_context_table(path, provider)
                    state["tables"][provider] = table
                if stat_fingerprint(path.stat()) != table["file_state"]:
                    raise RuntimeError
                rows = table["rows"].get(gene, ())
            except RuntimeError:
                outcomes.append(_lookup(provider, gene, None, "query_failed", "context_resource_changed"))
                continue
            except (ValueError, EOFError, gzip.BadGzipFile, zlib.error):
                outcomes.append(_lookup(provider, gene, None, "schema_mismatch", "context_resource_invalid"))
                continue
            except OSError:
                outcomes.append(
                    _lookup(provider, gene, None, "resource_unavailable", "context_resource_unavailable")
                )
                continue
            outcomes.append(_lookup(provider, gene, None, "found" if rows else "not_found"))
            for number, raw in enumerate(dict.fromkeys(rows), 1):
                source = Source(
                    provider=provider,
                    kind="annotation" if provider == "mane" else "prediction",
                    locator=f"{_CONTEXT_URL[provider]}#gene/{gene}/row/{number}",
                    version=f"{release};compressed-sha256={table['digest']}",
                    genes=(gene,),
                    content=json.dumps(
                        {"header": table["header"], "row": raw}, ensure_ascii=False, sort_keys=True
                    ),
                )
                documents.append(
                    _document(
                        source,
                        f"{provider.upper()}:{gene}:{number}",
                        f"{provider} gene/transcript context",
                        f"{_CONTEXT_SCOPE[provider]}\n{source.content}",
                        "gene_context",
                    )
                )
    return _knowledge_result(gene, documents=documents, lookups=outcomes)


def _load_context_table(path, provider):
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NONBLOCK", 0) | getattr(os, "O_BINARY", 0))
    with os.fdopen(descriptor, "rb") as stream:
        before = os.fstat(stream.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise OSError
        compressed = stream.read(_MAX_COMPRESSED + 1)
        after = os.fstat(stream.fileno())
    if len(compressed) > _MAX_COMPRESSED or stat_fingerprint(before) != stat_fingerprint(after):
        raise RuntimeError
    with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as stream:
        content = stream.read(_MAX_TEXT + 1)
    if len(content) > _MAX_TEXT:
        raise ValueError
    lines = content.decode().splitlines()
    if not lines:
        raise ValueError
    header = tuple(lines[0].lstrip("#").split("\t"))
    if header != _CONTEXT_HEADERS[provider]:
        raise ValueError
    column, grouped = header.index(_CONTEXT_FIELD[provider]), {}
    for line in lines[1:]:
        if not line:
            continue
        cells = line.split("\t")
        if len(cells) != len(header):
            raise ValueError
        gene = cells[column]
        if gene and gene not in {"-", ".", "NA"}:
            grouped.setdefault(gene, []).append(line)
    return {
        "header": header,
        "rows": {gene: tuple(rows) for gene, rows in grouped.items()},
        "digest": hashlib.sha256(compressed).hexdigest(),
        "file_state": stat_fingerprint(after),
    }


def protein_context(settings, gene, *, variant_key=None, transcript=None, protein_hgvs=None):
    validate_gene(gene)
    if variant_key is not None:
        variant_key = _key(variant_key, "protein_context.variant_key")
    if any(value is not None and not isinstance(value, str) for value in (transcript, protein_hgvs)):
        raise MendelPrimeError("annotation_text_or_none_required", field="protein_context")
    match = _CODON.search(protein_hgvs or "")
    if variant_key is None or not transcript or match is None:
        return _knowledge_result(
            gene, lookups=(_lookup("clinvar-protein", gene, variant_key, "not_applicable"),)
        )
    documents, status, error = [], "not_found", ""
    try:
        database = settings["clinvar_protein_disease_db"]
        if not database.is_file():
            raise OSError
        with closing(sqlite3.connect(database.as_uri() + "?mode=ro&immutable=1", uri=True)) as connection:
            connection.row_factory = sqlite3.Row
            names = tuple(row["name"] for row in connection.execute("PRAGMA table_info(protein_variants)"))
            if names != _PROTEIN_COLUMNS:
                return _knowledge_result(
                    gene,
                    lookups=(
                        _lookup(
                            "clinvar-protein",
                            gene,
                            variant_key,
                            "schema_mismatch",
                            "protein_context_schema_invalid",
                        ),
                    ),
                )
            rows = connection.execute(
                "SELECT * FROM protein_variants WHERE gene=? AND transcript=? AND codon=? LIMIT 20",
                (gene, transcript, int(match[1])),
            ).fetchall()
        status = "found" if rows else "not_found"
        for number, row in enumerate(rows, 1):
            if (
                row["gene"] != gene
                or row["transcript"] != transcript
                or row["codon"] != int(match[1])
                or not str(row["variation_id"]).isdigit()
            ):
                return _knowledge_result(
                    gene,
                    lookups=(
                        _lookup(
                            "clinvar-protein",
                            gene,
                            variant_key,
                            "schema_mismatch",
                            "protein_context_record_invalid",
                        ),
                    ),
                )
            content = json.dumps(dict(row), ensure_ascii=False, sort_keys=True)
            record_id = hashlib.sha256(content.encode()).hexdigest()[:16]
            source = Source(
                provider="clinvar-protein",
                kind="annotation",
                locator=f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{row['variation_id']}/#protein-row/{record_id}",
                version=settings["clinvar_protein_release"],
                genes=(gene,),
                content=content,
            )
            documents.append(
                _document(
                    source,
                    f"PCODON:{gene}:{variant_text(variant_key)}:{number}",
                    "ClinVar same-transcript/codon comparator context",
                    "This record concerns the comparator variant, not necessarily the target allele. Same codon alone establishes neither PS1 nor PM5; disease and allele independence remain to be reviewed.\n"
                    + source.content,
                    "protein_comparison",
                )
            )
    except OSError:
        status, error = "resource_unavailable", "protein_context_resource_unavailable"
    except (sqlite3.Error, ValueError, TypeError, MendelPrimeError):
        status, error = "query_failed", "protein_context_query_failed"
    return _knowledge_result(
        gene,
        documents=documents if not error else (),
        lookups=(_lookup("clinvar-protein", gene, variant_key, status, error),),
    )


def _document(source, prompt_id, title, prompt_text, section):
    pmids = (
        literature.context_pmids(source)
        if source.provider in {"clingen", "genereviews"} and source.kind == "gene_disease"
        else ()
    )
    return literature.evidence_document(
        source, prompt_id, title, prompt_text, cited_pmids=pmids, section=section
    )


def _lookup(provider, gene, key, status, error_code=""):
    return {
        "provider": provider,
        "gene": gene,
        "variant_key": key,
        "status": status,
        "error_code": error_code,
    }


def _knowledge_result(gene, *, documents=(), lookups=()):
    return {
        "gene": gene,
        "variant_evidence": (),
        "gene_context": tuple(documents),
        "citation_seeds": (),
        "lookups": tuple(lookups),
    }


def _key(value, field):
    if type(value) is not tuple or len(value) != 5:
        raise MendelPrimeError("variant_key_required", field=field)
    return variant_key(*value)


def _load_protocol(settings, *, report=True):
    paths = {
        name: settings[f"ai_cura_{name}"]
        for name in (("protocol", "workbook", "knowledgebase", "manifest") if report else ("knowledgebase", "manifest"))
    }
    raw = {name: _read_protocol(path, name) for name, path in paths.items()}
    digests = {name: hashlib.sha256(value).hexdigest() for name, value in raw.items()}
    if report and digests["protocol"] != _EXPECTED_PROTOCOL:
        raise MendelPrimeError("hash_mismatch", field="protocol")
    if digests["manifest"] not in (_EXPECTED_MANIFEST, _PUBLIC_MANIFEST):
        raise MendelPrimeError("hash_mismatch", field="manifest")
    manifest = _json_object(raw["manifest"], "manifest")
    knowledgebase = _json_object(raw["knowledgebase"], "knowledgebase")
    release = _verify_manifest(manifest, paths, raw, digests)
    workbook = manifest["primary_method"]["tables_workbook"]
    workbook_path = paths.get("workbook") or Path(PureWindowsPath(workbook["path"]).name)
    workbook_digest = digests.get("workbook", workbook["sha256"])
    _verify_knowledgebase(knowledgebase, workbook_path, workbook_digest)
    version = (
        f"ai-cura:{release};doi={_DOI};workbook-sha256={workbook_digest};"
        f"derived-sha256={digests['knowledgebase']};manifest-sha256={digests['manifest']}"
    )
    protocol_source = Source(
        provider="ai-cura",
        kind="guideline",
        locator=f"{_DOI_URL}#retained-protocol",
        version=f"ai-cura:{release};doi={_DOI};protocol-sha256={digests['protocol']};manifest-sha256={digests['manifest']}",
        content=_decode(raw["protocol"], "protocol"),
    ) if report else None
    consistency = _table_sources(
        knowledgebase,
        version,
        "gene_disease_consistency",
        8,
        ("specificity", "genes", "disease_phenotypes"),
        gene_field="genes",
    )
    segregation = _table_sources(
        knowledgebase,
        version,
        "pp1_recommendations",
        9,
        ("segregation", "genotype", "phenotype", "points_per_individual"),
        optional=("remark",),
    ) if report else ()
    phenotype = _table_sources(
        knowledgebase,
        version,
        "pp4_knowledgebase",
        10,
        (
            "disease",
            "gene",
            "clinical_observable_phenotype",
            "regular_genetic_tests",
            "specific_tests",
            "pp4_evidence_level_to_be_assigned",
        ),
        gene_field="gene",
    )
    generic = {source.source_id for source in segregation}
    generic.update(
        source.source_id
        for source, row in zip(consistency, knowledgebase["gene_disease_consistency"], strict=True)
        if row["genes"] == _GENERIC_S8
    )
    generic.update(
        source.source_id
        for source, row in zip(phenotype, knowledgebase["pp4_knowledgebase"], strict=True)
        if row["gene"] in _GENERIC_S10
    )
    return {
        "protocol_source": protocol_source,
        "rows": (*consistency, *segregation, *phenotype),
        "ranking_rows": (*consistency, *phenotype),
        "generic_source_ids": frozenset(generic),
    }


def _read_protocol(path, field):
    try:
        with path.open("rb") as handle:
            value = handle.read(_MAX_MANIFEST + 1)
    except OSError:
        raise MendelPrimeError("resource_unavailable", field=field) from None
    if len(value) > _MAX_MANIFEST:
        raise MendelPrimeError("resource_too_large", field=field)
    return value


def _decode(value, field):
    try:
        text = value.decode()
    except UnicodeError:
        raise MendelPrimeError("utf8_required", field=field) from None
    if not text:
        raise MendelPrimeError("nonempty_resource_required", field=field)
    return text


def _json_object(value, field):
    def pairs(items):
        result = {}
        for key, item in items:
            if key in result:
                raise ValueError
            result[key] = item
        return result

    try:
        parsed = json.loads(
            _decode(value, field),
            object_pairs_hook=pairs,
            parse_constant=lambda _: (_ for _ in ()).throw(ValueError()),
        )
    except (json.JSONDecodeError, ValueError):
        raise MendelPrimeError("malformed_json", field=field) from None
    if not isinstance(parsed, dict):
        raise MendelPrimeError("object_required", field=field)
    return parsed


def _verify_manifest(value, paths, raw, digests):
    if value.get("schema_version") != _MANIFEST_SCHEMA:
        raise MendelPrimeError("schema_mismatch", field="manifest.schema_version")
    release = value.get("release")
    try:
        if not isinstance(release, str) or date.fromisoformat(release).isoformat() != release:
            raise ValueError
    except ValueError:
        raise MendelPrimeError("iso_date_required", field="manifest.release") from None
    primary = value.get("primary_method")
    if not isinstance(primary, dict):
        raise MendelPrimeError("object_required", field="manifest.primary_method")
    if primary.get("doi") != _DOI:
        raise MendelPrimeError("source_mismatch", field="manifest.primary_method.doi")
    if "workbook" in raw:
        _verify_artifact(
            primary.get("tables_workbook"),
            "manifest.primary_method.tables_workbook",
            paths["workbook"],
            len(raw["workbook"]),
            digests["workbook"],
        )
    _verify_artifact(
        primary.get("derived_knowledgebase"),
        "manifest.primary_method.derived_knowledgebase",
        paths["knowledgebase"],
        len(raw["knowledgebase"]),
        digests["knowledgebase"],
    )
    return release


def _verify_artifact(value, field, path, size, digest):
    if not isinstance(value, dict):
        raise MendelPrimeError("object_required", field=field)
    declared = value.get("path")
    name = (
        PureWindowsPath(declared).name
        if isinstance(declared, str) and "\\" in declared
        else PurePosixPath(declared).name
        if isinstance(declared, str)
        else ""
    )
    if name != path.name:
        raise MendelPrimeError("source_mismatch", field=f"{field}.path")
    if type(value.get("size_bytes")) is not int or value["size_bytes"] < 0:
        raise MendelPrimeError("nonnegative_integer_required", field=f"{field}.size_bytes")
    if not isinstance(value.get("sha256"), str) or not _SHA256.fullmatch(value["sha256"]):
        raise MendelPrimeError("sha256_required", field=f"{field}.sha256")
    if value["size_bytes"] != size:
        raise MendelPrimeError("size_mismatch", field=field)
    if value["sha256"] != digest:
        raise MendelPrimeError("hash_mismatch", field=field)


def _verify_knowledgebase(value, workbook, workbook_digest):
    if set(value) - _KNOWLEDGEBASE_KEYS:
        raise MendelPrimeError("unexpected_top_level_field", field="knowledgebase")
    required = _KNOWLEDGEBASE_KEYS - {"ps4_thresholds", "ps2_pm6_strength_thresholds"}
    if not required.issubset(value):
        raise MendelPrimeError("required_field_missing", field="knowledgebase")
    if value.get("schema_version") != _KNOWLEDGEBASE_SCHEMA:
        raise MendelPrimeError("schema_mismatch", field="knowledgebase.schema_version")
    source = value.get("source")
    if not isinstance(source, dict) or set(source) != {"doi", "file", "sha256"}:
        raise MendelPrimeError("schema_mismatch", field="knowledgebase.source")
    expected = {"doi": _DOI, "file": workbook.name, "sha256": workbook_digest}
    for name, expected_value in expected.items():
        if source.get(name) != expected_value:
            raise MendelPrimeError(
                "hash_mismatch" if name == "sha256" else "source_mismatch",
                field=f"knowledgebase.source.{name}",
            )


def _table_sources(value, version, name, table, required, *, optional=(), gene_field=None):
    rows = value.get(name)
    if not isinstance(rows, list):
        raise MendelPrimeError("array_required", field=name)
    if not rows:
        raise MendelPrimeError("nonempty_array_required", field=name)
    if any(not isinstance(row, dict) for row in rows):
        raise MendelPrimeError("object_rows_required", field=name)
    sources = []
    for index, row in enumerate(rows, 1):
        if not set(required).issubset(row) or set(row) - set(required) - set(optional):
            raise MendelPrimeError("schema_mismatch", field=f"{name}[{index - 1}]")
        for field in (*required, *(field for field in optional if field in row)):
            if not isinstance(row.get(field), str):
                raise MendelPrimeError("text_required", field=f"{name}[{index - 1}].{field}")
        genes = (
            tuple(token for part in row[gene_field].split(",") if _GENE.fullmatch(token := part.strip()))
            if gene_field
            else ()
        )
        content = json.dumps(row, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        sources.append(
            Source(
                provider="ai-cura",
                kind="guideline",
                locator=f"{_DOI_URL}#supplementary-table-{table}/row-{index}",
                version=version,
                content=content,
                genes=genes,
            )
        )
    return tuple(sources)


__all__ = ["open_resources", "close_resources", "gene_context", "protein_context"]
