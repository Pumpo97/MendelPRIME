"""Offline VEP, gnomAD, REVEL and SpliceAI readers with source-bound records."""

import gzip
import hashlib
import http.client
import json
import math
import os
import re
import socket
import ssl
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pysam

from mendelprime.domain import MendelPrimeError, Source, validate_gene, variant_key, variant_text
from mendelprime.knowledge.literature import _RETRYABLE_HTTP, _http_retry_delay
from mendelprime.runtime.http import urlopen
from mendelprime.runtime.process import batch_outcome, complete_batch, fail_batch, run_process
from mendelprime.runtime.storage import atomic_replace_bytes, private_directory, read_regular_bounded

MINIMUM_OVER20_FRACTION = 0.9
_MAX_VEP_OUTPUT = 64 * 1024 * 1024
_MAX_VEP_LINE = 4 * 1024 * 1024
_MAX_COVERAGE_HEADER = 16 * 1024
_SEQUENCE = re.compile(r"[ACGTN]+", re.ASCII)
_MISSING = {"", ".", "NA"}
_REVEL_COLUMNS = (
    "chr",
    "hg19_pos",
    "grch38_pos",
    "ref",
    "alt",
    "aaref",
    "aaalt",
    "REVEL",
    "Ensembl_transcriptid",
)
_VCF_COLUMNS = ("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")


def public_resource(url, cache):
    """Cache bounded public resource bytes, not model judgments or patient labels."""
    allowed = {"grch37.rest.ensembl.org", "search.clinicalgenome.org", "www.ncbi.nlm.nih.gov"}
    parsed = urllib.parse.urlsplit(url)
    if (
        parsed.scheme != "https"
        or parsed.hostname not in allowed
        or (
            parsed.hostname == "www.ncbi.nlm.nih.gov"
            and not re.fullmatch(r"/(?:sites/)?books/NBK\d+/", parsed.path)
        )
    ):
        raise MendelPrimeError("public_resource_host_not_allowed", category="resource")
    private_directory(cache)
    path = cache / (hashlib.sha256(url.encode()).hexdigest() + ".data")
    try:
        return read_regular_bounded(path, 20_000_000, follow_symlinks=False)[0]
    except FileNotFoundError:
        request = urllib.request.Request(url, headers={"User-Agent": "MendelPRIME-evidence/1.0"})
        for attempt in range(2):
            try:
                with urlopen(request, timeout=45) as response:
                    redirected = urllib.parse.urlsplit(response.geturl())
                    if (
                        redirected.scheme != "https"
                        or redirected.hostname not in allowed
                        or (
                            redirected.hostname == "www.ncbi.nlm.nih.gov"
                            and not re.fullmatch(r"/(?:sites/)?books/NBK\d+/", redirected.path)
                        )
                    ):
                        raise MendelPrimeError("public_resource_redirect_not_allowed", category="resource")
                    data = response.read(20_000_001)
                if len(data) > 20_000_000:
                    raise MendelPrimeError("public_resource_too_large", category="resource")
                break
            except MendelPrimeError:
                raise
            except Exception as error:
                failure = _public_resource_failure(error, attempt + 1, parsed.hostname)
                if isinstance(error, urllib.error.HTTPError):
                    error.close()
                if failure.retryable and attempt == 0:
                    time.sleep(_http_retry_delay(error, attempt))
                    continue
                raise failure from error
        atomic_replace_bytes(path, data, sync=False)
        return data


def _public_resource_failure(error, attempts, host):
    if isinstance(error, urllib.error.HTTPError):
        return MendelPrimeError(
            "public_resource_http_error",
            category="resource",
            attempts=attempts,
            http_status=error.code,
            retryable=error.code in _RETRYABLE_HTTP,
            source=host,
        )
    reason = error.reason if isinstance(error, urllib.error.URLError) else error
    if isinstance(reason, ssl.SSLCertVerificationError):
        code, retryable = "public_resource_tls_error", False
    elif isinstance(reason, (TimeoutError, socket.timeout)):
        code, retryable = "public_resource_timeout", True
    elif isinstance(reason, ssl.SSLError):
        code, retryable = "public_resource_tls_error", False
    else:
        retryable = isinstance(
            reason,
            (ConnectionError, http.client.IncompleteRead, http.client.RemoteDisconnected, OSError),
        ) or isinstance(error, urllib.error.URLError)
        code = "public_resource_transport_error"
    return MendelPrimeError(code, category="resource", attempts=attempts, retryable=retryable, source=host)


def mane_transcripts(sources, gene):
    """MANE supplies versioned RNA candidates, never assembly coordinates or HGVS."""
    candidates = {}
    for source in sources:
        if source.provider != "mane":
            continue
        if source.kind != "annotation" or source.genes != (gene,) or source.variant_keys or source.case_key:
            raise MendelPrimeError("mane_scope_invalid", field="transcript")
        fact = json.loads(source.content)
        header = fact["header"]
        if len(set(header)) != len(header) or not {"symbol", "RefSeq_nuc", "MANE_status"} <= set(header):
            raise MendelPrimeError("mane_schema_invalid", field="transcript")
        row = dict(zip(header, fact["row"].replace("\\t", "\t").split("\t"), strict=True))
        transcript = row["RefSeq_nuc"]
        if row["symbol"] != gene or not re.fullmatch(r"NM_\d+\.\d+", transcript):
            raise MendelPrimeError("mane_identity_invalid", field="transcript")
        if row["MANE_status"] in {"MANE Select", "MANE Plus Clinical"}:
            candidates.setdefault(transcript, source)
    return candidates


def merged_refseq_annotation(settings, work_root, target, transcripts, *, gene_aliases=None, preferred=()):
    """Reannotate an allele on an exact resource version; never transplant HGVS."""
    key = _key(target["key"])
    aliases, alias_digest = gene_aliases or ({}, "")
    candidates = target.get("mane_candidates", {})
    if candidates and candidates != mane_transcripts(candidates.values(), target["gene"]):
        raise MendelPrimeError("mane_candidate_binding_invalid", field="transcript")
    if key[0] != "GRCh37":
        raise MendelPrimeError("post_annotation_assembly_unsupported", category="resource")
    representations = (target, *target.get("aliases", ()))
    wanted = tuple(
        dict.fromkeys(
            transcript
            for row in (tuple({"transcript": value} for value in preferred) if preferred else representations)
            for transcript in transcripts.get(
                (target["gene"], str(row.get("transcript", "")).split(".")[0], key[1]), ()
            )
            if not preferred or "." not in row["transcript"] or transcript == row["transcript"]
        )
    )
    if len(wanted) != 1:
        raise MendelPrimeError("post_refseq_missing_or_ambiguous", category="resource")
    position, reference, alternate = key[2:]
    while reference and alternate and reference[0] == alternate[0]:
        position, reference, alternate = position + 1, reference[1:], alternate[1:]
    while reference and alternate and reference[-1] == alternate[-1]:
        reference, alternate = reference[:-1], alternate[:-1]
    region = f"{key[1]}:{position}-{position + len(reference) - 1}:1/{alternate or '-'}"
    url = (
        "https://grch37.rest.ensembl.org/vep/human/region/"
        + urllib.parse.quote(region, safe=":/-")
        + "?content-type=application/json&hgvs=1&numbers=1&merged=1"
    )
    data = public_resource(url, work_root / "public")
    rows = json.loads(data)
    if not isinstance(rows, list) or len(rows) != 1:
        raise MendelPrimeError("post_vep_response_identity_invalid", category="resource")
    raw = rows[0]
    if (
        raw.get("assembly_name"),
        str(raw.get("seq_region_name")),
        raw.get("start"),
        raw.get("end"),
        raw.get("strand"),
    ) != (key[0], key[1], position, position + len(reference) - 1, 1):
        raise MendelPrimeError("post_vep_allele_mismatch", category="resource")
    if raw.get("allele_string") != f"{reference or '-'}/{alternate or '-'}":
        raise MendelPrimeError("post_vep_reference_mismatch", category="resource")
    selected = [
        row
        for row in raw.get("transcript_consequences", ())
        if aliases.get(row.get("gene_symbol"), row.get("gene_symbol")) == target["gene"]
        and row.get("transcript_id") == wanted[0]
        and str(row.get("hgvsc", "")).split(":")[0] == wanted[0]
    ]
    if len(selected) != 1:
        raise MendelPrimeError("post_vep_exact_refseq_unavailable", category="resource")
    row = selected[0]
    canonical_splice = bool(
        {"splice_acceptor_variant", "splice_donor_variant"}.intersection(row.get("consequence_terms", ()))
    )
    if row.get("biotype") != "protein_coding" or (not row.get("hgvsp") and not canonical_splice):
        raise MendelPrimeError("post_vep_coding_effect_unavailable", category="resource")
    record = {
        "record_type": "verified_refseq_allele_reannotation_v1",
        "assembly": key[0],
        "gene": target["gene"],
        "variant_key": variant_text(key),
        "observed_transcripts": list(dict.fromkeys(r.get("transcript") for r in representations)),
        "assessed_transcript": wanted[0],
        "selected_transcript": row,
        "resolution": "independent GRCh37 allele annotation on exact AutoPVS1 resource transcript; HGVS not copied across versions",
    }
    if row["gene_symbol"] != target["gene"]:
        record["gene_alias_binding"] = {
            "reported_symbol": row["gene_symbol"],
            "canonical_symbol": target["gene"],
            "resource_sha256": alias_digest,
        }
    if preferred:
        record["requested_transcript_candidates"] = list(preferred)
        record["selection_scope"] = (
            "Transcript identity candidate; not proof of disease relevance or allele effect."
        )
        record["transcript_selection_source_ids"] = (
            [target["mane_candidates"][wanted[0]].source_id]
            if wanted[0] in target.get("mane_candidates", {})
            else []
        )
    source = Source(
        provider="ensembl-vep",
        kind="annotation",
        locator=url,
        version="GRCh37 REST;sha256="
        + hashlib.sha256(data).hexdigest()
        + (";transcript=" + wanted[0] if preferred else ""),
        genes=(target["gene"],),
        variant_keys=(key,),
        content=json.dumps(record, ensure_ascii=False, sort_keys=True, separators=(",", ":")),
    )
    return {
        "transcript_id": wanted[0],
        "hgvsc": row["hgvsc"],
        "hgvsp": row.get("hgvsp") or "",
        "consequence_terms": tuple(row.get("consequence_terms", ())),
        "source": source,
        "protein_id": row.get("protein_id"),
        "exon": row.get("exon"),
        "intron": row.get("intron"),
        "canonical": row.get("canonical") == 1,
    }


def annotate(settings, work_root, targets, *, enabled=True, gene_aliases=None):
    """Run one private VEP batch and return one explicit outcome per target."""
    if enabled:
        _validate_vep_settings(settings, work_root)
    if type(targets) is not tuple:
        raise MendelPrimeError("annotation_target_tuple_required", field="vep.targets")
    normalized = tuple(_annotation_target(value) for value in targets)
    if not enabled:
        return tuple(
            batch_outcome(target, "not_applicable", "vep_disabled_by_configuration") for target in normalized
        )
    if not normalized:
        return ()
    outcomes = [None] * len(normalized)
    eligible = []
    for index, target in enumerate(normalized):
        if target["variant_key"][0] != "GRCh37":
            outcomes[index] = batch_outcome(target, "unknown", "vep_unsupported_assembly")
        else:
            eligible.append((index, target, f"MP{len(eligible)}"))
    if not eligible:
        return complete_batch(normalized, outcomes, "vep_internal_outcome_missing")
    resource_error = _vep_resource_error(settings, work_root)
    if resource_error:
        fail_batch(outcomes, eligible, resource_error)
        return complete_batch(normalized, outcomes, "vep_internal_outcome_missing")
    try:
        with tempfile.TemporaryDirectory(prefix="mendelprime-vep-", dir=work_root) as temporary:
            root = Path(temporary)
            root.chmod(0o700)
            input_vcf, output_jsonl = root / "input.vcf", root / "output.jsonl"
            indexed = tuple((identifier, target) for _, target, identifier in eligible)
            _write_vep_vcf(input_vcf, indexed)
            receipt = run_process(
                command=_vep_command(settings, input_vcf, output_jsonl, root),
                root=root,
                timeout_seconds=settings["vep_timeout_seconds"],
                artifacts=(output_jsonl,),
            )
            if receipt["error_code"]:
                fail_batch(outcomes, eligible, receipt["error_code"])
            else:
                try:
                    rows = _read_vep_rows(output_jsonl, frozenset(identifier for identifier, _ in indexed))
                except MendelPrimeError as exc:
                    fail_batch(outcomes, eligible, exc.code)
                else:
                    for index, target, identifier in eligible:
                        raw = rows.get(identifier)
                        if raw is None:
                            outcomes[index] = batch_outcome(target, "unknown", "vep_result_missing")
                            continue
                        issue = _vep_identity_error(target, identifier, raw)
                        if issue:
                            outcomes[index] = batch_outcome(target, "unknown", issue)
                            continue
                        selected, status, issue = _select_annotation(target, raw, gene_aliases)
                        if selected is None:
                            outcomes[index] = batch_outcome(target, status, issue)
                            continue
                        try:
                            annotation = _build_annotation(target, selected, raw, settings)
                        except (MendelPrimeError, TypeError, ValueError):
                            outcomes[index] = batch_outcome(
                                target, "failed", "vep_annotation_schema_mismatch"
                            )
                        else:
                            outcomes[index] = batch_outcome(target, "annotated", annotation=annotation)
    except OSError:
        fail_batch(outcomes, eligible, "vep_private_workdir_failed")
    return complete_batch(normalized, outcomes, "vep_internal_outcome_missing")


def _validate_vep_settings(settings, work_root):
    command = settings.get("vep_command")
    if (
        type(command) not in {tuple, list}
        or not command
        or any(not isinstance(part, str) or not part or "\x00" in part for part in command)
    ):
        raise MendelPrimeError("nonempty_argv_required", field="vep.command")
    for field in ("vep_cache", "vep_fasta"):
        if not isinstance(settings.get(field), Path) or not settings[field].is_absolute():
            raise MendelPrimeError("absolute_path_required", field=f"vep.{field[4:]}")
    if settings.get("assembly") != "GRCh37":
        raise MendelPrimeError("grch37_required", field="vep.assembly")
    for field in ("vep_cache_version", "vep_tool_version"):
        value = settings.get(field)
        if not isinstance(value, str) or not value.strip() or value != value.strip():
            raise MendelPrimeError("version_required", field=f"vep.{field[4:]}")
    timeout = settings.get("vep_timeout_seconds")
    if type(timeout) not in {int, float} or not math.isfinite(timeout) or timeout <= 0:
        raise MendelPrimeError("positive_finite_timeout_required", field="vep.timeout_seconds")
    if not isinstance(work_root, Path) or not work_root.is_absolute():
        raise MendelPrimeError("absolute_path_required", field="vep.work_root")


def _annotation_target(value):
    if type(value) is not dict or set(value) != {"gene", "variant_key", "transcript_id"}:
        raise MendelPrimeError("annotation_target_required", field="vep.targets")
    validate_gene(value["gene"])
    key = _key(value["variant_key"], "annotation_target.variant_key")
    transcript = value["transcript_id"]
    if transcript is not None and (
        not isinstance(transcript, str)
        or not transcript
        or transcript != transcript.strip()
        or ":" in transcript
        or any(char.isspace() for char in transcript)
    ):
        raise MendelPrimeError("transcript_identifier_required", field="annotation_target.transcript_id")
    return {"gene": value["gene"], "variant_key": key, "transcript_id": transcript}


def _vep_command(settings, input_vcf, output_jsonl, workdir):
    prefix = [part.replace("{workdir}", str(workdir)) for part in settings["vep_command"]]
    return [
        *prefix,
        "--offline",
        "--cache",
        "--cache_version",
        settings["vep_cache_version"],
        "--assembly",
        settings["assembly"],
        "--dir_cache",
        str(settings["vep_cache"]),
        "--fasta",
        str(settings["vep_fasta"]),
        "--input_file",
        str(input_vcf),
        "--output_file",
        str(output_jsonl),
        "--json",
        "--hgvs",
        "--symbol",
        "--numbers",
        "--protein",
        "--transcript_version",
        "--canonical",
        "--force_overwrite",
        "--no_stats",
    ]


def _vep_resource_error(settings, work_root):
    try:
        if not work_root.is_dir() or not os.access(work_root, os.W_OK | os.X_OK):
            return "vep_work_root_unavailable"
        if not settings["vep_cache"].is_dir() or not os.access(settings["vep_cache"], os.R_OK | os.X_OK):
            return "vep_cache_unavailable"
        fasta = settings["vep_fasta"]
        if not fasta.is_file() or fasta.stat().st_size == 0:
            return "vep_fasta_unavailable"
        index = Path(f"{fasta}.fai")
        if not index.is_file() or index.stat().st_size == 0:
            return "vep_fasta_index_unavailable"
        if not os.access(fasta, os.R_OK) or not os.access(index, os.R_OK):
            return "vep_fasta_unavailable"
    except OSError:
        return "vep_resource_unavailable"
    return ""


def _write_vep_vcf(path, indexed):
    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for identifier, target in indexed:
        _, chromosome, position, reference, alternate = target["variant_key"]
        lines.append("\t".join((chromosome, str(position), identifier, reference, alternate, ".", ".", ".")))
    with path.open("x", encoding="ascii", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")
    path.chmod(0o600)


def _read_vep_rows(path, expected):
    try:
        if not path.is_file() or path.stat().st_size == 0:
            raise MendelPrimeError("vep_output_missing")
        if path.stat().st_size > _MAX_VEP_OUTPUT:
            raise MendelPrimeError("vep_output_too_large")
        rows, line_count = {}, 0
        with path.open("rb") as handle:
            for encoded in handle:
                if not encoded.strip():
                    continue
                line_count += 1
                if len(encoded) > _MAX_VEP_LINE:
                    raise MendelPrimeError("vep_output_bounds_exceeded")
                try:
                    raw = json.loads(encoded.decode("utf-8"))
                except (UnicodeDecodeError, json.JSONDecodeError):
                    raise MendelPrimeError("vep_output_parse_failed") from None
                if not isinstance(raw, dict) or not isinstance(raw.get("input", ""), str):
                    raise MendelPrimeError("vep_output_schema_mismatch")
                tokens, raw_id = raw["input"].split(), raw.get("id", "")
                if (
                    len(tokens) < 3
                    or not tokens[2]
                    or (raw_id not in (None, "") and not isinstance(raw_id, str))
                    or raw_id not in (None, "", tokens[2])
                ):
                    raise MendelPrimeError("vep_malformed_id")
                identifier = raw_id or tokens[2]
                if not identifier.startswith("MP") or not identifier[2:].isdigit():
                    raise MendelPrimeError("vep_malformed_id")
                if identifier not in expected:
                    raise MendelPrimeError("vep_unknown_id")
                if identifier in rows:
                    raise MendelPrimeError("vep_duplicate_id")
                rows[identifier] = raw
                if line_count > len(expected):
                    raise MendelPrimeError("vep_output_bounds_exceeded")
        return rows
    except MendelPrimeError:
        raise
    except OSError:
        raise MendelPrimeError("vep_output_read_failed") from None


def _vep_identity_error(target, identifier, raw):
    key = target["variant_key"]
    if raw.get("assembly_name") != key[0]:
        return "vep_assembly_mismatch"
    tokens = raw.get("input", "").split()
    if len(tokens) < 5:
        return "vep_input_missing"
    if tokens[2] != identifier or raw.get("id") not in (None, "", identifier):
        return "vep_input_id_mismatch"
    try:
        position = int(tokens[1])
    except (TypeError, ValueError):
        return "vep_input_allele_mismatch"
    if (tokens[0], position, tokens[3], tokens[4]) != key[1:]:
        return "vep_input_allele_mismatch"
    return ""


def _select_annotation(target, raw, gene_aliases=None):
    consequences = raw.get("transcript_consequences")
    if not isinstance(consequences, list):
        return None, "unknown", "vep_transcript_consequences_missing"
    rows = [
        item for item in consequences if isinstance(item, dict)
        and (gene_aliases or {}).get(item.get("gene_symbol"), item.get("gene_symbol")) == target["gene"]
    ]
    if not rows:
        return None, "unknown", "vep_gene_not_found"
    candidates = [item for item in rows if isinstance(item.get("hgvsc"), str) and item["hgvsc"]]
    if target["transcript_id"]:
        selected = [item for item in candidates if item["hgvsc"].split(":", 1)[0] == target["transcript_id"]]
        missing, ambiguous = "vep_exact_transcript_not_found", "vep_exact_transcript_ambiguous"
    else:
        selected = [item for item in candidates if item.get("canonical") == 1]
        missing, ambiguous = "vep_canonical_transcript_not_found", "vep_canonical_transcript_ambiguous"
    if not selected:
        return None, "unknown", missing
    if len(selected) != 1:
        return None, "ambiguous", ambiguous
    return selected[0], "annotated", ""


def _build_annotation(target, selected, raw, settings):
    hgvsc = selected.get("hgvsc")
    if not isinstance(hgvsc, str) or not hgvsc:
        raise MendelPrimeError("text_required", field="vep.hgvsc")
    hgvsp, protein_id = selected.get("hgvsp", "") or "", selected.get("protein_id", "") or ""
    terms, exon, intron = selected.get("consequence_terms", []), selected.get("exon"), selected.get("intron")
    if not isinstance(hgvsp, str) or not isinstance(protein_id, str):
        raise MendelPrimeError("text_fields_required", field="vep.annotation")
    if not isinstance(terms, list) or any(not isinstance(term, str) or not term for term in terms):
        raise MendelPrimeError("text_list_required", field="vep.consequence_terms")
    if any(value is not None and (not isinstance(value, str) or not value) for value in (exon, intron)):
        raise MendelPrimeError("optional_text_fields_required", field="vep.annotation")
    tokens = raw["input"].split()
    transcript = hgvsc.split(":", 1)[0]
    content = json.dumps(
        {
            "assembly": raw["assembly_name"],
            "gene": target["gene"],
            "input_allele": {
                "chromosome": tokens[0],
                "position": int(tokens[1]),
                "reference": tokens[3],
                "alternate": tokens[4],
            },
            "selected_transcript": selected,
        },
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    key = target["variant_key"]
    scope = ":".join((*map(str, key), target["gene"], transcript))
    source = Source(
        provider="ensembl-vep",
        kind="annotation",
        locator=f"method:ensembl-vep/annotation#{scope}",
        version=(
            f"configured_tool={settings['vep_tool_version']};"
            f"configured_cache={settings['vep_cache_version']};configured_assembly={settings['assembly']}"
        ),
        content=content,
        genes=(target["gene"],),
        variant_keys=(key,),
    )
    return {
        "transcript_id": transcript,
        "hgvsc": hgvsc,
        "hgvsp": hgvsp,
        "protein_id": protein_id,
        "consequence_terms": tuple(terms),
        "exon": exon,
        "intron": intron,
        "canonical": selected.get("canonical") == 1,
        "source": source,
    }


def _key(value, field="variant"):
    if type(value) is not tuple or len(value) != 5:
        raise MendelPrimeError("variant_key_required", field=field)
    return variant_key(*value)


def _key_text(key):
    return variant_text(key)


def population(settings, reference, gene, key, *, details=False):
    """Return the exact gnomAD observation or a coverage-qualified absence."""
    validate_gene(gene)
    key = _key(key, "population.key")
    if key[0] != settings["assembly"]:
        return _population_failure(key, "identity_mismatch", "population_assembly_mismatch")
    sites = settings["gnomad_sites"]
    if not _has_tabix(sites):
        return _population_failure(key, "resource_unavailable", "population_sites_unavailable")
    start, end = key[2], key[2] + len(key[3]) - 1
    if reference is not None:
        resolved = reference.resolve(key)
        if resolved["status"] == "sequence_verified":
            start, end = resolved["equivalent_region"]
    rows, _, failure = _tabix_fetch(sites, key[1], start, end, "population_sites")
    if failure:
        return _population_failure(key, *failure)
    try:
        match = _population_site_match(rows, key, reference)
    except ValueError:
        return _population_failure(key, "schema_mismatch", "population_frequency_invalid")
    if match is not None:
        status = "found" if match["record_filter"] in {"PASS", "."} else "found_filtered"
        if status == "found" and match["allele_frequency"] is None:
            status = "frequency_unavailable"
        measurements = {
            name: match[name]
            for name in ("record_filter", "allele_count", "allele_number", "allele_frequency")
        }
        measurements["coverage_status"] = "not_assessed"
        coverage_rows, extra = (), {}
        if details:
            coverage = _population_coverage(settings, gene, key)
            measurements.update(
                {name: value for name, value in coverage.items() if name.startswith("coverage_")}
            )
            if coverage.get("source"):
                coverage_rows = tuple(json.loads(coverage["source"].content)["raw_records"]["coverage"])
            columns = match["raw"].split("\t")
            index = columns[4].split(",").index(match["key"][4])
            for name, value in _info(columns[7]).items():
                if not re.fullmatch(
                    r"(?:AC|AN|AF|nhomalt|faf95|faf99)(?:_(?:afr|amr|asj|eas|fin|nfe|oth|sas|popmax))?", name
                ):
                    continue
                values = value.split(",")
                try:
                    number = float(values[index] if len(values) > 1 else values[0])
                except (ValueError, IndexError):
                    continue
                if math.isfinite(number):
                    extra[name] = number
            measurements["population_details"] = extra
        source = _population_source(
            settings,
            gene,
            key,
            status=status,
            record_key=match["key"],
            match_method=match["method"],
            raw_sites=(match["raw"],),
            raw_coverage=coverage_rows,
            measurements=measurements,
        )
        result = _population_record(
            key,
            status,
            record_key=match["key"],
            match_method=match["method"],
            record_filter=match["record_filter"],
            allele_count=match["allele_count"],
            allele_number=match["allele_number"],
            allele_frequency=match["allele_frequency"],
            source=source,
            error_code="population_frequency_unavailable" if status == "frequency_unavailable" else "",
        )
        return (
            {
                **result,
                **{
                    k: v
                    for k, v in measurements.items()
                    if k.startswith("coverage_") or k == "population_details"
                },
            }
            if details
            else result
        )
    return _population_coverage(settings, gene, key)


def _population_record(
    key,
    status,
    *,
    record_key=None,
    match_method="none",
    record_filter="",
    allele_count=None,
    allele_number=None,
    allele_frequency=None,
    coverage_status="not_assessed",
    coverage_adequate=None,
    coverage_position_count=None,
    coverage_expected_count=None,
    coverage_over20_min=None,
    source=None,
    error_code="",
):
    return {
        "original_key": key,
        "status": status,
        "record_key": record_key,
        "match_method": match_method,
        "record_filter": record_filter,
        "allele_count": allele_count,
        "allele_number": allele_number,
        "allele_frequency": allele_frequency,
        "coverage_status": coverage_status,
        "coverage_adequate": coverage_adequate,
        "coverage_position_count": coverage_position_count,
        "coverage_expected_count": coverage_expected_count,
        "coverage_over20_min": coverage_over20_min,
        "source": source,
        "error_code": error_code,
    }


def _population_failure(key, status, error_code, *, coverage_status="not_assessed"):
    return _population_record(key, status, coverage_status=coverage_status, error_code=error_code)


def _has_tabix(path):
    try:
        return path.is_file() and any(Path(f"{path}{suffix}").is_file() for suffix in (".tbi", ".csi"))
    except OSError:
        return False


def _tabix_fetch(path, chromosome, start, end, resource):
    source = None
    try:
        source = pysam.TabixFile(str(path))
        headers = tuple(source.header)
        contig = next((value for value in (chromosome, f"chr{chromosome}") if value in source.contigs), None)
        if contig is None:
            return (), headers, ("contig_unavailable", f"{resource}_contig_unavailable")
        return tuple(source.fetch(contig, max(0, start - 1), end)), headers, None
    except (OSError, ValueError):
        return (), (), ("query_failed", f"{resource}_query_failed")
    finally:
        if source is not None:
            source.close()


def _info(raw):
    values = {}
    for item in raw.split(";"):
        if not item or "=" not in item:
            continue
        name, value = item.split("=", 1)
        if not name or name in values:
            raise ValueError
        values[name] = value
    return values


def _frequency_values(info, index):
    try:
        counts, numbers = info["AC"].split(","), info["AN"].split(",")
        raw_frequency = info.get("AF", ".")
        frequencies = raw_frequency.split(",")
        if len(counts) <= index or len(numbers) != 1 or (raw_frequency != "." and len(frequencies) <= index):
            raise ValueError
        value = raw_frequency if raw_frequency == "." else frequencies[index]
        count, number = int(counts[index]), int(numbers[0])
        frequency = None if value == "." else float(value)
    except (KeyError, ValueError, OverflowError):
        raise ValueError from None
    if (
        number < 0
        or count < 0
        or count > number
        or (frequency is not None and (not math.isfinite(frequency) or not 0 <= frequency <= 1))
    ):
        raise ValueError
    if number == 0:
        if frequency not in (None, 0):
            raise ValueError
        # No sampled denominator means no measured frequency, including a reported zero.
        frequency = None
    elif frequency is None:
        raise ValueError
    return count, number, frequency


def _population_site_match(rows, key, reference):
    for raw in rows:
        fields = raw.split("\t")
        if len(fields) < 8:
            raise ValueError
        chromosome = fields[0][3:] if fields[0].startswith("chr") else fields[0]
        if chromosome != key[1]:
            raise ValueError
        try:
            position = int(fields[1])
        except ValueError:
            raise ValueError from None
        if position < 1 or not fields[3] or not fields[4] or not fields[6]:
            raise ValueError
        for index, alternate in enumerate(fields[4].split(",")):
            exact = (position, fields[3], alternate) == (key[2], key[3], key[4])
            if not exact and (
                reference is None or not _SEQUENCE.fullmatch(fields[3]) or not _SEQUENCE.fullmatch(alternate)
            ):
                continue
            try:
                record_key = variant_key(key[0], chromosome, position, fields[3], alternate)
            except MendelPrimeError:
                if exact:
                    raise ValueError from None
                continue
            method = (
                "exact" if exact else "sequence_equivalent" if reference.equivalent(key, record_key) else None
            )
            if method:
                count, number, frequency = _frequency_values(_info(fields[7]), index)
                return {
                    "raw": raw,
                    "key": record_key,
                    "method": method,
                    "record_filter": fields[6],
                    "allele_count": count,
                    "allele_number": number,
                    "allele_frequency": frequency,
                }
    return None


def _affected_positions(key):
    # Only a shared leading base is padding; an MNV may change its first base.
    padded = len(key[3]) > 1 and key[3][0] == key[4][0]
    start = key[2] + int(padded)
    end = key[2] + len(key[3]) - 1 if len(key[4]) <= len(key[3]) else key[2] + len(key[3])
    return tuple(range(start, end + 1))


def _coverage_columns(path, headers):
    if headers:
        raw = headers[-1]
    else:
        try:
            with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
                raw = handle.readline(_MAX_COVERAGE_HEADER + 1)
        except (OSError, EOFError, UnicodeError):
            raise OSError from None
    if not raw or len(raw) > _MAX_COVERAGE_HEADER:
        raise ValueError
    return tuple(raw.strip().lstrip("#").split())


def _population_coverage(settings, gene, key):
    path = settings["gnomad_coverage"]
    if not _has_tabix(path):
        return _population_failure(
            key, "resource_unavailable", "population_coverage_unavailable", coverage_status="unavailable"
        )
    expected = _affected_positions(key)
    rows, headers, failure = _tabix_fetch(path, key[1], expected[0], expected[-1], "population_coverage")
    if failure:
        status, code = failure
        coverage_status = "contig_unavailable" if status == "contig_unavailable" else "query_failed"
        return _population_failure(key, status, code, coverage_status=coverage_status)
    try:
        names = _coverage_columns(path, headers)
    except OSError:
        return _population_failure(
            key, "query_failed", "population_coverage_header_failed", coverage_status="query_failed"
        )
    except ValueError:
        return _population_failure(
            key,
            "schema_mismatch",
            "population_coverage_schema_unrecognized",
            coverage_status="schema_mismatch",
        )
    if "over_20" not in names:
        return _population_failure(
            key,
            "schema_mismatch",
            "population_coverage_schema_unrecognized",
            coverage_status="schema_mismatch",
        )
    column, observations, expected_set = names.index("over_20"), {}, set(expected)
    try:
        for raw in rows:
            fields = raw.split("\t")
            if len(fields) <= column:
                raise ValueError
            chromosome = fields[0][3:] if fields[0].startswith("chr") else fields[0]
            position, value = int(fields[1]), float(fields[column])
            if (
                chromosome != key[1]
                or position not in expected_set
                or position in observations
                or not math.isfinite(value)
                or not 0 <= value <= 1
            ):
                raise ValueError
            observations[position] = value
    except (ValueError, OverflowError):
        return _population_failure(
            key, "schema_mismatch", "population_coverage_invalid", coverage_status="schema_mismatch"
        )
    measured = expected_set <= observations.keys()
    minimum = min(observations.values()) if observations else None
    adequate = measured and all(observations[position] >= MINIMUM_OVER20_FRACTION for position in expected)
    status = "not_found_covered" if adequate else "coverage_insufficient"
    coverage_status = "measured" if measured else "incomplete"
    source = _population_source(
        settings,
        gene,
        key,
        status=status,
        record_key=None,
        match_method="none",
        raw_sites=(),
        raw_coverage=rows,
        measurements={
            "coverage_status": coverage_status,
            "coverage_adequate": adequate,
            "coverage_position_count": len(observations),
            "coverage_expected_count": len(expected),
            "coverage_over20_min": minimum,
        },
    )
    return _population_record(
        key,
        status,
        coverage_status=coverage_status,
        coverage_adequate=adequate,
        coverage_position_count=len(observations),
        coverage_expected_count=len(expected),
        coverage_over20_min=minimum,
        source=source,
    )


def _population_source(
    settings,
    gene,
    key,
    *,
    status,
    record_key,
    match_method,
    raw_sites,
    raw_coverage,
    measurements,
):
    content = json.dumps(
        {
            "record_type": "gnomad_population_observation",
            "original_key": _key_text(key),
            "record_key": _key_text(record_key) if record_key else None,
            "match_method": match_method,
            "status": status,
            "raw_records": {"sites": list(raw_sites), "coverage": list(raw_coverage)},
            "measurements": measurements,
        },
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    fragment = hashlib.sha256(content.encode()).hexdigest()[:16]
    return Source(
        provider="gnomad",
        kind="population",
        locator=f"resource:gnomad/population#{fragment}",
        version=settings["population_release"],
        content=content,
        genes=(gene,),
        variant_keys=(key,),
    )


def predictions(settings, gene, key):
    """Query exact-allele REVEL and exact-gene SpliceAI resources."""
    validate_gene(gene)
    key = _key(key, "predictions.key")
    if key[0] != settings["assembly"]:
        return {
            "gene": gene,
            "variant_key": key,
            "revel": _prediction_failure(key, "identity_mismatch", "revel_assembly_mismatch"),
            "spliceai": _prediction_failure(key, "identity_mismatch", "spliceai_assembly_mismatch"),
        }
    return {
        "gene": gene,
        "variant_key": key,
        "revel": _lookup_revel(settings, gene, key),
        "spliceai": _lookup_spliceai(settings, gene, key),
    }


def _prediction_lookup(key, status, records=(), source=None, error_code=""):
    return {
        "variant_key": key,
        "status": status,
        "records": tuple(records),
        "source": source,
        "error_code": error_code,
    }


def _prediction_failure(key, status, error_code):
    return _prediction_lookup(key, status, error_code=error_code)


def _prediction_fetch(path, key, resource):
    if not _has_tabix(path):
        return (), (), ("resource_unavailable", f"{resource}_resource_unavailable")
    rows, headers, failure = _tabix_fetch(path, key[1], key[2], key[2], resource)
    return rows, headers, failure


def _header_columns(headers, expected):
    if not any(tuple(raw.lstrip("#").split("\t")) == expected for raw in reversed(headers)):
        raise ValueError


def _unit_score(value):
    try:
        score = float(value)
    except (ValueError, OverflowError):
        raise ValueError from None
    if not math.isfinite(score) or not 0 <= score <= 1:
        raise ValueError
    return score


def _offset(value):
    if not re.fullmatch(r"[+-]?\d+", value, re.ASCII):
        raise ValueError
    return int(value)


def _lookup_revel(settings, gene, key):
    rows, headers, failure = _prediction_fetch(settings["revel"], key, "revel")
    if failure:
        return _prediction_failure(key, *failure)
    try:
        _header_columns(headers, _REVEL_COLUMNS)
        records, matching, missing = [], [], False
        for raw in rows:
            fields = raw.split("\t")
            if len(fields) != len(_REVEL_COLUMNS):
                raise ValueError
            position = int(fields[1])
            chromosome = fields[0][3:] if fields[0].startswith("chr") else fields[0]
            if chromosome != key[1] or position != key[2]:
                raise ValueError
            if fields[3] != key[3] or fields[4] != key[4]:
                continue
            matching.append(raw)
            if fields[7] in _MISSING:
                missing = True
                continue
            if any(
                not isinstance(value, str)
                or not value
                or any(ord(character) < 32 or ord(character) == 127 for character in value)
                for value in (fields[5], fields[6], fields[8])
            ):
                raise ValueError
            records.append(
                {
                    "variant_key": key,
                    "score": _unit_score(fields[7]),
                    "transcript": fields[8],
                    "amino_acid_reference": fields[5],
                    "amino_acid_alternate": fields[6],
                }
            )
    except (ValueError, OverflowError):
        return _prediction_failure(key, "schema_mismatch", "revel_record_invalid")
    if records or missing:
        source = _prediction_source(
            "revel",
            settings["revel_release"],
            "revel_precomputed_prediction",
            gene,
            key,
            headers,
            {"raw_records": tuple(matching)},
            {"gene_match_basis": "unavailable", "prediction_missing": missing},
        )
        return _prediction_lookup(
            key,
            "found" if records else "missing_prediction",
            records,
            source,
            "" if records else "revel_prediction_missing",
        )
    return _prediction_lookup(key, "not_found")


def _lookup_spliceai(settings, gene, key):
    resource_class = "snv" if len(key[3]) == len(key[4]) == 1 else "indel"
    path = settings["spliceai_snv"] if resource_class == "snv" else settings["spliceai_indel"]
    rows, headers, failure = _prediction_fetch(path, key, f"spliceai_{resource_class}")
    if failure:
        return _prediction_failure(key, *failure)
    try:
        _header_columns(headers, _VCF_COLUMNS)
        records, loci, selected, missing = [], [], [], False
        for raw in rows:
            fields = raw.split("\t")
            if len(fields) < len(_VCF_COLUMNS):
                raise ValueError
            position = int(fields[1])
            chromosome = fields[0][3:] if fields[0].startswith("chr") else fields[0]
            if chromosome != key[1]:
                return _prediction_failure(key, "schema_mismatch", "spliceai_contig_mismatch")
            if position != key[2] or fields[3] != key[3] or key[4] not in fields[4].split(","):
                continue
            annotations = _info(fields[7]).get("SpliceAI")
            if annotations is None or annotations in _MISSING:
                loci.append("\t".join(fields[:7]))
                missing = True
                continue
            matched = False
            for annotation in annotations.split(","):
                values = annotation.split("|")
                if len(values) < 2 or values[:2] != [key[4], gene]:
                    continue
                matched = True
                selected.append(annotation)
                if len(values) != 10:
                    raise ValueError
                if any(value in _MISSING for value in values[2:]):
                    missing = True
                    continue
                records.append(
                    {
                        "variant_key": key,
                        "gene": gene,
                        "ds_ag": _unit_score(values[2]),
                        "ds_al": _unit_score(values[3]),
                        "ds_dg": _unit_score(values[4]),
                        "ds_dl": _unit_score(values[5]),
                        "dp_ag": _offset(values[6]),
                        "dp_al": _offset(values[7]),
                        "dp_dg": _offset(values[8]),
                        "dp_dl": _offset(values[9]),
                    }
                )
            if matched:
                loci.append("\t".join(fields[:7]))
    except (ValueError, OverflowError):
        return _prediction_failure(key, "schema_mismatch", "spliceai_record_invalid")
    if records or missing:
        source = _prediction_source(
            "spliceai",
            settings["spliceai_release"],
            "spliceai_precomputed_prediction",
            gene,
            key,
            headers,
            {"raw_loci": tuple(loci), "selected_annotations": tuple(selected)},
            {
                "gene_match_basis": "exact_symbol_and_alternate",
                "resource_class": resource_class,
                "prediction_missing": missing,
            },
        )
        return _prediction_lookup(
            key,
            "found" if records else "missing_prediction",
            records,
            source,
            "" if records else "spliceai_prediction_missing",
        )
    return _prediction_lookup(key, "not_found")


def _prediction_source(provider, release, record_type, gene, key, headers, evidence, extra):
    content = json.dumps(
        {
            "record_type": record_type,
            "assembly": key[0],
            "query_gene": gene,
            "variant_key": _key_text(key),
            "header": list(headers),
            **evidence,
            **extra,
        },
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    fragment = hashlib.sha256(content.encode()).hexdigest()[:16]
    return Source(
        provider=provider,
        kind="prediction",
        locator=f"resource:{provider}/prediction#{fragment}",
        version=release,
        content=content,
        genes=(gene,),
        variant_keys=(key,),
    )


__all__ = ["annotate", "population", "predictions"]
