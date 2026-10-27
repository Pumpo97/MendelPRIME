"""Private AutoPVS1 child protocol and source-bound parent adapter."""

import argparse
import configparser
import hashlib
import importlib
import importlib.util
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path

from mendelprime.domain import MendelPrimeError, Source, validate_gene, variant_key, variant_text
from mendelprime.runtime.process import batch_outcome, complete_batch, fail_batch, run_process

PROTOCOL = "mendelprime.pvs1.v1"
_CHILD_BOOTSTRAP = (
    "import os,runpy,sys; "
    "source=sys.argv.pop(1); "
    "os.environ['PYTHONPATH']=source; "
    "sys.path.insert(0,source); "
    "runpy.run_module('mendelprime.tools.pvs1_driver',run_name='__main__')"
)
_MAX_REQUEST = 16 * 1024 * 1024
_MAX_OUTPUT = 32 * 1024 * 1024
_MAX_TARGETS = 10_000
_STRENGTHS = {"VeryStrong", "Strong", "Moderate", "Supporting", "Unmet", "Unset"}
_POSITIVE = {"VeryStrong", "Strong", "Moderate", "Supporting"}
_CONTEXT_NOTE = "Native algorithm context, not an independent assay or additional criterion; a missing population record does not establish adequate coverage."
_REQUEST_ID = re.compile(r"PVS1_[0-9]+", re.ASCII)
_ERROR_CODE = re.compile(r"pvs1_[a-z0-9_]+", re.ASCII)
_PACKAGE_FILES = (
    "__init__.py",
    "pvs1.py",
    "read_data.py",
    "strength.py",
    "utils.py",
    "splicing.py",
    "autoPVS1.py",
    "cnv.py",
    "pyhgvs/__init__.py",
    "pyhgvs/models.py",
    "pyhgvs/utils.py",
    "pyhgvs/variants.py",
    "maxentpy/__init__.py",
    "maxentpy/maxent.py",
)
_MATRIX_FILES = ("maxentpy/data/score3_matrix.txt", "maxentpy/data/score5_matrix.txt")
_DATA_FIELDS = (
    ("DEFAULT", "pvs1levels"),
    ("DEFAULT", "gene_alias"),
    ("DEFAULT", "gene_trans"),
    *(
        (build, name)
        for build in ("HG19", "HG38")
        for name in (
            "transcript",
            "domain",
            "hotspot",
            "curated_region",
            "exon_lof_popmax",
            "pathogenic_site",
        )
    ),
)
_ANNOTATION_FIELDS = (
    "criterion",
    "strength_raw",
    "strength",
    "met",
    "decision_tree_supported",
    "gene_adjustment_status",
    "transcript",
    "consequence",
    "actual_gene",
    "gene_match",
    "gene_bindings",
    "method",
    "research_use_only",
)
_PROVENANCE_FIELDS = (
    "source_commit",
    "git_diff_binary_sha256",
    "read_data_sha256",
    "config_sha256",
    "code_tree_sha256",
    "algorithm_assets_sha256",
)


def preflight_source(root):
    """Validate every import-time file and the usable hg19 FASTA index."""
    try:
        if not root.is_dir() or not os.access(root, os.R_OK | os.X_OK):
            return "pvs1_source_root_unavailable"
        package_files = tuple(root / name for name in (*_PACKAGE_FILES, "config.ini"))
        if any(
            not path.is_file()
            or not os.access(path, os.R_OK)
            or (path.name != "__init__.py" and path.stat().st_size == 0)
            for path in package_files
        ):
            return "pvs1_source_files_unavailable"
        config = configparser.ConfigParser()
        with (root / "config.ini").open(encoding="utf-8") as handle:
            config.read_file(handle)
        resources = _required_resource_paths(root, config) + tuple(root / name for name in _MATRIX_FILES)
        if any(
            not path.is_file() or path.stat().st_size == 0 or not os.access(path, os.R_OK)
            for path in resources
        ):
            return "pvs1_resource_unavailable"
        for section in ("HG19", "HG38"):
            fasta = _configured_path(root, config[section]["genome"])
            if section == "HG38" and not fasta.exists():
                continue
            index = Path(f"{fasta}.fai")
            if not fasta.is_file() or fasta.stat().st_size == 0 or not os.access(fasta, os.R_OK):
                return "pvs1_fasta_unavailable"
            if not index.is_file() or index.stat().st_size == 0 or not os.access(index, os.R_OK):
                return "pvs1_fasta_index_unavailable"
            if index.stat().st_mtime_ns < fasta.stat().st_mtime_ns:
                return "pvs1_fasta_index_stale"
    except (KeyError, configparser.Error, OSError, UnicodeError):
        return "pvs1_resource_preflight_failed"
    return ""


def run_pvs1(settings, work_root, variants):
    """Assess eligible native variant dictionaries in an isolated child process."""
    root, timeout = settings["autopvs1_root"], settings["pvs1_timeout_seconds"]
    for field, path in (("root", root), ("workroot", work_root)):
        if not isinstance(path, Path) or not path.is_absolute():
            raise MendelPrimeError("absolute_path_required", field=f"pvs1.{field}")
    if type(timeout) not in {int, float} or not math.isfinite(timeout) or timeout <= 0:
        raise MendelPrimeError("positive_finite_timeout_required", field="pvs1.timeout")
    if work_root.resolve().is_relative_to(root.resolve()):
        raise MendelPrimeError("workroot_overlaps_source", field="pvs1.workroot")
    if type(variants) is not tuple:
        raise MendelPrimeError("native_variant_tuple_required", field="pvs1.variants")
    if not variants:
        return ()
    outcomes, eligible = [None] * len(variants), []
    for index, target in enumerate(variants):
        target = _native_target(target)
        issue = _eligibility_error(target, root=root)
        if issue:
            outcomes[index] = batch_outcome(target, *issue)
            continue
        identifier = f"PVS1_{len(eligible)}"
        eligible.append((index, target, identifier, _request_target(identifier, target)))
    if not eligible:
        return complete_batch(variants, outcomes, "pvs1_internal_outcome_missing")
    error = _resource_error(root, work_root)
    if error:
        fail_batch(outcomes, eligible, error)
        return complete_batch(variants, outcomes, "pvs1_internal_outcome_missing")
    try:
        with tempfile.TemporaryDirectory(prefix="mendelprime-pvs1-", dir=work_root) as temporary:
            work = Path(temporary)
            work.chmod(0o700)
            request, response = work / "request.json", work / "response.json"
            _write_json(request, {"protocol": PROTOCOL, "targets": [item[3] for item in eligible]})
            receipt = run_process(
                command=[
                    sys.executable,
                    "-I",
                    "-B",
                    "-c",
                    _CHILD_BOOTSTRAP,
                    str(Path(__file__).resolve().parents[2]),
                    "--source-root",
                    str(root),
                    "--input",
                    str(request),
                    "--output",
                    str(response),
                ],
                root=work,
                timeout_seconds=timeout,
                artifacts=(response,),
            )
            if receipt["error_code"]:
                fail_batch(outcomes, eligible, receipt["error_code"])
            else:
                try:
                    rows = _read_response(response, frozenset(item[2] for item in eligible))
                except MendelPrimeError as exc:
                    fail_batch(outcomes, eligible, exc.code)
                else:
                    for index, target, identifier, request_row in eligible:
                        row = rows.get(identifier)
                        outcomes[index] = (
                            batch_outcome(target, "failed", "pvs1_result_missing")
                            if row is None
                            else _row_outcome(target, request_row, row)
                        )
    except OSError:
        fail_batch(outcomes, eligible, "pvs1_private_workdir_failed")
    return complete_batch(variants, outcomes, "pvs1_internal_outcome_missing")


def _native_target(target):
    required = {"gene", "key", "transcript", "hgvs_c", "hgvs_p", "details", "source"}
    if not isinstance(target, Mapping) or not required.issubset(target):
        raise MendelPrimeError("native_variant_required", field="pvs1.variants")
    target = dict(target)
    validate_gene(target["gene"])
    if target["key"] is not None:
        target["key"] = variant_key(*target["key"])
    if not isinstance(target["source"], Source):
        raise MendelPrimeError("source_required", field="pvs1.variant.source")
    return target


def resource_gene_aliases(root):
    config = configparser.ConfigParser()
    config.read(root / "config.ini")
    data = _configured_path(root, config["DEFAULT"]["gene_alias"]).read_bytes()
    aliases = {}
    for line in data.decode().splitlines():
        fields = line.split("\t")
        if len(fields) == 2:
            canonical, alias = fields
            if alias in aliases and aliases[alias] != canonical:
                raise ValueError("ambiguous resource gene alias")
            aliases[alias] = canonical
    return aliases, hashlib.sha256(data).hexdigest()


def _eligibility_error(target, *, root=None):
    if (
        target["source"].provider not in {"exomiser", "ensembl-vep"}
        and not (target["source"].provider == "lirical" and target.get("verify_native_snv") is True)
    ) or not isinstance(target["details"], Mapping):
        return "not_applicable", "pvs1_native_source_not_applicable"
    if target["source"].provider == "ensembl-vep":
        try:
            record = json.loads(target["source"].content)
            row = record["selected_transcript"]
            selection = record.get("transcript_selection_source_ids", ())
            if selection:
                from mendelprime.knowledge.annotation import mane_transcripts

                candidates = mane_transcripts(target.get("mane_candidates", {}).values(), target["gene"])
                candidate = candidates.get(target["transcript"])
                if candidate is None or list(selection) != [candidate.source_id]:
                    return "unknown", "pvs1_verified_annotation_required"
            gene_matches = row.get("gene_symbol") == target["gene"]
            if not gene_matches and root is not None:
                aliases, digest = resource_gene_aliases(root)
                gene_matches = aliases.get(row.get("gene_symbol")) == target["gene"] and record.get(
                    "gene_alias_binding"
                ) == {
                    "reported_symbol": row.get("gene_symbol"),
                    "canonical_symbol": target["gene"],
                    "resource_sha256": digest,
                }
            verified = (
                record.get("record_type") == "verified_refseq_allele_reannotation_v1"
                and record.get("gene") == target["gene"]
                and record.get("variant_key") == variant_text(target["key"])
                and record.get("assembly") == target["key"][0]
                and target["source"].genes == (target["gene"],)
                and target["source"].variant_keys == (target["key"],)
                and record.get("assessed_transcript") == row.get("transcript_id") == target["transcript"]
                and gene_matches
                and row.get("hgvsc") == target["transcript"] + ":" + target["hgvs_c"]
                and (row.get("hgvsp") or "").split(":", 1)[-1] == target["hgvs_p"]
                and "&".join(row.get("consequence_terms", ())) == target["details"]["consequence"]
            )
        except (MendelPrimeError, ValueError, KeyError, TypeError, OSError, configparser.Error):
            verified = False
        if not verified:
            return "unknown", "pvs1_verified_annotation_required"
    consequence = _consequence(target["details"].get("consequence") or "")
    if not consequence:
        return "not_applicable", "pvs1_consequence_not_applicable"
    if target["key"] is None:
        return "unknown", "pvs1_variant_identity_unavailable"
    if target["key"][0] != "GRCh37":
        return "unknown", "pvs1_unsupported_assembly"
    if not target["transcript"]:
        return "unknown", "pvs1_transcript_unavailable"
    if not target["hgvs_c"]:
        return "unknown", "pvs1_hgvsc_unavailable"
    if consequence in {"frameshift", "nonsense"} and not target["hgvs_p"]:
        return "unknown", "pvs1_hgvsp_unavailable"
    return None


def _request_target(identifier, target):
    key = target["key"]
    request = {
        "request_id": identifier,
        "allele": dict(
            zip(("assembly", "chromosome", "position", "reference", "alternate"), key, strict=True)
        ),
        "gene": target["gene"],
        "transcript": target["transcript"],
        "hgvs_c": target["hgvs_c"],
        "hgvs_p": target["hgvs_p"],
        "consequence": target["details"]["consequence"],
    }
    if target.get("verify_native_snv") is True:
        request["verify_native_snv"] = True
    if target.get("retain_context") is True:
        request["retain_context"] = True
    return request


def _resource_error(root, work_root):
    try:
        if not work_root.is_dir() or not os.access(work_root, os.W_OK | os.X_OK):
            return "pvs1_work_root_unavailable"
    except OSError:
        return "pvs1_work_root_unavailable"
    return preflight_source(root)


def _write_json(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=True, sort_keys=True, separators=(",", ":"))
    path.chmod(0o600)


def _read_response(path, expected):
    try:
        if not path.is_file() or path.stat().st_size == 0:
            raise MendelPrimeError("pvs1_output_missing")
        if path.stat().st_size > _MAX_OUTPUT:
            raise MendelPrimeError("pvs1_output_too_large")
        payload = json.loads(path.read_text(encoding="utf-8"))
    except MendelPrimeError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError):
        raise MendelPrimeError("pvs1_output_parse_failed") from None
    if not isinstance(payload, dict) or payload.get("protocol") != PROTOCOL:
        raise MendelPrimeError("pvs1_output_protocol_mismatch")
    results = payload.get("results")
    if not isinstance(results, list) or len(results) > len(expected):
        raise MendelPrimeError("pvs1_output_bounds_exceeded")
    rows = {}
    for row in results:
        if not isinstance(row, dict) or not isinstance(row.get("request_id"), str):
            raise MendelPrimeError("pvs1_output_schema_mismatch")
        identifier = row["request_id"]
        if identifier not in expected:
            raise MendelPrimeError("pvs1_unknown_request_id")
        if identifier in rows:
            raise MendelPrimeError("pvs1_duplicate_request_id")
        rows[identifier] = row
    return rows


def _row_outcome(target, request, row):
    if not _identity_matches(row, request):
        return batch_outcome(target, "failed", "pvs1_result_identity_mismatch")
    if (
        request.get("verify_native_snv")
        and row.get("status") == "assessed"
        and not row.get("native_hgvs_verification")
    ):
        return batch_outcome(target, "failed", "pvs1_native_verification_missing")
    status = row.get("status")
    if status in {"unknown", "failed"}:
        code, available = row.get("error_code"), row.get("annotation_available")
        if not isinstance(code, str) or not _ERROR_CODE.fullmatch(code) or type(available) is not bool:
            return batch_outcome(target, "failed", "pvs1_result_schema_mismatch")
        if not available:
            return batch_outcome(target, status, code)
        if status != "unknown" or code != "pvs1_strength_unset":
            return batch_outcome(target, "failed", "pvs1_result_schema_mismatch")
        try:
            return batch_outcome(target, "unknown", code, _annotation(target, row))
        except (MendelPrimeError, KeyError, TypeError, ValueError):
            return batch_outcome(target, "failed", "pvs1_result_schema_mismatch")
    if status != "assessed" or row.get("error_code") not in (None, ""):
        return batch_outcome(target, "failed", "pvs1_result_schema_mismatch")
    try:
        return batch_outcome(target, "assessed", "", _annotation(target, row))
    except (MendelPrimeError, KeyError, TypeError, ValueError):
        return batch_outcome(target, "failed", "pvs1_result_schema_mismatch")


def _identity_matches(row, request):
    return (
        row.get("allele") == request["allele"]
        and row.get("requested_gene") == request["gene"]
        and row.get("requested_transcript") == request["transcript"]
        and row.get("hgvs_c") == request["hgvs_c"]
        and row.get("hgvs_p") == request["hgvs_p"]
        and row.get("requested_consequence") == request["consequence"]
    )


def _annotation(target, row):
    provenance = row.get("provenance")
    if not isinstance(provenance, dict):
        raise MendelPrimeError("mapping_required", field="pvs1.provenance")
    fields = {name: row.get(name) for name in _ANNOTATION_FIELDS}
    if row.get("native_hgvs_verification"):
        fields["native_hgvs_verification"] = row["native_hgvs_verification"]
    if row.get("mechanism_context"):
        context = row["mechanism_context"]
        required = {"termination_aa", "cds_length_nt", "interpretation"}
        flags = {"is_nmd_target", "LoF_removes_more_than_10_percent_of_protein"}
        text = {"func_desc", "exon_lof_popmax_desc"}
        if (
            not isinstance(context, dict)
            or not required <= set(context)
            or set(context) - required - flags - text
            or context["interpretation"] != _CONTEXT_NOTE
            or type(context["termination_aa"]) is not int
            or context["termination_aa"] < 1
            or (
                context["cds_length_nt"] is not None
                and (type(context["cds_length_nt"]) is not int or context["cds_length_nt"] < 1)
            )
            or any(type(context[k]) is not bool for k in flags & set(context))
            or any(not isinstance(context[k], str) for k in text & set(context))
        ):
            raise MendelPrimeError("pvs1_context_schema_invalid", field="pvs1.mechanism_context")
        fields["mechanism_context"] = context
    bindings = fields["gene_bindings"]
    if not isinstance(bindings, list) or any(
        not isinstance(value, list) or len(value) != 2 for value in bindings
    ):
        raise MendelPrimeError("symbol_canonical_pairs_required", field="pvs1.gene_bindings")
    fields["gene_bindings"] = tuple((value[0], value[1]) for value in bindings)
    retained_provenance = {name: provenance.get(name) for name in _PROVENANCE_FIELDS}
    if row.get("transcript") != target["transcript"]:
        raise MendelPrimeError("requested_transcript_mismatch", field="pvs1.transcript")
    retained = {
        "allele": row["allele"],
        "gene": target["gene"],
        "hgvs_c": row.get("hgvs_c"),
        "hgvs_p": row.get("hgvs_p"),
        "native_consequence": row.get("requested_consequence"),
        "tool_result": fields,
        "provenance": retained_provenance,
    }
    content = json.dumps(retained, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    key = target["key"]
    version = "AutoPVS1;" + ";".join(
        f"{name}={value}" for name, value in retained_provenance.items() if value is not None
    )
    source = Source(
        provider="autopvs1",
        kind="prediction",
        locator=f"method:autopvs1/annotation#{variant_text(key)}:{target['gene']}:{target['transcript']}",
        version=version,
        content=content,
        genes=(target["gene"],),
        variant_keys=(key,),
        case_key=target["source"].case_key,
    )
    return {**fields, **retained_provenance, "source": source}


def _configured_path(root, value):
    path = Path(value)
    return path if path.is_absolute() else root / path


def _required_resource_paths(root, config):
    try:
        return tuple(_configured_path(root, config[section][name]) for section, name in _DATA_FIELDS)
    except (KeyError, configparser.Error):
        raise ValueError("pvs1_config_invalid") from None


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_commit(root):
    try:
        head = (root / ".git" / "HEAD").read_text(encoding="utf-8").strip()
        if not head.startswith("ref: "):
            return head if re.fullmatch(r"[0-9a-f]{40,64}", head) else "unknown"
        reference = head[5:]
        loose = root / ".git" / reference
        if loose.is_file():
            value = loose.read_text(encoding="utf-8").strip()
            return value if re.fullmatch(r"[0-9a-f]{40,64}", value) else "unknown"
        for line in (root / ".git" / "packed-refs").read_text(encoding="utf-8").splitlines():
            if line.endswith(f" {reference}"):
                value = line.split(" ", 1)[0]
                return value if re.fullmatch(r"[0-9a-f]{40,64}", value) else "unknown"
    except (OSError, UnicodeError):
        pass
    return "unknown"


def _provenance(root):
    commit, diff_digest = _git_commit(root), None
    if commit != "unknown":
        try:
            diff = subprocess.run(
                [
                    "git",
                    "--no-optional-locks",
                    "-C",
                    str(root),
                    "diff",
                    "--binary",
                    "--no-ext-diff",
                    "--no-textconv",
                    "HEAD",
                    "--",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                check=False,
                timeout=30,
            )
            if diff.returncode == 0:
                diff_digest = hashlib.sha256(diff.stdout).hexdigest()
        except (OSError, subprocess.SubprocessError):
            pass
    try:
        config = configparser.ConfigParser()
        with (root / "config.ini").open(encoding="utf-8") as handle:
            config.read_file(handle)
        assets = {
            f"{section}.{name}": _sha256(_configured_path(root, config[section][name]))
            for section, name in _DATA_FIELDS
        }
        assets.update({name: _sha256(root / name) for name in _MATRIX_FILES})
        return {
            "source_commit": None if commit == "unknown" else commit,
            "git_diff_binary_sha256": diff_digest,
            "read_data_sha256": _sha256(root / "read_data.py"),
            "config_sha256": _sha256(root / "config.ini"),
            "code_tree_sha256": hashlib.sha256(
                json.dumps(
                    {path.relative_to(root).as_posix(): _sha256(path) for path in sorted(root.rglob("*.py"))},
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode()
            ).hexdigest(),
            "algorithm_assets_sha256": hashlib.sha256(
                json.dumps(assets, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest(),
        }
    except (OSError, KeyError, configparser.Error, UnicodeError):
        return None


def _load_request(path):
    try:
        if not path.is_file() or path.stat().st_size == 0 or path.stat().st_size > _MAX_REQUEST:
            raise ValueError
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError):
        raise ValueError from None
    if not isinstance(payload, dict) or payload.get("protocol") != PROTOCOL:
        raise ValueError
    targets = payload.get("targets")
    if not isinstance(targets, list) or len(targets) > _MAX_TARGETS:
        raise ValueError
    seen = set()
    for target in targets:
        if not _valid_target(target) or target["request_id"] in seen:
            raise ValueError
        seen.add(target["request_id"])
    return targets


def _valid_target(target):
    if not isinstance(target, dict) or not _REQUEST_ID.fullmatch(str(target.get("request_id", ""))):
        return False
    if "verify_native_snv" in target and type(target["verify_native_snv"]) is not bool:
        return False
    if "retain_context" in target and type(target["retain_context"]) is not bool:
        return False
    allele = target.get("allele")
    if not isinstance(allele, dict) or set(allele) != {
        "assembly",
        "chromosome",
        "position",
        "reference",
        "alternate",
    }:
        return False
    if (
        allele.get("assembly") != "GRCh37"
        or type(allele.get("position")) is not int
        or allele["position"] < 1
    ):
        return False
    if any(
        not isinstance(allele.get(field), str) or not allele[field]
        for field in ("chromosome", "reference", "alternate")
    ):
        return False
    if any(
        not isinstance(target.get(field), str) or not target[field]
        for field in ("gene", "transcript", "hgvs_c", "consequence")
    ):
        return False
    return target.get("hgvs_p") is None or isinstance(target.get("hgvs_p"), str)


def _identity(target):
    return {
        "request_id": target["request_id"],
        "allele": target["allele"],
        "requested_gene": target["gene"],
        "requested_transcript": target["transcript"],
        "hgvs_c": target["hgvs_c"],
        "hgvs_p": target["hgvs_p"],
        "requested_consequence": target["consequence"],
    }


def _error(target, status, code, annotation=None):
    return {
        **_identity(target),
        "status": status,
        "error_code": code,
        "annotation_available": annotation is not None,
        **({} if annotation is None else annotation),
    }


def _consequence(value):
    normalized = value.strip().lower().replace("_variant", "")
    for marker, result in (
        ("frameshift", "frameshift"),
        ("stop_gained", "nonsense"),
        ("nonsense", "nonsense"),
        ("splice_donor", "splice-5"),
        ("splice_acceptor", "splice-3"),
        ("start_lost", "init-loss"),
    ):
        if marker in normalized:
            return result
    return ""


def _verify_native_nonsense(target, transcript):
    """Validate exact-version native SNV HGVS against genomic and codon sequence."""
    from bioutils.sequences import aa3_to_aa1, reverse_complement, translate_cds

    mapper = importlib.import_module("autopvs1.pyhgvs")
    data = importlib.import_module("autopvs1.read_data")
    genome = data.genome_hg19
    actual_gene = transcript.gene.name
    if transcript.full_name != target["transcript"] or data.gene_alias.get(
        actual_gene, actual_gene
    ) != data.gene_alias.get(target["gene"], target["gene"]):
        raise ValueError("exact transcript and gene required")
    coding = re.fullmatch(r"c\.(\d+)([ACGT])>([ACGT])", target["hgvs_c"])
    protein = re.fullmatch(r"p\.\(?([A-Z](?:[a-z]{2})?)(\d+)(?:Ter|\*)\)?", target["hgvs_p"] or "")
    if coding is None or protein is None:
        raise ValueError("simple nonsense SNV required")
    allele = target["allele"]
    chrom, pos, ref, alt = mapper.parse_hgvs_name(
        target["transcript"] + ":" + target["hgvs_c"], genome, transcript=transcript
    )
    if (chrom.removeprefix("chr"), pos, ref, alt) != (
        allele["chromosome"].removeprefix("chr"),
        allele["position"],
        allele["reference"],
        allele["alternate"],
    ):
        raise ValueError("HGVS genomic mismatch")
    codon_start = ((int(coding[1]) - 1) // 3) * 3 + 1
    if (codon_start + 2) // 3 != int(protein[2]):
        raise ValueError("protein coordinate mismatch")
    bases = []
    for coordinate in range(codon_start, codon_start + 3):
        position = mapper.cdna_to_genomic_coord(transcript, mapper.CDNACoord(coordinate))
        base = mapper.get_genomic_sequence(genome, transcript.tx_position.chrom, position, position)
        bases.append(base if transcript.tx_position.is_forward_strand else reverse_complement(base))
    codon = "".join(bases)
    offset = int(coding[1]) - codon_start
    changed = codon[:offset] + coding[3] + codon[offset + 1 :]
    amino = aa3_to_aa1(protein[1]) if len(protein[1]) == 3 else protein[1]
    if codon[offset] != coding[2] or translate_cds(codon) != amino or translate_cds(changed) != "*":
        raise ValueError("protein consequence mismatch")
    return {
        "method": "exact-version pyhgvs genomic roundtrip and local reference codon translation",
        "transcript": target["transcript"],
        "hgvs_c": target["hgvs_c"],
        "hgvs_p": target["hgvs_p"],
        "reference_codon": codon,
        "alternate_codon": changed,
        "allele": allele,
    }


def _evaluate(target, *, pvs1_type, strength_type, vcf_record_type, transcripts, aliases, provenance):
    try:
        transcript = transcripts.get(target["transcript"])
        actual_transcript = str(getattr(transcript, "full_name", "") or getattr(transcript, "name", ""))
        if transcript is None or actual_transcript != target["transcript"]:
            return _error(target, "unknown", "pvs1_exact_transcript_unavailable")
        actual_gene = str(getattr(getattr(transcript, "gene", None), "name", ""))
        if not actual_gene:
            return _error(target, "unknown", "pvs1_transcript_gene_unavailable")
        requested_current, actual_current = (
            aliases.get(target["gene"], target["gene"]),
            aliases.get(actual_gene, actual_gene),
        )
        if not all(
            isinstance(value, str) and value and value == value.strip()
            for value in (requested_current, actual_current)
        ):
            return _error(target, "failed", "pvs1_gene_alias_invalid")
        if actual_gene == target["gene"]:
            gene_match, bindings = "exact", ((target["gene"], requested_current),)
        elif requested_current == actual_current and (target["gene"] in aliases or actual_gene in aliases):
            gene_match, bindings = (
                "resource_alias",
                ((target["gene"], requested_current), (actual_gene, actual_current)),
            )
        else:
            return _error(target, "unknown", "pvs1_transcript_gene_mismatch")
        verification = None
        if target.get("verify_native_snv"):
            try:
                verification = _verify_native_nonsense(target, transcript)
            except Exception:
                return _error(target, "unknown", "pvs1_native_hgvs_verification_failed")
        consequence = _consequence(target["consequence"])
        if not consequence:
            return _error(target, "unknown", "pvs1_consequence_unavailable")
        allele = target["allele"]
        result = pvs1_type(
            vcf_record_type(
                allele["chromosome"].removeprefix("chr"),
                allele["position"],
                allele["reference"],
                allele["alternate"],
            ),
            consequence,
            target["hgvs_c"],
            target["hgvs_p"] or "",
            transcript,
            "GRCh37",
        )
        if consequence in {"frameshift", "nonsense"}:
            termination = getattr(result, "get_pHGVS_termination", None)
            if type(termination) is not int or termination < 1:
                return _error(target, "unknown", "pvs1_protein_termination_unavailable")
        strength_raw = _strength(result.strength_raw, strength_type)
        strength = _strength(result.strength, strength_type)
        if strength_raw not in _STRENGTHS or strength not in _STRENGTHS:
            return _error(target, "failed", "pvs1_strength_invalid")
        criterion = getattr(result, "criterion", None)
        if not isinstance(criterion, str) or not criterion:
            return _error(target, "failed", "pvs1_criterion_invalid")
        tree = None if strength_raw == "Unset" else strength_raw in _POSITIVE
        if strength == "Unset":
            met, adjustment = None, "unavailable"
        elif strength == "Unmet":
            met, adjustment = False, "evaluated"
        else:
            met, adjustment = True, "evaluated"
        annotation = {
            "actual_gene": actual_gene,
            "gene_match": gene_match,
            "gene_bindings": bindings,
            "transcript": target["transcript"],
            "consequence": consequence,
            "criterion": criterion,
            "strength_raw": strength_raw,
            "strength": strength,
            "met": met,
            "decision_tree_supported": tree,
            "gene_adjustment_status": adjustment,
            "method": "AutoPVS1 decision tree from supplied consequence and HGVS",
            "research_use_only": True,
            "provenance": provenance,
        }
        if verification:
            annotation["native_hgvs_verification"] = verification
        if target.get("retain_context") and consequence in {"frameshift", "nonsense"}:
            context = {
                "termination_aa": termination,
                "cds_length_nt": getattr(transcript, "cds_length", None),
            }
            for name in (
                "is_nmd_target",
                "LoF_removes_more_than_10_percent_of_protein",
                "func_desc",
                "exon_lof_popmax_desc",
            ):
                try:
                    value = getattr(result, name, None)
                except Exception:
                    value = None
                if type(value) in {bool, str}:
                    context[name] = value
            context["interpretation"] = _CONTEXT_NOTE
            annotation["mechanism_context"] = context
        if strength_raw == "Unset" or strength == "Unset":
            return _error(target, "unknown", "pvs1_strength_unset", annotation)
        return {**_identity(target), "status": "assessed", "error_code": "", **annotation}
    except Exception:
        return _error(target, "failed", "pvs1_evaluation_failed")


def _strength(value, strength_type):
    for name in _STRENGTHS:
        if value == getattr(strength_type, name):
            return name
    return str(value)


def _load_external(root):
    # Bind the configured directory, independent of its name or installed packages.
    for name in tuple(sys.modules):
        if name == "autopvs1" or name.startswith("autopvs1."):
            del sys.modules[name]
    spec = importlib.util.spec_from_file_location(
        "autopvs1", root / "__init__.py", submodule_search_locations=[str(root)]
    )
    package = importlib.util.module_from_spec(spec)
    sys.modules["autopvs1"] = package
    spec.loader.exec_module(package)
    if Path(str(getattr(package, "__file__", ""))).resolve().parent != root.resolve():
        raise ImportError
    pvs1_module = importlib.import_module("autopvs1.pvs1")
    read_data = importlib.import_module("autopvs1.read_data")
    strength_module = importlib.import_module("autopvs1.strength")
    utils = importlib.import_module("autopvs1.utils")
    return (
        pvs1_module.PVS1,
        strength_module.Strength,
        utils.VCFRecord,
        read_data.transcripts_hg19,
        read_data.gene_alias,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--check-source", type=Path)
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    try:
        arguments = parser.parse_args(argv)
        if arguments.check_source:
            root = arguments.check_source.resolve()
            if preflight_source(root):
                return 1
            try:
                _load_external(root)
            except Exception:
                return 1
            return 0
        if not all((arguments.source_root, arguments.input, arguments.output)):
            return 2
        targets = _load_request(arguments.input)
    except (SystemExit, ValueError):
        return 2
    root = arguments.source_root.resolve()
    code = preflight_source(root)
    provenance = None if code else _provenance(root)
    if not code and provenance is None:
        code = "pvs1_provenance_unavailable"
    if code:
        try:
            _write_json(
                arguments.output,
                {"protocol": PROTOCOL, "results": [_error(target, "failed", code) for target in targets]},
            )
        except OSError:
            return 3
        return 0
    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        try:
            components = _load_external(root)
        except Exception:
            results = [_error(target, "failed", "pvs1_dependency_import_failed") for target in targets]
        else:
            results = [
                _evaluate(
                    target,
                    pvs1_type=components[0],
                    strength_type=components[1],
                    vcf_record_type=components[2],
                    transcripts=components[3],
                    aliases=components[4],
                    provenance=provenance,
                )
                for target in targets
            ]
    finally:
        sys.dont_write_bytecode = previous
    if _provenance(root) != provenance:
        results = [_error(target, "failed", "pvs1_source_changed") for target in targets]
    try:
        _write_json(arguments.output, {"protocol": PROTOCOL, "results": results})
    except OSError:
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
