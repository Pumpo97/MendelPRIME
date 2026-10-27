"""Read-only ClinVar and ClinGen assertions, before independent ACMG assessment."""

import csv
import hashlib
import json
import re
import sqlite3
import threading
from contextlib import closing
from datetime import date, datetime
from pathlib import Path
from urllib.parse import quote

from mendelprime.domain import MendelPrimeError, Source, require, validate_gene, variant_key, variant_text
from mendelprime.runtime.storage import stat_fingerprint

from .literature import evidence_document

CLINVAR_SCHEMA = "mgagent_clinvar_index_v1"
CLINGEN_SCHEMA = "mgagent_clingen_resource_v1"
GENE_HEADERS = (
    "GENE SYMBOL",
    "GENE ID (HGNC)",
    "DISEASE LABEL",
    "DISEASE ID (MONDO)",
    "MOI",
    "SOP",
    "CLASSIFICATION",
    "ONLINE REPORT",
    "CLASSIFICATION DATE",
    "GCEP",
)
VCEP_HEADERS = (
    "Variation",
    "ClinVar Variation Id",
    "Allele Registry Id",
    "HGVS Expressions",
    "HGNC Gene Symbol",
    "Disease",
    "Mondo Id",
    "Mode of Inheritance",
    "Assertion",
    "Applied Evidence Codes (Met)",
    "Applied Evidence Codes (Not Met)",
    "Summary of interpretation",
    "PubMed Articles",
    "Expert Panel",
    "Guideline",
    "Approval Date",
    "Published Date",
    "Retracted",
    "Evidence Repo Link",
    "Uuid",
)
TABLE_COLUMNS = {
    "metadata": ("key", "value"),
    "variants": (
        "variation_id",
        "allele_id",
        "assembly",
        "chromosome",
        "position",
        "reference",
        "alternate",
        "gene_symbol",
        "rs_id",
        "clinical_significance",
        "review_status",
        "name",
    ),
    "citations": ("variation_id", "allele_id", "pmid"),
    "identifiers": ("identifier", "variation_id", "allele_id", "assembly"),
}
PRIMARY_KEYS = {**TABLE_COLUMNS, "metadata": ("key",), "variants": TABLE_COLUMNS["variants"][:7]}
INDEX_PREFIXES = {"variants": TABLE_COLUMNS["variants"][2:7], "citations": TABLE_COLUMNS["citations"][:2]}
CSV_SPECS = {
    "gene_validity": (
        "clingen_validity",
        "gene_disease_validity.csv",
        GENE_HEADERS,
        ",",
        "GENE SYMBOL",
        True,
    ),
    "vcep": ("clingen_vcep", "variant_pathogenicity.tsv", VCEP_HEADERS, "\t", "ClinVar Variation Id", False),
}
SUMMARY_PMID = re.compile(
    r"(?<![A-Za-z0-9])PMID(?P<separator>\s*:\s*|\s+)(?P<identifier>[^\s,;\)\]\}]+)", re.I
)
FIELD_PMID = re.compile(r"^PMID(?P<separator>\s*:\s*|\s+)(?P<identifier>\S+)$", re.I)
SHA256 = re.compile(r"[0-9a-f]{64}")


def new_variant_store(resources):
    paths = {
        name: resources[name]
        for name in ("clinvar_db", "clingen_validity", "clingen_vcep", "clingen_manifest")
    }
    for name, value in paths.items():
        require(isinstance(value, Path) and value.is_absolute(), name, "absolute_path_required")
    return {**paths, "manifest": None, "rows": {}, "clinvar": None, "lock": threading.RLock()}


def _json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _fingerprint(path):
    try:
        value = path.stat()
        return stat_fingerprint(value) if path.is_file() else None
    except OSError:
        return None


def _file_hash(path):
    try:
        with path.open("rb") as stream:
            return hashlib.file_digest(stream, "sha256").hexdigest()
    except OSError:
        return None


def _iso_date(value):
    try:
        return value if isinstance(value, str) and date.fromisoformat(value).isoformat() == value else None
    except ValueError:
        return None


def _manifest(state):
    with state["lock"]:
        path, cached = state["clingen_manifest"], state["manifest"]
        current = _fingerprint(path)
        if cached is not None:
            return (
                (cached, None)
                if current == cached["fingerprint"]
                else (None, ("schema_mismatch", "clingen_manifest_changed"))
            )
        if current is None:
            return None, ("resource_unavailable", "clingen_manifest_unavailable")
        try:
            with path.open("rb") as stream:
                raw = stream.read(1_000_001)
            if len(raw) > 1_000_000:
                raise ValueError
            value = json.loads(raw.decode("utf-8"))
        except OSError:
            return None, ("resource_unavailable", "clingen_manifest_unavailable")
        except (UnicodeError, ValueError):
            return None, ("schema_mismatch", "clingen_manifest_invalid")
        if not isinstance(value, dict) or value.get("schema_version") != CLINGEN_SCHEMA:
            return None, ("schema_mismatch", "clingen_manifest_schema_mismatch")
        release = _iso_date(value.get("release_date"))
        if release is None:
            return None, ("schema_mismatch", "clingen_manifest_release_invalid")
        if not isinstance(value.get("files"), dict):
            return None, ("schema_mismatch", "clingen_manifest_fields_missing")
        hashes = {}
        for _, name, *_ in CSV_SPECS.values():
            item = value["files"].get(name)
            digest = item.get("sha256") if isinstance(item, dict) else None
            if not isinstance(digest, str) or not SHA256.fullmatch(digest):
                return None, ("schema_mismatch", "clingen_manifest_hash_missing")
            hashes[name] = digest
        if _fingerprint(path) != current:
            return None, ("schema_mismatch", "clingen_manifest_changed")
        state["manifest"] = {"release": release, "hashes": hashes, "fingerprint": current}
        return state["manifest"], None


def _read_rows(path, headers, delimiter, key_field, locate):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter=delimiter)
        header = (
            next((row for row in reader if row and row[0] == key_field), None)
            if locate
            else next(reader, None)
        )
        require(
            header is not None and len(header) == len(set(header)) and set(headers) <= set(header),
            "clingen",
            "header_mismatch",
        )
        indexed = {}
        for values in reader:
            if not values or not any(values) or (locate and values[0].startswith("+++")):
                continue
            require(len(values) == len(header), "clingen", "row_mismatch")
            row = dict(zip(header, values, strict=True))
            if key := row[key_field].strip():
                indexed.setdefault(key, []).append(row)
        return {key: tuple(rows) for key, rows in indexed.items()}


def _rows(state, kind):
    field, name, *formatting = CSV_SPECS[kind]
    path, code = state[field], "clingen_gene_validity" if kind == "gene_validity" else "clingen_vcep"
    with state["lock"]:
        manifest, failure = _manifest(state)
        if failure:
            return None, failure
        current, cached = _fingerprint(path), state["rows"].get(kind)
        if cached is not None:
            return (
                (cached, None)
                if current == cached["fingerprint"]
                else (None, ("schema_mismatch", code + "_changed"))
            )
        if current is None:
            return None, ("resource_unavailable", code + "_unavailable")
        digest = _file_hash(path)
        if digest is None:
            return None, ("query_failed", code + "_read_failed")
        if digest != manifest["hashes"][name]:
            return None, ("schema_mismatch", code + "_hash_mismatch")
        try:
            records = _read_rows(path, *formatting)
        except MendelPrimeError as error:
            return None, ("schema_mismatch", f"{code}_{error.code}")
        except (OSError, UnicodeError, csv.Error):
            return None, ("query_failed", code + "_read_failed")
        if _fingerprint(path) != current:
            return None, ("schema_mismatch", code + "_changed")
        state["rows"][kind] = {
            "by_key": records,
            "version": f"clingen:{manifest['release']};sha256={digest}",
            "fingerprint": current,
        }
        return state["rows"][kind], None


def _connect(state):
    path = quote(state["clinvar_db"].as_posix(), safe="/:")
    return sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)


def _quoted(value):
    return '"' + value.replace('"', '""') + '"'


def _indexed(connection, table, prefix):
    for row in connection.execute(f"PRAGMA index_list({_quoted(table)})"):
        if len(row) > 4 and row[4]:
            continue
        columns = tuple(item[2] for item in connection.execute(f"PRAGMA index_info({_quoted(str(row[1]))})"))
        if columns[: len(prefix)] == prefix:
            return True
    return False


def _clinvar_info(state):
    current, cached = _fingerprint(state["clinvar_db"]), state["clinvar"]
    if cached is not None:
        return (
            (cached, None)
            if current == cached["fingerprint"]
            else (None, ("schema_mismatch", "clinvar_resource_changed"))
        )
    if current is None:
        return None, ("resource_unavailable", "clinvar_snapshot_unavailable")
    try:
        with closing(_connect(state)) as connection:
            for table, columns in TABLE_COLUMNS.items():
                rows = tuple(connection.execute(f"PRAGMA table_info({table})"))
                primary = tuple(row[1] for row in sorted(rows, key=lambda row: row[5]) if row[5])
                if tuple(row[1] for row in rows) != columns or primary != PRIMARY_KEYS[table]:
                    return None, ("schema_mismatch", f"clinvar_{table}_schema_mismatch")
            for table, prefix in INDEX_PREFIXES.items():
                if not _indexed(connection, table, prefix):
                    return None, ("schema_mismatch", f"clinvar_{table}_index_mismatch")
            metadata = dict(connection.execute("SELECT key,value FROM metadata"))
    except sqlite3.Error:
        return None, ("query_failed", "clinvar_open_failed")
    if metadata.get("schema_version") != CLINVAR_SCHEMA:
        return None, ("schema_mismatch", "clinvar_schema_version_mismatch")
    release = _iso_date(metadata.get("clinvar_release"))
    if release is None:
        return None, ("schema_mismatch", "clinvar_release_invalid")
    text_fields = ("built_at", "variant_summary_path", "var_citations_path")
    integer_fields = (
        "variant_summary_size",
        "var_citations_size",
        "variant_row_count",
        "citation_row_count",
        "identifier_row_count",
    )
    valid = all(isinstance(metadata.get(field), str) and metadata[field].strip() for field in text_fields)
    valid = valid and all(
        isinstance(metadata.get(field), str) and metadata[field].isdigit() for field in integer_fields
    )
    try:
        if valid:
            datetime.fromisoformat(metadata["built_at"].replace("Z", "+00:00"))
    except ValueError:
        valid = False
    if not valid:
        return None, ("schema_mismatch", "clinvar_metadata_invalid")
    if _fingerprint(state["clinvar_db"]) != current:
        return None, ("schema_mismatch", "clinvar_resource_changed")
    state["clinvar"] = {"version": f"clinvar:{release};schema={CLINVAR_SCHEMA}", "fingerprint": current}
    return state["clinvar"], None


def _citation_data(items):
    sources, ignored, matches = {}, [], []
    for pmid, field, value, kind in items:
        if pmid:
            fields = sources.setdefault(pmid, [])
            if field not in fields:
                fields.append(field)
            matches.append((pmid, field, kind))
        else:
            ignored.append((field, value))
    return {
        "pmids": tuple(sources),
        "sources": {key: tuple(value) for key, value in sources.items()},
        "ignored": tuple(ignored),
        "matches": tuple(dict.fromkeys(matches)),
    }


def _clinvar_citations(rows):
    items = []
    for row in rows:
        value = str(row.get("pmid", "")).strip()
        items.append((value if value.isdigit() else "", "citations.pmid", value, "numeric_field_value"))
    return _citation_data(items)


def _vcep_citations(row):
    items = []
    for raw in row.get("PubMed Articles", "").split(","):
        value = raw.strip()
        if not value:
            continue
        match = FIELD_PMID.fullmatch(value)
        identifier = value if value.isdigit() else match.group("identifier").rstrip(".") if match else ""
        kind = (
            "numeric_field_value"
            if value.isdigit()
            else ("pmid_label_colon" if match and ":" in match.group("separator") else "pmid_label_space")
        )
        items.append((identifier if identifier.isdigit() else "", "PubMed Articles", value, kind))
    for match in SUMMARY_PMID.finditer(row.get("Summary of interpretation", "")):
        value = match.group("identifier").rstrip(".")
        kind = "pmid_label_colon" if ":" in match.group("separator") else "pmid_label_space"
        items.append((value if value.isdigit() else "", "Summary of interpretation", match.group(0), kind))
    return _citation_data(items)


def _citation_content(data):
    return {
        "pmid_sources": {key: list(fields) for key, fields in data["sources"].items()},
        "pmid_matches": [
            {"pmid": pmid, "field": field, "format": kind} for pmid, field, kind in data["matches"]
        ],
        "ignored_citation_identifiers": [
            {"field": field, "value": value} for field, value in data["ignored"]
        ],
    }


def _document(
    provider, kind, locator, version, content, gene, keys, prompt_id, title, provenance, pmids=(), section=""
):
    source = Source(
        provider=provider,
        kind=kind,
        locator=locator,
        version=version,
        content=content,
        genes=(gene,),
        variant_keys=keys,
    )
    return evidence_document(
        source, prompt_id, title, content, provenance=provenance, cited_pmids=pmids, section=section
    )


def _assertion(gene, key, document, basis, *, normalized="", linked="", retracted=False):
    return {
        "gene": gene,
        "variant_key": key,
        "match_basis": basis,
        "assertion_scope": "disease_specific_expert" if linked else "aggregate_allele",
        "document": document,
        "cited_pmids": document["cited_pmids"],
        "normalized_hgvs": normalized,
        "linked_source_id": linked,
        "retracted": retracted,
    }


def _clinvar_match(gene, target, row, citations, version, source_key=None, resolved=None):
    row, citations = dict(row), [dict(item) for item in citations]
    data = _clinvar_citations(citations)
    equivalent = source_key is not None
    basis = "sequence_equivalent" if equivalent else "exact_coordinate"
    record = {
        "record_type": "clinvar_sql_sequence_equivalent_match"
        if equivalent
        else "clinvar_sql_normalized_coordinate_match",
        "record": row,
    }
    if equivalent:
        record["equivalence_proof"] = {
            "match_basis": basis,
            "binding_target_key": variant_text(target),
            "source_database_key": variant_text(source_key),
            "normalized_hgvs": resolved["normalized_hgvs"],
            "reference_status": resolved["status"],
        }
    content = _json({**record, "citation_rows": citations, **_citation_content(data)})
    provenance = (
        ("clinvar_sequence_equivalent", "sql_normalized_record", "local_reference_verified")
        if equivalent
        else ("clinvar_exact_coordinate", "sql_normalized_record")
    )
    document = _document(
        "clinvar",
        "annotation",
        f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{row['variation_id']}/",
        version,
        content,
        gene,
        (source_key or target,),
        f"CLINVAR:{row['variation_id']}:{row['allele_id']}:{hashlib.sha256(content.encode()).hexdigest()[:12]}",
        f"ClinVar: {row['name']}",
        provenance,
        data["pmids"],
        "aggregate_allele_record",
    )
    return {
        "variation_id": str(row["variation_id"]),
        "allele_id": str(row["allele_id"]),
        "citation_sources": data["sources"],
        "evidence": _assertion(
            gene, target, document, basis, normalized=resolved["normalized_hgvs"] if equivalent else ""
        ),
    }


def _associated(raw, gene):
    require(isinstance(raw, str), "clinvar.gene_symbol", "text_required")
    return any(value.strip().upper() == gene for value in re.split(r"[;,|]", raw))


def scope_native_transcripts(state, tools):
    """Keep regional candidates, but do not assign another gene's coding HGVS to them."""
    info, failure = _clinvar_info(state)
    if failure:
        return tools
    scoped, cache = [], {}
    try:
        with closing(_connect(state)) as connection:
            connection.row_factory = sqlite3.Row
            for tool in tools:
                if tool["tool_id"] != "lirical":
                    scoped.append(tool)
                    continue
                genes = []
                for gene in tool["genes"]:
                    aliases = {gene["gene"]}
                    for record in gene.get("records", ()):
                        mapping = record.get("gene_mapping_source")
                        if mapping is not None and mapping.provider == "hgnc":
                            row = json.loads(mapping.content)
                            names = {
                                value.strip() for field in ("symbol", "prev_symbol", "alias_symbol")
                                for value in (row.get(field) or "").split("|") if value.strip()
                            }
                            if gene["gene"] in names:
                                aliases.update(names)
                    variants = []
                    for variant in gene["variants"]:
                        key = variant.get("key")
                        transcript, hgvsc = variant.get("transcript"), variant.get("hgvs_c")
                        owners = {}
                        if key and transcript and hgvsc:
                            if key not in cache:
                                cache[key] = connection.execute(
                                    "SELECT variation_id,allele_id,gene_symbol,name FROM variants WHERE assembly=? AND chromosome=? AND position=? AND reference=? AND alternate=?",
                                    key,
                                ).fetchall()
                            pattern = rf"^{re.escape(transcript.split('.')[0])}(?:\.\d+)?\(([^()]+)\):{re.escape(hgvsc)}(?:\s|$)"
                            for row in cache[key]:
                                match = re.match(pattern, row["name"] or "")
                                if match and _associated(row["gene_symbol"], match[1]):
                                    owners.setdefault(match[1], []).append(dict(row))
                        if len(owners) == 1 and not aliases.intersection(owners):
                            owner = next(iter(owners))
                            validate_gene(owner)
                            scope = {
                                "scope": "regional_context",
                                "candidate_gene": gene["gene"],
                                "transcript_gene": owner,
                                "reported_identity": {
                                    name: variant.get(name) for name in ("transcript", "hgvs_c", "hgvs_p")
                                },
                                "clinvar_records": owners[owner],
                                "database_version": info["version"],
                                "basis": "same_allele_transcript_accession_and_coding_hgvs",
                            }
                            variant = {
                                **variant, "transcript": None, "hgvs_c": None, "hgvs_p": None,
                                "details": {**variant["details"], "gene_scope": scope},
                            }
                        variants.append(variant)
                    genes.append({**gene, "variants": tuple(variants)})
                scoped.append({**tool, "genes": tuple(genes)})
    except (sqlite3.Error, MendelPrimeError):
        return tools
    return tuple(scoped) if _fingerprint(state["clinvar_db"]) == info["fingerprint"] else tools


def _read_matches(connection, gene, target, rows, version, resolved=None):
    output = []
    for row in rows:
        citations = connection.execute(
            "SELECT variation_id,allele_id,pmid FROM citations WHERE variation_id=? AND allele_id=? ORDER BY pmid",
            (row["variation_id"], row["allele_id"]),
        ).fetchall()
        source_key = (
            variant_key(target[0], row["chromosome"], row["position"], row["reference"], row["alternate"])
            if resolved
            else None
        )
        output.append(_clinvar_match(gene, target, row, citations, version, source_key, resolved))
    return tuple(output)


def _clinvar_matches(state, gene, target, reference, exclusions=None):
    if target[1] == "MT":
        if reference is None:
            return (), ("not_applicable", "mitochondrial_reference_unverified")
        try:
            resolved = reference["resolve"](target)
        except Exception:
            return (), ("query_failed", "clinvar_reference_failed")
        if not isinstance(resolved, dict) or resolved.get("raw_key") != target:
            return (), ("query_failed", "clinvar_reference_result_invalid")
        if "mitochondrial_reference_unsupported" in resolved.get("issues", ()):
            return (), ("not_applicable", "mitochondrial_reference_unsupported")
        if resolved.get("status") != "sequence_verified" or resolved.get("reference_accession") != "NC_012920.1":
            return (), ("query_failed", "clinvar_reference_unresolved")
    info, failure = _clinvar_info(state)
    if failure:
        return (), failure
    try:
        with closing(_connect(state)) as connection:
            connection.row_factory = sqlite3.Row
            rows = connection.execute(
                "SELECT * FROM variants WHERE assembly=? AND chromosome=? AND position=? AND reference=? AND alternate=? ORDER BY variation_id,allele_id",
                target,
            ).fetchall()
            scoped = [row for row in rows if _associated(row["gene_symbol"], gene)]
            if rows and not scoped:
                for row in rows:
                    symbols = tuple(value.strip() for value in re.split(r"[;,|]", row["gene_symbol"]))
                    for symbol in symbols:
                        validate_gene(symbol)
                    if exclusions is not None:
                        exclusions.append({
                            "gene_symbols": symbols, "name": row["name"],
                            "variation_id": str(row["variation_id"]), "allele_id": str(row["allele_id"]),
                            "variant_key": target, "database_version": info["version"],
                            "locator": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{row['variation_id']}/",
                        })
            matches = _read_matches(connection, gene, target, scoped, info["version"])
    except sqlite3.Error:
        return (), ("query_failed", "clinvar_query_failed")
    except MendelPrimeError:
        return (), ("schema_mismatch", "clinvar_record_invalid")
    if _fingerprint(state["clinvar_db"]) != info["fingerprint"]:
        return (), ("schema_mismatch", "clinvar_resource_changed")
    if matches:
        return matches, None
    failure = (
        ("not_applicable", "clinvar_assertion_outside_gene") if rows and not scoped else ("not_found", "")
    )
    return _equivalent_matches(state, gene, target, reference, info, failure) if reference else ((), failure)


def _equivalent_matches(state, gene, target, reference, info, exact_failure):
    try:
        resolved = reference["resolve"](target)
    except Exception:
        return (), ("query_failed", "clinvar_reference_failed")
    if not isinstance(resolved, dict) or resolved.get("raw_key") != target:
        return (), ("query_failed", "clinvar_reference_result_invalid")
    if resolved.get("status") != "sequence_verified":
        return (), ("query_failed", "clinvar_reference_unresolved")
    normalized, region = resolved.get("normalized_hgvs"), resolved.get("equivalent_region")
    if (
        not isinstance(normalized, str)
        or not normalized
        or normalized != normalized.strip()
        or any(c.isspace() for c in normalized)
    ):
        return (), ("query_failed", "clinvar_reference_result_invalid")
    if (
        type(region) is not tuple
        or len(region) != 2
        or any(type(value) is not int for value in region)
        or region[0] < 0
        or region[1] < region[0]
        or not region[0] <= target[2] <= region[1]
    ):
        return (), ("query_failed", "clinvar_equivalent_region_invalid")
    start, end = max(1, region[0]), region[1]
    if end - start + 1 > 1_000_000:
        return (), ("query_failed", "clinvar_equivalent_region_too_large")
    try:
        with closing(_connect(state)) as connection:
            connection.row_factory = sqlite3.Row
            rows = connection.execute(
                "SELECT * FROM variants WHERE assembly=? AND chromosome=? AND position BETWEEN ? AND ? ORDER BY variation_id,allele_id LIMIT ?",
                (*target[:2], start, end, 10_001),
            ).fetchall()
            if len(rows) > 10_000:
                return (), ("query_failed", "clinvar_equivalent_query_too_broad")
            scoped = [row for row in rows if _associated(row["gene_symbol"], gene)]
            candidates = [
                (
                    row,
                    variant_key(
                        target[0], row["chromosome"], row["position"], row["reference"], row["alternate"]
                    ),
                )
                for row in scoped
            ]
            selected = []
            for row, source_key in candidates:
                try:
                    equivalent = reference["equivalent"](target, source_key)
                except Exception:
                    return (), ("query_failed", "clinvar_reference_failed")
                if type(equivalent) is not bool:
                    return (), ("query_failed", "clinvar_reference_result_invalid")
                if equivalent:
                    selected.append(row)
            matches = _read_matches(connection, gene, target, selected, info["version"], resolved)
    except sqlite3.Error:
        return (), ("query_failed", "clinvar_query_failed")
    except MendelPrimeError:
        return (), ("schema_mismatch", "clinvar_record_invalid")
    if _fingerprint(state["clinvar_db"]) != info["fingerprint"]:
        return (), ("schema_mismatch", "clinvar_resource_changed")
    if matches:
        return matches, None
    # Unrelated records in a normalization window are not gene-scope failures.
    return (), exact_failure


def _locator(value, fallback):
    return value.strip() if value.strip().startswith("https://") else fallback


def _gene_context(state, gene):
    rows, failure = _rows(state, "gene_validity")
    if failure:
        return (), failure
    output = []
    try:
        for row in rows["by_key"].get(gene, ()):
            content = _json({"record_type": "clingen_gene_validity", "row": row})
            digest = hashlib.sha256(content.encode()).hexdigest()
            output.append(
                _document(
                    "clingen",
                    "gene_disease",
                    _locator(
                        row.get("ONLINE REPORT", ""),
                        f"resource:clingen/gene_disease_validity.csv#record/{digest[:16]}",
                    ),
                    rows["version"],
                    content,
                    gene,
                    (),
                    f"CLINGEN:GENE_VALIDITY:{gene}:{digest[:12]}",
                    f"ClinGen gene validity: {gene} - {row.get('DISEASE LABEL', '')}",
                    ("clingen_gene_validity",),
                    section="gene_disease_validity",
                )
            )
    except MendelPrimeError:
        return (), ("schema_mismatch", "clingen_gene_record_invalid")
    return (tuple(output), None) if output else ((), ("not_found", ""))


def _vcep_evidence(gene, target, match, row, version, citations):
    anchor = match["evidence"]
    linked = anchor["document"]["source"].source_id
    join = {
        "match_basis": "clinvar_variation_id_join",
        "variation_id": match["variation_id"],
        "allele_id": match["allele_id"],
        "clinvar_source_id": linked,
        "variant_key": variant_text(target),
    }
    provenance = ["clingen_vcep", "clinvar_variation_id_join"]
    if anchor["match_basis"] == "sequence_equivalent":
        join.update(
            clinvar_anchor_match_basis="sequence_equivalent",
            clinvar_source_database_key=variant_text(anchor["document"]["source"].variant_keys[0]),
            normalized_hgvs=anchor["normalized_hgvs"],
        )
        provenance.append("sequence_equivalent_clinvar_anchor")
    retracted = row.get("Retracted", "").strip().casefold() in {"1", "true", "yes"}
    content = _json(
        {
            "record_type": "clingen_vcep_assertion",
            "row": row,
            "verified_clinvar_join": join,
            "retraction_status": "retracted" if retracted else "not_retracted",
            **_citation_content(citations),
        }
    )
    if retracted:
        provenance.append("retracted_source_record")
    digest = hashlib.sha256(content.encode()).hexdigest()
    document = _document(
        "clingen",
        "annotation",
        _locator(
            row.get("Evidence Repo Link", ""),
            f"resource:clingen/variant_pathogenicity.tsv#record/{digest[:16]}",
        ),
        version,
        content,
        gene,
        (target,),
        f"CLINGEN:VCEP:{row.get('Uuid') or match['variation_id']}:{digest[:12]}",
        f"ClinGen VCEP: {gene} - {row.get('Disease', '')}",
        tuple(provenance),
        citations["pmids"],
        "disease_specific_expert_assertion",
    )
    return _assertion(gene, target, document, "clinvar_variation_id_join", linked=linked, retracted=retracted)


def _seeds(evidence, sources, route):
    return [
        {
            "gene": evidence["gene"],
            "pmid": pmid,
            "variant_key": evidence["variant_key"],
            "discovered_from_source_id": evidence["document"]["source"].source_id,
            "route": route,
            "source_fields": sources[pmid],
        }
        for pmid in evidence["cited_pmids"]
    ]


def _vcep_matches(state, gene, target, matches):
    rows, failure = _rows(state, "vcep")
    if failure:
        return (), (), failure
    joined, seeds, other = [], [], False
    try:
        for match in matches:
            records = rows["by_key"].get(match["variation_id"], ())
            other = other or any(row.get("HGNC Gene Symbol") != gene for row in records)
            for row in records:
                if row.get("HGNC Gene Symbol") != gene:
                    continue
                data = _vcep_citations(row)
                evidence = _vcep_evidence(gene, target, match, row, rows["version"], data)
                joined.append(evidence)
                if not evidence["retracted"]:
                    seeds.extend(_seeds(evidence, data["sources"], "clingen_citation"))
    except MendelPrimeError:
        return (), (), ("schema_mismatch", "clingen_vcep_record_invalid")
    if joined:
        return tuple(joined), tuple(seeds), None
    return (), (), ("gene_scope_mismatch", "clingen_vcep_gene_scope_mismatch") if other else ("not_found", "")


def _outcome(provider, gene, key, values, failure=None):
    status, code = ("found", "") if values else failure or ("not_found", "")
    return {"provider": provider, "gene": gene, "variant_key": key, "status": status, "error_code": code}


def lookup(state, gene, variants, reference=None):
    validate_gene(gene)
    require(
        type(variants) is tuple and all(type(key) is tuple and len(key) == 5 for key in variants),
        "variants",
        "variant_key_tuple_required",
    )
    for key in variants:
        variant_key(*key)
    require(len(variants) == len(set(variants)), "variants", "duplicate_variant")
    context, failure = _gene_context(state, gene)
    outcomes, evidence, seeds = [_outcome("clingen_gene_validity", gene, None, context, failure)], [], []
    for key in variants:
        exclusions = []
        matches, failure = _clinvar_matches(state, gene, key, reference, exclusions)
        outcome = _outcome("clinvar", gene, key, matches, failure)
        if exclusions:
            outcome["out_of_scope_records"] = tuple(exclusions)
        outcomes.append(outcome)
        for match in matches:
            evidence.append(match["evidence"])
            seeds.extend(_seeds(match["evidence"], match["citation_sources"], "clinvar_citation"))
        if matches:
            joined, citations, failure = _vcep_matches(state, gene, key, matches)
            evidence.extend(joined)
            seeds.extend(citations)
            outcomes.append(_outcome("clingen_vcep", gene, key, joined, failure))
        else:
            status = failure[0] if failure else "not_found"
            code = {
                "not_found": "no_exact_clinvar_join",
                "not_applicable": "no_applicable_clinvar_join",
                "gene_scope_mismatch": "no_gene_scoped_clinvar_join",
            }.get(status, "clinvar_join_unavailable")
            outcomes.append(_outcome("clingen_vcep", gene, key, (), (status, code)))
    unique = {(row["document"]["source"].source_id, row["variant_key"]): row for row in evidence}
    unique_seeds = {tuple(row.items()): row for row in seeds}
    result = {
        "gene": gene,
        "variant_evidence": tuple(unique.values()),
        "gene_context": context,
        "citation_seeds": tuple(unique_seeds.values()),
        "lookups": tuple(outcomes),
    }
    validate_result(result)
    return result


def validate_result(result):
    gene = result["gene"]
    validate_gene(gene)
    evidence = {row["document"]["source"].source_id: row for row in result["variant_evidence"]}
    for row in result["variant_evidence"]:
        source, key = row["document"]["source"], row["variant_key"]
        require(row["gene"] == gene, "variant_evidence", "gene_mismatch")
        require(source.genes == (gene,), "variant_evidence.document", "source_gene_scope_mismatch")
        basis, scope = row["match_basis"], row["assertion_scope"]
        require(
            basis in {"exact_coordinate", "sequence_equivalent", "clinvar_variation_id_join"},
            "variant_evidence.match_basis",
            "unsupported_match_basis",
        )
        require(
            scope in {"aggregate_allele", "disease_specific_expert"},
            "variant_evidence.assertion_scope",
            "unsupported_assertion_scope",
        )
        if basis == "sequence_equivalent":
            require(
                len(source.variant_keys) == 1
                and source.variant_keys[0] != key
                and source.variant_keys[0][:2] == key[:2],
                "variant_evidence.document",
                "equivalent_source_variant_scope_mismatch",
            )
            normalized = row["normalized_hgvs"]
            require(
                isinstance(normalized, str)
                and normalized
                and normalized == normalized.strip()
                and not any(c.isspace() for c in normalized),
                "variant_evidence.normalized_hgvs",
                "equivalence_proof_required",
            )
        else:
            require(
                source.variant_keys == (key,), "variant_evidence.document", "source_variant_scope_mismatch"
            )
            require(
                not row["normalized_hgvs"],
                "variant_evidence.normalized_hgvs",
                "equivalence_proof_not_allowed",
            )
        pmids = row["cited_pmids"]
        require(
            type(pmids) is tuple and all(isinstance(item, str) for item in pmids),
            "variant_evidence.cited_pmids",
            "pmid_tuple_required",
        )
        require(
            all(item.isdigit() for item in pmids), "variant_evidence.cited_pmids", "numeric_pmids_required"
        )
        require(len(pmids) == len(set(pmids)), "variant_evidence.cited_pmids", "duplicate_pmid")
        require(
            row["cited_pmids"] == row["document"]["cited_pmids"],
            "variant_evidence.cited_pmids",
            "document_pmids_mismatch",
        )
        linked = row["linked_source_id"]
        require(
            isinstance(linked, str)
            and (not linked or linked == linked.strip() and not any(c.isspace() for c in linked)),
            "variant_evidence.linked_source_id",
            "canonical_source_id_required",
        )
        require(
            bool(linked) if basis == "clinvar_variation_id_join" else not linked,
            "variant_evidence.linked_source_id",
            "join_source_required" if basis == "clinvar_variation_id_join" else "join_source_not_allowed",
        )
        require(type(row["retracted"]) is bool, "variant_evidence.retracted", "boolean_required")
        if linked:
            require(linked in evidence, "variant_evidence", "linked_source_missing")
            anchor = evidence[linked]
            require(
                anchor["match_basis"] in {"exact_coordinate", "sequence_equivalent"}
                and anchor["gene"] == gene
                and anchor["variant_key"] == key,
                "variant_evidence",
                "linked_exact_allele_mismatch",
            )
    for document in result["gene_context"]:
        require(
            document["source"].genes == (gene,) and not document["source"].variant_keys,
            "gene_context",
            "gene_only_source_required",
        )
    for field in ("citation_seeds", "lookups"):
        require(all(row["gene"] == gene for row in result[field]), field, "gene_mismatch")
    for row in result["citation_seeds"]:
        require(
            isinstance(row["pmid"], str) and row["pmid"].isdigit(),
            "citation_seed.pmid",
            "numeric_pmid_required",
        )
        identifier = row["discovered_from_source_id"]
        require(
            isinstance(identifier, str) and identifier and not any(c.isspace() for c in identifier),
            "citation_seed.discovered_from_source_id",
            "canonical_source_id_required",
        )
        require(
            row["route"] in {"clinvar_citation", "clingen_citation"},
            "citation_seed.route",
            "unsupported_citation_route",
        )
        fields = row["source_fields"]
        require(
            type(fields) is tuple
            and fields
            and len(fields) == len(set(fields))
            and set(fields) <= {"citations.pmid", "PubMed Articles", "Summary of interpretation"},
            "citation_seed.source_fields",
            "citation_source_fields_required",
        )
    for row in result["lookups"]:
        provider, code, status = row["provider"], row["error_code"], row["status"]
        require(
            isinstance(provider, str) and provider and not any(c.isspace() for c in provider),
            "lookup.provider",
            "canonical_provider_required",
        )
        require(
            status
            in {
                "found",
                "not_found",
                "not_applicable",
                "resource_unavailable",
                "query_failed",
                "schema_mismatch",
                "gene_scope_mismatch",
            },
            "lookup.status",
            "unsupported_lookup_status",
        )
        require(
            isinstance(code, str) and not any(c.isspace() for c in code),
            "lookup.error_code",
            "safe_error_code_required",
        )
        require(status != "found" or not code, "lookup.error_code", "found_must_not_have_error")


def sources(result):
    documents = (*[row["document"] for row in result["variant_evidence"]], *result["gene_context"])
    return tuple({document["source"].source_id: document["source"] for document in documents}.values())
