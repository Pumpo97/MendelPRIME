"""Source-bound claim selection, English composition, and direct Markdown output."""

from __future__ import annotations

import html
import json
import re
from pathlib import PurePosixPath, PureWindowsPath

from mendelprime import prompting
from mendelprime.domain import (
    MendelPrimeError,
    Source,
    require,
    validate_case_key,
    validate_gene,
    variant_key,
    variant_text,
)
from mendelprime.harness.memory import OTEM, bound_alleles, validate_claim_shape

_MAX_TOKENS = 6000
_INPUT_CHAR_LIMIT = 250_000
_PACKET_CHAR_TARGET = 200_000
_POPULATION_SUBSETS = ("controls_", "non_neuro_", "non_cancer_", "non_topmed_")
_SCOPES = {"case", "gene", "allele"}
_POLARITIES = {"support", "limitation", "background"}
_ROLES = {"summary", "support", "limitation", "ranking"}
_MEDICAL = {
    "guideline",
    "gene_disease",
    "literature",
    "population",
    "prediction",
    "annotation",
    "similar_case",
}
_VARIANT_ALIAS = re.compile(r"(?<![A-Za-z0-9_])V[0-9]{3,}(?![A-Za-z0-9_])")
_CJK = re.compile(
    r"[\u3040-\u30ff\u31f0-\u31ff\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff\uac00-\ud7af"
    r"\U00020000-\U0002a6df\U0002a700-\U0002b73f\U0002b740-\U0002b81f"
    r"\U0002b820-\U0002ceaf\U0002ceb0-\U0002ebef\U0002ebf0-\U0002ee5f"
    r"\U0002f800-\U0002fa1f\U00030000-\U0003134f\U00031350-\U000323af]"
)
_ORDERED_LIST = re.compile(r"(?m)^(\s*\d+)\.(?=\s)")
_LEADING_LIST = re.compile(r"(?m)^(\s*)([-+])(?=\s)")
_HORIZONTAL_RULE = re.compile(r"(?m)^(\s*)(-{3,})(\s*)$")
_LONG_HEX = re.compile(r"(?i)(?<![0-9a-f])[0-9a-f]{32,}(?![0-9a-f])")
_ABSOLUTE_PATH = re.compile(
    r"(?i)(?:(?<![a-z0-9])(?:[a-z]:[\\/]|\\\\[^\\/\s]+[\\/])"
    r"|(?<![a-z0-9/])/(?!/)|(?<![a-z0-9:/])//[^/\s]+/)"
)


def _error(code):
    raise MendelPrimeError(code, category="response")


def _check(condition, code):
    if not condition:
        _error(code)


def _object(value, code):
    _check(type(value) is dict, code)
    return value


def _rows(value, code, *, nonempty=True):
    _check(type(value) is list and (value or not nonempty), code)
    return value


def _text(value, code):
    _check(isinstance(value, str) and bool(value.strip()), code)
    return value


def _strings(value, code, *, nonempty=True):
    rows = _rows(value, code, nonempty=nonempty)
    _check(all(isinstance(item, str) and item for item in rows), code)
    result = tuple(rows)
    _check(len(set(result)) == len(result), code)
    return result


def _identity(record):
    row = record.get("identity") if type(record) is dict else None
    require(
        type(row) is dict and set(row) == {"gene", "variant_key"}, "assessment", "assessment_identity_invalid"
    )
    validate_gene(row["gene"], "assessment.identity.gene")
    value = row["variant_key"]
    if type(value) is tuple:
        key = variant_key(*value)
    elif isinstance(value, str):
        parts = value.split(":")
        require(len(parts) == 5, "assessment", "variant_key_invalid")
        key = variant_key(parts[0], parts[1], int(parts[2]), parts[3], parts[4])
    elif type(value) is dict and set(value) == {
        "assembly",
        "chromosome",
        "position",
        "reference",
        "alternate",
    }:
        key = variant_key(**value)
    else:
        raise MendelPrimeError("variant_key_invalid", field="assessment")
    return row["gene"], key


def _sequence(value, field):
    require(type(value) in {list, tuple}, field, "assessment_record_invalid")
    return tuple(value)


def _source_scope(source):
    if source.variant_keys:
        return "patient_exact_allele" if source.kind == "patient_observation" else "exact_allele"
    if source.genes:
        return "gene_or_region"
    return "case_or_general"


def _source_uses(source):
    if source.kind == "patient_observation":
        return ["patient_qc", "case_observation"] if source.variant_keys else ["case_observation"]
    if source.kind == "derived":
        return ["assessment_narration"]
    if source.kind == "tool_observation":
        if source.provider == "mendelprime" and source.locator == "method:ranking-observations/frozen-order":
            return ["ranking_context"]
        return ["tool_reported_observation", "identity_context"]
    if source.kind in _MEDICAL:
        return ["exact_allele"] if source.variant_keys else ["gene_background", "regional_context"]
    return ["descriptive_context"]


def _pack_info(tokens):
    passages, current = [], ""
    for token in tokens:
        if len(token) > 600:
            if current:
                passages.append(current)
                current = ""
            key, separator, value = token.partition("=")
            parts = value.split(",") if key.lower() == "vep" and separator else [token]
            passages.extend(passage for part in parts for passage in _prose_passages(part))
            continue
        candidate = f"{current};{token}" if current else token
        if len(candidate) <= 600:
            current = candidate
        else:
            passages.append(current)
            current = token
    if current:
        passages.append(current)
    return passages


def _population_passages(raw, omitted_info=(), *, alternate=None):
    if not isinstance(raw, dict):
        return []
    passages = [passage for record in raw.get("coverage", ()) for passage in _prose_passages(str(record))]
    for record in raw.get("sites", ()):
        fields = record.rstrip("\n").split("\t", 7)
        if len(fields) != 8:
            passages.extend(_prose_passages(record))
            continue
        passages.append("\t".join(fields[:7]))
        tokens = fields[7].split(";")
        indexed = {token.partition("=")[0]: token for token in tokens}
        alternates = fields[4].split(",")
        allele_index = (
            alternates.index(alternate) if alternate in alternates else 0 if len(alternates) == 1 else None
        )
        # Keep subgroup aggregates and extrema with their denominators, not every
        # cross-product of study subset, ancestry and sex in the raw VCF matrix.
        selected = set()
        for prefix in _POPULATION_SUBSETS:
            for metric in ("AF", "nhomalt", "nhemialt", "faf95", "faf99"):
                values = []
                for key, token in indexed.items():
                    if key.startswith(prefix + metric):
                        try:
                            entries = token.partition("=")[2].split(",")
                            number = float(
                                entries[allele_index]
                                if len(entries) > 1 and allele_index is not None
                                else entries[0]
                            )
                            values.append((number, key))
                        except (ValueError, IndexError):
                            selected.add(key)
                if values:
                    for _, key in (min(values), max(values)):
                        suffix = key[len(prefix + metric) :]
                        selected.update(
                            prefix + field + suffix
                            for field in ("AC", "AN", "AF", "nhomalt", "nhemialt", "faf95", "faf99")
                        )
        visible = []
        for token in tokens:
            key = token.partition("=")[0]
            subset = next((prefix for prefix in _POPULATION_SUBSETS if key.startswith(prefix)), None)
            auxiliary = "hist" in key or (
                subset is not None
                and allele_index is not None
                and any(
                    key[len(subset) :].startswith(field + "_")
                    for field in ("AC", "AN", "AF", "nhomalt", "nhemialt", "faf95", "faf99")
                )
                and key not in selected
                and key[len(subset) :] not in {"AC", "AN", "AF", "nhomalt", "nhemialt", "popmax"}
                and not key.endswith("_popmax")
            )
            if key.lower() in omitted_info or auxiliary:
                passages.extend(_pack_info(visible))
                visible = []
            else:
                visible.append(token)
        passages.extend(_pack_info(visible))
    return [json.dumps(passage, ensure_ascii=False)[1:-1] for passage in passages]


def _field_span(content, key, expected):
    decoder = json.JSONDecoder()
    token = json.dumps(key, ensure_ascii=False)
    spans = []
    start = 0
    while True:
        index = content.find(token, start)
        if index < 0:
            break
        cursor = index + len(token)
        while cursor < len(content) and content[cursor].isspace():
            cursor += 1
        if cursor < len(content) and content[cursor] == ":":
            cursor += 1
            while cursor < len(content) and content[cursor].isspace():
                cursor += 1
            try:
                value, end = decoder.raw_decode(content, cursor)
            except json.JSONDecodeError:
                value = None
            else:
                if value == expected:
                    spans.append(content[index:end])
        start = index + len(token)
    return spans[0] if spans else None


def _exact_passages(content, value, *, source_kind):
    limit = None if source_kind == "patient_observation" else 600
    passages = []

    def visit(item):
        if isinstance(item, dict):
            for key, child in item.items():
                span = _field_span(content, key, child)
                if span is not None and (limit is None or len(span) <= limit):
                    passages.append(span)
                else:
                    visit(child)
        elif isinstance(item, list):
            for child in item:
                visit(child)
        elif isinstance(item, str) and item:
            encoded = json.dumps(item, ensure_ascii=False)[1:-1]
            if encoded in content:
                passages.extend([encoded] if limit is None else _prose_passages(encoded))

    visit(value)
    return list(dict.fromkeys(passages))


def _prose_passages(content):
    passages = []
    for block in (part for part in re.split(r"\n+", content) if part):
        while len(block) > 600:
            boundary = block.rfind(" ", 0, 601)
            boundary = boundary if boundary > 0 else 600
            passages.append(block[:boundary])
            block = block[boundary:].lstrip()
        if block:
            passages.append(block)
    return passages


def _reliable_annotation_refs(source, prepared):
    refs = []
    for alias, other in prepared["source_by_alias"].items():
        if other.kind != "annotation" or not set(source.variant_keys).intersection(other.variant_keys):
            continue
        try:
            record = json.loads(other.content)
        except (json.JSONDecodeError, TypeError):
            continue
        if (
            isinstance(record, dict)
            and record.get("status") in (None, "annotated")
            and record.get("selected_transcript")
            and set(source.genes).intersection(other.genes)
            and record.get("gene", source.genes[0]) in source.genes
        ):
            refs.append(alias)
    return refs


def _selected_json_passages(content, document, keys, *, source_kind):
    passages = []
    for key in keys:
        if key not in document:
            continue
        span = _field_span(content, key, document[key])
        if span is not None and (source_kind == "patient_observation" or len(span) <= 600):
            passages.append(span)
        else:
            passages.extend(_exact_passages(content, document[key], source_kind=source_kind))
    return passages


def _source_view(alias, source, prepared):
    document = None
    try:
        document = json.loads(source.content)
    except (json.JSONDecodeError, RecursionError):
        passages = _prose_passages(source.content)
    else:
        factored_annotations = (
            _reliable_annotation_refs(source, prepared) if source.kind == "population" else []
        )
        if source.kind == "population" and isinstance(document, dict):
            passages = _selected_json_passages(
                source.content,
                document,
                tuple(key for key in document if key != "raw_records"),
                source_kind=source.kind,
            )
            passages.extend(
                _population_passages(
                    document.get("raw_records"),
                    {"vep"} if factored_annotations else set(),
                    alternate=source.variant_keys[0][4] if source.variant_keys else None,
                )
            )
        elif source.kind == "patient_observation" and isinstance(document, dict):
            passages = _selected_json_passages(
                source.content,
                document,
                (
                    "alternate_index",
                    "filter_text",
                    "quality",
                    "samples",
                    "target_mappings",
                    "vcf_variant_key",
                    "phenotypes",
                    "clinical_description",
                ),
                source_kind=source.kind,
            )
        elif source.provider == "spliceai" and isinstance(document, dict):
            passages = _selected_json_passages(
                source.content,
                document,
                tuple(key for key in document if key != "header"),
                source_kind=source.kind,
            )
            passages.extend(
                json.dumps(line, ensure_ascii=False)[1:-1]
                for line in document.get("header", ())
                if line.startswith("##INFO=<ID=SpliceAI,")
            )
        elif source.kind == "derived" and source.locator.startswith("method:post-ranking-acmg/"):
            decision = document.get("decision", {}) if isinstance(document, dict) else {}
            passages = _selected_json_passages(
                source.content,
                decision,
                (
                    "classification",
                    "confidence",
                    "disease",
                    "identity",
                    "patient_fit",
                    "rationale",
                    "summary",
                ),
                source_kind=source.kind,
            )
        else:
            passages = _exact_passages(source.content, document, source_kind=source.kind)
    if not passages:
        passages = [source.content] if len(source.content) <= 600 else _prose_passages(source.content)
    factoring = {}
    if isinstance(document, dict):
        if source.kind == "population" and "raw_records" in document:
            factoring["population_raw_records"] = (
                "broad populations, subset aggregates/extrema, counts, FILTER and coverage; "
                "full subgroup matrix and histograms retained in immutable provenance"
            )
            if factored_annotations:
                factoring["vep"] = {
                    "duplicate_representation": "explicit selected-transcript annotation source",
                    "exact_passage_source_refs": factored_annotations,
                }
        if source.kind == "patient_observation" and any(
            isinstance(sample, dict) and "raw_fields" in sample for sample in document.get("samples", ())
        ):
            factoring["raw_sample_fields"] = "duplicated_by_named_sample_fields"
    publication_binding = {}
    if source.kind == "literature" and ";scope=verified-exact-paper-binding-v1;" in source.version:
        parents = [prepared["source_by_id"][source_id] for source_id in source.parent_source_ids]
        original = next(
            (
                parent
                for parent in parents
                if parent.kind == "literature"
                and parent.provider == source.provider
                and parent.content == source.content
            ),
            None,
        )
        receipt = next(
            (parent for parent in parents if parent.version == "verified-exact-paper-binding-v1"),
            None,
        )
        if original is not None:
            ref = prepared["source_alias_by_id"].get(original.source_id)
            if ref is not None:
                publication_binding["same_publication_as"] = ref
        if receipt is not None:
            ref = prepared["source_alias_by_id"].get(receipt.source_id)
            if ref is not None:
                publication_binding["binding_receipt_ref"] = ref
    return {
        "source_ref": alias,
        "provider": source.provider,
        "kind": source.kind,
        "locator": source.locator,
        "version": source.version,
        "genes": list(source.genes),
        "variant_refs": [prepared["variant_alias_by_key"][key] for key in source.variant_keys],
        "evidence_scope": _source_scope(source),
        "eligible_uses": _source_uses(source),
        "claim_scopes": ["allele"]
        if source.variant_keys
        else ["gene", "allele"]
        if source.genes
        else ["case"]
        if source.provider == "plink2-king" and source.locator == "method:family/king-qc"
        else ["case", "gene", "allele"],
        "exact_passages": list(dict.fromkeys(passages)),
        "parent_source_refs": [
            prepared["source_alias_by_id"][parent]
            for parent in source.parent_source_ids
            if parent in prepared["source_alias_by_id"]
        ],
        **publication_binding,
        **(
            {
                "assessment_context_ref": {
                    "gene": source.genes[0],
                    "variant_ref": prepared["variant_alias_by_key"][source.variant_keys[0]],
                }
            }
            if source.kind == "derived"
            and source.locator.startswith("method:post-ranking-acmg/")
            and len(source.genes) == len(source.variant_keys) == 1
            else {}
        ),
        **({"factoring": factoring} if factoring else {}),
    }


def _context(value):
    required = {"case_key", "fixed_genes", "sources", "available_source_ids", "assessment"}
    require(
        type(value) is dict
        and required <= set(value)
        and not set(value) - required - {"assessment_availability"},
        "report_context",
        "schema_mismatch",
    )
    validate_case_key(value["case_key"])
    genes = value["fixed_genes"]
    require(
        type(genes) is tuple and 1 <= len(genes) <= 5 and len(set(genes)) == len(genes),
        "fixed_genes",
        "top1_to_top5_tuple_required",
    )
    for gene in genes:
        validate_gene(gene)
    sources = value["sources"]
    require(
        type(sources) is tuple and all(isinstance(source, Source) for source in sources),
        "sources",
        "source_tuple_required",
    )
    source_by_id = {source.source_id: source for source in sources}
    require(len(source_by_id) == len(sources), "sources", "duplicate_source")
    memory = OTEM(value["case_key"])
    memory.add(sources)
    available = value["available_source_ids"]
    require(
        type(available) is tuple and available and len(set(available)) == len(available),
        "available_source_ids",
        "nonempty_tuple_required",
    )
    for source_id in available:
        require(
            isinstance(source_id, str) and source_id in source_by_id,
            "available_source_ids",
            "unknown_source",
        )
        source = memory.get_source(source_id)
        require(
            not source.genes or not set(source.genes).isdisjoint(genes),
            "available_source_ids",
            "unrelated_gene_source",
        )
    assessment, records = value["assessment"], None
    if assessment is not None:
        if type(assessment) is dict and assessment.get("workflow") == "clinical_profile_v1":
            require(
                set(assessment) == {"workflow", "profiles", "final", "sources"},
                "assessment",
                "assessment_record_invalid",
            )
            from mendelprime.tools.acmg import validate_result

            validate_result(assessment)
            require(
                all(source_by_id.get(source.source_id) == source for source in assessment["sources"]),
                "assessment.sources",
                "assessment_source_unavailable",
            )
            profiles = _sequence(assessment["profiles"], "assessment.profiles")
            final = _sequence(assessment["final"], "assessment.final")
            profile_ids = tuple(_identity(row) for row in profiles)
            final_ids = tuple(_identity(row) for row in final)
            require(
                bool(profile_ids) and profile_ids == final_ids,
                "assessment",
                "assessment_coverage_or_order_invalid",
            )
            require(all(gene in genes for gene, _ in final_ids), "assessment", "unrelated_gene_assessment")
            records = {
                "workflow": "clinical_profile_v1",
                "profiles": dict(zip(profile_ids, profiles, strict=True)),
                "final": tuple(zip(final_ids, final, strict=True)),
            }
        else:
            require(
                type(assessment) is dict and {"independent", "partners", "final"} <= set(assessment),
                "assessment",
                "assessment_record_invalid",
            )
            groups = tuple(
                _sequence(assessment[name], f"assessment.{name}")
                for name in ("independent", "partners", "final")
            )
            require(
                bool(groups[0]) and all(len(group) == len(groups[0]) for group in groups),
                "assessment",
                "assessment_coverage_or_order_invalid",
            )
            identities = tuple(tuple(_identity(row) for row in group) for group in groups)
            require(
                identities[0] == identities[1] == identities[2],
                "assessment",
                "assessment_coverage_or_order_invalid",
            )
            require(
                all(gene in genes for gene, _ in identities[2]),
                "assessment",
                "unrelated_gene_assessment",
            )
            records = {
                "independent": dict(zip(identities[0], groups[0], strict=True)),
                "partners": dict(zip(identities[1], groups[1], strict=True)),
                "final": tuple(zip(identities[2], groups[2], strict=True)),
            }
    availability = value.get("assessment_availability")
    if availability is not None:
        completed = tuple(_identity(row) for row in assessment["final"]) if assessment is not None else ()
        require(
            type(availability) is dict
            and set(availability) == {"expected_alleles", "classified_alleles", "unclassified_alleles"},
            "assessment_availability",
            "schema_mismatch",
        )
        missing = tuple(_identity({"identity": row}) for row in availability["unclassified_alleles"])
        require(
            type(availability["expected_alleles"]) is int
            and type(availability["classified_alleles"]) is int
            and availability["expected_alleles"] > 0
            and availability["classified_alleles"] == len(completed)
            and availability["expected_alleles"] == len(completed) + len(missing)
            and missing
            and len(set(missing)) == len(missing)
            and not set(missing).intersection(completed)
            and all(gene in genes for gene, _ in missing),
            "assessment_availability",
            "classification_coverage_invalid",
        )
    return memory, records, source_by_id


def _prepare(context):
    memory, assessment, source_by_id = _context(context)
    sources, variants = [], {}
    for number, source_id in enumerate(context["available_source_ids"], 1):
        source = memory.get_source(source_id)
        sources.append((f"S{number:03d}", source))
        for key in source.variant_keys:
            variants.setdefault(key, None)
    if assessment is not None:
        for (_, key), _ in assessment["final"]:
            variants.setdefault(key, None)
    variant_by_alias = {f"V{number:03d}": key for number, key in enumerate(variants, 1)}
    prepared = {
        "source_by_id": source_by_id,
        "source_by_alias": dict(sources),
        "source_alias_by_id": {source.source_id: alias for alias, source in sources},
        "variant_by_alias": variant_by_alias,
        "variant_alias_by_key": {key: alias for alias, key in variant_by_alias.items()},
        "source_views": {},
        "assessment_sources": {},
        "assessment_by_gene": {},
    }
    for alias, source in sources:
        prepared["source_views"][alias] = _source_view(alias, source, prepared)
        if source.kind == "derived" and len(source.genes) == 1:
            gene = source.genes[0]
            for key in source.variant_keys:
                if source.locator == f"method:post-ranking-acmg/{gene}/{variant_text(key)}":
                    prepared["assessment_sources"].setdefault((gene, key), []).append(source.source_id)
    if assessment is not None:
        for (gene, key), _ in assessment["final"]:
            prepared["assessment_by_gene"].setdefault(gene, []).append(key)
    base = {
        "fixed_genes": list(context["fixed_genes"]),
        "source_budget": 50 if assessment is not None or context.get("assessment_availability") else 10,
        "available_sources": list(prepared["source_views"].values()),
        "variant_aliases": [
            {"variant_ref": alias, "canonical_variant_key": variant_text(key)}
            for alias, key in variant_by_alias.items()
        ],
    }
    if context.get("assessment_availability") is not None:
        base["assessment_availability"] = context["assessment_availability"]
    assertions = _assessment_payload(assessment, prepared)
    if assertions is not None:
        base["assessment_context"] = assertions
        prepared["assessment_by_identity"] = {
            (row["gene"], row["variant_ref"]): row for row in assertions["classification_assertions"]
        }
    else:
        prepared["assessment_by_identity"] = {}
    pending = [assertions]
    while pending:
        value = pending.pop()
        if isinstance(value, dict):
            if set(value) >= {"source_ref", "quote"}:
                alias, quote = value["source_ref"], value["quote"]
                source = prepared["source_by_alias"][alias]
                _check(quote in source.content, "report_assessment_quote_not_in_source")
                passages = prepared["source_views"][alias]["exact_passages"]
                if quote not in passages:
                    passages.append(quote)
            pending.extend(value.values())
        elif isinstance(value, (list, tuple)):
            pending.extend(value)
    prepared["input"] = base
    return memory, prepared


def _wire_packet(base, views, assertions):
    """Intern literal passages without merging their source or allele bindings."""
    passage_ids, catalog, sources = {}, {}, {}
    for view in views:
        refs = []
        passages = [
            part
            for literal in view["exact_passages"]
            for part in ([literal] if view["kind"] == "patient_observation" else _prose_passages(literal))
        ]
        for passage in passages:
            if passage not in passage_ids:
                ref = f"P{len(passage_ids) + 1:04d}"
                passage_ids[passage] = ref
                catalog[ref] = passage
            refs.append(passage_ids[passage])
        row = sources.setdefault(
            view["source_ref"],
            {
                **{key: value for key, value in view.items() if key != "exact_passages"},
                "passage_refs": [],
            },
        )
        row["passage_refs"] = list(dict.fromkeys((*row["passage_refs"], *refs)))
    packet = {
        key: value for key, value in base.items() if key not in {"available_sources", "assessment_context"}
    }
    packet.update(available_sources=list(sources.values()), passages=catalog)
    if assertions:
        packet["assessment_context"] = {
            **{
                key: value
                for key, value in base["assessment_context"].items()
                if key != "classification_assertions"
            },
            # Final narrative fields are already literal passages of the bound
            # assessment source. Keep the distinct profile facts and limitations.
            "classification_assertions": [
                {
                    key: value
                    for key, value in row.items()
                    if key not in {"rationale", "summary", "patient_fit", "input_source_refs"}
                }
                for row in assertions
            ],
        }
    return packet


def _packet_size(packet, system=""):
    return len(system) + len(json.dumps(packet, ensure_ascii=False, sort_keys=True, separators=(",", ":")))


def _selection_packets(prepared, system):
    base = prepared["input"]
    assertions = base.get("assessment_context", {}).get("classification_assertions", [])
    views = list(prepared["source_views"].values())
    whole = _wire_packet(base, views, assertions)
    if _packet_size(whole, system) <= _PACKET_CHAR_TARGET:
        return [whole]
    packets = []

    for index, gene in enumerate(base["fixed_genes"]):
        units, covered = [], set()
        for decision in (row for row in assertions if row["gene"] == gene):
            # An assessment and its cited literal premises travel together.
            refs, quotes = set(decision["assessment_source_refs"]), {}
            pending = [value for key, value in decision.items() if key != "input_source_refs"]
            while pending:
                value = pending.pop()
                if isinstance(value, dict):
                    if set(value) >= {"source_ref", "quote"}:
                        quotes.setdefault(value["source_ref"], []).append(value["quote"])
                    pending.extend(value.values())
                elif isinstance(value, (list, tuple)):
                    pending.extend(value)
                elif isinstance(value, str) and value in prepared["source_views"]:
                    refs.add(value)
            bound = []
            for alias in sorted(refs):
                view = prepared["source_views"][alias]
                passages = quotes.get(alias, view["exact_passages"])
                bound.append({**view, "exact_passages": passages})
                covered.update((alias, passage) for passage in passages)
            units.append((bound, [decision]))
        for view in views:
            if gene not in view["genes"] and (view["genes"] or index and view["kind"] == "guideline"):
                continue
            remaining = [p for p in view["exact_passages"] if (view["source_ref"], p) not in covered]
            if remaining:
                units.append(([{**view, "exact_passages": remaining}], []))
        rows, decisions = [], []
        while units:
            added_rows, added_decisions = units.pop(0)
            packet = _wire_packet(base, rows + added_rows, decisions + added_decisions)
            if _packet_size(packet, system) <= _INPUT_CHAR_LIMIT - 8192:
                rows.extend(added_rows)
                decisions.extend(added_decisions)
            elif rows:
                packets.append({**_wire_packet(base, rows, decisions), "focus_genes": [gene]})
                rows, decisions = [], []
                units.insert(0, (added_rows, added_decisions))
            elif not added_decisions and len(added_rows[0]["exact_passages"]) > 1:
                passages = added_rows[0]["exact_passages"]
                middle = len(passages) // 2
                units[:0] = [
                    ([{**added_rows[0], "exact_passages": part}], [])
                    for part in (passages[:middle], passages[middle:])
                ]
            else:
                _error("report_atomic_record_budget_exceeded")
        if rows:
            packets.append({**_wire_packet(base, rows, decisions), "focus_genes": [gene]})
    return packets


def _select_claims(model, context, memory, prepared, system):
    packets = _selection_packets(prepared, system)
    selected, seen = [], set()
    for number, packet in enumerate(packets, 1):
        packet = {**packet, "partial_input": len(packets) > 1}
        aliases = {row["source_ref"] for row in packet["available_sources"]}
        scoped = {
            **prepared,
            "source_by_alias": {
                key: value for key, value in prepared["source_by_alias"].items() if key in aliases
            },
            "wire_packet": packet,
        }
        stage = "report_claim_selection" + (f":part{number:03d}" if len(packets) > 1 else "")
        claims = prompting.call_with_repair(
            model,
            context["case_key"],
            stage,
            system,
            packet,
            _MAX_TOKENS,
            lambda raw: _decode_claims(raw, context, memory, scoped, coverage=len(packets) == 1),
            single_retry=True,
            max_input_characters=_INPUT_CHAR_LIMIT,
        )
        for claim in claims:
            identity = json.dumps(
                {key: value for key, value in claim.items() if key != "claim_id"}, sort_keys=True
            )
            if identity not in seen:
                seen.add(identity)
                selected.append({**claim, "claim_id": f"C{len(selected) + 1:03d}"})
    claims = tuple(selected)
    # Do not discard distinct validated support or limitations to meet a
    # readability target. Shards are merged locally, without another LLM call.
    _validate_claims(
        claims, memory, max(prepared["input"]["source_budget"], len(_claim_sources(claims))), prepared
    )
    _claim_coverage(claims, context, prepared)
    return claims


def _compose(model, context, claims, prepared, system):
    base = {
        "fixed_genes": list(context["fixed_genes"]),
        **(
            {"assessment_availability": context["assessment_availability"]}
            if context.get("assessment_availability")
            else {}
        ),
    }
    pending, packets = [list(claims)], []
    while pending:
        group = pending.pop(0)
        payload = {**base, "claims": [_claim_payload(claim, prepared) for claim in group]}
        if _packet_size(payload, system) > _PACKET_CHAR_TARGET and len(group) > 1:
            middle = len(group) // 2
            pending[:0] = [group[:middle], group[middle:]]
        elif _packet_size(payload, system) > _PACKET_CHAR_TARGET and len(group[0]["citations"]) > 1:
            citations = group[0]["citations"]
            middle = len(citations) // 2
            pending[:0] = [
                [{**group[0], "citations": part}] for part in (citations[:middle], citations[middle:])
            ]
        else:
            packets.append((group, payload))
    parts = []
    for number, (group, payload) in enumerate(packets, 1):
        partial = len(packets) > 1
        instruction = (
            "\nThis is a bounded part of the selected evidence. Describe only claims supplied here. "
            "Keep the full fixed gene order; include sections only for genes supported in this part. "
            "Do not treat unsupplied evidence as absent. Only compare genes when both have claims here. "
            "Sections will be joined in the original frozen order without reranking."
            if partial
            else ""
        )
        parts.append(
            prompting.call_with_repair(
                model,
                context["case_key"],
                "report_composition" + (f":part{number:03d}" if partial else ""),
                system + instruction,
                payload,
                _MAX_TOKENS,
                lambda raw: _decode_report(raw, context, tuple(group), prepared),
                single_retry=True,
                max_input_characters=_INPUT_CHAR_LIMIT,
            )
        )
    if len(parts) == 1:
        return parts[0]
    return {
        "case_key": context["case_key"],
        "fixed_genes": context["fixed_genes"],
        "claims": claims,
        "summary": tuple(row for part in parts for row in part["summary"]),
        "gene_sections": tuple(
            {
                "gene": gene,
                "paragraphs": tuple(
                    row
                    for part in parts
                    for section in part["gene_sections"]
                    if section["gene"] == gene
                    for row in section["paragraphs"]
                ),
            }
            for gene in context["fixed_genes"]
        ),
        "comparisons": tuple(row for part in parts for row in part["comparisons"]),
    }


def _assessment_payload(assessment, prepared):
    if assessment is None:
        return None
    if assessment.get("workflow") == "clinical_profile_v1":
        return _clinical_profile_payload(assessment, prepared)
    source_aliases = prepared["source_alias_by_id"]
    variants = prepared["variant_alias_by_key"]
    assertions = []
    for identity, final in assessment["final"]:
        review, partner = assessment["independent"][identity], assessment["partners"][identity]
        criteria = {item["code"]: item for item in _sequence(review["criteria"], "assessment.criteria")}
        criteria[partner["criterion"]["code"]] = partner["criterion"]
        decisive = tuple(criteria[code] for code in final["decisive_codes"])
        required = dict.fromkeys(
            (
                *final["decisive_source_ids"],
                *final["contextual_source_ids"],
                *(citation["source_id"] for criterion in decisive for citation in criterion["citations"]),
            )
        )
        _check(
            all(source_id in source_aliases for source_id in required), "report_assessment_source_unavailable"
        )
        assessment_sources = [source_aliases[key] for key in prepared["assessment_sources"].get(identity, ())]
        _check(bool(assessment_sources), "report_assessment_source_unavailable")
        gene, key = identity
        assertions.append(
            {
                "variant_ref": variants[key],
                "gene": gene,
                "disease": dict(final["disease"]),
                "classification": final["classification"],
                "assessment_source_refs": assessment_sources,
                "decisive_criteria": [
                    {
                        "code": code,
                        **_fields(criteria[code], "status", "strength", "rationale"),
                        "source_refs": [
                            source_aliases[citation["source_id"]] for citation in criteria[code]["citations"]
                        ],
                    }
                    for code in final["decisive_codes"]
                ],
                "limitation_codes": list(final["limitation_codes"]),
                "decisive_source_refs": [
                    source_aliases[source_id] for source_id in final["decisive_source_ids"]
                ],
                "contextual_source_refs": [
                    source_aliases[source_id] for source_id in final["contextual_source_ids"]
                ],
                **_fields(final, "rationale", "patient_fit"),
            }
        )
    return {
        "instruction": (
            "These are immutable post-ranking assessment assertions. Explain them only when "
            "supported by selected claims; do not generate, revise, or recompute a classification."
        ),
        "classification_assertions": assertions,
    }


def _profile_citations(citations, source_aliases):
    output = []
    for citation in citations:
        source_id = citation["source_id"]
        _check(source_id in source_aliases, "report_assessment_source_unavailable")
        output.append({"source_ref": source_aliases[source_id], "quote": citation["quote"]})
    return output


def _profile_evidence(item, source_aliases):
    return {
        key: (_profile_citations(value, source_aliases) if key == "citations" else value)
        for key, value in item.items()
    }


def _clinical_profile_payload(assessment, prepared):
    variants = prepared["variant_alias_by_key"]
    assertions = []
    for identity, final in assessment["final"]:
        source_aliases = dict(prepared["source_alias_by_id"])
        for alias, view in prepared["source_views"].items():
            if (
                view.get("same_publication_as")
                and view["genes"] == [identity[0]]
                and view["variant_refs"] == [variants[identity[1]]]
            ):
                original = prepared["source_by_alias"][view["same_publication_as"]]
                source_aliases[original.source_id] = alias
        profile = assessment["profiles"][identity]
        support = {item["id"]: item for item in profile["key_support"]}
        limits = {item["id"]: item for item in profile["key_limits"]}
        _check(
            all(item_id in support for item_id in final["key_evidence_ids"]),
            "report_assessment_evidence_unavailable",
        )
        _check(
            all(item_id in limits for item_id in final["limitation_ids"]),
            "report_assessment_limitation_unavailable",
        )
        assessment_sources = [source_aliases[key] for key in prepared["assessment_sources"].get(identity, ())]
        _check(bool(assessment_sources), "report_assessment_source_unavailable")
        for source_id in (*final["decisive_source_ids"], *final["contextual_source_ids"]):
            _check(source_id in source_aliases, "report_assessment_source_unavailable")
        allele_disease = profile["allele_disease"]
        patient_fit = profile["patient_fit"]
        anchor = profile["validated_pvs1_anchor"]
        for source_id in profile["input_source_ids"]:
            _check(source_id in source_aliases, "report_assessment_source_unavailable")
        assertions.append(
            {
                "variant_ref": variants[identity[1]],
                "gene": identity[0],
                "disease": dict(final["disease"]),
                "classification": final["classification"],
                "confidence": final["confidence"],
                "allele_disease": {
                    **_fields(allele_disease, "label", "identifier", "inheritance", "reason"),
                    "citations": _profile_citations(allele_disease["citations"], source_aliases),
                },
                "profile_patient_fit": {
                    **_fields(patient_fit, "status", "reason"),
                    "citations": _profile_citations(patient_fit["citations"], source_aliases),
                },
                "decisive_evidence": [
                    _profile_evidence(support[item_id], source_aliases)
                    for item_id in final["key_evidence_ids"]
                ],
                "limiting_evidence": [
                    _profile_evidence(limits[item_id], source_aliases) for item_id in final["limitation_ids"]
                ],
                "unknowns": list(profile["unknowns"]),
                "profile_summary": profile["summary"],
                **(
                    {
                        "retained_observations": [
                            _profile_evidence(row, source_aliases) for row in profile["retained_observations"]
                        ],
                        "observation_review": final["observation_review"],
                    }
                    if profile.get("retained_observations")
                    else {}
                ),
                "validated_pvs1_anchor": (
                    _profile_evidence(anchor, source_aliases) if anchor is not None else None
                ),
                "input_source_refs": [source_aliases[source_id] for source_id in profile["input_source_ids"]],
                "assessment_source_refs": assessment_sources,
                "decisive_source_refs": [
                    source_aliases[source_id] for source_id in final["decisive_source_ids"]
                ],
                "contextual_source_refs": [
                    source_aliases[source_id] for source_id in final["contextual_source_ids"]
                ],
                **_fields(final, "rationale", "patient_fit", "summary"),
            }
        )
    return {
        "workflow": "clinical_profile_v1",
        "instruction": (
            "These are immutable post-ranking ungraded clinical evidence profiles and final "
            "five-tier assertions. Explain the supplied profile and decision only; do not "
            "generate independent classifications, partner grades, or criterion rows."
        ),
        "classification_assertions": assertions,
    }


def _assessment_for_claims(prepared, claims):
    rows = prepared["assessment_by_identity"]
    identities = dict.fromkeys(
        (gene, prepared["variant_alias_by_key"][key])
        for claim in claims
        for gene in claim["genes"]
        for key in claim["variant_keys"]
        if key in prepared["variant_alias_by_key"]
    )
    bound = {
        prepared["source_alias_by_id"][citation["source_id"]]
        for claim in claims
        for citation in claim["citations"]
        if citation["source_id"] in prepared["source_alias_by_id"]
    }
    selected = []
    for identity in identities:
        if identity not in rows:
            continue
        row = rows[identity]
        assessment_refs = [ref for ref in row["assessment_source_refs"] if ref in bound]
        if not assessment_refs:
            continue
        selected.append(
            {
                **_fields(row, "gene", "variant_ref", "disease", "classification"),
                **({"confidence": row["confidence"]} if "confidence" in row else {}),
                "assessment_source_refs": assessment_refs,
            }
        )
    if not selected:
        return None
    original = prepared["input"]["assessment_context"]
    return {key: value for key, value in original.items() if key != "classification_assertions"} | {
        "classification_assertions": selected
    }


def _claim_payload(claim, prepared):
    payload = {
        **_fields(claim, "claim_id", "scope"),
        "genes": list(claim["genes"]),
        "variant_refs": [prepared["variant_alias_by_key"][key] for key in claim["variant_keys"]],
        **_fields(claim, "polarity", "text"),
        "citations": [
            {"source_ref": prepared["source_alias_by_id"][row["source_id"]], "quote": row["quote"]}
            for row in claim["citations"]
        ],
    }
    assessment = _assessment_for_claims(prepared, (claim,))
    if assessment is not None:
        payload["relevant_assessment"] = assessment
    return payload


def _quote_pattern(text):
    parts = re.split(r"\s+", text.strip())
    pattern = re.escape(parts[0])
    for previous, current in zip(parts, parts[1:]):
        separator = (
            r"\s*" if previous[-1:] in "=:,;[{" or current[:1] in "=:,;]}" else r"(?:\s|\[PMID:\d+\])+"
        )
        pattern += separator + re.escape(current)
    return pattern


def recover_quote(source, quote):
    if len(quote) > 600 or not quote.strip():
        _error("quote_not_in_original_source")
    if quote in source:
        return quote
    try:
        document = json.loads(source)
    except (json.JSONDecodeError, RecursionError):
        document = None
    encoded = json.dumps(quote, ensure_ascii=False)[1:-1]
    if isinstance(document, (dict, list)) and encoded != quote and encoded in source:
        if source.count(encoded) != 1 or len(encoded) > 600:
            _error("quote_not_in_original_source")
        return encoded
    assigned = re.fullmatch(r"([A-Za-z_][A-Za-z_0-9]*)=([^\n]+)", quote.strip())
    if assigned and isinstance(document, (dict, list)):
        literal = (
            re.escape(json.dumps(assigned[1]))
            + r"\s*:\s*"
            + re.escape(json.dumps(assigned[2], ensure_ascii=False))
        )
        fields = list(re.finditer(literal, source))
        if len(fields) == 1:
            return source[fields[0].start() : fields[0].end()]
    matches = list(re.finditer(_quote_pattern(quote), source, flags=re.IGNORECASE))
    if not matches and quote[-1:] in {".", "。", ";", "；"}:
        matches = list(
            re.finditer(_quote_pattern(quote.strip()[:-1]) + r"[.;；。](?!\w)", source, flags=re.IGNORECASE)
        )
    if len(matches) != 1:
        _error("quote_not_in_original_source")
    original = source[matches[0].start() : matches[0].end()]
    if len(original) > 600:
        _error("quote_not_in_original_source")
    return original


def _expand_aliases(text, prepared, allowed, unknown, unbound):
    def expand(match):
        key = prepared["variant_by_alias"].get(match.group())
        _check(key is not None, unknown)
        _check(key in allowed, unbound)
        return variant_text(key)

    return _VARIANT_ALIAS.sub(expand, text)


def _medical_source(source, claim):
    if source.kind not in _MEDICAL:
        return False
    if claim["scope"] == "gene":
        return claim["genes"][0] in source.genes
    if claim["scope"] == "allele":
        return claim["genes"][0] in source.genes and bool(
            set(claim["variant_keys"]) & set(source.variant_keys)
        )
    if claim["scope"] == "comparison":
        return bool(set(claim["genes"]) & set(source.genes))
    return False


def _validate_claims(claims, memory, max_sources, prepared=None):
    require(type(claims) is tuple and bool(claims), "report.claims", "nonempty_tuple_required")
    for claim in claims:
        try:
            memory.validate_claim(claim)
        except MendelPrimeError as error:
            if prepared is None:
                raise
            refs = tuple(
                prepared["source_alias_by_id"].get(citation["source_id"], "unknown")
                for citation in claim["citations"]
            )
            variants = tuple(
                prepared["variant_alias_by_key"].get(key, "unknown") for key in claim["variant_keys"]
            )
            locator = f"claim={claim['claim_id']};sources={','.join(refs) or 'none'}"
            if variants:
                locator += f";variants={','.join(variants)}"
            raise MendelPrimeError(f"{error.code};{locator}", category="response") from error
        require(
            not _CJK.search(claim["text"]),
            f"claims[{claim['claim_id']}].text",
            "cjk_generated_prose_forbidden",
        )
        if claim["scope"] != "case" and claim["polarity"] == "support":
            sources = tuple(memory.get_source(citation["source_id"]) for citation in claim["citations"])
            if claim["scope"] == "allele":
                valid = all(
                    any(_medical_source(source, {**claim, "variant_keys": (key,)}) for source in sources)
                    for key in claim["variant_keys"]
                )
            elif claim["scope"] == "comparison":
                valid = all(
                    any(gene in source.genes and _medical_source(source, claim) for source in sources)
                    for gene in claim["genes"]
                )
            else:
                valid = any(_medical_source(source, claim) for source in sources)
            require(valid, f"claims[{claim['claim_id']}]", "causal_support_requires_medical_source")
    count = len({citation["source_id"] for claim in claims for citation in claim["citations"]})
    require(count <= max_sources, "report.claims", "source_budget_exceeded")


def _decode_claims(raw, context, memory, prepared, *, coverage=True):
    _check(set(raw) == {"claims"}, "report_claim_selection_schema_invalid")
    claims, binding_errors = [], []
    for number, item in enumerate(_rows(raw["claims"], "report_claim_selection_empty", nonempty=coverage), 1):
        row = _object(item, "report_claim_invalid")
        _check(
            set(row) == {"scope", "genes", "variant_refs", "polarity", "text", "citations"},
            "report_claim_schema_invalid",
        )
        _check(row["scope"] in _SCOPES and row["polarity"] in _POLARITIES, "report_claim_kind_invalid")
        genes = _strings(row["genes"], "report_claim_genes_invalid", nonempty=False)
        variant_refs = _strings(row["variant_refs"], "report_claim_variants_invalid", nonempty=False)
        try:
            variants = tuple(prepared["variant_by_alias"][alias] for alias in variant_refs)
        except KeyError:
            _error("report_claim_variant_unknown")
        text = _expand_aliases(
            _text(row["text"], "report_claim_text_invalid"),
            prepared,
            set(variants),
            "report_claim_variant_alias_unknown",
            "report_claim_variant_alias_unbound",
        )
        citations, declared_sources, claim_errors = [], [], []
        for citation_number, item in enumerate(_rows(row["citations"], "report_claim_citations_missing"), 1):
            citation = _object(item, "report_claim_citation_invalid")
            _check(
                set(citation) in ({"source_ref", "quote"}, {"source_ref", "passage_ref"}),
                "report_claim_citation_invalid",
            )
            ref = _text(citation["source_ref"], "report_claim_source_unknown")
            source = prepared["source_by_alias"].get(ref)
            _check(source is not None, "report_claim_source_unknown")
            declared_sources.append(source)
            if "passage_ref" in citation:
                passage_ref = _text(citation["passage_ref"], "report_claim_quote_invalid")
                packet = prepared.get("wire_packet", {})
                allowed = next(
                    (
                        row["passage_refs"]
                        for row in packet.get("available_sources", ())
                        if row["source_ref"] == ref
                    ),
                    (),
                )
                if passage_ref not in allowed:
                    claim_errors.append(
                        MendelPrimeError(
                            "report_passage_source_mismatch",
                            category="response",
                            field=f"claims[{number}].citations[{citation_number}]",
                            source={
                                "claim": row,
                                "invalid_citation": citation,
                                "passage_owner_refs": [
                                    entry["source_ref"]
                                    for entry in packet.get("available_sources", ())
                                    if passage_ref in entry["passage_refs"]
                                ],
                                "instruction": "Choose a source/passage pair valid for this claim. Passage ownership does not establish gene or allele eligibility.",
                            },
                        )
                    )
                    continue
                offered = packet["passages"][passage_ref]
            else:
                offered = _text(citation["quote"], "report_claim_quote_invalid")
            if source.kind != "patient_observation" and len(offered) > 600:
                _error(f"report_quote_length_exceeded_{ref}_limit_600")
            quote = (
                offered
                if source.kind == "patient_observation" and offered in source.content
                else recover_quote(source.content, offered)
            )
            if prepared.get("wire_packet", {}).get("partial_input"):
                packet = prepared["wire_packet"]
                refs = next(
                    (row["passage_refs"] for row in packet["available_sources"] if row["source_ref"] == ref),
                    (),
                )
                _check(any(quote in packet["passages"][key] for key in refs), "report_quote_not_in_packet")
            citations.append({"source_id": source.source_id, "quote": quote})
        support = (
            declared_sources
            if row["polarity"] == "support"
            else memory.closure(tuple(s.source_id for s in declared_sources))
        )
        anchored = bound_alleles(support, genes)
        missing = [alias for alias, key in zip(variant_refs, variants) if key not in anchored]
        if missing:
            offered = {
                entry["source_ref"]: entry["passage_refs"]
                for entry in prepared.get("wire_packet", {}).get("available_sources", ())
            }
            options = {
                alias: [
                    {
                        "source_ref": name,
                        "kind": source.kind,
                        "provider": source.provider,
                        "passage_refs": offered[name],
                    }
                    for name, source in prepared["source_by_alias"].items()
                    if name in offered
                    and prepared["variant_by_alias"][alias] in bound_alleles((source,), genes)
                    and (row["polarity"] != "support" or source.kind in _MEDICAL)
                ]
                for alias in missing
            }
            if row["polarity"] != "support":
                for alias, candidates in options.items():
                    identity = [
                        item
                        for item in candidates
                        if item["kind"] in {"patient_observation", "tool_observation", "annotation"}
                    ]
                    options[alias] = identity or candidates
            claim_errors.append(
                MendelPrimeError(
                    f"report_claim_selection_invalid:allele_source_required;claim=C{number:03d}",
                    category="response",
                    field=f"claims[{number}]",
                    source={
                        "claim": row,
                        "missing_variant_refs": missing,
                        "identity_source_options": options,
                        "instruction": "Repair only the identified binding; preserve every other already-valid source/passage pair. Choose a listed passage of its own identity source that actually establishes the claimed identity, together with the paper. Patient/native identity observations are not pathogenicity evidence. Causal support still requires exact-allele medical evidence. If identity cannot be established, retain only genuinely supported gene-level context; never add an unrelated citation to satisfy the check.",
                    },
                )
            )
        if claim_errors:
            binding_errors.extend(claim_errors)
            continue
        claim = {
            "case_key": context["case_key"],
            "claim_id": f"C{number:03d}",
            "scope": row["scope"],
            "genes": genes,
            "text": text,
            "polarity": row["polarity"],
            "citations": tuple(citations),
            "variant_keys": variants,
        }
        try:
            validate_claim_shape(claim, context["case_key"])
        except MendelPrimeError as error:
            raise MendelPrimeError("report_claim_scope_invalid", category="response") from error
        claims.append(claim)
    if binding_errors:
        first = binding_errors[0]
        if len(binding_errors) > 1:
            first.source = {
                **first.source,
                "additional_binding_errors": [
                    {"error": error.code, "location": error.field, **error.source}
                    for error in binding_errors[1:]
                ],
            }
        raise first
    selected = tuple(claims)
    if not selected and not coverage:
        return ()
    _check(
        all(
            set(claim["genes"])
            <= set(prepared.get("wire_packet", {}).get("focus_genes", context["fixed_genes"]))
            for claim in selected
        ),
        "report_claim_gene_outside_frozen_ranking",
    )
    try:
        _validate_claims(selected, memory, prepared["input"]["source_budget"], prepared)
    except MendelPrimeError as error:
        raise MendelPrimeError(f"report_claim_selection_invalid:{error.code}", category="response") from error
    if coverage:
        _claim_coverage(selected, context, prepared)
    return selected


def _claim_coverage(selected, context, prepared):
    covered = {
        claim["genes"][0]
        for claim in selected
        if claim["scope"] in {"gene", "allele"} and len(claim["genes"]) == 1
    }
    _check(set(context["fixed_genes"]) <= covered, "report_gene_premise_incomplete")
    for gene in context["fixed_genes"]:
        keys = prepared["assessment_by_gene"].get(gene, ())
        if not keys:
            continue
        valid = any(
            claim["scope"] == "allele"
            and claim["genes"] == (gene,)
            and any(
                key in claim["variant_keys"]
                and {citation["source_id"] for citation in claim["citations"]}.intersection(
                    prepared["assessment_sources"].get((gene, key), ())
                )
                for key in keys
            )
            for claim in selected
        )
        _check(valid, "report_assessment_claim_coverage_incomplete")


def _paragraph(raw, prepared, claims, *, summary=False):
    row = _object(raw, "report_paragraph_invalid")
    _check(set(row) == {"text", "claim_ids", "role"}, "report_paragraph_schema_invalid")
    role = row["role"]
    _check(role in ({"summary"} if summary else _ROLES - {"summary"}), "report_paragraph_role_invalid")
    claim_ids = _strings(row["claim_ids"], "report_paragraph_claims_invalid")
    _check(all(key in claims for key in claim_ids), "report_paragraph_claim_unknown")
    allowed = {
        key for claim_id in claim_ids if claim_id in claims for key in claims[claim_id]["variant_keys"]
    }
    return {
        "text": _expand_aliases(
            _text(row["text"], "report_paragraph_text_invalid"),
            prepared,
            allowed,
            "report_paragraph_variant_alias_unknown",
            "report_paragraph_variant_alias_unbound",
        ),
        "claim_ids": claim_ids,
        "role": role,
    }


def _decode_report(raw, context, claims, prepared):
    _check(set(raw) == {"summary", "gene_sections", "comparisons"}, "report_composition_schema_invalid")
    claim_map = {claim["claim_id"]: claim for claim in claims}
    summary = tuple(
        _paragraph(item, prepared, claim_map, summary=True)
        for item in _rows(raw["summary"], "report_summary_missing")
    )
    sections = []
    for item in _rows(raw["gene_sections"], "report_gene_sections_missing"):
        row = _object(item, "report_gene_section_invalid")
        _check(set(row) == {"gene", "paragraphs"}, "report_gene_section_schema_invalid")
        gene = _text(row["gene"], "report_gene_invalid")
        try:
            validate_gene(gene)
        except MendelPrimeError as error:
            raise MendelPrimeError("report_gene_section_invalid", category="response") from error
        sections.append(
            {
                "gene": gene,
                "paragraphs": tuple(
                    _paragraph(item, prepared, claim_map)
                    for item in _rows(row["paragraphs"], "report_gene_paragraphs_missing")
                ),
            }
        )
    comparisons = []
    for item in _rows(raw["comparisons"], "report_comparisons_invalid", nonempty=False):
        row = _object(item, "report_comparison_invalid")
        _check(
            set(row) == {"first_gene", "second_gene", "text", "evidence_claim_ids", "favored_gene"},
            "report_comparison_schema_invalid",
        )
        favored = row["favored_gene"]
        _check(favored is None or isinstance(favored, str), "report_comparison_favored_gene_invalid")
        claim_ids = _strings(row["evidence_claim_ids"], "report_comparison_claims_invalid")
        _check(all(key in claim_map for key in claim_ids), "report_comparison_claim_unknown")
        allowed = {
            key
            for claim_id in claim_ids
            if claim_id in claim_map
            for key in claim_map[claim_id]["variant_keys"]
        }
        comparisons.append(
            {
                "first_gene": _text(row["first_gene"], "report_comparison_gene_invalid"),
                "second_gene": _text(row["second_gene"], "report_comparison_gene_invalid"),
                "text": _expand_aliases(
                    _text(row["text"], "report_comparison_text_invalid"),
                    prepared,
                    allowed,
                    "report_comparison_variant_alias_unknown",
                    "report_comparison_variant_alias_unbound",
                ),
                "evidence_claim_ids": claim_ids,
                "favored_gene": favored,
            }
        )
    report = {
        "case_key": context["case_key"],
        "fixed_genes": context["fixed_genes"],
        "claims": claims,
        "summary": summary,
        "gene_sections": tuple(sections),
        "comparisons": tuple(comparisons),
    }
    return report


def _claim_sources(claims):
    return tuple(dict.fromkeys(citation["source_id"] for claim in claims for citation in claim["citations"]))


def _fields(record, *names):
    return {name: record[name] for name in names}


def _escape(value):
    escaped = html.escape(value, quote=False).replace("\\", "\\\\")
    for character in "`*_{}[]#+!|":
        escaped = escaped.replace(character, f"\\{character}")
    escaped = _ORDERED_LIST.sub(r"\1\\.", escaped)
    escaped = _LEADING_LIST.sub(r"\1\\\2", escaped)
    return _HORIZONTAL_RULE.sub(r"\1\\\2\3", escaped)


def _version_label(version):
    compact = _LONG_HEX.sub(lambda match: f"{match.group()[:12]}...", version)
    label = compact if len(compact) <= 160 else compact[:157].rstrip() + "..."
    require(
        not (
            PureWindowsPath(label).is_absolute()
            or PurePosixPath(label).is_absolute()
            or label.startswith(("\\\\", "//"))
            or _ABSOLUTE_PATH.search(label)
        ),
        "source.version",
        "absolute_resource_path_forbidden",
    )
    return label


def render_markdown(report, sources, *, availability=None, assessment=None):
    """Render parsed report content without post-generation review."""
    source_by_id = {source.source_id: source for source in sources}
    claims = {claim["claim_id"]: claim for claim in report["claims"]}
    groups = (
        *(row["claim_ids"] for row in report["summary"]),
        *(row["claim_ids"] for section in report["gene_sections"] for row in section["paragraphs"]),
        *(row["evidence_claim_ids"] for row in report["comparisons"]),
    )
    used_claims = tuple(dict.fromkeys(claim_id for ids in groups for claim_id in ids))
    source_order = _claim_sources(claims[claim_id] for claim_id in used_claims)
    numbers = {source_id: number for number, source_id in enumerate(source_order, 1)}

    def paragraph(text, ids):
        source_ids = _claim_sources(claims[claim_id] for claim_id in ids)
        return f"{_escape(text)} {' '.join(f'[{numbers[source_id]}]' for source_id in source_ids)}"

    lines = ["# Traceable Gene Ranking Explanation", ""]
    if availability:
        genes = ", ".join(dict.fromkeys(row["gene"] for row in availability["unclassified_alleles"]))
        lines.extend(
            [
                f"> ACMG assessment is incomplete: {availability['classified_alleles']} of "
                f"{availability['expected_alleles']} candidate alleles have a completed classification. "
                f"Unclassified alleles occur in {_escape(genes)}.",
                "",
            ]
        )
    lines.extend(["## Summary", ""])
    for row in report["summary"]:
        lines.extend([paragraph(row["text"], row["claim_ids"]), ""])
    lines.extend(["## Candidate Overview", ""])
    if assessment is not None or availability:
        labels = {
            "pathogenic": "Pathogenic",
            "likely_pathogenic": "Likely pathogenic",
            "uncertain_significance": "VUS",
            "likely_benign": "Likely benign",
            "benign": "Benign",
        }
        entries = [
            (*_identity(row), labels[row["classification"]]) for row in (assessment or {}).get("final", ())
        ]
        entries.extend(
            (*_identity({"identity": row}), "Not assessed")
            for row in (availability or {}).get("unclassified_alleles", ())
        )
        lines.extend(["| Rank | Gene | Assessed allele | Classification |", "|---:|---|---|---|"])
        for rank, gene in enumerate(report["fixed_genes"], 1):
            alleles = [(variant_text(key), label) for name, key, label in entries if name == gene]
            for allele, label in alleles or [("Not available", "Not assessed")]:
                lines.append(f"| {rank} | {gene} | {_escape(allele)} | {label} |")
        lines.append("")
    else:
        lines.extend(["| Rank | Gene |", "|---:|---|"])
        lines.extend(f"| {rank} | {gene} |" for rank, gene in enumerate(report["fixed_genes"], 1))
        lines.append("")
    lines.extend(["## Ranked Genes", ""])
    for rank, section in enumerate(report["gene_sections"], 1):
        lines.extend([f"### {rank}. {section['gene']}", ""])
        for row in section["paragraphs"]:
            lines.extend([paragraph(row["text"], row["claim_ids"]), ""])
    if report["comparisons"]:
        lines.extend(["## Comparisons With Rank 1", ""])
        for number, row in enumerate(report["comparisons"], 1):
            lines.extend(
                [
                    f"{number}. **{_escape(row['first_gene'])} and {_escape(row['second_gene'])}.** {paragraph(row['text'], row['evidence_claim_ids'])}",
                    "",
                ]
            )
    quotes = {source_id: [] for source_id in source_order}
    for claim in (claims[claim_id] for claim_id in used_claims):
        for citation in claim["citations"]:
            if citation["source_id"] in quotes and citation["quote"] not in quotes[citation["source_id"]]:
                quotes[citation["source_id"]].append(citation["quote"])
    lines.extend(["## References", ""])
    for source_id in source_order:
        source = source_by_id[source_id]
        metadata = f"{_escape(source.provider)}; {_escape(source.kind)}; locator: {_escape(source.locator)}; version label: {_escape(_version_label(source.version))}"
        lines.append(f"{numbers[source_id]}. {metadata}")
        lines.extend(["", "<details><summary>Original source excerpts</summary>", ""])
        for quote in quotes[source_id]:
            escaped = _escape(quote).splitlines() or [""]
            if len(escaped) == 1:
                lines.append(f'   > "{escaped[0]}"')
            else:
                lines.append(f'   > "{escaped[0]}')
                lines.extend(f"   > {line}" for line in escaped[1:-1])
                lines.append(f'   > {escaped[-1]}"')
        lines.extend(["", "</details>", ""])
    return "\n".join(lines).rstrip() + "\n"


def run(model, context):
    """Generate a source-bound report; no post-generation audit is run."""
    sources = context.get("sources", ()) if type(context) is dict else ()
    report, markdown = None, ""
    try:
        memory, prepared = _prepare(context)
        partial_instruction = (
            "\n" + prompting.text("report_assessment_availability")
            if context.get("assessment_availability")
            else ""
        )
        claims = _select_claims(
            model, context, memory, prepared, prompting.text("report_claim_selection") + partial_instruction
        )
        parsed = _compose(
            model, context, claims, prepared, prompting.text("report_composition") + partial_instruction
        )
        markdown = render_markdown(
            parsed,
            sources,
            availability=context.get("assessment_availability"),
            assessment=context["assessment"],
        )
        report = parsed
        status, issues = "ok", ()
    except (MendelPrimeError, OSError) as error:
        issue = "explanation_resource_io_error" if isinstance(error, OSError) else error.code
        status, issues = "degraded", (issue,)
    return {
        "status": status,
        "report": report,
        "review_status": "not_performed",
        "sources": sources,
        "markdown": markdown,
        "issues": issues,
    }


__all__ = ["recover_quote", "render_markdown", "run"]
