"""Ungraded evidence profiles, profile-only classification, and saved-ledger validation."""

import hashlib
import json
import math
import re
from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path

from mendelprime import prompting
from mendelprime.domain import (
    MendelPrimeError,
    Source,
    case_tuple,
    is_evidence_limitation,
    variant_text,
    variant_tuple,
)
from mendelprime.tools.native import variant_is_candidate

try:
    import hgvs.parser
except ImportError:  # pragma: no cover - deployment dependency is checked by doctor
    hgvs = None


_CLASSES = {
    "pathogenic",
    "likely_pathogenic",
    "uncertain_significance",
    "likely_benign",
    "benign",
}
_FIT_STATUSES = {"compatible", "uncertain", "not_supported"}
_MEASUREMENT_ROLES = {"patient_observation", "observation_qc"}
_POSITIVE_CLASSES = {"pathogenic", "likely_pathogenic"}
_PROFILE_INHERITANCE = {"AR", "AD", "XL", "unknown"}
_PROFILE_SUPPORT_KINDS = {
    "pvs1",
    "family",
    "population",
    "functional",
    "segregation",
    "case_series",
    "molecular",
    "other",
}
_PROFILE_RELIABILITY = {"high", "moderate", "limited"}
_PROFILE_RELEVANCE = {"high", "moderate", "supporting", "uncertain"}
_PVS1_STRENGTHS = {"very_strong", "strong", "moderate", "supporting"}
_OBSERVATION_LABELS = {
    "revel": "Original missense-effect prediction; no criterion or classification assigned.",
    "spliceai": "Original splice-effect prediction; not a test of missense protein function.",
    "publication-family-ledger": "Original within-study family observations; not extra cases or assigned ACMG strength.",
    "clingen": "Original disease-bound gene mechanism context; not proof of this allele's effect or an independent criterion.",
    "genereviews": "Original disease-bound gene mechanism context; not proof of this allele's effect or an independent criterion.",
    "autopvs1": "Native transcript/mechanism context; not an independent assay, validated PVS1 or additional criterion.",
}
_BACKGROUND_PAPERS = 2
_RNA_TERMS = re.compile(r"\bRNA\b|splic(?:e|ing)|transcript", re.IGNORECASE)
_CLINICAL_TITLE = re.compile(
    r"syndrome|encephalopath|neurodevelopment|germline|Mendelian|patients|"
    r"phenotyp|case report|mutation spectrum",
    re.IGNORECASE,
)
_MECHANISM_TEXT = re.compile(
    r"loss.of.function|haploinsufficien|truncat|nonsense|frameshift",
    re.IGNORECASE,
)


def _fail(code):
    raise MendelPrimeError(code, category="response")


def _text(value, code, *, allow_empty=False):
    if not isinstance(value, str) or (not allow_empty and not value.strip()):
        _fail(code)
    return value


def _list(value, code):
    if not isinstance(value, list):
        _fail(code)
    return value


def _unique_strings(value, code):
    rows = _list(value, code)
    if any(not isinstance(item, str) or not item for item in rows) or len(rows) != len(set(rows)):
        _fail(code)
    return tuple(rows)


def _canonical_json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _key(value):
    return variant_tuple(value, MendelPrimeError("variant_key_invalid", category="response"))


def _key_text(value):
    return variant_text(_key(value))


def _identity(target):
    return target["primary"]["gene"], _key(target["primary"]["key"])


def _criterion(code, status, strength, rationale, citations=(), independence_group=None, origin=""):
    return {
        "code": code,
        "status": status,
        "strength": strength,
        "rationale": rationale,
        "citations": tuple(citations),
        "independence_group": independence_group or code,
        "origin": origin,
    }


def _citation(source_id, quote):
    return {"source_id": source_id, "quote": quote}


@lru_cache(maxsize=512)
def _protein_forms(value):
    try:
        parsed = hgvs.parser.Parser().parse_hgvs_variant("NP_000001.1:" + value)
        return tuple(parsed.format({"p_3_letter": flag}).split(":", 1)[1] for flag in (True, False))
    except Exception:
        return (value,)


def _allele_tokens(target):
    values = [target["reference"].get("normalized_hgvs")]
    for variant in (target["primary"], *target.get("aliases", ())):
        values.extend((variant.get("hgvs_c"), variant.get("hgvs_p")))
    output = set()
    for value in values:
        token = str(value or "").split(":")[-1].strip()
        if token in {"", "p.?", "p.(?)", "p.=", "p.(=)"}:
            continue
        output.add(token)
        if token.startswith("p."):
            for form in _protein_forms(token):
                output.update((form, form.removeprefix("p.").strip("()")))
    return tuple(sorted(alias for alias in output if len(alias) >= 5))


def _contains(text, value):
    pattern = re.escape(value).replace(">", r"\s*>\s*")
    return bool(re.search(r"(?<![A-Za-z0-9])" + pattern + r"(?![A-Za-z0-9])", text, re.I))


def _mentions_allele(target, text):
    return any(_contains(text, alias) for alias in _allele_tokens(target))


def _publication_title(source):
    try:
        value = json.loads(source.content)
    except (json.JSONDecodeError, RecursionError):
        value = None
    if isinstance(value, dict):
        for key in ("publication_title", "title", "article_title"):
            title = value.get(key)
            if isinstance(title, str) and title.strip():
                return title.strip()[:300]
    return source.content.split("..", 1)[0].split("\n", 1)[0].split(". ", 1)[0][:300]


def _publication_key(source):
    locator = source.locator.split("#", 1)[0].casefold().rstrip("/")
    if source.kind == "literature" and (
        source.provider in {"pubmed", "pmc"} or "pubmed" in locator or "pmc" in locator
    ):
        return locator
    title = _publication_title(source).casefold().strip(" .")
    return title or locator or source.source_id


def _disease_terms(address):
    values = list(address["target"].get("disease_labels", ()))
    for item in address["sources"]:
        if item["binding"]["role"] != "gene_disease_context":
            continue
        try:
            raw = json.loads(item["source"].content)
        except (json.JSONDecodeError, RecursionError):
            raw = None
        if isinstance(raw, dict):
            values.extend(raw.get(key, "") for key in ("disease_label", "hi_disease_label", "disease_name"))
    return tuple(sorted({str(value).strip() for value in values if len(str(value).strip()) > 5}))


def _relevant_background(address, item):
    source = item["source"]
    text = source.content + " " + _publication_title(source)
    return (
        item["binding"]["match_scope"] == "context"
        or _contains(text, _identity(address["target"])[0])
        or any(_contains(text, label) for label in _disease_terms(address))
    )


def _background_priority(address, item):
    source = item["source"]
    title, text = _publication_title(source), source.content
    traceable = source.kind == "literature" and (
        source.provider in {"pubmed", "pmc"}
        or "pubmed" in source.locator.casefold()
        or "pmc" in source.locator.casefold()
    )
    return (
        not any(_contains(title, label) for label in _disease_terms(address)),
        not bool(_CLINICAL_TITLE.search(title)),
        not traceable,
        not _contains(title, _identity(address["target"])[0]),
        not bool(_MECHANISM_TEXT.search(text)),
        _publication_key(source),
        source.source_id,
    )


def _selected_sources(address):
    anchors = {
        citation["source_id"]
        for criterion in address["target"].get("preserved_criteria", ())
        for citation in criterion["citations"]
    }
    retained, background = [], []
    for item in address["sources"]:
        binding, source = item["binding"], item["source"]
        if source.source_id in anchors or binding["role"] != "literature":
            retained.append(item)
        elif binding["match_scope"] in {"exact_allele", "complex_allele"}:
            retained.append(item)
        elif _mentions_allele(address["target"], source.content) and _relevant_background(address, item):
            retained.append(item)
        elif _relevant_background(address, item):
            background.append(item)
    publications = {}
    groups = {
        item["alias"]: item["binding"]["dependency_group"]
        if item["binding"]["dependency_group"].startswith("publication:")
        else _publication_key(item["source"])
        for item in background
    }
    for item in sorted(background, key=lambda value: _background_priority(address, value)):
        publications.setdefault(groups[item["alias"]], item)
    chosen = set(tuple(publications)[:_BACKGROUND_PAPERS])
    retained.extend(item for item in background if groups[item["alias"]] in chosen)
    allowed = {item["alias"] for item in retained}
    return tuple(item for item in address["sources"] if item["alias"] in allowed)


def _view(packet):
    sources = {f"SRC_{index:03d}": source for index, source in enumerate(packet["sources"], 1)}
    aliases = {source.source_id: alias for alias, source in sources.items()}
    targets = []
    for index, target in enumerate(packet["targets"], 1):
        bound = tuple(
            {
                "alias": aliases[binding["source_id"]],
                "source": sources[aliases[binding["source_id"]]],
                "binding": binding,
            }
            for binding in target["bindings"]
        )
        targets.append({"alias": f"A{index:03d}", "target": target, "sources": bound})
    return {"input": packet, "targets": tuple(targets), "sources": sources, "aliases": aliases}


def _source_alias(view, source_id):
    return view["aliases"][source_id]


def _call_value(call, name, default=None):
    original = call.get("original", call)
    if name == "ploidy":
        return len(original.get("alleles", ()))
    return call.get(name, original.get(name, default))


def _call_projection(call):
    return {
        "role": _call_value(call, "role"),
        "sample_id": _call_value(call, "sample_id"),
        "gt": _call_value(call, "gt"),
        "alleles": list(_call_value(call, "alleles", ())),
        "ad": None if _call_value(call, "ad") is None else list(_call_value(call, "ad")),
        "dp": _call_value(call, "dp"),
        "gq": _call_value(call, "gq"),
        "ploidy": _call_value(call, "ploidy"),
        "phase_set": _call_value(call, "phase_set"),
        "filter": _call_value(call, "filter_text", _call_value(call, "filter")),
        "genotype_state": call["genotype_state"],
        "quality_state": call["quality_state"],
        "usable": call["usable"],
        "carried": call["carried"],
        "selected_alt_dosage": call["selected_alt_dosage"],
        "warnings": list(call.get("warnings", ())),
    }


def _criterion_projection(view, criterion):
    projected = {
        "code": criterion["code"],
        "status": criterion["status"],
        "strength": criterion["strength"],
        "rationale": criterion["rationale"],
        "source_refs": [_source_alias(view, item["source_id"]) for item in criterion["citations"]],
        "independence_group": criterion["independence_group"],
        "origin": criterion["origin"],
    }
    if criterion.get("disease_binding"):
        projected["disease_binding"] = criterion["disease_binding"]
    return projected


def _lookup_projection(value, view, kind):
    if value is None:
        return None
    output = {"status": value["status"], "error_code": value.get("error_code", "")}
    item = value.get("annotation")
    if item is None:
        return output
    if kind == "annotation":
        fields = ("transcript_id", "hgvsc", "hgvsp", "protein_id", "exon", "intron", "canonical")
        record = {name: item.get(name) for name in fields}
        record["consequence_terms"] = list(item.get("consequence_terms", ()))
    else:
        fields = (
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
            "method",
            "research_use_only",
        )
        record = {name: item.get(name) for name in fields}
    if item.get("mechanism_context"):
        record["mechanism_context"] = item["mechanism_context"]
    record["source_ref"] = _source_alias(view, item["source"].source_id)
    if kind == "pvs1" and item.get("observed_transcript"):
        record["observed_transcript"] = item["observed_transcript"]
        resolved = item.get("resolved_annotation")
        if resolved:
            record["independent_annotation_source_ref"] = _source_alias(view, resolved["source"].source_id)
            record["assessed_hgvsc"], record["assessed_hgvsp"] = resolved["hgvsc"], resolved["hgvsp"]
    output["annotation"] = record
    return output


def _population_projection(value, view):
    if value is None:
        return None
    names = (
        "status",
        "match_method",
        "record_filter",
        "allele_count",
        "allele_number",
        "allele_frequency",
        "coverage_status",
        "coverage_adequate",
        "coverage_position_count",
        "coverage_expected_count",
        "coverage_over20_min",
        "error_code",
    )
    output = {name: value.get(name) for name in names}
    if "population_details" in value:
        output["population_details"] = value["population_details"]
    output["record_key"] = None if value.get("record_key") is None else _key_text(value["record_key"])
    output["source_ref"] = (
        None if value.get("source") is None else _source_alias(view, value["source"].source_id)
    )
    return output


def _predictions_projection(value, view):
    if value is None:
        return None
    output = {}
    fields = {
        "revel": ("score", "transcript", "amino_acid_reference", "amino_acid_alternate"),
        "spliceai": ("gene", "ds_ag", "ds_al", "ds_dg", "ds_dl", "dp_ag", "dp_al", "dp_dg", "dp_dl"),
    }
    for name, record_fields in fields.items():
        lookup = value[name]
        output[name] = {
            "status": lookup["status"],
            "records": [
                {field: row.get(field) for field in record_fields} for row in lookup.get("records", ())
            ],
            "source_ref": None
            if lookup.get("source") is None
            else _source_alias(view, lookup["source"].source_id),
            "error_code": lookup.get("error_code", ""),
        }
    return output


def _target_projection(view, address):
    target, primary = address["target"], address["target"]["primary"]
    return {
        "allele_id": address["alias"],
        "gene": _identity(target)[0],
        "variant_key": _key_text(_identity(target)[1]),
        "identity": {
            "transcript": primary.get("transcript"),
            "hgvs_c": primary.get("hgvs_c"),
            "hgvs_p": primary.get("hgvs_p"),
            "verified_aliases": [
                {
                    "transcript": item.get("transcript"),
                    "hgvs_c": item.get("hgvs_c"),
                    "hgvs_p": item.get("hgvs_p"),
                }
                for item in target.get("aliases", ())
            ],
        },
        "reference": {
            "status": target["reference"]["status"],
            "equivalent_region": list(target["reference"]["equivalent_region"]),
            "reference_accession": target["reference"].get("reference_accession"),
            "reference_observed": target["reference"].get("reference_observed"),
            "normalized_hgvs": target["reference"].get("normalized_hgvs"),
            "normalization_method": target["reference"].get("normalization_method"),
            "issues": list(target["reference"].get("issues", ())),
        },
        "observation": {
            "status": target["observation"]["status"],
            "origin": target["observation"]["origin"],
            "parentage_status": target["observation"]["parentage_status"],
            "family_conflict": target.get("family_conflict", False),
            "vcf_variant_key": _key_text(target["observation"]["original"]["vcf_variant_key"]),
            "quality": target["observation"]["original"].get("quality"),
            "filter": target["observation"]["original"].get("filter_text"),
            "calls": [_call_projection(item) for item in target["observation"]["calls"]],
        },
        "population": _population_projection(target.get("population"), view),
        "annotation": _lookup_projection(target.get("annotation"), view, "annotation"),
        "predictions": _predictions_projection(target.get("predictions"), view),
        "pvs1": _lookup_projection(target.get("pvs1"), view, "pvs1"),
        "disease_hypotheses": list(target.get("disease_labels", ())),
        "disease_ids": list(target.get("disease_ids", ())),
        "inheritance_hypotheses": list(target.get("inheritance_hypotheses", ())),
        "disease_contexts": [
            {
                **{name: row.get(name) for name in ("disease_id", "disease_name", "moi")},
                "source_refs": [_source_alias(view, source_id) for source_id in row["source_ids"]],
            }
            for row in target.get("disease_contexts", ())
        ],
    }


def _source_catalog(view, targets, *, include_aliases=None):
    targets = tuple(targets)
    roles, scopes, dependencies = {}, {}, {}
    for address in targets:
        for item in _selected_sources(address):
            alias, binding = item["alias"], item["binding"]
            roles.setdefault(alias, set()).add(binding["role"])
            group = dependencies.setdefault(binding["dependency_group"], f"D{len(dependencies) + 1:03d}")
            scopes.setdefault(alias, []).append(
                {
                    "allele_id": address["alias"],
                    "role": binding["role"],
                    "match_scope": binding["match_scope"],
                    "dependency_group": group,
                }
            )
    selected = set(scopes)
    if include_aliases is not None:
        selected &= include_aliases
    pending = list(selected)
    while pending:
        alias = pending.pop()
        if "native_annotation" in roles.get(alias, ()) or view["sources"][alias].locator.startswith(
            "method:assessment/disease-context/"
        ):
            continue
        for parent_id in view["sources"][alias].parent_source_ids:
            parent = _source_alias(view, parent_id)
            if parent not in selected:
                selected.add(parent)
                pending.append(parent)
    output = []
    for alias, source in view["sources"].items():
        if alias not in selected:
            continue
        output.append(
            {
                "source_ref": alias,
                "provider": source.provider,
                "kind": source.kind,
                "locator": source.locator,
                "genes": list(source.genes),
                "variant_keys": [_key_text(item) for item in source.variant_keys],
                "bindings": scopes.get(alias, []),
                "original_excerpt": source.content,
            }
        )
    return output


def _patient_projection(packet):
    profile = packet["profile"]
    key = case_tuple(profile["case_key"])
    terms = profile["terms"]
    settings = {
        json.loads(s.content).get("assay_type")
        for s in packet.get("sources", ())
        if s.provider == "patient" and s.locator == "input:patient/genomic-settings" and s.case_key == key
    }
    if len(settings) > 1:
        _fail("conflicting_assay_records")
    family = packet["family"]
    members = {r["sample_id"]: r for r in (family.get("pedigree") or {}).get("members", ())}
    return {
        "case": {"cohort": key[0], "case_id": key[1]},
        "sex": packet["family"].get("sex", "UNKNOWN"),
        "family_mode": packet["family"]["mode"],
        "assembly": family["assembly"],
        "assay_type": next(iter(settings), None),
        "pedigree_records": [
            {
                "role": role,
                "reported_sex": members[sample]["sex"],
                "reported_affected_code": members[sample]["affected"],
            }
            for role, sample in family["roles"].items()
            if sample in members
        ],
        "pedigree_interpretation": "PED codes: 2 affected, 1 unaffected, 0/-9 unknown. Recorded context only, not comprehensive clinical assessment, karyotype, independent segregation, or benign evidence.",
        "present_phenotypes": [
            {"hpo_id": item["hpo_id"], "label": item["label"]} for item in terms if item["present"]
        ],
        "negated_phenotypes": [
            {"hpo_id": item["hpo_id"], "label": item["label"]} for item in terms if not item["present"]
        ],
        "clinical_description": profile.get("clinical_description", ""),
    }


def _calls(address):
    return {_call_value(item, "role"): item for item in address["target"]["observation"]["calls"]}


def _prediction_source_ids(address):
    value = address["target"].get("predictions")
    if value is None:
        return set()
    return {
        lookup["source"].source_id
        for lookup in (value["revel"], value["spliceai"])
        if lookup["status"] == "found" and lookup.get("records") and lookup.get("source") is not None
    }


def _relationship_qc(address):
    text = "\n".join(
        item["source"].content for item in address["sources"] if item["binding"]["role"] == "observation_qc"
    )
    return all(
        token in text
        for token in (
            "status=kinship_supports_recorded_trio",
            "first_degree_both=true",
            "informative_sites_sufficient=true",
        )
    )


def _eligibility(code, address, refs):
    target = address["target"]
    if code == "PVS1" and target.get("local_haplotype_unresolved"):
        _fail("pvs1_local_haplotype_unresolved")
    source_ids = {item["source"].source_id for item in refs}
    if code in {"PS2", "PM6"}:
        calls = _calls(address)
        valid = (
            target["observation"]["origin"] == "candidate_de_novo"
            and not target.get("family_conflict", False)
            and all(role in calls and calls[role]["usable"] for role in ("proband", "father", "mother"))
            and _call_value(calls["proband"], "selected_alt_dosage") > 0
            and _call_value(calls["father"], "selected_alt_dosage") == 0
            and _call_value(calls["mother"], "selected_alt_dosage") == 0
            and any(item["binding"]["role"] == "patient_observation" for item in refs)
            and not any(item["binding"]["role"] in {"literature", "family_series_ledger"} for item in refs)
        )
        if not valid:
            _fail("de_novo_unverified_calls")
        if code == "PS2" and not _relationship_qc(address):
            _fail("ps2_parentage_not_verified")
        if code == "PM6" and _relationship_qc(address):
            _fail("same_case_verified_parentage_requires_ps2_scope")
    if code in {"PS3", "BS3", "PS4"}:
        eligible = []
        for item in refs:
            binding = item["binding"]
            if binding["role"] not in {"literature", "protein_comparator"}:
                continue
            if code in {"PS3", "BS3"} and binding["match_scope"] == "complex_allele":
                continue
            if (
                binding["match_scope"] == "exact_allele"
                or (binding["match_scope"] == "complex_allele" and code == "PS4")
                or _mentions_allele(target, item["source"].content)
            ):
                eligible.append(item)
        if not eligible:
            _fail("variant_specific_experiment_or_cases_missing_" + code.lower())
    if code in {"PP3", "BP4"} and not source_ids.intersection(_prediction_source_ids(address)):
        _fail("variant_prediction_missing_" + code.lower())
    if code == "PVS1" and not any(
        item["code"] == "PVS1" and item["status"] == "met" for item in target.get("preserved_criteria", ())
    ):
        mechanism = any(
            item["binding"]["role"] in {"gene_disease_context", "gene_context"}
            or (item["binding"]["role"] == "literature" and item["binding"]["match_scope"] == "context")
            for item in refs
        )
        annotation = target.get("annotation")
        transcript = bool(
            annotation
            and annotation["status"] == "annotated"
            and annotation.get("annotation")
            and annotation["annotation"]["source"].source_id in source_ids
            and set(annotation["annotation"].get("consequence_terms", ())).intersection(
                {
                    "frameshift_variant",
                    "stop_gained",
                    "splice_acceptor_variant",
                    "splice_donor_variant",
                }
            )
        )
        pvs1 = target.get("pvs1")
        transcript = transcript or bool(
            pvs1
            and pvs1["status"] == "assessed"
            and pvs1.get("annotation")
            and pvs1["annotation"].get("met") is True
            and pvs1["annotation"]["source"].source_id in source_ids
        )
        rna = any(
            item["binding"]["role"] in {"pvs1", "literature"}
            and item["binding"]["match_scope"] != "complex_allele"
            and _RNA_TERMS.search(item["source"].content)
            and (
                item["binding"]["match_scope"] == "exact_allele"
                or _mentions_allele(target, item["source"].content)
            )
            for item in refs
        )
        if not rna and not (mechanism and transcript):
            _fail("pvs1_mechanism_and_exact_transcript_sources_required")
    if code in {"PM2", "BA1", "BS1"}:
        value = target.get("population")
        valid = bool(value and value.get("source") and value["source"].source_id in source_ids)
        if code == "PM2":
            valid = (
                valid and value["status"] == "not_found_covered" and value.get("coverage_adequate") is True
            )
        else:
            valid = (
                valid
                and value["status"] == "found"
                and all(
                    value.get(name) is not None
                    for name in ("allele_count", "allele_number", "allele_frequency")
                )
            )
        if not valid:
            _fail("population_measurement_ineligible_" + code.lower())
    if (
        code in {"PP1", "BS2", "BS4"}
        and refs
        and all(item["binding"]["role"] in _MEASUREMENT_ROLES for item in refs)
        and not any("clinical_family_phenotypes" in item["source"].content for item in refs)
    ):
        _fail("family_clinical_phenotypes_missing_" + code.lower())


def _identity_key(value):
    return value["gene"], _key(value["variant_key"])


def _phase(left, right):
    if left.get("family_conflict", False) or right.get("family_conflict", False):
        return "conflicting_family_calls"
    left_calls = {_call_value(item, "role"): item for item in left["observation"]["calls"]}
    right_calls = {_call_value(item, "role"): item for item in right["observation"]["calls"]}
    for role in ("proband", "father", "mother"):
        first, second = left_calls.get(role), right_calls.get(role)
        if (
            not first
            or not second
            or not first["usable"]
            or not second["usable"]
            or _call_value(first, "sample_id") != _call_value(second, "sample_id")
        ):
            return "phase_unknown"
    if any(
        _call_value(calls["proband"], "selected_alt_dosage") != 1
        or _call_value(calls["proband"], "ploidy") != 2
        for calls in (left_calls, right_calls)
    ):
        return "phase_unknown"
    origins = {left["observation"]["origin"], right["observation"]["origin"]}
    if origins == {"paternal", "maternal"}:
        return "confirmed_in_trans"
    if left["observation"]["origin"] == right["observation"]["origin"] and left["observation"]["origin"] in {
        "paternal",
        "maternal",
    }:
        return "same_parent_origin_phase_unknown"
    return "phase_unknown"


def validate_input(packet):
    if not isinstance(packet, Mapping) or set(packet) != {"profile", "family", "targets", "sources"}:
        raise MendelPrimeError("schema_mismatch", field="assessment_input")
    if not packet["targets"] or type(packet["targets"]) is not tuple:
        raise MendelPrimeError("nonempty_tuple_required", field="assessment_input.targets")
    if type(packet["sources"]) is not tuple or any(
        not isinstance(item, Source) for item in packet["sources"]
    ):
        raise MendelPrimeError("source_tuple_required", field="assessment_input.sources")
    by_id = {source.source_id: source for source in packet["sources"]}
    if len(by_id) != len(packet["sources"]):
        raise MendelPrimeError("source_registered_more_than_once", field="assessment_input.sources")
    identities = [_identity(item) for item in packet["targets"]]
    targets_by_identity = dict(zip(identities, packet["targets"], strict=True))
    if len(identities) != len(set(identities)):
        raise MendelPrimeError("duplicate_allele_identity", field="assessment_input.targets")
    case_key = case_tuple(packet["profile"]["case_key"])
    if any(identity[1][0] != packet["family"]["assembly"] for identity in identities):
        raise MendelPrimeError("family_assembly_mismatch", field="assessment_input.targets")
    if any(source.case_key is not None and source.case_key != case_key for source in packet["sources"]):
        raise MendelPrimeError("case_scope_mismatch", field="assessment_input.sources")
    if any(parent not in by_id for source in packet["sources"] for parent in source.parent_source_ids):
        raise MendelPrimeError("parent_source_missing", field="assessment_input.sources")
    typed_kinds = {
        "population": "population",
        "transcript_annotation": "annotation",
        "prediction": "prediction",
        "pvs1": "prediction",
    }
    for target in packet["targets"]:
        gene, key = _identity(target)
        if _key(target["observation"]["original"]["variant_key"]) != key:
            raise MendelPrimeError("target_allele_mismatch", field="assessment_input.targets")
        if target["observation"]["status"] not in {"carried", "carried_unverified"}:
            raise MendelPrimeError("carried_observation_required", field="assessment_input.targets")
        if _key(target["reference"]["raw_key"]) != key:
            raise MendelPrimeError("target_allele_mismatch", field="assessment_input.targets")
        if target["reference"]["status"] == "reference_mismatch":
            raise MendelPrimeError("reference_mismatch_not_assessable", field="assessment_input.targets")
        if type(target.get("family_conflict", False)) is not bool:
            raise MendelPrimeError("boolean_required", field="assessment_input.targets.family_conflict")
        patient = [item for item in target["bindings"] if item["role"] == "patient_observation"]
        if not patient:
            raise MendelPrimeError(
                "patient_observation_source_binding_required", field="assessment_input.targets"
            )
        for binding in target["bindings"]:
            source = by_id.get(binding["source_id"])
            if source is None:
                raise MendelPrimeError("bound_source_missing", field="assessment_input.sources")
            if source.case_key is not None and source.case_key != case_key:
                raise MendelPrimeError("binding_case_scope_mismatch", field="assessment_input.sources")
            if source.genes and gene not in source.genes:
                raise MendelPrimeError("binding_gene_scope_mismatch", field="assessment_input.sources")
            if source.variant_keys and key not in source.variant_keys:
                raise MendelPrimeError("binding_allele_scope_mismatch", field="assessment_input.sources")
            if binding["match_scope"] == "exact_allele" and (
                source.genes != (gene,) or key not in source.variant_keys
            ):
                raise MendelPrimeError(
                    "exact_allele_binding_scope_required", field="assessment_input.sources"
                )
            if binding["match_scope"] in {"candidate_allele", "complex_allele"} and (
                (
                    binding["role"] != "literature"
                    and not (
                        binding["role"] == "partner_evidence" and binding["match_scope"] == "complex_allele"
                    )
                )
                or gene not in source.genes
            ):
                raise MendelPrimeError("literature_gene_scope_required", field="assessment_input.sources")
            if binding["role"] == "patient_observation" and (
                source.kind != "patient_observation"
                or source.case_key != case_key
                or source.genes != (gene,)
                or key not in source.variant_keys
            ):
                raise MendelPrimeError("patient_observation_scope_required", field="assessment_input.sources")
            expected_kind = typed_kinds.get(binding["role"])
            if expected_kind is not None and source.kind != expected_kind:
                raise MendelPrimeError("typed_fact_source_kind_mismatch", field="assessment_input.sources")
        typed = {
            "population": (
                ()
                if target.get("population") is None or target["population"].get("source") is None
                else (target["population"]["source"],)
            ),
            "transcript_annotation": (
                ()
                if target.get("annotation") is None or target["annotation"].get("annotation") is None
                else (target["annotation"]["annotation"]["source"],)
            ),
            "prediction": tuple(
                lookup["source"]
                for lookup in (
                    target.get("predictions", {}).get("revel", {}),
                    target.get("predictions", {}).get("spliceai", {}),
                )
                if lookup.get("source") is not None
            ),
            "pvs1": (
                ()
                if target.get("pvs1") is None or target["pvs1"].get("annotation") is None
                else (target["pvs1"]["annotation"]["source"],)
            ),
        }
        for role, sources in typed.items():
            if not sources:
                continue
            actual = tuple(item for item in target["bindings"] if item["role"] == role)
            if {item["source_id"] for item in actual} != {source.source_id for source in sources} or any(
                item["match_scope"] != "exact_allele" for item in actual
            ):
                raise MendelPrimeError("typed_fact_binding_mismatch", field="assessment_input.targets")
        for criterion in target.get("preserved_criteria", ()):
            bound = {item["source_id"] for item in target["bindings"]}
            for citation in criterion["citations"]:
                source = by_id.get(citation["source_id"])
                if source is None or citation["source_id"] not in bound:
                    raise MendelPrimeError("bound_source_required", field="criterion.citations")
                if citation["quote"] not in source.content:
                    raise MendelPrimeError("quote_not_in_source", field="criterion.citations")
        _validate_post_facts(target, targets_by_identity, by_id, case_key)
    return packet


def _validate_post_facts(target, targets, sources, case_key):
    """Keep recovered hypotheses and independent partner facts at the input boundary."""
    gene, key = _identity(target)
    bound = {b["source_id"] for b in target["bindings"]}
    for criterion in target.get("prior_criteria", ()):
        if (
            set(criterion)
            != {"code", "status", "strength", "rationale", "citations", "independence_group", "counted"}
            or type(criterion["counted"]) is not bool
        ):
            raise MendelPrimeError("prior_criterion_schema_invalid", field="assessment_input")
        if not re.fullmatch(
            r"PVS1(?:_RNA)?|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7]", str(criterion["code"])
        ):
            raise MendelPrimeError("prior_criterion_code_invalid", field="assessment_input")
        for citation in criterion["citations"]:
            source = sources.get(citation.get("source_id"))
            if (
                source is None
                or source.source_id not in bound
                or not citation.get("quote")
                or citation["quote"] not in source.content
            ):
                raise MendelPrimeError("prior_criterion_quote_invalid", field="assessment_input")
    for other in target.get("local_haplotype_unresolved", ()):
        partner = targets.get((gene, _key(other)))
        if partner is None or _key(other) == key or _phase(target, partner) == "confirmed_in_trans":
            raise MendelPrimeError("haplotype_identity_invalid", field="assessment_input")
    facts = target.get("partner_facts", ())
    actual = {b["source_id"] for b in target["bindings"] if b["role"] == "partner_evidence"}
    if actual != {f["source_id"] for f in facts} or len(actual) != len(facts):
        raise MendelPrimeError("partner_fact_binding_invalid", field="assessment_input")
    for fact in facts:
        other = _key(fact["partner_key"])
        partner, source = targets.get((gene, other)), sources.get(fact["source_id"])
        if (
            partner is None
            or other == key
            or source is None
            or source.provider != "mendelprime"
            or source.kind != "derived"
            or source.case_key != case_key
            or source.genes != (gene,)
            or set(source.variant_keys) != {key, other}
            or _phase(target, partner) != "confirmed_in_trans"
            or fact["classification"] not in _POSITIVE_CLASSES
        ):
            raise MendelPrimeError("independent_partner_scope_invalid", field="assessment_input")
        shared = {
            (d["disease_id"], _ranking_mode(d.get("mode_of_inheritance", d.get("moi"))))
            for d in target.get("disease_contexts", ())
            if d.get("disease_id")
        }
        paired = {
            (d["disease_id"], _ranking_mode(d.get("mode_of_inheritance", d.get("moi"))))
            for d in partner.get("disease_contexts", ())
            if d.get("disease_id")
        }
        if len(shared) != 1 or shared != paired or shared != {(fact["disease"]["disease_id"], "AR")}:
            raise MendelPrimeError("independent_partner_disease_invalid", field="assessment_input")
        record = json.loads(source.content)
        if (
            record.get("PM3_excluded") is not True
            or record.get("phase") != "confirmed_in_trans"
            or record.get("target_key") != variant_text(key)
            or record.get("partner_key") != variant_text(other)
            or record.get("independent_partner_pathogenicity") != fact["classification"]
            or _canonical_json(record.get("disease")) != _canonical_json(fact["disease"])
        ):
            raise MendelPrimeError("independent_partner_record_invalid", field="assessment_input")
        criteria = []
        for item in record["basis"]:
            if item["code"] in {"PM3", "PP5", "BP6"}:
                raise MendelPrimeError("circular_partner_criterion", field="assessment_input")
            for citation in item["citations"]:
                parent = sources.get(citation["source_id"])
                if (
                    parent is None
                    or parent.source_id not in source.parent_source_ids
                    or not citation["quote"]
                    or citation["quote"] not in parent.content
                ):
                    raise MendelPrimeError("partner_parent_quote_invalid", field="assessment_input")
            criteria.append(
                {
                    **item,
                    "status": "met",
                    "direction": "benign" if item["code"].startswith("B") else "pathogenic",
                }
            )
        if (
            len({c["code"] for c in criteria}) != len(criteria)
            or _ranking_classify(criteria) != fact["classification"]
        ):
            raise MendelPrimeError("independent_partner_class_invalid", field="assessment_input")


def validate_result(result):
    if isinstance(result, Mapping) and result.get("workflow") == "clinical_profile_v1":
        return _validate_profiles(result)
    if not isinstance(result, Mapping) or set(result) != {"independent", "partners", "final", "sources"}:
        raise MendelPrimeError("schema_mismatch", field="assessment_result")
    expected = tuple(_identity_key(item["identity"]) for item in result["independent"])
    if not expected or len(expected) != len(set(expected)):
        raise MendelPrimeError("duplicate_or_empty_identity", field="assessment_result.independent")
    if tuple(_identity_key(item["identity"]) for item in result["partners"]) != expected:
        raise MendelPrimeError("independent_coverage_or_order_invalid", field="assessment_result.partners")
    if tuple(_identity_key(item["identity"]) for item in result["final"]) != expected:
        raise MendelPrimeError("independent_coverage_or_order_invalid", field="assessment_result.final")
    by_source = {source.source_id: source for source in result["sources"]}
    if len(by_source) != len(result["sources"]):
        raise MendelPrimeError("source_registered_more_than_once", field="assessment_result.sources")
    reviews = {_identity_key(item["identity"]): item for item in result["independent"]}
    partners = {_identity_key(item["identity"]): item for item in result["partners"]}
    for review in result["independent"]:
        by_code = {item["code"]: item for item in review["criteria"]}
        if any(code not in by_code or by_code[code]["status"] != "met" for code in review["decisive_codes"]):
            raise MendelPrimeError("met_criterion_required", field="independent.decisive_codes")
        for criterion in review["criteria"]:
            for citation in criterion["citations"]:
                source = by_source.get(citation["source_id"])
                if source is None or citation["quote"] not in source.content:
                    raise MendelPrimeError("citation_source_or_quote_invalid", field="result.independent")
        gene, key = _identity_key(review["identity"])
        for source_id in review["decisive_source_ids"]:
            source = by_source.get(source_id)
            if source is None:
                raise MendelPrimeError("decisive_source_missing", field="result.independent")
            if (source.genes and gene not in source.genes) or (
                source.variant_keys and key not in source.variant_keys
            ):
                raise MendelPrimeError("decisive_source_scope_mismatch", field="result.independent")
    for partner in result["partners"]:
        for citation in partner["criterion"]["citations"]:
            source = by_source.get(citation["source_id"])
            if source is None or citation["quote"] not in source.content:
                raise MendelPrimeError("citation_source_or_quote_invalid", field="result.partners")
    for item in result["final"]:
        identity = _identity_key(item["identity"])
        review, partner = reviews[identity], partners[identity]
        if item["disease"] != review["disease"]:
            raise MendelPrimeError("disease_context_changed", field="result.final")
        criteria = {criterion["code"]: criterion for criterion in review["criteria"]}
        criteria[partner["criterion"]["code"]] = partner["criterion"]
        if any(code not in criteria or criteria[code]["status"] != "met" for code in item["decisive_codes"]):
            raise MendelPrimeError("decisive_reviewed_criterion_required", field="result.final")
        if any(
            code not in criteria or criteria[code]["status"] == "met" for code in item["limitation_codes"]
        ):
            raise MendelPrimeError("limitation_reviewed_criterion_required", field="result.final")
        partner_sources = (
            {citation["source_id"] for citation in partner["criterion"]["citations"]}
            if partner["criterion"]["status"] == "met"
            else set()
        )
        gene, key = identity
        for source_id in (*item["decisive_source_ids"], *item["contextual_source_ids"]):
            source = by_source.get(source_id)
            if source is None:
                raise MendelPrimeError("final_source_missing", field="result.final")
            if (
                (source.genes and gene not in source.genes)
                or (source.variant_keys and key not in source.variant_keys)
            ) and source_id not in partner_sources:
                raise MendelPrimeError("source_scope_mismatch", field="result.final")
    return result


def _progress_value(value):
    if isinstance(value, Source):
        return {"source_id": value.source_id, "content_sha256": value.content_sha256}
    if isinstance(value, Path):
        return {"path": str(value)}
    if isinstance(value, Mapping):
        return {name: _progress_value(item) for name, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_progress_value(item) for item in value]
    if type(value) is float:
        return value if math.isfinite(value) else None
    if value is None or type(value) in {str, int, bool}:
        return value
    raise MendelPrimeError("assessment_progress_input_not_serializable", field="assessment_progress")


def _progress_digest(value):
    return hashlib.sha256(_canonical_json(_progress_value(value)).encode("utf-8")).hexdigest()


def _progress_identities(packet):
    return tuple(
        {"gene": _identity(target)[0], "variant_key": _identity(target)[1]} for target in packet["targets"]
    )


def _progress_aliases(identities, sources):
    return {
        "alleles": tuple(
            {"alias": f"A{index:03d}", "identity": identity} for index, identity in enumerate(identities, 1)
        ),
        "sources": tuple(
            {"alias": f"SRC_{index:03d}", "source_id": source.source_id}
            for index, source in enumerate(sources, 1)
        ),
    }


def _progress_guard(packet):
    """Only run() computes current code/prompt constraints; loading is pure."""
    constraints = {
        "policy": "clinical-profile-resume-v1",
        "acmg_code": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "prompting_code": hashlib.sha256(Path(prompting.__file__).read_bytes()).hexdigest(),
        "profile_prompt": prompting.text("acmg_evidence_profile"),
        "classification_prompt": prompting.text("acmg_profile_classification"),
        "qualitative_policy": prompting.text("acmg_qualitative_policy"),
        "profile_tokens": 2800,
        "classification_tokens": 1500,
        "empty_support_contract": "key_evidence_ids=[]; no S source IDs",
    }
    return {
        "input_sha256": _progress_digest(packet),
        "constraints_sha256": _progress_digest(constraints),
        "aliases": _progress_aliases(_progress_identities(packet), packet["sources"]),
    }


def validate_progress(progress):
    """Validate a hydrated partial record without I/O, inference or source loading.

    Runtime codecs encode sources as references and hydrate them from the single
    authoritative provenance registry before calling this function.
    """
    try:
        return _validate_progress_record(progress)
    except MendelPrimeError:
        raise
    except (KeyError, TypeError, IndexError, AttributeError) as error:
        raise MendelPrimeError("assessment_progress_invalid", field="assessment_progress") from error


def _validate_progress_record(progress):
    fields = {
        "workflow",
        "status",
        "case_key",
        "profiles",
        "final",
        "sources",
        "expected_identities",
        "guard",
        "issues",
    }
    if (
        not isinstance(progress, Mapping)
        or set(progress) != fields
        or progress["workflow"] not in {"clinical_profile_partial_v1", "clinical_profile_partial_v2"}
        or progress["status"] != "degraded"
    ):
        raise MendelPrimeError("assessment_progress_schema_invalid", field="assessment_progress")
    if any(
        type(progress[name]) is not tuple
        for name in ("profiles", "final", "sources", "expected_identities", "issues")
    ):
        raise MendelPrimeError("assessment_progress_tuple_required", field="assessment_progress")
    case_key = case_tuple(progress["case_key"])
    identities = tuple(_identity_key(row) for row in progress["expected_identities"])
    if not identities or len(identities) != len(set(identities)) or len(progress["final"]) >= len(identities):
        raise MendelPrimeError("assessment_progress_target_coverage_invalid", field="assessment_progress")
    if any(set(row) != {"gene", "variant_key"} for row in progress["expected_identities"]):
        raise MendelPrimeError("assessment_progress_identity_invalid", field="assessment_progress")
    profile_ids = tuple(_identity_key(row["identity"]) for row in progress["profiles"])
    final_ids = tuple(_identity_key(row["identity"]) for row in progress["final"])
    if profile_ids != tuple(key for key in identities if key in profile_ids):
        raise MendelPrimeError("assessment_progress_profile_order_invalid", field="assessment_progress")
    if any(
        not isinstance(source, Source) or source.case_key is not None and source.case_key != case_key
        for source in progress["sources"]
    ):
        raise MendelPrimeError("assessment_progress_source_case_mismatch", field="assessment_progress")
    guard = progress["guard"]
    if (
        not isinstance(guard, Mapping)
        or set(guard) != {"input_sha256", "constraints_sha256", "aliases"}
        or any(
            not isinstance(guard[name], str) or not re.fullmatch(r"[0-9a-f]{64}", guard[name])
            for name in ("input_sha256", "constraints_sha256")
        )
    ):
        raise MendelPrimeError("assessment_progress_guard_invalid", field="assessment_progress")
    if _canonical_json(guard["aliases"]) != _canonical_json(
        _progress_aliases(progress["expected_identities"], progress["sources"])
    ):
        raise MendelPrimeError("assessment_progress_alias_mismatch", field="assessment_progress")
    if any(
        set(row) != {"stage", "code"} or not isinstance(row["code"], str) or not row["code"]
        for row in progress["issues"]
    ):
        raise MendelPrimeError("assessment_progress_issue_invalid", field="assessment_progress")
    if progress["workflow"] == "clinical_profile_partial_v1":
        if (
            profile_ids != identities[: len(profile_ids)]
            or final_ids != profile_ids[: len(final_ids)]
            or not len(final_ids) <= len(profile_ids) <= len(final_ids) + 1
        ):
            raise MendelPrimeError("assessment_progress_profile_prefix_invalid", field="assessment_progress")
        pending = (len(final_ids),)
    else:
        pending = tuple(index for index, key in enumerate(identities) if key not in final_ids)
    expected_stages = tuple(
        f"{'vcf_profile_classification' if identities[index] in profile_ids else 'vcf_evidence_profile'}:A{index + 1:03d}"
        for index in pending
    )
    if tuple(row["stage"] for row in progress["issues"]) != expected_stages:
        raise MendelPrimeError("assessment_progress_stage_invalid", field="assessment_progress")
    _validate_profiles(
        {
            "workflow": "clinical_profile_v1",
            "profiles": progress["profiles"],
            "final": progress["final"],
            "sources": progress["sources"],
        },
        allow_pending=True,
    )
    return progress


def _resume_profile(profile, address, visible, anchor):
    aliases = {source.source_id: alias for alias, source in visible.items()}

    def citations(rows):
        if any(row["source_id"] not in aliases for row in rows):
            raise MendelPrimeError("assessment_progress_profile_source_changed", field="assessment_progress")
        return [{"source_id": aliases[row["source_id"]], "quote": row["quote"]} for row in rows]

    raw = {
        "allele_disease": {
            **profile["allele_disease"],
            "citations": citations(profile["allele_disease"]["citations"]),
        },
        "patient_fit": {
            **profile["patient_fit"],
            "citations": citations(profile["patient_fit"]["citations"]),
        },
        "key_support": [
            {
                **{
                    name: row[name]
                    for name in ("kind", "finding", "pathogenic_relevance", "independence_group")
                },
                "reliability": row.get("reported_reliability", row["reliability"]),
                "citations": citations(row["citations"]),
            }
            for row in profile["key_support"]
        ],
        "key_limits": [
            {"finding": row["finding"], "effect": row["effect"], "citations": citations(row["citations"])}
            for row in profile["key_limits"]
        ],
        "unknowns": list(profile["unknowns"]),
        "summary": profile["summary"],
    }
    decoded = _decode_profile(raw, address, visible, anchor)
    if _canonical_json(decoded) != _canonical_json(profile):
        raise MendelPrimeError("assessment_progress_profile_changed", field="assessment_progress")
    return decoded


def _validated_resume(progress, packet, view, guard):
    validate_progress(progress)
    if (
        progress["case_key"] != case_tuple(packet["profile"]["case_key"])
        or _canonical_json(progress["expected_identities"]) != _canonical_json(_progress_identities(packet))
        or progress["guard"]["input_sha256"] != guard["input_sha256"]
        or _canonical_json(progress["guard"]["aliases"]) != _canonical_json(guard["aliases"])
        or progress["guard"]["constraints_sha256"] != guard["constraints_sha256"]
    ):
        raise MendelPrimeError(
            "assessment_progress_input_or_constraints_changed", field="assessment_progress"
        )
    if tuple(source.source_id for source in progress["sources"]) != tuple(
        source.source_id for source in packet["sources"]
    ):
        raise MendelPrimeError("assessment_progress_sources_changed", field="assessment_progress")
    profiles = []
    addresses = {_identity(row["target"]): row for row in view["targets"]}
    for profile in progress["profiles"]:
        address = addresses[_identity_key(profile["identity"])]
        _, visible, anchor = _profile_input(view, address)
        profiles.append(_resume_profile(profile, address, visible, anchor))
    return tuple(profiles), progress["final"]


def completed_assessment(progress):
    """Expose only complete profile/classification pairs, never a substitute grade."""
    validate_progress(progress)
    if not progress["final"]:
        return None
    identities = {_identity_key(row["identity"]) for row in progress["final"]}
    return validate_result(
        {
            "workflow": "clinical_profile_v1",
            "profiles": tuple(
                row for row in progress["profiles"] if _identity_key(row["identity"]) in identities
            ),
            "final": progress["final"],
            "sources": progress["sources"],
        }
    )


def run(model, assessment_input, *, progress=None, workers=1, target_identity=None):
    """Review ungraded evidence, then classify each allele from that profile only."""
    if not callable(getattr(model, "complete", None)):
        raise MendelPrimeError("model_required", field="acmg.model")
    if type(workers) is not int or workers < 1:
        raise MendelPrimeError("acmg_workers_invalid", field="acmg.workers")
    packet = validate_input(assessment_input)
    patient = Source(
        provider="patient",
        kind="patient_observation",
        locator="input:patient/acmg-profile",
        version="normalized-profile-v2",
        case_key=case_tuple(packet["profile"]["case_key"]),
        content=_canonical_json(_patient_projection(packet)),
    )
    packet = {
        **packet,
        "sources": packet["sources"]
        if any(s.source_id == patient.source_id for s in packet["sources"])
        else (*packet["sources"], patient),
        "targets": tuple(
            {
                **target,
                "bindings": (
                    *target["bindings"],
                    {
                        "source_id": patient.source_id,
                        "role": "patient_phenotype",
                        "match_scope": "context",
                        "dependency_group": "patient_phenotype",
                    },
                ),
            }
            for target in packet["targets"]
        ),
    }
    view = _view(packet)
    selected = view["targets"]
    if target_identity is not None:
        if progress is not None:
            raise MendelPrimeError("selected_assessment_cannot_resume", field="acmg.progress")
        selected = tuple(row for row in selected if _identity(row["target"]) == target_identity)
        if len(selected) != 1:
            raise MendelPrimeError("assessment_target_not_found", field="acmg.target_identity")
    guard = _progress_guard(packet)
    saved_profiles, saved_final = (
        _validated_resume(progress, packet, view, guard) if progress is not None else ((), ())
    )
    profiles, final, issues = [], [], []
    profile_system = prompting.text("acmg_evidence_profile")
    classification_template = prompting.text("acmg_profile_classification")
    saved_profiles = {_identity_key(row["identity"]): row for row in saved_profiles}
    saved_final = {_identity_key(row["identity"]): row for row in saved_final}
    case_key = case_tuple(packet["profile"]["case_key"])

    def assess(address):
        identity_key = _identity(address["target"])
        if identity_key in saved_final:
            return saved_profiles[identity_key], saved_final[identity_key], None
        payload, visible, anchor = _profile_input(view, address)
        stage = f"vcf_evidence_profile:{address['alias']}"
        profile = saved_profiles.get(identity_key)
        try:
            if profile is None:
                profile = prompting.call_with_repair(
                    model,
                    case_key,
                    stage,
                    profile_system,
                    payload,
                    2800,
                    lambda raw: _decode_profile(raw, address, visible, anchor),
                )
            identity = profile["identity"]
            grade_input = {
                "case_id": case_key[1],
                "variant_key": _key_text(identity["variant_key"]),
                "allele": {"gene": identity["gene"], "variant_key": _key_text(identity["variant_key"])},
                **{
                    name: profile[name]
                    for name in ("allele_disease", "patient_fit", "key_support", "key_limits", "unknowns")
                },
                "evidence_summary": profile["summary"],
                **(
                    {"retained_observations": profile["retained_observations"]}
                    if profile.get("retained_observations")
                    else {}
                ),
            }
            classification_system = classification_template
            if profile.get("retained_observations"):
                classification_system += (
                    "\nretained_observations are immutable neutral source facts, not extra support criteria. "
                    "Return observation_review with exactly one {id,status,reason} for each O-ID. "
                    "status is considered or not_applicable; give a nonempty reason. Never place O-IDs in "
                    "key_evidence_ids. Consider REVEL as a computational prediction and SpliceAI as a separate "
                    "splicing prediction; a low splice score does not refute a missense protein effect. "
                    "Do not double-count a prediction and its already-cited source or a functional assay as "
                    "independent demonstrations of the same effect. A ledger and its paper are one study, not "
                    "two cohorts. No PP3/PS4 strength or positive class is implied by these observations. "
                    "Gene mechanism observations are disease-bound background for interpreting an existing "
                    "allele-level fact, not proof that this allele produces that effect. They cannot add "
                    "PVS1, functional evidence, or a separate independent criterion."
                    " Native transcript/mechanism observations are descriptive algorithm context, not "
                    "an assay or a PVS1 anchor; missing population records do not prove adequate coverage."
                )
            if not profile["key_support"]:
                classification_system += (
                    "\nkey_support is empty: key_evidence_ids must be []. Do not return S source IDs "
                    "or invent evidence IDs. Source citations are not profile support IDs."
                )
            stage = f"vcf_profile_classification:{address['alias']}"
            classification_system += "\n\n" + prompting.text("acmg_qualitative_policy").strip()
            decision = prompting.call_with_repair(
                model,
                case_key,
                stage,
                classification_system,
                grade_input,
                1500,
                lambda raw: _decode_profile_classification(raw, profile),
            )
            return profile, decision, None
        except MendelPrimeError as error:
            if error.category not in {"model", "response"}:
                raise
            return profile, None, {"stage": stage, "code": error.code}

    with ThreadPoolExecutor(max_workers=min(workers, len(selected)) or 1) as pool:
        for profile, decision, issue in pool.map(assess, selected):
            if profile is not None:
                profiles.append(profile)
            if decision is not None:
                final.append(decision)
            if issue is not None:
                issues.append(issue)
    if issues:
        if target_identity is not None:
            raise MendelPrimeError(issues[0]["code"], category="response")
        failure = MendelPrimeError("assessment_incomplete", category="execution")
        failure.partial_result = validate_progress(
            {
                "workflow": "clinical_profile_partial_v2",
                "status": "degraded",
                "case_key": case_key,
                "profiles": tuple(profiles),
                "final": tuple(final),
                "sources": packet["sources"],
                "expected_identities": _progress_identities(packet),
                "guard": guard,
                "issues": tuple(issues),
            }
        )
        raise failure
    return validate_result(
        {
            "workflow": "clinical_profile_v1",
            "profiles": tuple(profiles),
            "final": tuple(final),
            "sources": packet["sources"],
        }
    )


def _profile_anchor(address):
    target = address["target"]
    if target.get("local_haplotype_unresolved"):
        return None
    for criterion in target.get("preserved_criteria", ()):
        if criterion["code"] == "PVS1" and criterion["status"] == "met":
            return criterion
    if target.get("post_evidence_recovered"):
        return None
    lookup = target.get("pvs1") or {}
    annotation = lookup.get("annotation") or {}
    strength = {
        "VeryStrong": "very_strong",
        "Strong": "strong",
        "Moderate": "moderate",
        "Supporting": "supporting",
    }.get(annotation.get("strength"), annotation.get("strength"))
    if (
        lookup.get("status") == "assessed"
        and annotation.get("met") is True
        and strength in _PVS1_STRENGTHS
        and annotation.get("gene_match") in {True, "exact", "resource_alias"}
        and annotation.get("decision_tree_supported") is True
    ):
        source = annotation["source"]
        return _criterion(
            "PVS1",
            "met",
            strength,
            "Validated source-bound PVS1 resource assessment.",
            (_citation(source.source_id, source.content),),
            "validated_pvs1_mechanism_and_allele_effect",
            "validated_resource",
        )
    return None


def _profile_input(view, address):
    target, packet = address["target"], view["input"]
    identity = _identity(target)
    projected = _target_projection(view, address)
    catalog = _source_catalog(view, (address,))
    visible = {row["source_ref"]: view["sources"][row["source_ref"]] for row in catalog}
    phase = []
    for partner in view["targets"]:
        if _identity(partner["target"])[0] != identity[0] or _identity(partner["target"]) == identity:
            continue
        partner_refs = [
            item["alias"] for item in partner["sources"] if item["binding"]["role"] == "patient_observation"
        ]
        for row in _source_catalog(view, (partner,), include_aliases=set(partner_refs)):
            if row["source_ref"] not in visible:
                catalog.append(row)
                visible[row["source_ref"]] = view["sources"][row["source_ref"]]
        fact = next(
            (
                row
                for row in target.get("partner_facts", ())
                if row["partner_key"] == _identity(partner["target"])[1]
            ),
            None,
        )
        phase.append(
            {
                "configuration": _phase(target, partner["target"]),
                "partner_key": _key_text(_identity(partner["target"])[1]),
                "partner_observation_source_refs": partner_refs,
                "independent_partner_pathogenicity": {
                    "classification": fact["classification"],
                    "disease": fact["disease"],
                    "source_ref": _source_alias(view, fact["source_id"]),
                    "PM3_excluded": True,
                }
                if fact
                else None,
                "interpretation": "Same-disease AR, confirmed-trans, independently computed PM3-free partner evidence."
                if fact
                else "Measured family context only; no independent partner grade is supplied.",
            }
        )
    anchor = _profile_anchor(address)
    patient = _patient_projection(packet)
    published = _publication_families(address)
    relationship = _relationship_qc(address)
    if target.get("post_evidence_recovered"):
        projected["observation"] = {
            **projected["observation"],
            "parentage_status": "genetically_supported" if relationship else "unconfirmed",
            "clinical_identity_documentation": "not_supplied",
        }
    evidence = {
        "case_id": patient["case"]["case_id"],
        "patient": patient,
        "target": {
            name: projected[name]
            for name in (
                "gene",
                "variant_key",
                "identity",
                "reference",
                "annotation",
                "predictions",
                "pvs1",
                *(
                    ("disease_contexts",)
                    if "disease_contexts" in target
                    else ("disease_hypotheses", "disease_ids", "inheritance_hypotheses")
                ),
            )
        },
        "pvs1_validated_anchor": _criterion_projection(view, anchor)
        if anchor
        else {"status": "not_assessed"},
        "patient_de_novo_review": {
            "current_family_observations": [projected["observation"]],
            "counted_external_families": [
                {
                    "publication": value["publication"],
                    "families": value["families"],
                    "within_publication_family_count": value["within_publication_family_count"],
                    "ledger_source_ref": _source_alias(view, source.source_id),
                }
                for source, value in published
            ],
            "reported_external_de_novo_groups": [
                value["reported_de_novo"] for _, value in published if value["reported_de_novo"] is not None
            ],
            "literature_limitations": [
                "Within-study recurrence (PS4 context) is separate from PS2/PM6. Aggregate de novo counts have no assigned family IDs or confirmed relationships. No cross-publication total or evidence strength is assigned."
                if published
                else "No independent external-family count ledger is supplied."
            ],
        },
        "population_observations": projected["population"],
        "independent_source_ids": sorted({row["alias"] for row in address["sources"] if _support_source(row, target)}),
        "recessive_phase_context": {"observations": phase},
        "sources": [
            {
                "id": row["source_ref"],
                "text": row["original_excerpt"],
                **{
                    name: row[name]
                    for name in ("provider", "kind", "locator", "genes", "variant_keys", "bindings")
                },
            }
            for row in catalog
        ],
        "input_contract": "clinical_profile_v1: source facts; no Gold, rank, prior model class or partner grade",
    }
    if "prior_criteria" in target:
        evidence["prior_criterion_proposals"] = [
            {
                **{
                    name: criterion.get(name)
                    for name in ("code", "status", "strength", "rationale", "counted")
                },
                "interpretation": "Frozen ranking hypothesis, not an independent pathogenicity assertion; verify the cited source facts.",
                "citations": [
                    {"source_id": _source_alias(view, c["source_id"]), "quote": c["quote"]}
                    for c in criterion["citations"]
                ],
            }
            for criterion in target["prior_criteria"]
        ]
        evidence["patient_de_novo_review"]["genetic_relationships_supported"] = relationship
        evidence["patient_de_novo_review"]["relationship_evidence"] = {
            "genetic_relationship_support": "supported" if relationship else "not_established",
            "clinical_identity_documentation": "not_supplied",
            "scope": "current_patient_only; no inference about families in publications",
        }
        evidence["local_haplotype_uncertainty"] = [
            variant_text(key) for key in target.get("local_haplotype_unresolved", ())
        ]
        evidence["curated_allele_assertions"] = []
        for row in catalog:
            if not any(b["role"] == "allele_assertion" for b in row.get("bindings", ())):
                continue
            raw = _ranking_json(visible[row["source_ref"]])
            raw = raw.get("record", raw.get("row", raw))
            evidence["curated_allele_assertions"].append(
                {
                    "source_id": row["source_ref"],
                    "record": {
                        k: v
                        for k, v in raw.items()
                        if k.lower()
                        in {
                            "name",
                            "gene_symbol",
                            "classification",
                            "clinical_significance",
                            "review_status",
                            "variation_id",
                            "disease",
                            "disease_id",
                            "disease_name",
                            "criteria_met",
                        }
                    },
                    "interpretation": "Curated context only; not an independent functional experiment or PP5 criterion.",
                }
            )
        evidence["input_contract"] = (
            "clinical_profile_v1: no Gold or rank; no previous model grade; partner evidence is independently deterministic and PM3-free"
        )
    return (
        {
            "case_id": patient["case"]["case_id"],
            "variant_key": _key_text(identity[1]),
            "evidence_packet": evidence,
        },
        visible,
        anchor,
    )


def validate_family_series_source(source, sources):
    """Validate the same publication facts on resource load, inference and saved-result read."""
    value = json.loads(source.content)
    fields = {
        "schema_version",
        "gene",
        "variant_key",
        "article_source_id",
        "identity_source_id",
        "publication",
        "variant_table_original",
        "phenotype_header_original",
        "families",
        "within_publication_family_count",
        "recruitment_original",
        "methods_original",
        "phenotype_limitations_original",
        "reported_de_novo",
        "cross_publication_deduplication",
        "assigned_ps2_pm6_points",
    }
    if (
        set(value) != fields
        or source.provider != "publication-family-ledger"
        or source.kind != "annotation"
        or len(source.genes) != 1
        or len(source.variant_keys) != 1
        or value["schema_version"] != "publication_family_series_v1"
        or value["gene"] != source.genes[0]
        or value["variant_key"] != variant_text(source.variant_keys[0])
    ):
        _fail("publication_family_ledger_scope_invalid")
    paper, identity = sources.get(value["article_source_id"]), sources.get(value["identity_source_id"])
    if (
        paper is None
        or identity is None
        or paper.kind != "literature"
        or paper.provider not in {"pubmed", "pmc"}
        or identity.provider != "clinvar"
        or identity.kind != "annotation"
        or paper.source_id == identity.source_id
        or paper.locator != value["publication"]
        or source.locator != paper.locator.rstrip("/") + "/#family-series"
        or source.version
        != "publication-family-series-v1;sha256=" + hashlib.sha256(source.content.encode()).hexdigest()
        or any(s.genes != source.genes or s.variant_keys != source.variant_keys for s in (paper, identity))
    ):
        _fail("publication_family_ledger_scope_invalid")
    rows, dn = value["families"], value["reported_de_novo"]
    ids = [r["family_id"] for r in rows]
    match = re.search(r"(\d+(?:,\s*\d+)*)$", value["variant_table_original"])
    if (
        not rows
        or len(set(ids)) != len(rows)
        or value["within_publication_family_count"] != len(rows)
        or match is None
        or match[1].replace(" ", "").split(",") != ids
        or value["variant_table_original"].split()[0] not in identity.content
        or value["cross_publication_deduplication"] != "not_performed"
        or value["assigned_ps2_pm6_points"] is not None
    ):
        _fail("publication_family_count_invalid")
    for row in rows:
        if set(row) != {
            "family_id",
            "affected_members",
            "sporadic_marker",
            "phenotypes",
            "original_table_row",
        }:
            _fail("publication_family_row_invalid")
        raw = row["original_table_row"].split()
        sporadic = raw[2:3] == ["(S)"]
        header = re.sub(r"[^a-z]", "", value["phenotype_header_original"].lower())
        positions = {name: header.find(re.sub(r"[^a-z]", "", name.lower())) for name in row["phenotypes"]}
        cells = [row["phenotypes"][name] for name in sorted(positions, key=positions.get)]
        if (
            row["original_table_row"] not in paper.content
            or raw[:2] != [row["family_id"], str(row["affected_members"])]
            or type(row["affected_members"]) is not int
            or row["affected_members"] < 1
            or not positions
            or min(positions.values()) < 0
            or len(set(positions.values())) != len(positions)
            or row["sporadic_marker"] is not sporadic
            or cells != (raw[3:] if sporadic else raw[2:])
        ):
            _fail("publication_family_row_invalid")
    if dn is not None:
        if (
            set(dn)
            != {
                "count",
                "family_ids",
                "family_id_mapping",
                "parent_variant_analysis",
                "parental_relationships",
                "original",
            }
            or type(dn["count"]) is not int
            or not 0 <= dn["count"] <= len(rows)
            or dn["family_ids"] is not None
            or dn["family_id_mapping"] != "not_reported"
            or dn["parent_variant_analysis"] != "reported"
            or dn["parental_relationships"] != "not_documented"
        ):
            _fail("publication_de_novo_group_invalid")
        words = "zero one two three four five six seven eight nine ten".split()
        counts = (str(dn["count"]), words[dn["count"]]) if dn["count"] < len(words) else (str(dn["count"]),)
        totals = (str(len(rows)), words[len(rows)]) if len(rows) < len(words) else (str(len(rows)),)
        if not any(
            f"{n} of {total} unrelated families" in " ".join(dn["original"]).lower()
            for n in counts
            for total in totals
        ):
            _fail("publication_de_novo_count_quote_mismatch")
    for quote in (*(dn["original"] if dn else ()), *(value[k] for k in fields if k.endswith("_original"))):
        if not isinstance(quote, str) or not quote or quote not in paper.content:
            _fail("publication_family_quote_invalid")
    return value


def _publication_families(address):
    sources = {row["source"].source_id: row["source"] for row in address["sources"]}
    return tuple(
        (row["source"], validate_family_series_source(row["source"], sources))
        for row in address["sources"]
        if row["binding"]["role"] == "family_series_ledger"
    )


def _mechanism_context(source, gene, disease):
    if (
        source.provider not in {"clingen", "genereviews"}
        or source.kind != "gene_disease"
        or source.genes != (gene,)
        or source.variant_keys
        or not disease.get("identifier")
        or disease.get("inheritance") == "unknown"
    ):
        return False
    fact = _ranking_json(source)
    if fact.get("gene") != gene:
        return False
    if fact.get("record_type") == "clingen_dosage":
        row = fact.get("row", {})
        chromosome = str(row.get("Genomic Location", "")).split(":", 1)[0].removeprefix("chr")
        mode = (
            "XL"
            if chromosome == "X"
            else "AD"
            if chromosome.isdigit() and 1 <= int(chromosome) <= 22
            else None
        )
        return (
            row.get("Gene Symbol") == gene
            and row.get("Haploinsufficiency Disease ID") == disease["identifier"]
            and row.get("Haploinsufficiency Score") == "3"
            and disease["inheritance"] == mode
        )
    return bool(
        (
            fact.get("disease_id") == disease["identifier"]
            or source.provider == "genereviews"
            and disease["identifier"] in {f"OMIM:{value}" for value in fact.get("disease_ids", ())}
        )
        and fact.get("inheritance") == disease["inheritance"]
        and re.search(
            r"haploinsufficien|loss[- ]of[- ]function|gain[- ]of[- ]function|dominant[- ]negative|mechanism of disease",
            fact.get("evidence_summary", ""),
            re.I,
        )
    )


def _observation_quotes(source, kind):
    if kind == "molecular_mechanism_context":
        context = _ranking_json(source).get("tool_result", {}).get("mechanism_context")
        text = _canonical_json(context) if isinstance(context, dict) else ""
        if not text or text not in source.content or len(text) > 6000:
            _fail("molecular_context_original_excerpt_invalid")
        return (text,)
    if kind != "gene_disease_mechanism" or len(source.content) <= 6000:
        return (source.content,)
    fact = _ranking_json(source)
    text = fact.get("mechanism_original") or fact.get("evidence_summary", "")
    sentences = re.split(r"(?<=[.!?])\s+", text)
    quotes = tuple(
        dict.fromkeys(
            json.dumps(s, ensure_ascii=False)[1:-1]
            for s in sentences
            if re.search(
                r"haploinsufficien|loss[- ]of[- ]function|gain[- ]of[- ]function|dominant[- ]negative|mechanism",
                s,
                re.I,
            )
        )
    )
    if not quotes or sum(map(len, quotes)) > 12000 or any(q not in source.content for q in quotes):
        _fail("mechanism_original_excerpt_unavailable")
    return quotes


def _retained_observations(address, disease=None):
    """Carry measured facts across LLM summarization without assigning their weight."""
    if not address["target"].get("post_evidence_recovered"):
        return ()
    output = []
    seen = set()
    for row in sorted(
        address["sources"], key=lambda row: row["binding"]["role"] in {"gene_context", "gene_disease_context"}
    ):
        source, role = row["source"], row["binding"]["role"]
        mechanism = (
            disease
            and role in {"gene_context", "gene_disease_context"}
            and _mechanism_context(source, address["target"]["primary"]["gene"], disease)
        )
        molecular = (
            role == "pvs1"
            and source.provider == "autopvs1"
            and bool(_ranking_json(source).get("tool_result", {}).get("mechanism_context"))
        )
        if (
            not (
                role == "prediction"
                and source.provider in {"revel", "spliceai"}
                or role == "family_series_ledger"
                or mechanism
                or molecular
            )
            or source.source_id in seen
        ):
            continue
        seen.add(source.source_id)
        kind = (
            "gene_disease_mechanism"
            if mechanism
            else "molecular_mechanism_context"
            if molecular
            else "publication_family_series"
            if role == "family_series_ledger"
            else "computational_prediction"
        )
        output.append(
            {
                "id": f"O{len(output) + 1}",
                "kind": kind,
                "finding": _OBSERVATION_LABELS[source.provider],
                "citations": tuple(
                    _citation(source.source_id, quote) for quote in _observation_quotes(source, kind)
                ),
            }
        )
    return tuple(output)


def _profile_source_scope(source, gene, key, *, allow_related=False):
    return (not source.genes or gene in source.genes) and (
        not source.variant_keys
        or key in source.variant_keys
        or allow_related
        and source.kind == "patient_observation"
        and source.genes == (gene,)
    )


def _support_source(row, target):
    if _identity(target)[1][1] == "MT" and not (
        target["reference"].get("status") == "sequence_verified"
        and target["reference"].get("reference_accession") == "NC_012920.1"
    ):
        return False
    binding = row["binding"]
    return binding["role"] not in {
        "gene_context", "gene_disease_context", "patient_phenotype", "guideline", "prior_criterion_source"
    } and not (binding["role"] == "literature" and binding["match_scope"] == "context")


def _profile_quotes(rows, sources, allowed, *, required=True):
    rows = _list(rows, "profile_citations_required")
    if required and not rows:
        _fail("profile_citation_count_invalid")
    citations = []
    for row in rows:
        if not isinstance(row, dict) or set(row) != {"source_id", "quote"}:
            _fail("profile_citation_schema_invalid")
        source = sources.get(row["source_id"])
        if source is None or source.source_id not in allowed:
            _fail("profile_source_scope_invalid")
        quote = _text(row["quote"], "profile_quote_required")
        if quote not in source.content:
            _fail(f"profile_quote_not_in_source:{row['source_id']}")
        citations.append(_citation(source.source_id, quote))
    return tuple(citations)


def _decode_profile(raw, address, visible, anchor):
    fields = {"allele_disease", "patient_fit", "key_support", "key_limits", "unknowns", "summary"}
    if set(raw) != fields:
        _fail("profile_must_not_classify_or_score")
    gene, key = _identity(address["target"])
    direct = {source.source_id for source in visible.values() if _profile_source_scope(source, gene, key)}
    related = {
        source.source_id
        for source in visible.values()
        if _profile_source_scope(source, gene, key, allow_related=True)
    }
    disease, fit = raw["allele_disease"], raw["patient_fit"]
    if not isinstance(disease, dict) or set(disease) != {
        "label",
        "identifier",
        "inheritance",
        "reason",
        "citations",
    }:
        _fail("profile_disease_schema_invalid")
    if disease["inheritance"] not in _PROFILE_INHERITANCE:
        _fail("profile_inheritance_invalid")
    disease = {
        **disease,
        "label": _text(disease["label"], "profile_disease_required"),
        "identifier": _text(disease["identifier"], "profile_disease_id_invalid", allow_empty=True),
        "reason": _text(disease["reason"], "profile_disease_reason_required"),
        "citations": _profile_quotes(disease["citations"], visible, direct),
    }
    if (
        not isinstance(fit, dict)
        or set(fit) != {"status", "reason", "citations"}
        or fit["status"] not in _FIT_STATUSES
    ):
        _fail("profile_patient_fit_invalid")
    fit = {
        **fit,
        "reason": _text(fit["reason"], "profile_patient_fit_reason_required"),
        "citations": _profile_quotes(fit["citations"], visible, related, required=False),
    }
    support, limits = (
        _list(raw["key_support"], "profile_support_required"),
        _list(raw["key_limits"], "profile_limits_required"),
    )
    binding = anchor.get("disease_binding") if anchor is not None else None
    if binding and any(disease.get(name) != binding.get(name) for name in ("identifier", "inheritance")):
        if any(isinstance(item, dict) and item.get("kind") == "pvs1" for item in support):
            _fail("pvs1_disease_binding_mismatch")
        anchor = None
    if len(support) > 6 or len(limits) > 3:
        _fail("profile_size_invalid")
    normalized = []
    audited = {fact["article_source_id"] for _, fact in _publication_families(address)}
    for index, item in enumerate(support, 1):
        if not isinstance(item, dict) or set(item) != {
            "kind",
            "finding",
            "reliability",
            "pathogenic_relevance",
            "independence_group",
            "citations",
        }:
            _fail("profile_support_schema_invalid")
        if (
            item["kind"] not in _PROFILE_SUPPORT_KINDS
            or item["reliability"] not in _PROFILE_RELIABILITY
            or item["pathogenic_relevance"] not in _PROFILE_RELEVANCE
        ):
            _fail("profile_support_kind_invalid")
        if item["kind"] == "pvs1" and address["target"].get("post_evidence_recovered") and anchor is None:
            _fail("pvs1_requires_verified_disease_anchor")
        citations = _profile_quotes(item["citations"], visible, direct)
        cited = {row["source_id"] for row in citations}
        refs = tuple(row for row in address["sources"] if row["source"].source_id in cited)
        if (anchor is None or item["kind"] != "pvs1") and any(
            row["source"].provider == "autopvs1"
            and _ranking_json(row["source"]).get("tool_result", {}).get("mechanism_context")
            for row in refs
        ):
            _fail("unvalidated_native_context_is_not_independent_support")
        if not any(_support_source(row, address["target"]) for row in refs):
            _fail("gene_background_is_not_independent_allele_support")
        if item["kind"] == "pvs1" and anchor is not None:
            if not {citation["source_id"] for citation in citations}.intersection(
                citation["source_id"] for citation in anchor["citations"]
            ):
                _fail("validated_pvs1_anchor_source_required")
        elif item["kind"] in {"pvs1", "functional", "case_series"}:
            _eligibility(
                {"pvs1": "PVS1", "functional": "PS3", "case_series": "PS4"}[item["kind"]], address, refs
            )
        finding = _text(item["finding"], "profile_support_finding_required")
        papers = {r["source"].source_id for r in refs if r["binding"]["role"] == "literature"}
        family_de_novo = item["kind"] == "family" and re.search(r"\bde[- ]novo\b", finding, re.I)
        published_family = (
            family_de_novo
            and any(r["binding"]["role"] in {"literature", "family_series_ledger"} for r in refs)
            and all(
                r["binding"]["role"] in {"literature", "family_series_ledger", "allele_identity"}
                and (
                    r["binding"]["match_scope"] == "exact_allele"
                    or r["binding"]["role"] == "literature"
                    and r["binding"]["match_scope"] == "candidate_allele"
                    and any(
                        c["source_id"] == r["source"].source_id
                        and _mentions_allele(address["target"], c["quote"])
                        and re.search(r"\bde[- ]novo\b", c["quote"], re.I)
                        for c in citations
                    )
                )
                for r in refs
            )
        )
        if family_de_novo and not published_family:
            _eligibility("PS2" if _relationship_qc(address) else "PM6", address, refs)
        row = {
            **item,
            "id": f"E{index}",
            "finding": finding,
            "independence_group": _text(item["independence_group"], "profile_independence_required"),
            "citations": citations,
        }
        if published_family or (
            row["kind"] == "case_series"
            and row["reliability"] != "limited"
            and not (len(papers) == 1 and papers <= audited)
        ):
            row.update(
                reported_reliability=row["reliability"],
                reliability="limited",
                normalization_note=(
                    "Published-family observation only; limited reliability. This does not establish "
                    "current-patient de novo origin or verified PS2/PM6 evidence."
                    if published_family
                    else "No independently audited external-family ledger is supplied."
                ),
            )
        normalized.append(row)
    normalized_limits = []
    for index, item in enumerate(limits, 1):
        if not isinstance(item, dict) or set(item) != {"finding", "effect", "citations"}:
            _fail("profile_limit_schema_invalid")
        normalized_limits.append(
            {
                "id": f"L{index}",
                "finding": _text(item["finding"], "profile_limit_required"),
                "effect": _text(item["effect"], "profile_limit_effect_required"),
                "citations": _profile_quotes(item["citations"], visible, related),
            }
        )
    if anchor is not None and not any(item["kind"] == "pvs1" for item in normalized):
        _fail("validated_pvs1_anchor_omitted")
    observations = _retained_observations(address, disease)
    return {
        "identity": {"gene": gene, "variant_key": key},
        "allele_disease": disease,
        "patient_fit": fit,
        "key_support": tuple(normalized),
        "key_limits": tuple(normalized_limits),
        "unknowns": _unique_strings(raw["unknowns"], "profile_unknowns_invalid"),
        "summary": _text(raw["summary"], "profile_summary_required"),
        "validated_pvs1_anchor": anchor,
        "input_source_ids": tuple(source.source_id for source in visible.values()),
        **({"retained_observations": observations} if observations else {}),
    }


def _decode_profile_classification(raw, profile):
    fields = {
        "assessed_disease_label",
        "classification",
        "confidence",
        "key_evidence_ids",
        "limitation_ids",
        "rationale",
        "patient_fit",
        "summary",
    }
    observations = profile.get("retained_observations", ())
    if observations:
        fields.add("observation_review")
    if set(raw) != fields:
        _fail("profile_classification_schema_invalid")
    if raw["classification"] not in _CLASSES or raw["confidence"] not in {"high", "moderate", "low"}:
        _fail("profile_classification_invalid")
    if raw["assessed_disease_label"] != profile["allele_disease"]["label"]:
        _fail("profile_disease_changed")
    support = {item["id"]: item for item in profile["key_support"]}
    limits = {item["id"]: item for item in profile["key_limits"]}
    used = _unique_strings(raw["key_evidence_ids"], "profile_evidence_ids_invalid")
    limited = _unique_strings(raw["limitation_ids"], "profile_limit_ids_invalid")
    if not set(used) <= set(support) or not set(limited) <= set(limits):
        _fail("profile_unknown_evidence_id")
    if (support or raw["classification"] in _POSITIVE_CLASSES) and not used:
        _fail("profile_decisive_support_required")
    decisive = tuple(dict.fromkeys(cite["source_id"] for name in used for cite in support[name]["citations"]))
    context = (profile["allele_disease"], profile["patient_fit"], *(limits[name] for name in limited))
    contextual = tuple(dict.fromkeys(cite["source_id"] for item in context for cite in item["citations"]))
    review = ()
    if observations:
        review = tuple(_list(raw["observation_review"], "observation_review_required"))
        if (
            any(not isinstance(r, dict) for r in review)
            or len(review) != len(observations)
            or {r.get("id") for r in review} != {r["id"] for r in observations}
            or any(
                set(r) != {"id", "status", "reason"}
                or r["status"] not in {"considered", "not_applicable"}
                or not isinstance(r["reason"], str)
                or not r["reason"].strip()
                for r in review
            )
        ):
            _fail("observation_review_coverage_invalid")
        by_id = {r["id"]: r for r in review}
        review = tuple(by_id[r["id"]] for r in observations)
        contextual = tuple(
            dict.fromkeys((*contextual, *(c["source_id"] for r in observations for c in r["citations"])))
        )
    return {
        "identity": profile["identity"],
        "disease": {name: profile["allele_disease"][name] for name in ("label", "identifier", "inheritance")},
        **{name: raw[name] for name in ("classification", "confidence")},
        **{
            name: _text(raw[name], "profile_final_" + name + "_required")
            for name in ("rationale", "patient_fit", "summary")
        },
        "key_evidence_ids": used,
        "limitation_ids": limited,
        "decisive_source_ids": decisive,
        "contextual_source_ids": contextual,
        **({"observation_review": review} if observations else {}),
    }


def _validate_profiles(result, *, allow_pending=False):
    if set(result) != {"workflow", "profiles", "final", "sources"}:
        raise MendelPrimeError("schema_mismatch", field="assessment_result")
    if (
        any(type(result[name]) is not tuple for name in ("profiles", "final", "sources"))
        or not result["profiles"]
        and not allow_pending
    ):
        raise MendelPrimeError("nonempty_profile_tuple_required", field="assessment_result")
    registry = {source.source_id: source for source in result["sources"] if isinstance(source, Source)}
    if len(registry) != len(result["sources"]):
        raise MendelPrimeError("source_registered_more_than_once", field="assessment_result.sources")
    if any(parent not in registry for source in registry.values() for parent in source.parent_source_ids):
        raise MendelPrimeError("parent_source_missing", field="assessment_result.sources")
    if len({source.case_key for source in registry.values() if source.case_key is not None}) != 1:
        raise MendelPrimeError("case_scope_mismatch", field="assessment_result.sources")
    identities = tuple(_identity_key(profile["identity"]) for profile in result["profiles"])
    if (
        len(set(identities)) != len(identities)
        or tuple(_identity_key(row["identity"]) for row in result["final"])
        != tuple(
            key for key in identities if key in {_identity_key(row["identity"]) for row in result["final"]}
        )
        or (not allow_pending and len(result["final"]) != len(identities))
    ):
        raise MendelPrimeError("profile_coverage_or_order_invalid", field="assessment_result")
    profile_fields = {
        "identity",
        "allele_disease",
        "patient_fit",
        "key_support",
        "key_limits",
        "unknowns",
        "summary",
        "validated_pvs1_anchor",
        "input_source_ids",
    }
    decisions = {_identity_key(row["identity"]): row for row in result["final"]}
    for profile in result["profiles"]:
        decision = decisions.get(_identity_key(profile["identity"]))
        if set(profile) not in (profile_fields, profile_fields | {"retained_observations"}):
            raise MendelPrimeError("profile_must_not_classify_or_score", field="assessment_result.profiles")
        disease, fit = profile["allele_disease"], profile["patient_fit"]
        if (
            set(disease) != {"label", "identifier", "inheritance", "reason", "citations"}
            or disease["inheritance"] not in _PROFILE_INHERITANCE
            or set(fit) != {"status", "reason", "citations"}
            or fit["status"] not in _FIT_STATUSES
        ):
            raise MendelPrimeError("profile_context_invalid", field="assessment_result.profiles")
        for value in (disease["label"], disease["reason"], fit["reason"], profile["summary"]):
            _text(value, "profile_text_required")
        _text(disease["identifier"], "profile_disease_id_invalid", allow_empty=True)
        if type(profile["unknowns"]) is not tuple:
            _fail("profile_unknowns_invalid")
        _unique_strings(list(profile["unknowns"]), "profile_unknowns_invalid")
        gene, key = _identity_key(profile["identity"])
        allowed = profile["input_source_ids"]
        if (
            type(allowed) is not tuple
            or len(set(allowed)) != len(allowed)
            or not set(allowed) <= set(registry)
        ):
            raise MendelPrimeError("profile_source_missing", field="assessment_result.profiles")
        observations = profile.get("retained_observations", ())
        if type(observations) is not tuple:
            _fail("retained_observations_tuple_required")
        for index, row in enumerate(observations, 1):
            if (
                set(row) != {"id", "kind", "finding", "citations"}
                or row["id"] != f"O{index}"
                or row["kind"]
                not in {
                    "computational_prediction",
                    "publication_family_series",
                    "gene_disease_mechanism",
                    "molecular_mechanism_context",
                }
                or type(row["citations"]) is not tuple
                or not row["citations"]
            ):
                _fail("retained_observation_schema_invalid")
            source = registry.get(row["citations"][0]["source_id"])
            providers = (
                {"revel", "spliceai"}
                if row["kind"] == "computational_prediction"
                else {"clingen", "genereviews"}
                if row["kind"] == "gene_disease_mechanism"
                else {"autopvs1"}
                if row["kind"] == "molecular_mechanism_context"
                else {"publication-family-ledger"}
            )
            if (
                source is None
                or source.provider not in providers
                or row["finding"] != _OBSERVATION_LABELS[source.provider]
                or tuple(c["quote"] for c in row["citations"]) != _observation_quotes(source, row["kind"])
                or any(c["source_id"] != source.source_id for c in row["citations"])
                or source.source_id not in allowed
                or row["kind"] == "gene_disease_mechanism"
                and not _mechanism_context(source, gene, profile["allele_disease"])
                or row["kind"] == "molecular_mechanism_context"
                and (source.genes != (gene,) or key not in source.variant_keys)
            ):
                _fail("retained_observation_not_original")
            if source.provider == "publication-family-ledger":
                validate_family_series_source(source, registry)
        if (
            type(profile["key_support"]) is not tuple
            or type(profile["key_limits"]) is not tuple
            or len(profile["key_support"]) > 6
            or len(profile["key_limits"]) > 3
        ):
            raise MendelPrimeError("profile_size_invalid", field="assessment_result.profiles")
        for prefix, rows in (("E", profile["key_support"]), ("L", profile["key_limits"])):
            if tuple(row["id"] for row in rows) != tuple(
                f"{prefix}{index}" for index in range(1, len(rows) + 1)
            ):
                raise MendelPrimeError("profile_evidence_ids_invalid", field="assessment_result.profiles")
            for row in rows:
                fields = {"id", "finding", "citations"} | (
                    {"kind", "reliability", "pathogenic_relevance", "independence_group"}
                    if prefix == "E"
                    else {"effect"}
                )
                extras = {"reported_reliability", "normalization_note"} if prefix == "E" else set()
                if not fields <= set(row) or set(row) - fields - extras:
                    _fail("profile_fact_schema_invalid")
                _text(row["finding"], "profile_finding_required")
                if prefix == "E":
                    if (
                        row["kind"] not in _PROFILE_SUPPORT_KINDS
                        or row["reliability"] not in _PROFILE_RELIABILITY
                        or row["pathogenic_relevance"] not in _PROFILE_RELEVANCE
                    ):
                        _fail("profile_support_kind_invalid")
                    _text(row["independence_group"], "profile_independence_required")
                else:
                    _text(row["effect"], "profile_limit_effect_required")
        for item in (
            profile["allele_disease"],
            profile["patient_fit"],
            *profile["key_support"],
            *profile["key_limits"],
            *observations,
        ):
            if type(item["citations"]) is not tuple or (item is not fit and not item["citations"]):
                _fail("profile_citation_count_invalid")
            for citation in item["citations"]:
                if set(citation) != {"source_id", "quote"}:
                    _fail("profile_citation_schema_invalid")
                _text(citation["quote"], "profile_quote_required")
                source = registry.get(citation["source_id"])
                if (
                    source is None
                    or source.source_id not in allowed
                    or citation["quote"] not in source.content
                ):
                    raise MendelPrimeError(
                        "citation_source_or_quote_invalid", field="assessment_result.profiles"
                    )
                allow_related = item is fit or any(item is limit for limit in profile["key_limits"])
                if not _profile_source_scope(source, gene, key, allow_related=allow_related):
                    raise MendelPrimeError("source_scope_mismatch", field="assessment_result.profiles")
        anchor = profile["validated_pvs1_anchor"]
        if anchor is not None:
            binding = anchor.get("disease_binding")
            if binding and any(
                profile["allele_disease"].get(name) != binding.get(name)
                for name in ("identifier", "inheritance")
            ):
                _fail("pvs1_disease_binding_mismatch")
            if (
                anchor.get("code") != "PVS1"
                or anchor.get("status") != "met"
                or anchor.get("strength") not in _PVS1_STRENGTHS
                or not anchor.get("citations")
                or not any(item["kind"] == "pvs1" for item in profile["key_support"])
            ):
                _fail("validated_pvs1_anchor_invalid")
            for citation in anchor["citations"]:
                source = registry.get(citation["source_id"])
                if (
                    source is None
                    or source.source_id not in allowed
                    or not citation["quote"]
                    or citation["quote"] not in source.content
                    or (source.genes and gene not in source.genes)
                    or (source.variant_keys and key not in source.variant_keys)
                ):
                    _fail("validated_pvs1_source_invalid")
        if decision is None:
            continue
        raw = {
            "assessed_disease_label": decision["disease"]["label"],
            **{
                name: decision[name]
                for name in ("classification", "confidence", "rationale", "patient_fit", "summary")
            },
            "key_evidence_ids": list(decision["key_evidence_ids"]),
            "limitation_ids": list(decision["limitation_ids"]),
            **({"observation_review": list(decision["observation_review"])} if observations else {}),
        }
        if _decode_profile_classification(raw, profile) != decision:
            raise MendelPrimeError("profile_final_projection_mismatch", field="assessment_result.final")
    return result


_RANKING_RULES = "September-AI-CURA-2015-combinations+approved-quality-PM3+PM2-AF-20261022"
_RANKING_POINTS = {"supporting": 1, "moderate": 2, "strong": 4, "very_strong": 8, "stand_alone": 8}
_LITERATURE_CODES = {"PVS1_RNA", "PS2", "PM6", "PS3", "PS4", "PM3", "PP1", "PP4"}


def _ranking_packet(prepared):
    packet = prepared.get("input", prepared)
    packet = {name: packet[name] for name in ("profile", "family", "targets", "sources")}
    if packet["targets"]:
        validate_input(
            {
                **packet,
                "targets": tuple(
                    {**target, "predictions": target.get("predictions") or {}} for target in packet["targets"]
                ),
            }
        )
        registry = {source.source_id: source for source in packet["sources"]}
        for target in packet["targets"]:
            gene, key = _identity(target)
            for native in (target["primary"], *target.get("aliases", ())):
                source = native["source"]
                if (
                    registry.get(source.source_id) != source
                    or gene not in source.genes
                    or key not in source.variant_keys
                ):
                    raise MendelPrimeError("native_identity_source_mismatch", field="ranking_acmg.input")
    else:
        if type(packet["targets"]) is not tuple or type(packet["sources"]) is not tuple:
            raise MendelPrimeError("tuple_required", field="ranking_acmg.input")
        if any(not isinstance(source, Source) for source in packet["sources"]):
            raise MendelPrimeError("source_required", field="ranking_acmg.sources")
        registry = {source.source_id: source for source in packet["sources"]}
        if len(registry) != len(packet["sources"]) or any(
            source.case_key is not None and source.case_key != case_tuple(packet["profile"]["case_key"])
            for source in packet["sources"]
        ):
            raise MendelPrimeError("source_identity_or_case_mismatch", field="ranking_acmg.sources")
        if any(parent not in registry for source in packet["sources"] for parent in source.parent_source_ids):
            raise MendelPrimeError("parent_source_missing", field="ranking_acmg.sources")
    return packet


def _ranking_json(source):
    try:
        value = json.loads(source.content)
    except (ValueError, RecursionError):
        return {}
    return value if isinstance(value, dict) else {}


def _ranking_number(value):
    return value if type(value) in (int, float) and math.isfinite(value) else None


def _ranking_mode(value):
    normalized = re.sub(r"[^A-Z0-9]+", "_", str(value).upper()).strip("_")
    if normalized in {"AR", "AUTOSOMAL_RECESSIVE", "AUTOSOMAL_RECESSIVE_COMP_HET"}:
        return "AR"
    if normalized in {"AD", "AUTOSOMAL_DOMINANT"}:
        return "AD"
    if normalized in {"XL", "X_LINKED", "X_LINKED_RECESSIVE", "X_LINKED_DOMINANT"}:
        return "XL"
    return normalized or "unknown"


def _ranking_contexts(target, sources):
    rows = []
    for row in target.get("disease_contexts", ()):
        ids = tuple(row.get("source_ids", ()))
        if any(identifier not in sources for identifier in ids):
            raise MendelPrimeError("disease_context_source_missing", field="ranking_acmg.contexts")
        if any(
            sources[identifier].genes and _identity(target)[0] not in sources[identifier].genes
            for identifier in ids
        ):
            raise MendelPrimeError("disease_context_gene_mismatch", field="ranking_acmg.contexts")
        item = {
            "disease_id": row.get("disease_id") or "",
            "disease_label": row.get("disease_label", row.get("disease_name")) or "",
            "mode_of_inheritance": _ranking_mode(row.get("mode_of_inheritance", row.get("moi")) or "unknown"),
            "source_ids": ids,
        }
        same = next(
            (
                current
                for current in rows
                if all(
                    current[name] == item[name]
                    for name in ("disease_id", "disease_label", "mode_of_inheritance")
                )
            ),
            None,
        )
        if same is None:
            rows.append(item)
        else:
            same["source_ids"] = tuple(dict.fromkeys((*same["source_ids"], *ids)))
    return tuple(rows)


def _ranking_criterion(
    code, strength, refs, *, method, group, status="met", facts=None, limitations=(), origin="deterministic"
):
    direction = "benign" if code.startswith(("BA", "BS", "BP")) else "pathogenic"
    points = _RANKING_POINTS.get(strength, 0) * (-1 if direction == "benign" else 1) if status == "met" else 0
    citations = tuple(_citation(source.source_id, source.content) for source in dict.fromkeys(refs))
    return {
        "code": code,
        "strength": strength,
        "direction": direction,
        "status": status,
        "points": points,
        "source_refs": tuple(row["source_id"] for row in citations),
        "excerpts": tuple(row["quote"] for row in citations),
        "citations": citations,
        "method": method,
        "independence_group": group,
        "origin": origin,
        "raw_facts": dict(facts or {}),
        "limitations": tuple(limitations),
        "counted": status == "met",
    }


def _ranking_with_citations(criterion, citations):
    citations = tuple(citations)
    output = dict(criterion)
    output.update(
        source_refs=tuple(row["source_id"] for row in citations),
        excerpts=tuple(row["quote"] for row in citations),
        citations=citations,
    )
    return output


def _ranking_select(criteria):
    # Retain all source facts, but count one code and one dependent evidence group.
    output = [{**row, "counted": False} for row in criteria]
    selected_codes, selected_groups = set(), set()
    ordered = sorted(
        range(len(output)),
        key=lambda index: (
            -abs(output[index]["points"]),
            output[index]["origin"] != "expert_assertion",
            output[index]["origin"] != "literature_aggregate",
            index,
        ),
    )
    for index in ordered:
        item = output[index]
        group = (item["direction"], item["independence_group"])
        if item["status"] == "met" and item["code"] not in selected_codes and group not in selected_groups:
            item["counted"] = True
            selected_codes.add(item["code"])
            selected_groups.add(group)
    return tuple(output)


def _ranking_classify(criteria):
    """Classify met criteria by strength combinations, not audit-point totals."""
    met = [item for item in criteria if item["status"] == "met" and item.get("counted", True)]
    pathogenic = [item for item in met if item["direction"] == "pathogenic"]
    benign = [item for item in met if item["direction"] == "benign"]
    if pathogenic and benign:
        return "uncertain_significance"
    if any(item["code"] == "BA1" for item in benign):
        return "benign"
    bs = sum(item["strength"] in {"strong", "very_strong", "stand_alone"} for item in benign)
    bp = sum(item["strength"] == "supporting" for item in benign)
    if bs >= 2:
        return "benign"
    if (bs >= 1 and bp >= 1) or bp >= 2:
        return "likely_benign"
    vs, s, m, p = (
        sum(item["strength"] == strength for item in pathogenic)
        for strength in ("very_strong", "strong", "moderate", "supporting")
    )
    if (
        (vs >= 1 and s >= 1)
        or (vs >= 1 and m >= 2)
        or (vs >= 1 and m >= 1 and p >= 1)
        or (vs >= 1 and p >= 2)
        or s >= 2
        or (s >= 1 and m >= 3)
        or (s >= 1 and m >= 2 and p >= 2)
        or (s >= 1 and m >= 1 and p >= 4)
    ):
        return "pathogenic"
    if (
        (vs >= 1 and m >= 1)
        or (s >= 1 and 1 <= m <= 2)
        or (s >= 1 and p >= 2)
        or m >= 3
        or (m >= 2 and p >= 2)
        or (m >= 1 and p >= 4)
    ):
        return "likely_pathogenic"
    return "uncertain_significance"


def _ranking_summary(variant, criteria):
    criteria = _ranking_select(criteria)
    met = [item for item in criteria if item["counted"]]
    return {
        **variant,
        "criteria": criteria,
        "classification": _ranking_classify(criteria),
        "points_for_audit_only": sum(item["points"] for item in met),
        "conflicting_evidence": bool(
            any(item["direction"] == "pathogenic" for item in met)
            and any(item["direction"] == "benign" for item in met)
        ),
    }


def _ranking_independent(variant):
    without_pm3 = _ranking_summary(variant, tuple(row for row in variant["criteria"] if row["code"] != "PM3"))
    return {**variant, "independent_classification": without_pm3["classification"]}


def _ranking_bound(target, sources):
    return tuple({"source": sources[item["source_id"]], "binding": item} for item in target["bindings"])


def _ranking_consequence(target):
    annotation = target.get("annotation") or {}
    selected = annotation.get("annotation") or {}
    if annotation.get("status") == "annotated":
        return tuple(selected.get("consequence_terms", ())), selected.get("source")
    details = target["primary"].get("details", {})
    consequence = details.get("consequence") or details.get("variant_effect") or ""
    return tuple(value.strip().lower() for value in re.split(r"[,&|]", consequence)) if isinstance(
        consequence, str
    ) else (), target["primary"]["source"]


def _ranking_population(target):
    value = target.get("population") or {}
    source = value.get("source")
    if source is None:
        return []
    if value.get("status") == "not_found_covered" and value.get("coverage_adequate") is True:
        return [
            _ranking_criterion(
                "PM2",
                "supporting",
                (source,),
                method="coverage-qualified population absence",
                group="population_frequency",
            )
        ]
    if value.get("status") != "found" or value.get("record_filter") not in {"PASS", "."}:
        return []
    af = _ranking_number(value.get("allele_frequency"))
    if (
        af is None
        or not 0 <= af <= 1
        or any(_ranking_number(value.get(name)) is None for name in ("allele_count", "allele_number"))
        or value["allele_number"] <= 0
    ):
        return []
    values = [af]
    popmax = _ranking_number(value.get("popmax_allele_frequency"))
    if value.get("popmax_allele_frequency") is not None and (popmax is None or not 0 <= popmax <= 1):
        return []
    if popmax is not None:
        values.append(popmax)
    maximum = max(values)
    code, strength = (
        ("BA1", "stand_alone")
        if maximum > 0.05
        else ("BS1", "strong")
        if maximum > 0.01
        else ("PM2", "supporting")
        if maximum <= 0.005
        else ("", "")
    )
    return (
        [
            _ranking_criterion(
                code,
                strength,
                (source,),
                method="September AI-CURA population AF rule",
                group="population_frequency",
                facts={"maximum_allele_frequency": maximum},
                limitations=(
                    "Authorized 2026-10-22: measured AF <=0.005 restores PM2 supporting; covered absence remains separately qualified.",
                )
                if code == "PM2"
                else (),
            )
        ]
        if code
        else []
    )


def _ranking_predictions(target, terms):
    predictions = target.get("predictions") or {}
    revel = predictions.get("revel") or {}
    transcript = target["primary"].get("transcript")
    records = [
        row
        for row in revel.get("records", ())
        if row.get("transcript") == transcript and _ranking_number(row.get("score")) is not None
    ]
    if "missense_variant" in terms or "missense" in terms:
        scores = {row["score"] for row in records}
        if transcript and revel.get("status") == "found" and revel.get("source") and len(scores) == 1:
            score = next(iter(scores))
            intervals = (
                (score <= 0.003, "BP4", "very_strong"),
                (0.003 < score <= 0.016, "BP4", "strong"),
                (0.016 < score <= 0.183, "BP4", "moderate"),
                (0.183 < score <= 0.290, "BP4", "supporting"),
                (0.644 <= score < 0.73, "PP3", "supporting"),
                (0.773 <= score < 0.932, "PP3", "moderate"),
                (score >= 0.932, "PP3", "strong"),
            )
            for applies, code, strength in intervals:
                if applies:
                    return [
                        _ranking_criterion(
                            code,
                            strength,
                            (revel["source"],),
                            method="September calibrated REVEL intervals",
                            group="computational",
                            facts={"score": score},
                        )
                    ]
    splice = predictions.get("spliceai") or {}
    rows = [row for row in splice.get("records", ()) if row.get("gene") == _identity(target)[0]]
    scores = [_ranking_number(row.get(name)) for row in rows for name in ("ds_ag", "ds_al", "ds_dg", "ds_dl")]
    if (
        splice.get("status") != "found"
        or not splice.get("source")
        or not scores
        or any(score is None for score in scores)
    ):
        return []
    score = max(scores)
    if score >= 0.2:
        return [
            _ranking_criterion(
                "PP3",
                "supporting",
                (splice["source"],),
                method="September SpliceAI >=0.2",
                group="computational",
                facts={"score": score},
            )
        ]
    if "synonymous_variant" in terms or "synonymous" in terms:
        return [
            _ranking_criterion(
                "BP7",
                "supporting",
                (splice["source"],),
                method="September synonymous SpliceAI <0.2",
                group="computational",
                facts={"score": score},
            )
        ]
    return []


def _ranking_specificity(packet, target, refs):
    contexts = target.get("disease_contexts", ())
    label = contexts[0].get("disease_label", "") if len(contexts) == 1 else ""
    text = re.sub(
        r"[^a-z0-9]+",
        " ",
        " ".join((label, *(row["label"] for row in packet["profile"]["terms"] if row["present"]))).casefold(),
    ).strip()
    candidates = []
    for item in refs:
        source = item["source"]
        row = _ranking_json(source)
        row = row.get("row", row)
        if not isinstance(row, dict) or not {"specificity", "genes", "disease_phenotypes"} <= row.keys():
            continue
        genes = {value.strip().upper() for value in re.split(r"[,;|]", row["genes"])}
        if _identity(target)[0] not in genes and "ANY" not in genes:
            continue
        rule = re.sub(r"[^a-z0-9]+", " ", row["disease_phenotypes"].casefold()).strip()
        tokens = {word for word in rule.split() if len(word) > 3}
        if rule and (rule in text or text in rule or len(set(text.split()) & tokens) >= min(2, len(tokens))):
            return row["specificity"], (source,)
        candidates.append(source)
    if candidates:
        return "Consistent with gene but not highly specific and high genetic heterogeneity", tuple(
            candidates
        )
    return "", ()


def _ranking_denovo_weight(specificity, confirmed, unconfirmed):
    label = specificity.casefold().strip()
    if label == "highly specific":
        return confirmed * 2.0 + unconfirmed
    if label == "consistent with gene but not highly specific":
        return confirmed + unconfirmed * 0.5
    if "high genetic heterogeneity" in label:
        return confirmed * 0.5 + unconfirmed * 0.25
    return 0.0


def _ranking_strength(score):
    return next(
        (
            strength
            for threshold, strength in (
                (4, "very_strong"),
                (2, "strong"),
                (1, "moderate"),
                (0.5, "supporting"),
            )
            if score >= threshold
        ),
        "unassigned",
    )


def _ranking_family(packet, target, refs):
    address = {"target": target, "sources": refs}
    calls = _calls(address)
    if (
        target["observation"]["origin"] != "candidate_de_novo"
        or target.get("family_conflict")
        or not all(role in calls and calls[role]["usable"] for role in ("proband", "father", "mother"))
    ):
        return []
    code = "PS2" if _relationship_qc(address) else "PM6"
    try:
        _eligibility(code, address, refs)
    except MendelPrimeError:
        return []
    specificity, rules = _ranking_specificity(packet, target, refs)
    confirmed = int(code == "PS2")
    score = _ranking_denovo_weight(specificity, confirmed, 1 - confirmed)
    sources = tuple(item["source"] for item in refs if item["binding"]["role"] in _MEASUREMENT_ROLES) + rules
    return [
        _ranking_criterion(
            code,
            _ranking_strength(score),
            sources,
            method="September Table S8 + ClinGen PS2/PM6 points",
            group="de_novo",
            status="met" if score >= 0.5 else "unknown",
            facts={
                "de_novo_score": score,
                "de_novo_confirmed_count": confirmed,
                "de_novo_unconfirmed_count": 1 - confirmed,
                "phenotypic_consistency": specificity,
                "subject_scope": "current_case",
                "subject_group": ":".join(case_tuple(packet["profile"]["case_key"])),
            },
        )
    ]


def _ranking_assertions(target, refs, contexts):
    criteria, assertions = [], []
    context = contexts[0] if len(contexts) == 1 else None
    for item in refs:
        source, binding = item["source"], item["binding"]
        if binding["role"] not in {"allele_assertion", "protein_comparator"}:
            continue
        payload = _ranking_json(source)
        if payload.get("record_type") == "verified_sequence_equivalent_clinical_evidence":
            try:
                payload = json.loads(payload["original_content"])
            except (KeyError, ValueError):
                continue
        row = payload.get("row", payload.get("record", payload))
        if not isinstance(row, dict):
            continue
        classification = str(row.get("Assertion", row.get("clinical_significance", "")))
        assertion = {
            "source_id": source.source_id,
            "dependency_group": binding["dependency_group"],
            "classification": classification,
            "expert_panel": str(row.get("Expert Panel", "")),
            "disease_id": str(row.get("Mondo Id", "")),
            "mode_of_inheritance": _ranking_mode(row.get("Mode of Inheritance", "unknown")),
            "conflicting_evidence": "conflict" in classification.casefold(),
            "retracted": False,
        }
        if payload.get("retraction_status") == "retracted" or str(row.get("Retracted", "")).casefold() in {
            "true",
            "1",
            "yes",
        }:
            assertions.append({**assertion, "status": "retracted", "retracted": True})
            continue
        assertions.append({**assertion, "status": "external_assertion"})
        if source.provider == "clingen" and binding["match_scope"] == "exact_allele":
            matches = bool(
                context
                and context["disease_id"]
                and context["disease_id"].upper()
                in re.split(r"[,;|\s]+", str(row.get("Mondo Id", "")).upper())
                and _ranking_mode(row.get("Mode of Inheritance")) == context["mode_of_inheritance"]
            )
            for raw in re.split(r"[,;|\s]+", str(row.get("Applied Evidence Codes (Met)", ""))):
                match = re.fullmatch(
                    r"(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])(?:_(Supporting|Moderate|Strong|VeryStrong|Very_Strong))?",
                    raw,
                    re.I,
                )
                if not match or match[1].upper() in {"PP5", "BP6", "PM3"}:
                    continue
                code = match[1].upper()
                strength = (
                    (
                        match[2]
                        or (
                            "VeryStrong"
                            if code == "PVS1"
                            else "Strong"
                            if code.startswith(("PS", "BS"))
                            else "Moderate"
                            if code.startswith("PM")
                            else "stand_alone"
                            if code == "BA1"
                            else "Supporting"
                        )
                    )
                    .lower()
                    .replace("verystrong", "very_strong")
                )
                criteria.append(
                    _ranking_criterion(
                        code,
                        strength if matches else "unassigned",
                        (source,),
                        method="ClinGen VCEP expert assertion",
                        group=f"{binding['dependency_group']}:{code}",
                        status="met" if matches else "unknown",
                        origin="expert_assertion",
                        limitations=()
                        if matches
                        else ("Exact disease and inheritance context not selected or mismatched",),
                    )
                )
        if source.provider == "clinvar-protein":
            primary = target["primary"]
            target_protein = str(primary.get("hgvs_p") or "").replace("(", "").replace(")", "")
            other_protein = (
                str(row.get("protein_hgvs") or "").split(":")[-1].replace("(", "").replace(")", "")
            )
            same_transcript = row.get("transcript") == primary.get("transcript")
            positive = str(row.get("clinical_significance", "")).lower() in {
                "pathogenic",
                "likely pathogenic",
                "pathogenic/likely pathogenic",
            }
            disease_ids = re.split(r"[,;|\s]+", str(row.get("phenotype_ids", "")))
            same_disease = context is not None and context["disease_id"] in disease_ids
            distinct = bool(
                row.get("coding_hgvs")
                and row["coding_hgvs"] != primary.get("hgvs_c")
                and str(row["coding_hgvs"]).split(":")[-1] != str(primary.get("hgvs_c", "")).split(":")[-1]
            )
            target_match = re.fullmatch(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", target_protein)
            other_match = re.fullmatch(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", other_protein)
            if (
                same_transcript
                and positive
                and distinct
                and target_match
                and other_match
                and target_match.group(1, 2) == other_match.group(1, 2)
            ):
                code = "PS1" if target_protein == other_protein else "PM5"
                criteria.append(
                    _ranking_criterion(
                        code,
                        "strong" if code == "PS1" else "moderate",
                        (source,),
                        method="September ClinVar same-transcript protein comparator",
                        group="protein_comparison",
                        status="met" if same_disease else "unknown",
                        limitations=()
                        if same_disease
                        else ("Comparator disease identifier does not match selected context",),
                    )
                )
    return criteria, tuple(assertions)


def _ranking_variant(packet, target, sources, index):
    gene, key = _identity(target)
    refs = _ranking_bound(target, sources)
    contexts = _ranking_contexts(target, sources)
    criteria, assertions = _ranking_assertions(target, refs, contexts)
    terms, consequence_source = _ranking_consequence(target)
    criteria.extend(_ranking_population(target))
    criteria.extend(_ranking_predictions(target, terms))
    criteria.extend(_ranking_family(packet, target, refs))
    anchor = _profile_anchor({"target": target})
    if anchor:
        criteria.append(
            _ranking_criterion(
                "PVS1",
                anchor["strength"],
                tuple(sources[item["source_id"]] for item in anchor["citations"]),
                method="validated AutoPVS1 decision tree",
                group="loss_of_function",
            )
        )
    if any(term in {"inframe_insertion", "inframe_deletion", "stop_lost"} for term in terms):
        criteria.append(
            _ranking_criterion(
                "PM4",
                "moderate",
                (consequence_source,),
                method="September in-frame protein length change",
                group="protein_length",
            )
        )
    if set(terms).intersection({"missense", "missense_variant"}):
        for item in refs:
            source = item["source"]
            if source.provider != "gnomad-constraint":
                continue
            raw = _ranking_json(source)
            if isinstance(raw.get("header"), (list, tuple)) and isinstance(raw.get("row"), str):
                row = dict(zip(raw["header"], raw["row"].split("\t")))
                try:
                    score = float(row.get("mis_z", "nan"))
                except (ValueError, TypeError):
                    continue
                if math.isfinite(score) and score >= 3.09 and row.get("gene") == gene:
                    criteria.append(
                        _ranking_criterion(
                            "PP2",
                            "supporting",
                            (source,),
                            method="September missense constraint Z >=3.09",
                            group="gene_constraint",
                            facts={"missense_z": score},
                        )
                    )
    primary = target["primary"]
    variant = {
        "variant_ref": f"{gene}:V{index}",
        "identity": {
            "gene": gene,
            "variant_key": key,
            "transcript": primary.get("transcript"),
            "hgvs_c": primary.get("hgvs_c"),
            "hgvs_p": primary.get("hgvs_p"),
        },
        "disease_contexts": contexts,
        "unresolved": tuple(
            value
            for value, missing in (
                ("disease_context_unavailable", not contexts),
                ("disease_context_selection_unresolved", len(contexts) > 1),
                ("proband_call_quality_unverified", target["observation"]["status"] != "carried"),
            )
            if missing
        ),
        "external_assertions": assertions,
    }
    variant = _ranking_summary(variant, criteria)
    return {**variant, "independent_classification": variant["classification"]}


def _ranking_pm3(packet, targets, variants, sources):
    variants = tuple(
        _ranking_summary(
            row,
            tuple(
                criterion
                for criterion in row["criteria"]
                if not (criterion["code"] == "PM3" and criterion["origin"] == "deterministic")
            ),
        )
        for row in variants
    )
    reviews = {_identity(target): row for target, row in zip(targets, variants)}
    output = []
    for target, review in zip(targets, variants):
        contexts = review["disease_contexts"]
        if len(contexts) != 1 or contexts[0]["mode_of_inheritance"] != "AR":
            output.append(review)
            continue
        child = _calls({"target": target}).get("proband")
        observations = []
        rare = any(item["code"] == "PM2" and item["counted"] for item in review["criteria"])
        compatible = target.get("patient_fit_status") == "compatible"
        if child and child["usable"] and not target.get("family_conflict"):
            if _call_value(child, "ploidy") == 2 and _call_value(child, "selected_alt_dosage") == 2:
                observations.append((0.5, "homozygous", (target,)))
            elif _call_value(child, "ploidy") == 2 and _call_value(child, "selected_alt_dosage") == 1:
                for partner in targets:
                    other = reviews[_identity(partner)]
                    other_contexts = other["disease_contexts"]
                    same_disease = (
                        len(other_contexts) == 1
                        and all(
                            other_contexts[0][name] == contexts[0][name]
                            for name in ("disease_id", "mode_of_inheritance")
                        )
                        and bool(contexts[0]["disease_id"])
                    )
                    if partner is target or _identity(partner)[0] != _identity(target)[0] or not same_disease:
                        continue
                    other_call = _calls({"target": partner}).get("proband")
                    if (
                        not other_call
                        or not other_call["usable"]
                        or _call_value(other_call, "ploidy") != 2
                        or _call_value(other_call, "selected_alt_dosage") != 1
                    ):
                        continue
                    if _key(
                        partner["observation"]["original"].get("vcf_variant_key", _identity(partner)[1])
                    ) == _key(target["observation"]["original"].get("vcf_variant_key", _identity(target)[1])):
                        continue
                    if (
                        not any(item["code"] == "PM2" and item["counted"] for item in other["criteria"])
                        or partner.get("patient_fit_status") != "compatible"
                    ):
                        continue
                    phase = _phase(target, partner)
                    grade = other["independent_classification"]
                    weight = (
                        (
                            1.0
                            if grade in _POSITIVE_CLASSES
                            else 0.25
                            if grade == "uncertain_significance"
                            else 0
                        )
                        if phase == "confirmed_in_trans"
                        else (0.5 if grade == "pathogenic" else 0.25 if grade == "likely_pathogenic" else 0)
                        if phase == "phase_unknown"
                        else 0
                    )
                    observations.append(
                        (
                            weight,
                            "vus_partner" if grade == "uncertain_significance" else "p_lp_partner",
                            (target, partner),
                        )
                    )
        weight, group, observed = max(
            observations, key=lambda row: (row[0], row[1] == "p_lp_partner"), default=(0, "", ())
        )
        capped = (
            min(0.5, weight)
            if group == "vus_partner"
            else min(1.0, weight)
            if group == "homozygous"
            else weight
        )
        score = capped if rare and compatible else 0
        refs = tuple(
            sources[binding["source_id"]]
            for item in observed
            for binding in item["bindings"]
            if binding["role"] == "patient_observation"
        )
        criterion = _ranking_criterion(
            "PM3",
            _ranking_strength(score),
            refs,
            method="preserved ClinGen SVI PM3 v1.0 independent-partner ledger",
            group="allelic_in_trans",
            status="met" if score >= 0.5 else "unknown",
            facts={
                "score": score,
                "contribution_group": group,
                "vus_partner_aggregate_cap": 0.5,
                "homozygous_aggregate_cap": 1.0,
                "partner_grades_exclude_pm3": True,
            },
            limitations=()
            if score >= 0.5
            else ("No quality-qualified rare affected-proband observation reached PM3 threshold",),
        )
        output.append(_ranking_summary(review, (*review["criteria"], criterion)))
    return tuple(output)


def _ranking_quality_genotype(packet, targets, variants):
    refs = tuple(row["variant_ref"] for row in variants)
    context_keys = {
        (row["disease_id"], row["mode_of_inheritance"])
        for variant in variants
        for row in variant["disease_contexts"]
    }
    mode = next(iter(context_keys))[1] if len(context_keys) == 1 else "unknown"
    status, phase, selected = "uncertain", "phase_unknown", refs
    usable = [
        (target, row, _calls({"target": target}).get("proband")) for target, row in zip(targets, variants)
    ]
    usable = [
        (target, row, call)
        for target, row, call in usable
        if call and call["usable"] and _call_value(call, "selected_alt_dosage", 0) > 0
    ]
    unresolved = list(dict.fromkeys(value for row in variants for value in row["unresolved"]))
    if not usable:
        status = "missing"
    elif any(target.get("family_conflict") for target in targets):
        status = "contradicted"
    elif len(context_keys) != 1 or not next(iter(context_keys))[0]:
        unresolved.append("coherent_disease_inheritance_context_unavailable")
    elif mode == "AR":
        homozygous = [
            row["variant_ref"]
            for _, row, call in usable
            if _call_value(call, "ploidy") == 2 and _call_value(call, "selected_alt_dosage") == 2
        ]
        pairs = [
            (left[1]["variant_ref"], right[1]["variant_ref"])
            for index, left in enumerate(usable)
            for right in usable[index + 1 :]
            if _phase(left[0], right[0]) == "confirmed_in_trans"
        ]
        hypotheses = tuple((ref,) for ref in homozygous) + tuple(pairs)
        if len(hypotheses) == 1:
            selected, status = hypotheses[0], "complete"
            phase = "homozygous" if len(selected) == 1 else "confirmed_in_trans"
        elif not hypotheses and len(usable) == 1:
            status = "incomplete"
            unresolved.append("second_recessive_allele_unavailable")
        else:
            unresolved.append("recessive_hypothesis_or_phase_unresolved")
    elif mode == "AD":
        if (
            len(usable) == 1
            and _call_value(usable[0][2], "ploidy") == 2
            and _call_value(usable[0][2], "selected_alt_dosage") == 1
        ):
            selected, status = (usable[0][1]["variant_ref"],), "complete"
        else:
            unresolved.append("dominant_genotype_hypothesis_unresolved")
    elif mode == "XL":
        if (
            len(usable) == 1
            and packet["family"].get("sex") == "MALE"
            and _call_value(usable[0][2], "ploidy") == 1
            and _call_value(usable[0][2], "selected_alt_dosage") == 1
        ):
            selected, status, phase = (usable[0][1]["variant_ref"],), "complete", "hemizygous"
        else:
            unresolved.append("x_linked_sex_or_genotype_hypothesis_unresolved")
    return {
        "genotype_status": status,
        "phase": phase,
        "unresolved": tuple(dict.fromkeys(unresolved)),
    }


def _ranking_normalized_mode(value):
    value = re.sub(r"[^A-Z0-9]+", "_", str(value or "").upper()).strip("_")
    aliases = {
        "AD": "AUTOSOMAL_DOMINANT",
        "AR": "AUTOSOMAL_RECESSIVE",
        "AR_COMP_HET": "AUTOSOMAL_RECESSIVE",
        "AUTOSOMAL_RECESSIVE_COMP_HET": "AUTOSOMAL_RECESSIVE",
        "AUTOSOMAL_RECESSIVE_COMPOUND_HETEROZYGOUS": "AUTOSOMAL_RECESSIVE",
        "AUTOSOMAL_RECESSIVE_HOM_ALT": "AUTOSOMAL_RECESSIVE",
        "XD": "X_LINKED_DOMINANT",
        "XR": "X_LINKED_RECESSIVE",
        "MITO": "MITOCHONDRIAL",
        "MT": "MITOCHONDRIAL",
        "": "UNKNOWN",
        "UNRESOLVED": "UNKNOWN",
        "UNSPECIFIED": "UNKNOWN",
        "NA": "UNKNOWN",
        "N_A": "UNKNOWN",
    }
    return aliases.get(value, value)


def _ranking_heterozygous_pair(values, target_by_key):
    heterozygous = []
    for value in values:
        if value["key"] not in target_by_key:
            continue
        target = target_by_key[value["key"]]
        child = _calls({"target": target}).get("proband")
        if child and _call_value(child, "ploidy") == 2 and _call_value(child, "selected_alt_dosage") == 1:
            heterozygous.append(value)
    for index, left in enumerate(heterozygous):
        for right in heterozygous[index + 1 :]:
            if left["key"][1:4] != right["key"][1:4]:
                return (left, right)
    return ()


def _ranking_native_hypotheses(packet, targets):
    """Build source-specific genotype hypotheses from native variant evidence."""
    by_key = {_identity(target)[1]: target for target in targets}
    by_source = {"exomiser": [], "lirical": []}
    for target in targets:
        for native in (target["primary"], *target.get("aliases", ())):
            provider = native["source"].provider
            if provider in by_source and variant_is_candidate(native):
                row = {
                    "key": _key(native["key"]),
                    "details": native.get("details", {}),
                    "source_id": native["source"].source_id,
                }
                if row not in by_source[provider]:
                    by_source[provider].append(row)
    gene = _identity(targets[0])[0] if targets else ""
    known_exomiser = {row["key"] for row in by_source["exomiser"]}
    for source in packet["sources"]:
        if source.provider != "exomiser" or source.kind != "tool_observation" or gene not in source.genes:
            continue
        raw = _ranking_json(source)
        if str(raw.get("CONTRIBUTING_VARIANT", "")).lower() not in {"true", "1"}:
            continue
        for key in source.variant_keys:
            if key not in known_exomiser:
                by_source["exomiser"].append(
                    {
                        "key": key,
                        "details": {"moi": raw.get("MOI"), "contributing": True},
                        "source_id": source.source_id,
                    }
                )
                known_exomiser.add(key)
    hypotheses = []

    def add(mode, values, provider, extra_sources=()):
        hypotheses.append(
            {
                "inheritance_mode": mode or "UNKNOWN",
                "variant_keys": tuple(dict.fromkeys(row["key"] for row in values)),
                "compatibility": "tool_supported",
                "phase": "unknown",
                "source": provider,
                "source_ids": tuple(dict.fromkeys((*extra_sources, *(row["source_id"] for row in values)))),
            }
        )

    modes = {}
    for row in by_source["exomiser"]:
        if row["details"].get("contributing") is True:
            modes.setdefault(str(row["details"].get("moi") or "UNKNOWN").upper(), []).append(row)
    compound = ()
    for mode, values in modes.items():
        if _ranking_normalized_mode(mode) == "AUTOSOMAL_RECESSIVE":
            compound = _ranking_heterozygous_pair(values, by_key)
        if compound:
            add("AUTOSOMAL_RECESSIVE_COMP_HET", compound, "exomiser")
            break
    for mode, values in modes.items():
        if compound and _ranking_normalized_mode(mode) == "AUTOSOMAL_RECESSIVE":
            continue
        for row in values:
            add(mode, (row,), "exomiser")
    lirical = {}
    for row in by_source["lirical"]:
        lirical.setdefault(row["details"].get("disease_id") or "", []).append(row)
    for disease, values in lirical.items():
        description, json_sources = [], []
        for source in packet["sources"]:
            if source.provider != "lirical" or gene not in source.genes:
                continue
            raw = _ranking_json(source)
            if raw.get("diseaseId") != disease or not isinstance(raw.get("genotypeLR"), dict):
                continue
            genotype = raw["genotypeLR"]
            description.extend((str(genotype.get("matchType") or ""), str(genotype.get("explanation") or "")))
            json_sources.append(source.source_id)
        text = " ".join(description).upper()
        recessive = bool(re.search(r"\bAR\b|RECESSIVE", text))
        pair = _ranking_heterozygous_pair(values, by_key) if recessive else ()
        if pair:
            add("AUTOSOMAL_RECESSIVE_COMP_HET", pair, "lirical", json_sources)
            continue
        mode = "AUTOSOMAL_RECESSIVE" if recessive else "UNKNOWN"
        if recessive and any(
            _call_value(_calls({"target": by_key[row["key"]]}).get("proband", {}), "selected_alt_dosage") == 2
            for row in values
        ):
            mode = "AUTOSOMAL_RECESSIVE_HOM_ALT"
        for marker, candidate in (
            ("AD", "AUTOSOMAL_DOMINANT"),
            ("XD", "X_LINKED_DOMINANT"),
            ("XR", "X_LINKED_RECESSIVE"),
            ("MITO", "MITOCHONDRIAL"),
        ):
            if not recessive and marker in text:
                mode = candidate
                break
        add(mode, values, "lirical", json_sources)
    return tuple(hypotheses)


def _ranking_hypothesis_fact(hypothesis, variants):
    by_key = {row["identity"]["variant_key"]: row for row in variants}
    observed = hypothesis["variant_keys"]
    keys = tuple(key for key in by_key if key in observed)
    missing = tuple(key for key in observed if key not in by_key)
    modes = []
    for key in keys:
        row = by_key[key]
        contexts = [
            context
            for context in row["disease_contexts"]
            if _ranking_normalized_mode(context["mode_of_inheritance"]) != "UNKNOWN"
        ]
        known = {_ranking_normalized_mode(context["mode_of_inheritance"]) for context in contexts}
        normalized = next(iter(known)) if len(known) == 1 else "UNKNOWN"
        mode = contexts[0]["mode_of_inheritance"] if len(known) == 1 else "unknown"
        disease_ids = tuple(
            dict.fromkeys(context["disease_id"] for context in contexts if context["disease_id"])
        )
        modes.append(
            {
                "variant_ref": row["variant_ref"],
                "observed": mode,
                "normalized": normalized,
                "disease_id": disease_ids[0] if len(disease_ids) == 1 else "",
            }
        )
    normalized = _ranking_normalized_mode(hypothesis["inheritance_mode"])
    known_modes = {row["normalized"] for row in modes if row["normalized"] != "UNKNOWN"}
    conflicts = (
        ("hypothesis_mode_conflicts_with_variant_disease_mode",)
        if normalized != "UNKNOWN" and any(mode != normalized for mode in known_modes)
        else ()
    )
    consistent = (
        normalized != "UNKNOWN" and bool(modes) and all(row["normalized"] == normalized for row in modes)
    )
    disease_ids = tuple(dict.fromkeys(row["disease_id"] for row in modes if row["disease_id"]))
    return {
        **hypothesis,
        "normalized_inheritance_mode": normalized,
        "variant_keys": keys,
        "variant_refs": tuple(by_key[key]["variant_ref"] for key in keys),
        "observed_variant_keys": observed,
        "variant_disease_modes": tuple(modes),
        "disease_ids": disease_ids,
        "missing_variant_keys": missing,
        "status": "unusable"
        if missing or not observed
        else "conflicting"
        if conflicts
        else "consistent"
        if consistent
        else "unresolved",
        "conflicts": conflicts,
        "uncertainties": ("selected_variants_have_nonidentical_disease_ids",) if len(disease_ids) > 1 else (),
    }


def _ranking_hypothesis_sort(row):
    compound = (
        "COMP_HET" in row["inheritance_mode"].upper() or "COMPOUND_HET" in row["inheritance_mode"].upper()
    )
    return (
        row["normalized_inheritance_mode"],
        row["variant_refs"],
        row["phase"].lower(),
        0 if compound else 1,
        row["inheritance_mode"].upper(),
        row["source"],
    )


def _ranking_hypothesis_select(hypotheses):
    groups = {}
    for index, row in enumerate(hypotheses):
        if row["status"] == "consistent":
            contexts = tuple(
                sorted(
                    (item["variant_ref"], item["normalized"], item["disease_id"])
                    for item in row["variant_disease_modes"]
                )
            )
            semantic = (
                row["normalized_inheritance_mode"],
                row["phase"].lower(),
                frozenset(row["variant_refs"]),
                contexts,
            )
            groups.setdefault(semantic, []).append(index)
    if len(groups) == 1:
        indices = next(iter(groups.values()))
        return indices[0], "equivalent_disease_mode_consistent_hypotheses" if len(
            indices
        ) > 1 else "unique_disease_mode_consistent_hypothesis"
    if len(groups) > 1:
        return None, "multiple_distinct_disease_mode_consistent_hypotheses"
    usable = [index for index, row in enumerate(hypotheses) if row["status"] != "unusable"]
    if len(usable) == 1 and hypotheses[usable[0]]["status"] == "unresolved":
        return usable[0], "sole_usable_hypothesis_with_incomplete_disease_mode_context"
    return None, "no_unique_nonconflicting_hypothesis" if hypotheses else "no_tool_hypothesis"


def _ranking_call_state(call):
    if call.get("genotype_state") != "called":
        return call.get("genotype_state", "unknown")
    dosage = _call_value(call, "selected_alt_dosage", 0)
    ploidy = _call_value(call, "ploidy", len(_call_value(call, "alleles", ())))
    return (
        "hemizygous_alternate"
        if ploidy == 1 and dosage == 1
        else "homozygous_alternate"
        if ploidy == 2 and dosage == 2
        else "heterozygous"
        if ploidy == 2 and dosage == 1
        else "homozygous_reference"
        if dosage == 0
        else "unknown"
    )


def _ranking_genotype(packet, targets, variants):
    hypotheses = tuple(
        {**row, "hypothesis_index": index}
        for index, row in enumerate(
            sorted(
                (
                    _ranking_hypothesis_fact(hypothesis, variants)
                    for hypothesis in _ranking_native_hypotheses(packet, targets)
                ),
                key=_ranking_hypothesis_sort,
            )
        )
    )
    index, reason = _ranking_hypothesis_select(hypotheses)
    selected = hypotheses[index] if index is not None else None
    observed = tuple(row["variant_ref"] for row in variants)
    refs = (
        selected["variant_refs"]
        if selected is not None
        else tuple(
            dict.fromkeys(
                ref for row in hypotheses if row["status"] != "unusable" for ref in row["variant_refs"]
            )
        )
        or observed
    )
    selected_rows = tuple(row for row in variants if row["variant_ref"] in refs)
    by_key = {_identity(target)[1]: target for target in targets}
    selected_targets = tuple(by_key[row["identity"]["variant_key"]] for row in selected_rows)
    quality_rows = tuple(
        {
            **row,
            "disease_contexts": tuple(
                context
                for context in row["disease_contexts"]
                if _ranking_normalized_mode(context["mode_of_inheritance"]) != "UNKNOWN"
            ),
        }
        for row in selected_rows
    )
    quality = _ranking_quality_genotype(packet, selected_targets, quality_rows)
    status = (
        "missing"
        if quality["genotype_status"] == "missing"
        else quality["genotype_status"]
        if selected is not None and selected["status"] == "consistent"
        else "contradicted"
        if hypotheses and all(row["status"] == "conflicting" for row in hypotheses)
        else "uncertain"
    )
    phase = (
        "homozygous"
        if quality["phase"] == "homozygous"
        else "parental_origin_supported_trans"
        if quality["phase"] == "confirmed_in_trans" and selected is not None
        else selected["phase"]
        if selected is not None
        else "unknown"
    )
    unresolved = list(quality["unresolved"])
    if hypotheses and selected is None:
        unresolved.append("genotype_hypothesis_selection_unresolved")
    if any(row["status"] == "conflicting" for row in hypotheses):
        unresolved.append("genotype_hypothesis_conflict")
    family = {}
    states, origins = {}, {}
    for target, row in zip(targets, variants):
        calls = {}
        for call in target["observation"]["calls"]:
            role = _call_value(call, "role")
            calls[role] = {
                "genotype": _call_value(call, "gt"),
                "allele_depths": _call_value(call, "ad"),
                "depth": _call_value(call, "dp"),
                "genotype_quality": _call_value(call, "gq"),
                "call_status": call["genotype_state"],
                "quality_status": "pass" if call["usable"] else call["quality_state"],
                "usable_for_inheritance": call["usable"],
                "state": _ranking_call_state(call),
                "carries_alternate": _call_value(call, "selected_alt_dosage", 0) > 0,
                "sample_id": _call_value(call, "sample_id"),
                "role": role,
                "ploidy": _call_value(call, "ploidy"),
                "pl": _call_value(call, "pl"),
            }
        family[row["variant_ref"]] = calls
        states[row["variant_ref"]] = calls.get("proband", {}).get("state", "unknown")
        origins[row["variant_ref"]] = target["observation"]["origin"]
    return {
        "inheritance_mode": selected["inheritance_mode"] if selected is not None else "unknown",
        "normalized_inheritance_mode": selected["normalized_inheritance_mode"]
        if selected is not None
        else "UNKNOWN",
        "compatibility": selected["compatibility"] if selected is not None else "uncertain",
        "genotype_status": status,
        "phase": phase,
        "phase_evidence": {
            "kind": "parental_origin_supported_trans",
            "selected_variant_refs": refs,
            "source": "patient_family_vcf",
        }
        if phase == "parental_origin_supported_trans"
        else {},
        "selected_variant_refs": refs,
        "selected_variant_refs_role": "selected_hypothesis"
        if selected is not None
        else "observed_unresolved_candidates",
        "observed_variant_refs": observed,
        "selected_hypothesis": selected or {},
        "hypotheses": hypotheses,
        "alternative_hypotheses": tuple(row for offset, row in enumerate(hypotheses) if offset != index),
        "hypothesis_selection_reason": reason,
        "hypothesis_conflicts": tuple(
            {
                "hypothesis_index": row["hypothesis_index"],
                "source": row["source"],
                "conflicts": row["conflicts"],
            }
            for row in hypotheses
            if row["conflicts"]
        ),
        "parental_origins": origins,
        "proband_states": states,
        "family_observations": family,
        "unresolved": tuple(dict.fromkeys(unresolved)),
    }


def _ranking_degraded(issues):
    return any(not is_evidence_limitation(issue) for issue in issues)


def ranking_preliminary(prepared):
    """Compute source-bound pre-ranking rules without model calls or resource I/O."""
    packet = _ranking_packet(prepared)
    sources = {source.source_id: source for source in packet["sources"]}
    targets = packet["targets"]
    counts = dict.fromkeys(prepared.get("clinical", {}), 0)
    variants = []
    for target in targets:
        gene = _identity(target)[0]
        counts[gene] = counts.get(gene, 0) + 1
        variants.append(_ranking_variant(packet, target, sources, counts[gene]))
    variants = _ranking_pm3(packet, targets, tuple(variants), sources)
    by_gene = {}
    for gene in counts:
        rows = tuple(row for row in variants if row["identity"]["gene"] == gene)
        gene_targets = [target for target in targets if _identity(target)[0] == gene]
        configuration = _ranking_genotype(packet, gene_targets, rows)
        by_gene[gene] = {
            "gene": gene,
            "phase": "preliminary",
            "variants": rows,
            "selected_variant_refs": configuration["selected_variant_refs"],
            "genotype_configuration": configuration,
            "unresolved": configuration["unresolved"],
        }
    issues = list(prepared.get("issues", ()))
    failed = {"query_failed", "resource_unavailable", "schema_mismatch", "identity_mismatch", "failed"}
    for target in targets:
        lookups = [(name, target.get(name) or {}) for name in ("population", "annotation", "pvs1")]
        lookups.extend(
            (name, value)
            for name, value in (target.get("predictions") or {}).items()
            if isinstance(value, Mapping)
        )
        for name, lookup in lookups:
            if lookup.get("status") in failed:
                issues.append(
                    f"ranking_preliminary:{_identity(target)[0]}:{name}:{lookup.get('error_code') or lookup['status']}"
                )
    issues = tuple(dict.fromkeys(issues))
    ledger = {
        "schema_version": "ranking_acmg_v1",
        "case_key": case_tuple(packet["profile"]["case_key"]),
        "phase": "preliminary",
        "status": "degraded" if _ranking_degraded(issues) else "ok",
        "rules_version": _RANKING_RULES,
        "by_gene": by_gene,
        "coverage": (),
        "issues": issues,
    }
    ledger["source_ids"] = _ranking_source_ids(ledger, packet)
    return ledger


def _ranking_articles(evidence):
    groups, abstracts = {}, []
    for document in evidence.get("documents", ()):
        source = document["source"]
        if source.kind != "literature":
            continue
        if source.provider == "pmc":
            metadata = document.get("metadata", {})
            parent = metadata.get("parent_pmid", "")
            article = (
                parent
                if parent.startswith("PMID:")
                else "PMC:" + source.locator.split("/articles/")[-1].split("/", 1)[0]
            )
            groups.setdefault(article, []).append(document)
        elif source.provider == "pubmed" and document["prompt_id"].startswith("PMID:"):
            abstracts.append(document)
    packets = [(key, tuple(rows[:4])) for key, rows in groups.items()]
    seen = set(groups)
    for document in abstracts:
        article = document["prompt_id"].split("#", 1)[0]
        if article not in seen:
            packets.append((article, (document,)))
            seen.add(article)
    return tuple(packets[:2])


def _ranking_pp4(refs):
    rows, sources = [], []
    for item in refs:
        if item["binding"]["role"] != "guideline":
            continue
        source = item["source"]
        raw = _ranking_json(source)
        row = raw.get("row", raw)
        if isinstance(row, dict) and "pp4_evidence_level_to_be_assigned" in row:
            rows.append(row)
            sources.append(source)
    strengths = {
        value.lower()
        for row in rows
        for value in re.findall(
            r"PP4(?:_)?(Supporting|Moderate|Strong)", str(row["pp4_evidence_level_to_be_assigned"]), re.I
        )
    }
    return tuple(rows), tuple(sources), strengths


def _ranking_literature_payload(packet, target, variant, article, documents, sources):
    refs = _ranking_bound(target, sources)
    specificity, specificity_sources = _ranking_specificity(packet, target, refs)
    pp4, pp4_sources, strengths = _ranking_pp4(refs)
    context = variant["disease_contexts"][0] if len(variant["disease_contexts"]) == 1 else {}
    payload = {
        "task": "Extract explicit literature evidence for the allowed ACMG criteria.",
        "patient_phenotypes": [
            {"hpo_id": term["hpo_id"], "label": term["label"]}
            for term in packet["profile"]["terms"]
            if term["present"]
        ],
        "target": {
            **variant["identity"],
            "variant_key": _key_text(variant["identity"]["variant_key"]),
            "variant_ref": variant["variant_ref"],
            **context,
            "observed_consequence": _ranking_consequence(target)[0],
            "observed_family_genotypes": {
                "origin": target["observation"]["origin"],
                "calls": [_call_projection(call) for call in target["observation"]["calls"]],
            },
        },
        "publication": {
            "evidence_id": article,
            "title": documents[0]["title"],
            "text": "\n\n".join(
                document["prompt_text"][:1800]
                if document["source"].provider == "pmc"
                else document["prompt_text"]
                for document in documents
            ),
            "provenance": list(
                dict.fromkeys(route for document in documents for route in document["provenance"])
            ),
        },
        "rules": {
            "allowed_criteria": sorted(_LITERATURE_CODES),
            "PVS1_RNA": "explicit patient RNA evidence showing the predicted abnormal transcript",
            "PS2_PM6": "the target variant is explicitly de novo; use PS2 only with confirmed parentage",
            "PS3": "a well-established functional assay directly tests the target variant",
            "PS4": "multiple unrelated affected probands carry the target variant; provide count",
            "PM3": "the target variant is explicitly observed in trans with a pathogenic variant",
            "PP1": "the target variant cosegregates with disease; provide informative meioses",
            "PP4": "the current phenotype is highly specific for the supplied disease mechanism",
        },
        "ai_cura_knowledgebase": {
            "source_ref": tuple(source.source_id for source in (*specificity_sources, *pp4_sources)),
            "gene_disease_consistency": {"specificity": specificity},
            "pp4_rules_for_gene": pp4,
        },
        "required_json_schema": {
            "findings": [
                {
                    "criterion": "PS3",
                    "status": "met|unknown|not_met",
                    "strength": "leave empty; rule strength is assigned deterministically",
                    "source_excerpt": "complete verbatim sentence or compact table row",
                    "reason": "short criterion-specific explanation",
                    "affected_count": 0,
                    "cohort_key": "author/institution/location/time identifier or empty",
                    "authors": [],
                    "affiliations": [],
                    "region": "",
                    "recruitment_period": "",
                    "cohort_description": "",
                    "deleterious_functional_effect": False,
                    "normal_control_present": False,
                    "abnormal_control_present": False,
                    "multicellular_model_organism": False,
                    "mave_high_throughput_assay": False,
                    "confirmed_parental_relationships": False,
                    "de_novo_confirmed_count": 0,
                    "de_novo_unconfirmed_count": 0,
                    "phenotypic_consistency": "copy one supplied Table S8 specificity label or empty",
                    "confirmed_in_trans": False,
                    "pp4_requirements_met": False,
                    "pp4_strength": "supporting|moderate|strong|empty",
                }
            ]
        },
    }
    return payload, specificity, strengths, tuple(dict.fromkeys((*specificity_sources, *pp4_sources)))


def _ranking_int(value):
    if type(value) is not int or value < 0:
        raise MendelPrimeError("nonnegative_integer_required", category="response")
    return value


def _ranking_literature_decode(raw, target, article, documents, specificity, pp4_strengths, rule_sources):
    if not isinstance(raw.get("findings"), list):
        raise MendelPrimeError("findings_must_be_list", category="response")
    findings = []
    for item in raw["findings"]:
        if not isinstance(item, dict) or item.get("status") != "met":
            continue
        code = str(item.get("criterion", "")).upper()
        excerpt = str(item.get("source_excerpt", "")).strip()
        if code not in _LITERATURE_CODES or len(excerpt) < 20:
            continue
        refs, citations = [], []
        for document in documents:
            source = document["source"]
            if not _contains(source.content, _identity(target)[0]) or not _mentions_allele(
                target, source.content
            ):
                continue
            pattern = r"\s+".join(re.escape(word) for word in excerpt.split())
            match = re.search(pattern, source.content, re.I)
            if match:
                refs.append(source)
                citations.append(_citation(source.source_id, match[0]))
        if not refs:
            continue
        facts = {
            name: item.get(name, default)
            for name, default in (
                ("cohort_key", ""),
                ("authors", []),
                ("affiliations", []),
                ("region", ""),
                ("recruitment_period", ""),
                ("cohort_description", ""),
            )
        }
        if any(
            not isinstance(facts[name], str)
            for name in ("cohort_key", "region", "recruitment_period", "cohort_description")
        ) or any(
            not isinstance(facts[name], list) or any(not isinstance(value, str) for value in facts[name])
            for name in ("authors", "affiliations")
        ):
            raise MendelPrimeError("cohort_metadata_invalid", category="response")
        strength, status, limitations = (
            "unassigned",
            "unknown",
            ("September criterion-specific scoring unavailable; source fact retained without scoring",),
        )
        canonical = "PVS1" if code == "PVS1_RNA" else code
        facts.update(
            article_key=article, subject_scope="published_cases", subject_group=facts["cohort_key"] or article
        )
        if code == "PS4":
            count = _ranking_int(item.get("affected_count", 0))
            facts["affected_count"] = count
            if count >= 1:
                strength, status = (
                    "strong" if count >= 4 else "moderate" if count >= 2 else "supporting",
                    "met",
                )
                limitations = ("September AI-CURA sensitivity-prioritized PS4 counts",)
        elif code in {"PS2", "PM6"}:
            if (
                not specificity
                or item.get("phenotypic_consistency", "").casefold().strip() != specificity.casefold().strip()
            ):
                continue
            confirmed = _ranking_int(item.get("de_novo_confirmed_count", 0))
            unconfirmed = _ranking_int(item.get("de_novo_unconfirmed_count", 0))
            score = _ranking_denovo_weight(specificity, confirmed, unconfirmed)
            canonical = "PS2" if confirmed else "PM6"
            facts.update(
                de_novo_score=score,
                de_novo_confirmed_count=confirmed,
                de_novo_unconfirmed_count=unconfirmed,
                phenotypic_consistency=specificity,
            )
            if score >= 0.5:
                strength, status, limitations = _ranking_strength(score), "met", ()
        elif (
            code == "PS3"
            and item.get("deleterious_functional_effect") is True
            and (
                (item.get("normal_control_present") is True and item.get("abnormal_control_present") is True)
                or item.get("multicellular_model_organism") is True
            )
        ):
            strength, status = "supporting", "met"
            limitations = (
                ("MAVE is not automatically upgraded",)
                if item.get("mave_high_throughput_assay") is True
                else ()
            )
        elif (
            code == "PP4"
            and item.get("pp4_requirements_met") is True
            and item.get("pp4_strength") in pp4_strengths
        ):
            strength, status, limitations = item["pp4_strength"], "met", ()
        if code in {"PS2", "PM6", "PP4"}:
            refs.extend(rule_sources)
            citations.extend(_citation(source.source_id, source.content) for source in rule_sources)
        group = (
            "de_novo"
            if code in {"PS2", "PM6"}
            else "loss_of_function"
            if code == "PVS1_RNA"
            else "allelic_in_trans"
            if code == "PM3"
            else f"literature:{canonical}:{facts['subject_group']}"
        )
        criterion = _ranking_criterion(
            canonical,
            strength,
            tuple(refs),
            method=f"September single-publication extraction: {code}",
            group=group,
            status=status,
            facts=facts,
            limitations=limitations,
            origin="literature",
        )
        findings.append(_ranking_with_citations(criterion, citations))
    return tuple(findings)


def _ranking_request(packet, stage, system, payload, tokens, decode):
    request = {
        "case_key": case_tuple(packet["profile"]["case_key"]),
        "stage_id": stage,
        "system_prompt": prompting.text(system),
        "user_prompt": json.dumps(payload, ensure_ascii=False, separators=(",", ":")),
        "max_tokens": tokens,
        "timeout_seconds": 180.0,
        "cache_choice": False,
    }
    return request, decode


def _ranking_ps4(model, packet, variant, criteria):
    criteria = tuple(
        row for row in criteria if not (row["code"] == "PS4" and row["origin"] == "literature_aggregate")
    )
    rows = [
        row
        for row in criteria
        if row["code"] == "PS4"
        and row["status"] == "met"
        and row["origin"] == "literature"
        and row["raw_facts"].get("article_key")
    ]
    by_article = {}
    for row in rows:
        key = row["raw_facts"]["article_key"]
        if (
            key not in by_article
            or row["raw_facts"]["affected_count"] > by_article[key]["raw_facts"]["affected_count"]
        ):
            by_article[key] = row
    if len(by_article) < 2:
        return criteria
    summaries = [
        {
            "source_ref": key,
            **{
                name: row["raw_facts"].get(
                    name, [] if name in {"authors", "affiliations"} else 0 if name == "affected_count" else ""
                )
                for name in (
                    "affected_count",
                    "authors",
                    "affiliations",
                    "region",
                    "recruitment_period",
                    "cohort_description",
                )
            },
        }
        for key, row in by_article.items()
    ]
    payload = {
        "target": {"gene": variant["identity"]["gene"], "variant_ref": variant["variant_ref"]},
        "publication_summaries": summaries,
        "required_json_schema": {
            "included_source_refs": list(by_article),
            "excluded_as_potential_overlap": [],
            "total_unique_probands": 0,
            "regional_groups": [],
        },
    }

    def decode(raw):
        included = raw.get("included_source_refs")
        if (
            not isinstance(included, list)
            or not included
            or len(included) != len(set(included))
            or any(key not in by_article for key in included)
        ):
            raise MendelPrimeError("ps4_included_sources_invalid", category="response")
        total = _ranking_int(raw.get("total_unique_probands"))
        if not 1 <= total <= sum(by_article[key]["raw_facts"]["affected_count"] for key in included):
            raise MendelPrimeError("ps4_unique_count_invalid", category="response")
        chosen = [by_article[key] for key in included]
        strength = "strong" if total >= 4 else "moderate" if total >= 2 else "supporting"
        merged = _ranking_with_citations(
            {
                **chosen[0],
                "origin": "literature_aggregate",
                "strength": strength,
                "points": _RANKING_POINTS[strength],
                "independence_group": "literature:PS4:deduplicated",
                "method": "September PS4 publication extraction and cohort de-duplication",
                "raw_facts": {
                    "affected_count": total,
                    "included_article_keys": tuple(included),
                    "cohort_deduplicated": True,
                },
            },
            (citation for row in chosen for citation in row["citations"]),
        )
        return tuple(criteria) + (merged,)

    request, decoder = _ranking_request(
        packet,
        f"vcf_ranking_ps4_dedup:{variant['variant_ref']}",
        "vcf_ranking_ps4_cohort_deduplication",
        payload,
        1200,
        decode,
    )
    return prompting.call_json(model, request, decoder)


def _ranking_denovo_merge(criteria):
    criteria = tuple(
        row
        for row in criteria
        if not (row["code"] in {"PS2", "PM6"} and row["origin"] == "literature_aggregate")
    )
    candidates = [
        row
        for row in criteria
        if row["code"] in {"PS2", "PM6"}
        and row["status"] == "met"
        and row["raw_facts"].get("de_novo_score", 0) > 0
    ]
    if len(candidates) < 2:
        return criteria
    groups = {}
    for row in candidates:
        facts = row["raw_facts"]
        group = facts.get("subject_group") or facts.get("article_key") or row["source_refs"][0]
        if group not in groups or facts["de_novo_score"] > groups[group]["raw_facts"]["de_novo_score"]:
            groups[group] = row
    selected = tuple(groups.values())
    score = sum(row["raw_facts"]["de_novo_score"] for row in selected)
    confirmed = sum(row["raw_facts"].get("de_novo_confirmed_count", 0) for row in selected)
    unconfirmed = sum(row["raw_facts"].get("de_novo_unconfirmed_count", 0) for row in selected)
    strength = _ranking_strength(score)
    merged = _ranking_with_citations(
        {
            **selected[0],
            "origin": "literature_aggregate",
            "code": "PS2" if confirmed else "PM6",
            "strength": strength,
            "points": _RANKING_POINTS[strength],
            "raw_facts": {
                "de_novo_score": score,
                "de_novo_confirmed_count": confirmed,
                "de_novo_unconfirmed_count": unconfirmed,
                "subject_groups": tuple(groups),
            },
            "method": "September Table S8 independent-source de novo aggregation",
        },
        (citation for row in selected for citation in row["citations"]),
    )
    return tuple(criteria) + (merged,)


def _ranking_source_ids(ledger, packet):
    identifiers = list(ledger.get("source_ids", ()))
    for target in packet["targets"]:
        identifiers.append(target["primary"]["source"].source_id)
        identifiers.extend(row["source"].source_id for row in target.get("aliases", ()))
        identifiers.extend(binding["source_id"] for binding in target["bindings"])
    for gene in ledger["by_gene"].values():
        identifiers.extend(
            identifier
            for hypothesis in gene["genotype_configuration"].get("hypotheses", ())
            for identifier in hypothesis.get("source_ids", ())
        )
        for variant in gene["variants"]:
            identifiers.extend(
                source_id for context in variant["disease_contexts"] for source_id in context["source_ids"]
            )
            identifiers.extend(row["source_id"] for row in variant["external_assertions"])
            identifiers.extend(
                source_id for criterion in variant["criteria"] for source_id in criterion["source_refs"]
            )
    identifiers.extend(source_id for row in ledger["coverage"] for source_id in row.get("source_ids", ()))
    return tuple(dict.fromkeys(identifiers))


def ranking_literature(model, prepared, ledger, evidence, *, targets, workers=1):
    """Extract facts only for bounded uncovered allele/article pairs; keep all genes."""
    packet = _ranking_packet(prepared)
    if ledger.get("schema_version") != "ranking_acmg_v1" or ledger.get("case_key") != case_tuple(
        packet["profile"]["case_key"]
    ):
        raise MendelPrimeError("ranking_acmg_case_or_schema_mismatch", field="ranking_acmg")
    if (
        type(targets) is not tuple
        or any(not isinstance(gene, str) for gene in targets)
        or type(workers) is not int
        or workers < 1
    ):
        raise MendelPrimeError("ranking_literature_targets_or_workers_invalid", field="ranking_acmg")
    fresh = ranking_preliminary(prepared)
    merged_genes = dict(ledger["by_gene"])
    for gene, value in fresh["by_gene"].items():
        if gene not in merged_genes:
            merged_genes[gene] = value
        else:
            current = merged_genes[gene]
            refreshed = {row["identity"]["variant_key"]: row for row in value["variants"]}
            rows = []
            for row in current["variants"]:
                replacement = refreshed.get(row["identity"]["variant_key"])
                if replacement is None:
                    rows.append(row)
                    continue
                literature = tuple(
                    criterion
                    for criterion in row["criteria"]
                    if criterion["origin"] in {"literature", "literature_aggregate"}
                )
                replacement = _ranking_summary(
                    {**replacement, "variant_ref": row["variant_ref"]},
                    (*replacement["criteria"], *literature),
                )
                rows.append(_ranking_independent(replacement))
            current = {**current, "variants": tuple(rows)}
            known = {row["identity"]["variant_key"] for row in current["variants"]}
            added = tuple(
                {**row, "variant_ref": f"{gene}:V{len(current['variants']) + index}"}
                for index, row in enumerate(
                    (row for row in value["variants"] if row["identity"]["variant_key"] not in known), 1
                )
            )
            rows = (*current["variants"], *added)
            by_identity = {_identity(target): target for target in packet["targets"]}
            if all((gene, row["identity"]["variant_key"]) in by_identity for row in rows):
                configuration = _ranking_genotype(
                    packet, tuple(by_identity[(gene, row["identity"]["variant_key"])] for row in rows), rows
                )
                merged_genes[gene] = {
                    **current,
                    "variants": rows,
                    "selected_variant_refs": configuration["selected_variant_refs"],
                    "genotype_configuration": configuration,
                    "unresolved": configuration["unresolved"],
                }
            else:
                merged_genes[gene] = {**current, "variants": rows}
    ledger = {
        **ledger,
        "by_gene": merged_genes,
        "issues": tuple(dict.fromkeys((*ledger["issues"], *fresh["issues"]))),
    }
    sources = {source.source_id: source for source in packet["sources"]}
    for row in evidence.values():
        for document in (*row.get("documents", ()), *row.get("omim_context", ())):
            source = document["source"]
            if source.source_id in sources and sources[source.source_id] != source:
                raise MendelPrimeError("source_id_collision", field="ranking_acmg")
            if source.case_key is not None and source.case_key != ledger["case_key"]:
                raise MendelPrimeError("source_case_mismatch", field="ranking_acmg")
            sources[source.source_id] = source
    if any(parent not in sources for source in sources.values() for parent in source.parent_source_ids):
        raise MendelPrimeError("parent_source_missing", field="ranking_acmg")
    if any(identifier not in sources for identifier in _ranking_source_ids(ledger, packet)):
        raise MendelPrimeError("retained_ledger_source_missing", field="ranking_acmg")
    prepared_targets = {_identity(row): row for row in packet["targets"]}
    coverage = list(ledger["coverage"])
    findings, jobs, issues = {}, [], []
    reviewed = tuple(dict.fromkeys(targets))[:10]
    for gene in reviewed:
        if gene not in ledger["by_gene"]:
            record = {
                "identity": {"gene": gene, "variant_key": None},
                "article_key": "",
                "status": "no_eligible_allele",
                "source_ids": (),
            }
            if record not in coverage:
                coverage.append(record)
            continue
        gene_evidence = evidence.get(gene, {})
        issues.extend(
            f"ranking_literature:{gene}:{error}"
            for error in gene_evidence.get("retrieval_errors", ())
            if not str(error).endswith("article_or_abstract_missing")
        )
        articles = _ranking_articles(gene_evidence)
        assessment = ledger["by_gene"][gene]
        configuration = assessment["genotype_configuration"]
        by_ref = {row["variant_ref"]: row for row in assessment["variants"]}
        refs = configuration["selected_variant_refs"]
        if configuration["selected_variant_refs_role"] == "selected_hypothesis":
            selected_variants = tuple(by_ref[ref] for ref in refs if ref in by_ref)[:2]
        else:
            selected_variants = assessment["variants"][:2]
        for variant in selected_variants:
            identity = variant["identity"]
            scoped_identity = {"gene": gene, "variant_key": identity["variant_key"]}
            target = prepared_targets.get((gene, identity["variant_key"]))
            if target is None:
                issues.append(f"ranking_literature:{gene}:prepared_allele_missing")
                continue
            if variant["classification"] in {"pathogenic", "benign"} and not variant["conflicting_evidence"]:
                record = {
                    "identity": scoped_identity,
                    "article_key": "",
                    "status": "skipped_resolved",
                    "source_ids": (),
                }
                if record not in coverage:
                    coverage.append(record)
                continue
            if not articles:
                record = {
                    "identity": scoped_identity,
                    "article_key": "",
                    "status": "no_documents",
                    "source_ids": (),
                }
                if record not in coverage:
                    coverage.append(record)
            for article, documents in articles:
                if any(
                    source.genes and gene not in source.genes
                    for source in (document["source"] for document in documents)
                ):
                    raise MendelPrimeError("literature_gene_scope_mismatch", field="ranking_acmg")
                source_ids = tuple(document["source"].source_id for document in documents)
                if any(
                    row["identity"] == scoped_identity
                    and row["article_key"] == article
                    and row["status"] in {"processed", "no_findings"}
                    and row["source_ids"] == source_ids
                    for row in coverage
                ):
                    continue
                payload, specificity, strengths, rule_sources = _ranking_literature_payload(
                    packet, target, variant, article, documents, sources
                )
                request, decode = _ranking_request(
                    packet,
                    f"vcf_ranking_literature:{variant['variant_ref']}:{article}",
                    "vcf_ranking_literature_extraction",
                    payload,
                    1800,
                    lambda raw, target=target, article=article, documents=documents, specificity=specificity, strengths=strengths, rule_sources=rule_sources: (
                        _ranking_literature_decode(
                            raw, target, article, documents, specificity, strengths, rule_sources
                        )
                    ),
                )
                jobs.append((scoped_identity, variant["variant_ref"], article, documents, request, decode))

    def extract(job):
        try:
            return prompting.call_json(model, job[4], job[5]), None
        except Exception as error:
            return (), error.code if isinstance(error, MendelPrimeError) else "literature_extraction_failed"

    if jobs:
        with ThreadPoolExecutor(max_workers=min(workers, len(jobs))) as pool:
            outcomes = tuple(pool.map(extract, jobs))
        for (identity, ref, article, documents, _, _), (values, error) in zip(jobs, outcomes):
            coverage = [
                row for row in coverage if not (row["identity"] == identity and row["article_key"] == article)
            ]
            coverage.append(
                {
                    "identity": identity,
                    "article_key": article,
                    "status": "failed" if error else "processed" if values else "no_findings",
                    "source_ids": tuple(document["source"].source_id for document in documents),
                    "inspection_role": "bounded_allele_inspection_not_causal_selection",
                }
            )
            findings.setdefault((identity["gene"], ref), []).extend(values)
            if error:
                issues.append(f"ranking_literature:{identity['gene']}:{article}:{error}")
    by_gene = dict(ledger["by_gene"])
    for gene in reviewed:
        if gene not in by_gene:
            continue
        variants = []
        for variant in by_gene[gene]["variants"]:
            added = findings.get((gene, variant["variant_ref"]), ())
            replaced_articles = {row["raw_facts"].get("article_key") for row in added}
            criteria = tuple(
                row
                for row in variant["criteria"]
                if row["raw_facts"].get("article_key") not in replaced_articles
                or row["origin"] != "literature"
            ) + tuple(added)
            unresolved_ps4 = sum(
                row["code"] == "PS4" and row["origin"] == "literature" and row["status"] == "met"
                for row in criteria
            ) >= 2 and not any(
                row["code"] == "PS4" and row["origin"] == "literature_aggregate" for row in criteria
            )
            if added or unresolved_ps4:
                try:
                    criteria = _ranking_ps4(model, packet, variant, criteria)
                except Exception as error:
                    criteria = tuple(
                        row
                        for row in criteria
                        if not (row["code"] == "PS4" and row["origin"] == "literature_aggregate")
                    )
                    issues.append(
                        f"ranking_literature:{gene}:PS4_dedup:{error.code if isinstance(error, MendelPrimeError) else 'failed'}"
                    )
                criteria = _ranking_denovo_merge(criteria)
            variants.append(_ranking_summary(variant, criteria))
        by_gene[gene] = {**by_gene[gene], "phase": "literature_enhanced", "variants": tuple(variants)}
    # Freeze all PM3-free grades before updating any partner relationship.
    available = []
    independent = []
    for gene in by_gene.values():
        for variant in gene["variants"]:
            target = prepared_targets.get((gene["gene"], variant["identity"]["variant_key"]))
            if target is not None:
                available.append(target)
                independent.append(_ranking_independent(variant))
    updated_partners = {
        _identity(target): variant
        for target, variant in zip(available, _ranking_pm3(packet, available, tuple(independent), sources))
    }
    for gene, value in tuple(by_gene.items()):
        by_gene[gene] = {
            **value,
            "variants": tuple(
                updated_partners.get((gene, row["identity"]["variant_key"]), row) for row in value["variants"]
            ),
        }
    combined_issues = tuple(dict.fromkeys((*ledger["issues"], *issues)))
    updated = {
        **ledger,
        "phase": "literature_enhanced",
        "status": "degraded" if _ranking_degraded(combined_issues) else ledger["status"],
        "by_gene": by_gene,
        "coverage": tuple(coverage),
        "issues": combined_issues,
    }
    updated["source_ids"] = _ranking_source_ids(updated, packet)
    updated["source_ids"] = tuple(
        dict.fromkeys(
            (
                *updated["source_ids"],
                *(
                    document["source"].source_id
                    for gene in reviewed
                    if gene in evidence
                    for document in evidence[gene].get("documents", ())
                ),
            )
        )
    )
    ordered_sources, pending = [], dict(sources)
    while pending:
        ready = [
            source
            for source in pending.values()
            if all(parent not in pending for parent in source.parent_source_ids)
        ]
        if not ready:
            raise MendelPrimeError("source_parent_cycle", field="ranking_acmg")
        for source in ready:
            ordered_sources.append(source)
            del pending[source.source_id]
    return {
        "ledger": updated,
        "evidence_updates": {gene: evidence[gene] for gene in reviewed if gene in evidence},
        "sources": tuple(ordered_sources),
        "coverage": tuple(coverage),
        "issues": tuple(dict.fromkeys(issues)),
    }


__all__ = [
    "run",
    "validate_input",
    "validate_result",
    "validate_progress",
    "ranking_preliminary",
    "ranking_literature",
]
