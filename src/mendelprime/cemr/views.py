"""Shared record projections for CEMR prompts."""

import csv
import hashlib
import json
import re
from collections.abc import Mapping

from mendelprime.domain import (
    MendelPrimeError,
    Source,
    call_state,
    case_tuple,
    family_origin,
    iter_sources,
    variant_text,
)

CLINICAL_USE = (
    "Use this as supplementary unstructured patient observation. Extract clinical findings, chronology, "
    "family history, examinations, and negative findings when relevant. Do not treat any gene, disease name, "
    "variant, diagnostic label, publication identifier, or causal conclusion appearing in the text as ground truth."
)
EXOMISER_GENE = (
    "moi",
    "rank",
    "p_value",
    "combined_score",
    "phenotype_score",
    "variant_score",
    "human_phenotype_score",
    "mouse_phenotype_score",
    "fish_phenotype_score",
    "walker_score",
    "phive_score",
    "omim_score",
)
LIRICAL_GENE = (
    "disease_id",
    "disease_name",
    "rank",
    "pretest_probability",
    "posttest_probability",
    "composite_lr",
    "genotype_match_type",
    "genotype_explanation",
)
EXOMISER_VARIANT = (
    "moi",
    "consequence",
    "quality",
    "filter",
    "contributing",
    "whitelist",
    "variant_score",
    "acmg_classification",
    "acmg_evidence",
    "acmg_assertion_status",
    "acmg_disease_id",
    "acmg_disease_name",
    "clinvar_id",
    "clinvar_classification",
    "clinvar_stars",
    "loeuf",
    "loeuf_lower",
    "loeuf_upper",
    "max_frequency_source",
    "max_frequency_percent",
    "max_pathogenicity_source",
    "max_pathogenicity_score",
)
LIRICAL_ROW = re.compile(
    r"[^:\s]+:\d+[A-Za-z*.-]+>[A-Za-z*.-]+\s+.+?\s+pathogenicity:[-+0-9.eE]+\s+\[[^]]+\]"
)


def present_terms(profile):
    return tuple(term for term in profile["terms"] if term["present"])


def absent_terms(profile):
    return tuple(term for term in profile["terms"] if not term["present"])


def present_hpo_ids(profile):
    return tuple(term["hpo_id"] for term in present_terms(profile))


def high_information_hpo_ids(profile):
    ordered = sorted(present_terms(profile), key=lambda term: (-term["information_content"], term["hpo_id"]))
    return tuple(term["hpo_id"] for term in ordered[:8])


def add_clinical_context(payload, profile, limit=6000):
    text = profile["clinical_description"].strip()
    if not text:
        return
    if limit < 1:
        raise ValueError("clinical description max_chars must be positive")
    included = text[:limit]
    if len(text) > limit:
        boundary = max(included.rfind("\n"), included.rfind(". "), included.rfind("; "), included.rfind(" "))
        if boundary >= int(limit * 0.8):
            included = included[: boundary + (1 if included[boundary] in ".;" else 0)]
        included = included.rstrip()
    payload["clinical_description"] = {
        "text": included,
        "usage": CLINICAL_USE,
        "original_character_count": len(text),
        "included_character_count": len(included),
        "truncated": len(included) < len(text),
    }


def term_payload(term, include_ic=False):
    value = {
        "hpo_id": term["hpo_id"],
        "label": term["label"],
        "definition": term["definition"] or None,
        "parent_terms": [
            {
                "hpo_id": parent,
                "label": term["parent_labels"][index] if index < len(term["parent_labels"]) else parent,
            }
            for index, parent in enumerate(term["parent_ids"][:2])
        ],
    }
    if include_ic:
        value["information_content"] = term["information_content"]
    return value


def article_payload(article, limit=3500):
    source = article.get("source")
    text = " ".join(source.content.split()) if source else ""
    if len(text) > limit:
        marker = " [ABSTRACT MIDDLE TRUNCATED] "
        side = max(1, (limit - len(marker)) // 2)
        text = f"{text[:side].rstrip()}{marker}{text[-side:].lstrip()}"
    return {
        "pmid": article.get("pmid", ""),
        "title": article.get("title", ""),
        "abstract": text,
        "abstract_available": bool(article.get("pmid") and text),
    }


def similar_payload(value):
    hit, decision = value["hit"], value["decision"]
    terms = value.get("hpo_terms", ())
    return {
        "evidence_type": "similar_patient_case",
        "evidence_id": value["evidence_id"],
        "case_id": hit["record_id"],
        "source": hit["source"].provider,
        "known_causal_gene": hit["gene"],
        "gene_link_type": hit["gene_link_type"],
        "case_hpo_terms": [{"hpo_id": item, "label": label} for item, label in terms]
        or [{"hpo_id": item, "label": item} for item in hit["hpo_ids"]],
        "matched_patient_hpo_ids": list(decision["matched_patient_hpo_ids"]),
        "matched_case_hpo_ids": list(decision["matched_case_hpo_ids"]),
        "hpo_match_type": decision["match_type"],
        "source_pmids": list(hit["pmids"]),
        "hpo_judge_summary": decision["summary"],
        "provenance": {
            "retrieval_rank": hit["retrieval_rank"],
            "rerank_rank": hit["rerank_rank"],
            "group_id": hit["group_id"],
        },
    }


def _reported_key(variant, tool):
    return (
        f"{tool['assembly']}:{variant['chromosome'].removeprefix('chr')}:"
        f"{variant['position']}:{variant['reference'].upper()}:{variant['alternate'].upper()}"
    )


def _lirical_visible(variant):
    from mendelprime.tools.native import variant_is_candidate

    details = variant["details"]
    return bool(
        LIRICAL_ROW.fullmatch(str(details.get("raw_annotation") or "").strip())
    ) and variant_is_candidate(variant)


def _eligible_variants(gene, tool):
    if tool["tool_id"] == "exomiser":

        def priority(variant):
            details = variant["details"]
            pathogenic = str(details.get("acmg_classification") or "").upper() in {
                "PATHOGENIC",
                "LIKELY_PATHOGENIC",
            }
            return (
                details.get("contributing") is not True,
                not pathogenic,
                details.get("whitelist") is not True,
                -(details.get("variant_score") or 0.0),
                _reported_key(variant, tool),
            )

        return sorted(gene.get("variants", ()), key=priority)
    selected = {}
    for variant in gene.get("variants", ()):
        if not _lirical_visible(variant):
            continue
        key, score = _reported_key(variant, tool), variant["details"].get("pathogenicity") or 0
        if key not in selected or score > (selected[key]["details"].get("pathogenicity") or 0):
            selected[key] = variant
    return sorted(
        selected.values(),
        key=lambda item: (-(item["details"].get("pathogenicity") or 0), _reported_key(item, tool)),
    )


def _gene_record(record, tool_id):
    fields = LIRICAL_GENE if tool_id == "lirical" else EXOMISER_GENE
    result = {name: record.get(name) for name in fields}
    if tool_id == "exomiser":
        for name, maximum in (
            ("human_phenotype_evidence", 6000),
            ("mouse_phenotype_evidence", 3000),
            ("fish_phenotype_evidence", 3000),
            ("human_ppi_evidence", 2000),
        ):
            result[name] = str(record.get(name) or "")[:maximum]
    return result


def _family_view(observation, key):
    if observation is None:
        return None
    roles = {"proband": "", "mother": "", "father": ""}
    calls, by_role = {}, {call["role"]: call for call in observation["samples"]}
    for role, call in by_role.items():
        if role not in roles:
            continue
        called, carries, state = call_state(call, observation["alternate_index"])
        roles[role] = call["sample_id"]
        gq = call.get("gq")
        calls[role] = {
            "GT": call.get("gt"),
            "AD": list(call["ad"]) if call.get("ad") is not None else None,
            "DP": call.get("dp"),
            "GQ": int(gq) if gq is not None and float(gq).is_integer() else gq,
            "PL": list(call["pl"]) if call.get("pl") is not None else None,
            "state": call.get("state", state),
            "carries_alternate": call.get("carries_alternate", carries),
            "called": call.get("called", called),
        }
    origin = family_origin(observation["samples"], observation["alternate_index"])
    view = {
        "variant_key": key,
        "site_quality": observation.get("quality"),
        "filter": observation.get("filter_text"),
        "sample_roles": roles,
        "observations": calls,
        "origin": observation.get("origin", origin),
        "parentage_status": observation.get(
            "parentage_status", "unconfirmed" if {"father", "mother"} <= set(by_role) else "unavailable"
        ),
    }
    method = observation.get("match_method") or (
        "exact" if observation["variant_key"] == observation.get("vcf_variant_key") else "trimmed_padding"
    )
    if method != "exact":
        view["allele_match"] = {
            "method": method,
            "vcf_variant_key": variant_text(observation["vcf_variant_key"]),
            "alternate_index": observation["alternate_index"],
            "reference_checked": False,
        }
    return view


def _variant_view(variant, tool, observations):
    view = {
        "variant_key": _reported_key(variant, tool),
        "chromosome": variant["chromosome"].removeprefix("chr"),
        "position": variant["position"],
        "reference": variant["reference"].upper(),
        "alternate": variant["alternate"].upper(),
        "transcript": variant.get("transcript"),
        "hgvs_c": variant.get("hgvs_c"),
        "hgvs_p": variant.get("hgvs_p"),
        "genotype_raw": variant.get("genotype_raw"),
    }
    details = variant["details"]
    if tool["tool_id"] == "exomiser":
        for name in EXOMISER_VARIANT:
            view[name] = (
                "external_tool_unverified"
                if name == "acmg_assertion_status"
                else list(details.get(name, ()))
                if name == "acmg_evidence"
                else details.get(name)
            )
    else:
        view.update(
            pathogenicity=details.get("pathogenicity"),
            disease_id=details.get("disease_id"),
            raw_annotation=details.get("raw_annotation"),
        )
    view["source"] = tool["tool_id"]
    view["family_observation"] = _family_view(observations.get(variant.get("key")), view["variant_key"])
    return view


def _present(values):
    return {
        key: value for key, value in values.items() if value is not None and value not in ("", (), [], {})
    }


def document_view(document):
    source = document["source"]
    return {
        "source_id": source.source_id,
        "evidence_id": document["prompt_id"],
        "kind": source.kind,
        "provider": source.provider,
        "title": document["title"],
        "locator": source.locator,
        "text": document["prompt_text"],
        "version": source.version,
        "genes": list(source.genes),
        "variant_keys": [variant_text(key) for key in source.variant_keys],
        "parent_source_ids": list(source.parent_source_ids),
        "section": document.get("section", ""),
        "provenance": list(document.get("provenance", ())),
    }


def _recall_record(record, tool_id):
    fields = {
        key: value for key, value in _gene_record(record, tool_id).items() if not key.endswith("_evidence")
    }
    context = None
    if tool_id == "lirical":
        context = {
            "match_type": fields.pop("genotype_match_type"),
            "explanation": fields.pop("genotype_explanation"),
        }
    projected = {"source_id": record["source"].source_id, **_present(fields)}
    if tool_id == "lirical" and record.get("json_source") is not None:
        projected["genotype_context"] = {"source_id": record["json_source"].source_id, **_present(context)}
    if tool_id == "lirical" and record.get("gene_mapping_source") is not None:
        projected["gene_mapping_source_id"] = record["gene_mapping_source"].source_id
    return projected


def candidate_cards(native, genes, annotations=None, include_anchor=False, include_unprioritized=False):
    from mendelprime.tools.native import variant_is_candidate

    tools = native["tools"]
    by_tool = {tool["tool_id"]: {item["gene"]: item for item in tool["genes"]} for tool in tools}
    known = set().union(*(set(items) for items in by_tool.values()))
    if type(genes) is not tuple or len(set(genes)) != len(genes) or not set(genes) <= known:
        raise MendelPrimeError("unique_native_genes_required", field="genes")
    observed = {item["variant_key"]: item for item in native["observations"].get("observations", ())}
    family_refs = {
        (source.genes[0], key): source.source_id
        for source in native.get("family_sources", ())
        for key in source.variant_keys
    }
    cards = []
    for gene in genes:
        documents = (annotations or {}).get(gene)
        available = (*documents.get("documents", ()), *documents.get("omim_context", ())) if documents else ()
        annotated = {key for document in available for key in document["source"].variant_keys}
        records, alleles, unparsed = [], {}, []
        for tool in tools:
            native_gene = by_tool[tool["tool_id"]].get(gene)
            if native_gene is None:
                continue
            records.append(
                {
                    "tool": tool["tool_id"],
                    "native_rank": tool["ranking"].position(gene),
                    "records": [_recall_record(item, tool["tool_id"]) for item in native_gene["records"]],
                    **(
                        {
                            "raw_annotation_count": len(native_gene.get("variants", ())),
                            "computation_ineligible_annotation_count": sum(
                                not variant_is_candidate(item) for item in native_gene.get("variants", ())
                            ),
                        }
                        if include_unprioritized or annotations is not None
                        else {}
                    ),
                }
            )
            unparsed.extend(source.source_id for source in native_gene.get("unparsed_variants", ()))
            variants = (
                _eligible_variants(native_gene, tool)
                if tool["tool_id"] == "exomiser"
                else [
                    item
                    for item in native_gene.get("variants", ())
                    if variant_is_candidate(item)
                    and (include_unprioritized or _lirical_visible(item) or item.get("key") in annotated)
                ]
            )
            for variant in variants:
                observation = observed.get(variant.get("key"))
                raw_key = observation.get("vcf_variant_key") if observation else None
                identity = (
                    variant_text(raw_key or variant["key"])
                    if variant.get("key")
                    else f"unresolved:{variant['source'].source_id}"
                )
                family = _family_view(observation, variant_text(raw_key)) if observation else None
                if family is not None:
                    family.pop("sample_roles")
                    family["alternate_index"] = observation["alternate_index"]
                allele = alleles.setdefault(
                    identity,
                    {
                        "allele_identity": identity,
                        "vcf_allele_key": variant_text(raw_key) if raw_key else None,
                        "family_observation": family,
                        "family_source_ids": [],
                        "tool_annotations": {},
                    },
                )
                source_id = family_refs.get((gene, variant.get("key")))
                if source_id and source_id not in allele["family_source_ids"]:
                    allele["family_source_ids"].append(source_id)
                projected = _variant_view(variant, tool, {})
                compact = _present(
                    {
                        key: value
                        for key, value in projected.items()
                        if key
                        not in {
                            "chromosome",
                            "position",
                            "reference",
                            "alternate",
                            "family_observation",
                            "raw_annotation",
                        }
                    }
                )
                if tool["tool_id"] == "lirical" and not _lirical_visible(variant):
                    compact["source_ranking_eligible"] = False
                signature = json.dumps(compact, ensure_ascii=False, sort_keys=True)
                annotation = allele["tool_annotations"].setdefault(signature, {**compact, "source_ids": []})
                if variant["source"].source_id not in annotation["source_ids"]:
                    annotation["source_ids"].append(variant["source"].source_id)
        for allele in alleles.values():
            allele["tool_annotations"] = list(allele["tool_annotations"].values())
        papers = tuple(item for item in available if item["source"].kind == "literature")
        cards.append(
            {
                "gene": gene,
                "native_tool_records": records,
                "alleles": list(alleles.values()),
                "unparsed_variant_source_ids": unparsed,
                "additional_evidence": [
                    document_view(item) for item in available if item["source"].kind != "literature"
                ],
                "anchor_article": document_view(papers[0]) if include_anchor and papers else None,
                "retrieval_errors": list(documents.get("retrieval_errors", ())) if documents else [],
                "missing_information_is_uncertainty": True,
            }
        )
    return cards


def projection_source_ids(value):
    found = set()
    if isinstance(value, Mapping):
        for key, item in value.items():
            if key in {"source_id", "gene_mapping_source_id", "linked_source_id", "dependency_group_id"}:
                if not isinstance(item, str) or not item:
                    raise MendelPrimeError("source_id_required", field="molecular_context")
                found.add(item)
            elif key in {"source_ids", "family_source_ids", "source_refs"}:
                if not isinstance(item, (tuple, list)) or any(
                    not isinstance(ref, str) or not ref for ref in item
                ):
                    raise MendelPrimeError("source_id_sequence_required", field="molecular_context")
                found.update(item)
            elif isinstance(item, (Mapping, list, tuple)):
                found.update(projection_source_ids(item))
    elif isinstance(value, (list, tuple)):
        for item in value:
            found.update(projection_source_ids(item))
    return frozenset(found)


def molecular_data(document, case_key, gene):
    source = document.get("source")
    if source is None or source.kind != "derived":
        raise MendelPrimeError("derived_source_required", field="molecular_context")
    checks = (
        (source.case_key == case_key, "case_mismatch"),
        (source.genes == (gene,), "gene_scope_mismatch"),
        (document.get("prompt_id") == f"VCF:{gene}", "registered_context_required"),
        (document.get("prompt_text") == source.content, "context_content_mismatch"),
    )
    for valid, code in checks:
        if not valid:
            raise MendelPrimeError(code, field="molecular_context")
    try:
        data = json.loads(document["prompt_text"])
    except (ValueError, RecursionError):
        raise MendelPrimeError("context_json_required", field="molecular_context") from None
    if not isinstance(data, dict) or data.get("gene") != gene:
        raise MendelPrimeError("context_gene_mismatch", field="molecular_context")
    alleles, aliases = data.get("alleles", ()), []
    if not isinstance(alleles, (tuple, list)) or any(not isinstance(item, dict) for item in alleles):
        raise MendelPrimeError("allele_records_required", field="molecular_context")
    for allele in alleles:
        alias = allele.get("variant_ref")
        if alias is None:
            continue
        if not isinstance(alias, str) or re.fullmatch(rf"VAR:{re.escape(gene)}:[1-9]\d*", alias) is None:
            raise MendelPrimeError("variant_reference_scope_mismatch", field="molecular_context")
        if not projection_source_ids(allele):
            raise MendelPrimeError("variant_reference_source_missing", field="molecular_context")
        aliases.append(alias)
    if len(aliases) != len(set(aliases)):
        raise MendelPrimeError("duplicate_variant_reference", field="molecular_context")
    if projection_source_ids(data) != frozenset(source.parent_source_ids):
        raise MendelPrimeError("context_parent_mismatch", field="molecular_context")
    return data


def validate_molecular_context(context, case_key, genes):
    if any(gene not in context for gene in genes):
        raise MendelPrimeError("complete_gene_context_required", field="molecular_context")
    for gene in genes:
        molecular_data(context[gene], case_key, gene)


def _source_hypotheses(record, tool, variants):
    source = tool["tool_id"]
    hypotheses = []

    def add(mode, values):
        if values:
            hypotheses.append(
                {
                    "inheritance_mode": mode or "UNKNOWN",
                    "variant_keys": [item["variant_key"] for item in values],
                    "variant_refs": [item["variant_ref"] for item in values],
                    "compatibility": "tool_supported",
                    "phase": "unknown",
                    "source": source,
                    "source_ids": list(dict.fromkeys(item["source_id"] for item in values)),
                }
            )

    if source == "exomiser":
        groups = {}
        for variant in variants:
            if variant.get("contributing") is True:
                groups.setdefault(
                    str(variant.get("moi") or record["records"][0].get("moi") or "UNKNOWN").upper(), []
                ).append(variant)
        for mode, values in groups.items():
            heterozygous = [item for item in values if item.get("proband_call_status") == "heterozygous"]
            if ("AR" in mode or "RECESSIVE" in mode) and len(heterozygous) >= 2:
                add("AUTOSOMAL_RECESSIVE_COMP_HET", heterozygous[:2])
            else:
                for variant in values:
                    add(mode, [variant])
    else:
        first = record["records"][0]
        raw = (
            str(first.get("genotype_match_type") or "") + " " + str(first.get("genotype_explanation") or "")
        ).upper()
        recessive = bool(re.search(r"\bAR\b|RECESSIVE", raw))
        heterozygous = [item for item in variants if item.get("proband_call_status") == "heterozygous"]
        if recessive and len(heterozygous) >= 2:
            add("AUTOSOMAL_RECESSIVE_COMP_HET", heterozygous[:2])
        else:
            mode = (
                "AUTOSOMAL_RECESSIVE_HOM_ALT"
                if recessive
                and any(item.get("proband_call_status") == "homozygous_alternate" for item in variants)
                else "AUTOSOMAL_RECESSIVE"
                if recessive
                else "UNKNOWN"
            )
            for marker, label in (
                ("AD", "AUTOSOMAL_DOMINANT"),
                ("XD", "X_LINKED_DOMINANT"),
                ("XR", "X_LINKED_RECESSIVE"),
                ("MITO", "MITOCHONDRIAL"),
            ):
                if not recessive and marker in raw:
                    mode = label
                    break
            add(mode, variants)
    return hypotheses


def compact_source_cards(tool, observations):
    """Keep source annotations isolated; only the explicitly mapped proband supplies GT."""
    from mendelprime.tools.native import variant_is_candidate

    if observations.get("case_key") != tool["case_key"]:
        raise MendelPrimeError("case_mismatch", field="observations")
    by_key = {item["variant_key"]: item for item in observations.get("observations", ())}
    by_gene = {item["gene"]: item for item in tool["genes"]}
    cards = []
    for rank, gene in enumerate(tool["ranking"].genes[:100], 1):
        record = by_gene[gene]
        unique = {}
        for item in _eligible_variants(record, tool):
            unique.setdefault(item.get("key") or item["source"].source_id, item)
        variants = list(unique.values())
        selected = []
        for index, variant in enumerate(variants, 1):
            value = _variant_view(variant, tool, {})
            value.pop("family_observation", None)
            value.pop("raw_annotation", None)
            observation = by_key.get(variant.get("key"))
            calls = (
                [item for item in observation.get("samples", ()) if item.get("role") == "proband"]
                if observation
                else []
            )
            if len(calls) > 1:
                raise MendelPrimeError("ambiguous_proband", field="observations")
            proband = calls[0] if calls else None
            value["genotype_raw"] = proband.get("gt") if proband else None
            value["proband_call_status"] = (
                call_state(proband, observation["alternate_index"])[2] if proband else "unavailable"
            )
            value["quality_warning"] = (
                "call_quality_requires_downstream_assessment"
                if proband
                else "proband_observation_unavailable"
            )
            value["source_id"] = variant["source"].source_id
            value["variant_ref"] = f"{gene}:V{index}"
            value["phase"] = variant["details"].get("phase") or "unknown"
            selected.append(_present(value))
        hypotheses, selected_keys = [], set()
        for hypothesis in _source_hypotheses(record, tool, selected):
            keys = set(hypothesis["variant_keys"])
            if len(hypotheses) < 2 and len(selected_keys | keys) <= 2:
                hypotheses.append(hypothesis)
                selected_keys.update(keys)
        for variant in selected:
            if len(selected_keys) == 2:
                break
            selected_keys.add(variant["variant_key"])
        selected = [item for item in selected if item["variant_key"] in selected_keys]
        first = record["records"][0]
        cards.append(
            {
                "gene": gene,
                "raw_tool_rank": first.get("rank", rank),
                "raw_tool_score": first.get("combined_score")
                if tool["tool_id"] == "exomiser"
                else first.get("posttest_probability"),
                "tool_annotations": _present(
                    {
                        key: first.get(key)
                        for key in (
                            "moi",
                            "phenotype_score",
                            "variant_score",
                            "disease_id",
                            "disease_name",
                            "composite_lr",
                            "genotype_match_type",
                        )
                    }
                ),
                "variant_context": {
                    "status": "prioritized_variants" if selected else "not_assessed",
                    "negative_evidence_eligible": False,
                    "variants": selected,
                    "genotype_hypotheses": hypotheses,
                    "omitted_variant_count": max(0, len(variants) - len(selected)),
                    "raw_annotation_count": len(record.get("variants", ())),
                    "computation_ineligible_annotation_count": sum(
                        not variant_is_candidate(item) for item in record.get("variants", ())
                    ),
                    "limitations": [
                        "Only this source's variants are shown; missing information is uncertainty."
                    ],
                },
            }
        )
    return sorted(cards, key=lambda item: item["gene"])


def compact_prompt_payload(value, preserve_strings=()):
    """Factor repeated JSON structure without dropping or rewriting scientific facts."""
    candidate_arrays = {"candidate_cards_alphabetical", "candidates_in_admission_order"}
    reserved = {"$mp_s", "$mp_r", "$mp_t", "$mp_o"}
    normalized = _ranking_value(value)

    def collision(item):
        if isinstance(item, Mapping):
            return (
                set(item) == {"compact_representation"}
                or (len(item) == 1 and next(iter(item)) in reserved)
                or any(collision(child) for child in item.values())
            )
        return isinstance(item, (tuple, list)) and any(collision(child) for child in item)

    if collision(normalized):
        return normalized

    plain_ids, candidate_ids = set(), set()

    def mark_plain(item):
        if isinstance(item, Mapping):
            plain_ids.add(id(item))
            for child in item.values():
                mark_plain(child)
        elif isinstance(item, (tuple, list)):
            plain_ids.add(id(item))
            for child in item:
                mark_plain(child)

    def packed(item, path=()):
        if isinstance(item, Mapping):
            if len(path) == 1 and path[0] not in candidate_arrays:
                copied = _ranking_value(item)
                mark_plain(copied)
                return copied
            result = {key: packed(child, (*path, key)) for key, child in item.items()}
            if len(path) == 2 and path[0] in candidate_arrays:
                candidate_ids.add(id(result))
            return result
        if not isinstance(item, (tuple, list)):
            return item
        if len(path) == 1 and path[0] not in candidate_arrays:
            copied = _ranking_value(item)
            mark_plain(copied)
            return copied
        rows = [packed(child, (*path, index)) for index, child in enumerate(item)]
        if (
            len(path) == 1
            and path[0] in candidate_arrays
            or len(rows) < 2
            or not all(isinstance(row, dict) for row in rows)
        ):
            return rows
        columns = list(rows[0])
        if any(set(row) != set(columns) for row in rows[1:]):
            return rows
        table = {"$mp_t": {"columns": columns, "rows": [[row[column] for column in columns] for row in rows]}}
        return (
            table
            if len(json.dumps(table, ensure_ascii=False, separators=(",", ":")))
            < len(json.dumps(rows, ensure_ascii=False, separators=(",", ":")))
            else rows
        )

    encoded = packed(normalized)
    preserved = set(preserve_strings)
    for item in encoded.values():
        if isinstance(item, str):
            preserved.add(item)
    for name in candidate_arrays:
        for card in encoded.get(name, ()):
            if isinstance(card, Mapping) and isinstance(card.get("gene"), str):
                preserved.add(card["gene"])
    record_counts, record_order, record_values, signatures = {}, [], {}, {}

    def fingerprint(item):
        if isinstance(item, Mapping):
            digest = hashlib.sha256(b"mapping")
            for key, child in item.items():
                digest.update(json.dumps(key, ensure_ascii=False).encode())
                digest.update(b"\0")
                digest.update(fingerprint(child).encode())
                digest.update(b"\0")
            signature = digest.hexdigest()
            signatures[id(item)] = signature
            if id(item) not in plain_ids and id(item) not in candidate_ids:
                if signature not in record_counts:
                    record_counts[signature] = 0
                    record_order.append(signature)
                    record_values[signature] = item
                record_counts[signature] += 1
            return signature
        if isinstance(item, (tuple, list)):
            digest = hashlib.sha256(b"sequence")
            for child in item:
                digest.update(fingerprint(child).encode())
                digest.update(b"\0")
            return digest.hexdigest()
        return hashlib.sha256(
            (type(item).__name__ + ":" + json.dumps(item, ensure_ascii=False, sort_keys=True)).encode()
        ).hexdigest()

    fingerprint(encoded)
    record_aliases = {}
    for signature in record_order:
        record_size = len(
            json.dumps(record_values[signature], ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        )
        alias = f"R{len(record_aliases) + 1:03d}"
        reference_size = len(json.dumps({"$mp_r": alias}, separators=(",", ":")))
        catalog_size = len(
            json.dumps({alias: record_values[signature]}, ensure_ascii=False, separators=(",", ":"))
        )
        if (
            record_size >= 120
            and record_counts[signature] > 1
            and (
                catalog_size + record_counts[signature] * reference_size
                < record_counts[signature] * record_size
            )
        ):
            record_aliases[signature] = alias

    def record_references(item, skip=None):
        if isinstance(item, Mapping):
            if id(item) in plain_ids:
                return item
            signature = signatures[id(item)]
            if id(item) not in candidate_ids and signature != skip and signature in record_aliases:
                return {"$mp_r": record_aliases[signature]}
            return {key: record_references(child, skip) for key, child in item.items()}
        if isinstance(item, (tuple, list)):
            return [record_references(child, skip) for child in item]
        return item

    encoded = record_references(encoded)
    records = {
        record_aliases[signature]: record_references(record_values[signature], signature)
        for signature in record_order
        if signature in record_aliases
    }
    occurrences, order = {}, []

    def collect(item):
        if id(item) in plain_ids:
            return
        if isinstance(item, str) and len(item) >= 24:
            if item not in occurrences:
                occurrences[item] = 0
                order.append(item)
            occurrences[item] += 1
        elif isinstance(item, Mapping):
            for child in item.values():
                collect(child)
        elif isinstance(item, (tuple, list)):
            for child in item:
                collect(child)

    collect(encoded)
    collect(records)
    aliases = {}
    for text in order:
        alias = f"S{len(aliases) + 1:03d}"
        reference = {"$mp_s": alias}
        original_size = len(json.dumps(text, ensure_ascii=False)) * occurrences[text]
        compact_size = len(json.dumps(reference, separators=(",", ":"))) * occurrences[text] + len(
            json.dumps({alias: text}, ensure_ascii=False, separators=(",", ":"))
        )
        if text not in preserved and occurrences[text] > 1 and compact_size < original_size:
            aliases[text] = alias

    def referenced(item):
        if id(item) in plain_ids:
            return item
        if isinstance(item, str) and item in aliases:
            return {"$mp_s": aliases[item]}
        if isinstance(item, Mapping):
            return {key: referenced(child) for key, child in item.items()}
        if isinstance(item, (tuple, list)):
            return [referenced(child) for child in item]
        return item

    compact = referenced(encoded)
    contract = {
        "schema": "vcf_prompt_compact_v1",
        "instructions": (
            "Resolve each $mp_s through shared_strings, each $mp_r through shared_records, and each "
            "$mp_t row by pairing it with columns. "
            "These are lossless presentation references; every candidate and scientific fact remains supplied."
        ),
        "shared_strings": {alias: text for text, alias in aliases.items()},
        "shared_records": referenced(records),
    }
    compact["compact_representation"] = contract
    return compact


def compact_prompt_columns(value):
    """Share repeated object field lists when ordinary lossless packing is insufficient."""
    contract = value.get("compact_representation", {})
    if contract.get("schema") != "vcf_prompt_compact_v1":
        return value
    arrays = ("candidate_cards_alphabetical", "candidates_in_admission_order")
    counts, stack = {}, [value[name] for name in arrays if name in value]
    while stack:
        item = stack.pop()
        if isinstance(item, dict):
            if len(item) >= 3 and "gene" not in item:
                keys = tuple(item)
                counts[keys] = counts.get(keys, 0) + 1
            stack.extend(item.values())
        elif isinstance(item, list):
            stack.extend(item)
    aliases = {}
    for keys, count in counts.items():
        alias = f"C{len(aliases) + 1:03d}"
        before = json.dumps(dict.fromkeys(keys), separators=(",", ":"))
        after = json.dumps({"$mp_o": [alias, [None] * len(keys)]}, separators=(",", ":"))
        if count > 1 and count * (len(before) - len(after)) > len(json.dumps({alias: keys})):
            aliases[keys] = alias
    if not aliases:
        return value

    def packed(item):
        if isinstance(item, dict):
            values = {key: packed(child) for key, child in item.items()}
            alias = aliases.get(tuple(item)) if "gene" not in item else None
            return {"$mp_o": [alias, list(values.values())]} if alias else values
        return [packed(child) for child in item] if isinstance(item, list) else item

    return {
        **value,
        **{name: packed(value[name]) for name in arrays if name in value},
        "compact_representation": {
            **contract,
            "schema": "vcf_prompt_compact_v2",
            "shared_columns": {alias: list(keys) for keys, alias in aliases.items()},
            "instructions": contract["instructions"]
            + (
                " Resolve each $mp_o [column_set, values] by pairing shared_columns[column_set] "
                "with values to reconstruct one object. Omitted fields remain omitted; nothing is truncated."
            ),
        },
    }


def expand_prompt_payload(value):
    """Restore a compact VCF Prompt payload for internal reference validation."""
    contract = value.get("compact_representation") if isinstance(value, Mapping) else None
    if not isinstance(contract, Mapping) or contract.get("schema") not in {
        "vcf_prompt_compact_v1",
        "vcf_prompt_compact_v2",
    }:
        return value
    strings = contract.get("shared_strings", {})
    records = contract.get("shared_records", {})

    def expanded(item):
        if isinstance(item, Mapping) and set(item) == {"$mp_o"} and "shared_columns" in contract:
            alias, values = item["$mp_o"]
            return dict(zip(contract["shared_columns"][alias], expanded(values), strict=True))
        if isinstance(item, Mapping) and set(item) == {"$mp_s"}:
            return strings[item["$mp_s"]]
        if isinstance(item, Mapping) and set(item) == {"$mp_r"}:
            return expanded(records[item["$mp_r"]])
        if isinstance(item, Mapping) and set(item) == {"$mp_t"}:
            table = expanded(item["$mp_t"])
            return [dict(zip(table["columns"], row, strict=True)) for row in table["rows"]]
        if isinstance(item, Mapping):
            return {key: expanded(child) for key, child in item.items() if key != "compact_representation"}
        if isinstance(item, (tuple, list)):
            return [expanded(child) for child in item]
        return item

    return expanded(value)


def _ranking_value(value):
    if isinstance(value, Source):
        return {"source_id": value.source_id}
    if isinstance(value, Mapping):
        return {
            key: _ranking_value(item)
            for key, item in value.items()
            if key not in {"raw_record", "raw_fields", "sample_id"}
        }
    if isinstance(value, (tuple, list)):
        return [_ranking_value(item) for item in value]
    return value


def ranking_assessment(ledger, case_key, gene):
    if ledger is None:
        return None
    if ledger.get("schema_version") != "ranking_acmg_v1" or ledger.get("case_key") != case_key:
        raise MendelPrimeError("ranking_acmg_scope_mismatch", field="ranking_acmg")
    assessment = ledger.get("by_gene", {}).get(gene)
    if assessment is not None and assessment.get("gene") != gene:
        raise MendelPrimeError("gene_scope_mismatch", field="ranking_acmg")
    return assessment


def _structured_gene_sources(gene, targets, clinical, registry):
    roles = {}
    for target in targets:
        for binding in target.get("bindings", ()):
            if binding["role"] in {"gene_context", "gene_disease_context", "native_annotation"}:
                roles.setdefault(binding["source_id"], set()).add(binding["role"])
    for document in clinical.get("gene_context", ()):
        roles.setdefault(document["source"].source_id, set()).add(
            "gene_disease_context" if document["source"].kind == "gene_disease" else "gene_context"
        )
    pending = list(roles)
    for source_id in pending:
        for parent_id in registry[source_id].parent_source_ids:
            if parent_id not in roles:
                roles[parent_id] = set(roles[source_id])
                pending.append(parent_id)
    validity, constraint, features = [], [], []

    def constraint_fields(value):
        if not isinstance(value, Mapping):
            return {}
        found = {
            key: _ranking_value(item)
            for key, item in value.items()
            if str(key).lower()
            in {
                "loeuf",
                "loeuf_lower",
                "loeuf_upper",
                "pli",
                "oe_lof",
                "oe_lof_upper",
                "lof_oe",
                "lof.oe",
                "lof.oe_ci.upper",
                "lof.oe_ci.lower",
                "lof_oe_ci_upper",
                "lof_oe_ci_lower",
            }
            or str(key).lower().startswith("gene_constraint")
        }
        for item in value.values():
            if isinstance(item, Mapping):
                found.update(constraint_fields(item))
        return found

    for source_id, source_roles in roles.items():
        source = registry[source_id]
        if source.genes and gene not in source.genes:
            raise MendelPrimeError("gene_scope_mismatch", field="structured_gene_evidence")
        try:
            facts = json.loads(source.content)
        except ValueError:
            facts = {"text": source.content}
        if not isinstance(facts, Mapping):
            facts = {"record": facts}
        header, row = facts.get("header", facts.get("headers")), facts.get("row")
        if isinstance(header, (tuple, list)) and isinstance(row, str):
            values = next(csv.reader([row], delimiter="\t"))
            if len(header) == len(values):
                facts = {**facts, "decoded_row": dict(zip(header, values, strict=True))}
        record = {
            "source_id": source_id,
            "provider": source.provider,
            "version": source.version,
            "facts": _ranking_value(facts),
        }
        if "gene_disease_context" in source_roles:
            validity.append(record)
        if values := constraint_fields(facts):
            constraint.append({**record, "facts": values})
        if isinstance(facts.get("ranking_features"), Mapping):
            features.append({**record, "facts": _ranking_value(facts["ranking_features"])})
    return {"gene_disease_assertions": validity, "gene_constraint": constraint, "ranking_features": features}


def ranking_context_views(case_key, genes, context):
    """Project registered facts without creating or changing any Source identity."""
    from mendelprime.tools.native import variant_is_candidate

    prepared, ranking_acmg, evidence = (context[name] for name in ("prepared", "acmg", "evidence"))
    packet = prepared["input"]
    if case_tuple(packet["profile"]["case_key"]) != case_key:
        raise MendelPrimeError("case_mismatch", field="prepared")
    registry = {source.source_id: source for source in iter_sources((prepared, evidence))}
    output = {}
    for gene in genes:
        card = {"gene": gene}
        assessment = ranking_assessment(ranking_acmg, case_key, gene)
        variant_refs = {}
        for variant in (assessment or {}).get("variants", ()):
            identity = variant["identity"]
            identity_key = (identity["gene"], tuple(identity["variant_key"]))
            if identity_key in variant_refs and variant_refs[identity_key] != variant["variant_ref"]:
                raise MendelPrimeError("ambiguous_allele_reference", field="ranking_acmg")
            variant_refs[identity_key] = variant["variant_ref"]
        raw_targets = tuple(item for item in packet["targets"] if item["primary"]["gene"] == gene)
        targets = tuple(
            item
            for item in raw_targets
            if any(
                variant_is_candidate(native)
                for native in (item["primary"], *item.get("aliases", ()))
            )
        )
        eligible_keys = {tuple(target["primary"]["key"]) for target in targets}
        ineligible_keys = {tuple(target["primary"]["key"]) for target in raw_targets} - eligible_keys
        card["raw_prepared_variant_count"] = len(raw_targets)
        card["computation_ineligible_variant_count"] = len(raw_targets) - len(targets)
        pool = evidence.get(gene, {"gene": gene, "documents": (), "omim_context": (), "retrieval_errors": ()})
        if pool.get("gene") != gene:
            raise MendelPrimeError("gene_scope_mismatch", field="evidence")
        for document in (*pool.get("documents", ()), *pool.get("omim_context", ())):
            source = document["source"]
            if source.case_key not in {None, case_key} or source.genes and gene not in source.genes:
                raise MendelPrimeError("document_scope_mismatch", field="evidence")
        card["prepared_variants"] = [
            _ranking_value(
                {
                    "variant_ref": variant_refs.get((gene, tuple(target["primary"]["key"]))),
                    "variant_key": target["primary"]["key"],
                    "reference": target["reference"],
                    "transcript": target["primary"].get("transcript"),
                    "hgvs_c": target["primary"].get("hgvs_c"),
                    "hgvs_p": target["primary"].get("hgvs_p"),
                    **{
                        key: target.get(key)
                        for key in (
                            "observation",
                            "population",
                            "annotation",
                            "predictions",
                            "pvs1",
                            "bindings",
                            "inheritance_hypotheses",
                        )
                    },
                }
            )
            for target in targets
        ]
        card["ranking_acmg"] = _ranking_value(assessment)
        if assessment:
            card["ranking_acmg"]["variants"] = [
                variant
                for variant in card["ranking_acmg"]["variants"]
                if tuple(variant["identity"]["variant_key"]) not in ineligible_keys
            ]
            card["ranking_acmg"]["computation_ineligible_variant_count"] = len(assessment["variants"]) - len(
                card["ranking_acmg"]["variants"]
            )
            card["ranking_acmg"]["_source_variant_keys"] = {
                source_id: registry[source_id].variant_keys
                for source_id in projection_source_ids(assessment.get("genotype_configuration", {}))
                if source_id in registry
            }
        clinical = prepared.get("clinical", {}).get(gene, {})
        if clinical and clinical.get("gene") != gene:
            raise MendelPrimeError("gene_scope_mismatch", field="prepared.clinical")
        card["clinical_context"] = {
            "gene": gene,
            "gene_context": [document_view(document) for document in clinical.get("gene_context", ())],
            "variant_assertions": [
                {
                    **_ranking_value(
                        {
                            key: value
                            for key, value in item.items()
                            if key != "document"
                            and not (
                                key in {"linked_source_id", "gene_mapping_source_id", "dependency_group_id"}
                                and not value
                            )
                        }
                    ),
                    "document": document_view(item["document"]),
                    "patient_use": "eligible_observation_requires_quality_review"
                    if tuple(item["variant_key"]) in eligible_keys
                    else "database_context_only",
                    "eligible_variant_ref": variant_refs.get((gene, tuple(item["variant_key"])))
                    if tuple(item["variant_key"]) in eligible_keys
                    else None,
                    "database_assertion_does_not_establish_patient_carriage": True,
                }
                for item in clinical.get("variant_evidence", ())
            ],
            "lookups": _ranking_value(clinical.get("lookups", ())),
            "missing_information_is_not_benign": True,
        }
        card["structured_gene_evidence"] = _structured_gene_sources(gene, targets, clinical, registry)
        card["unusable_observations"] = [
            {
                "identity": _ranking_value(exclusion[0]),
                "reason": exclusion[1],
                "facts": _ranking_value(exclusion[2]) if len(exclusion) == 3 else None,
                "patient_use": "not_eligible_for_assessment",
                "missing_or_unusable_is_not_benign": True,
            }
            for exclusion in prepared.get("excluded", ())
            if exclusion[0]["gene"] == gene
        ]
        card["ranking_acmg_phase"] = ranking_acmg["phase"]
        card["ranking_literature"] = {
            "documents": [
                {
                    **document_view(item),
                    "metadata": _ranking_value(item.get("metadata", {})),
                    "publication_types": list(item.get("publication_types", ())),
                }
                for item in pool.get("documents", ())
            ],
            "omim_context": [
                {
                    **document_view(item),
                    "metadata": _ranking_value(item.get("metadata", {})),
                    "cited_pmids": list(item.get("cited_pmids", ())),
                }
                for item in pool.get("omim_context", ())
            ],
            "retrieval_errors": list(pool.get("retrieval_errors", ())),
        }
        parents = set(projection_source_ids(card))
        if assessment:
            parents.update(
                ref
                for variant in assessment["variants"]
                for criterion in variant["criteria"]
                for ref in criterion["source_refs"]
            )
        if parents - set(registry):
            raise MendelPrimeError("ranking_source_missing", field="molecular_context")
        for reference in parents:
            source = registry[reference]
            if source.case_key not in {None, case_key} or source.genes and gene not in source.genes:
                raise MendelPrimeError("ranking_source_scope_mismatch", field="molecular_context")
        output[gene] = card
    return output


def variant_projection(variant):
    """A bounded display record; the canonical variant and its criteria remain untouched."""
    value = {
        key: _ranking_value(variant[key])
        for key in (
            "variant_ref",
            "identity",
            "classification",
            "independent_classification",
            "points_for_audit_only",
            "conflicting_evidence",
        )
        if key in variant
    }
    omitted = {}
    for name, limit in (("disease_contexts", 2), ("external_assertions", 2), ("unresolved", 4)):
        rows = variant.get(name, ())
        value[name] = _ranking_value(rows[:limit])
        if len(rows) > limit:
            omitted[name] = len(rows) - limit
    value["omitted_counts"] = omitted
    return value


def clinical_document_projection(document):
    """Expose named clinical/resource facts without repeating the raw source document."""
    value = {
        name: document[name]
        for name in ("source_id", "evidence_id", "provider", "kind", "version", "variant_keys")
        if name in document
    }
    try:
        data = json.loads(document.get("text", ""))
    except (TypeError, ValueError):
        value.update(text=document.get("text", ""), facts={})
        return value
    if (
        isinstance(data, Mapping)
        and data.get("record_type") == "verified_sequence_equivalent_clinical_evidence"
    ):
        try:
            data = json.loads(data["original_content"])
        except (KeyError, TypeError, ValueError):
            pass
    facts = {}
    omitted = 0
    fields = {
        "classification",
        "clinicalsignificance",
        "assertion",
        "reviewstatus",
        "reviewstars",
        "expertpanel",
        "disease",
        "diseases",
        "diseaseid",
        "diseaseidmondo",
        "diseaselabel",
        "diseasename",
        "mondoid",
        "modeofinheritance",
        "moi",
        "clinvarclassification",
        "clinvarreviewstatus",
        "clinvarvariationid",
        "clinicalsignificancedescription",
        "clnsig",
        "clnrevstat",
        "clndn",
        "clndisdb",
        "clnhgvs",
        "alleleid",
        "variationid",
        "pmids",
        "pmid",
        "appliedevidencecodesmet",
        "appliedevidencecodesnotmet",
        "retracted",
        "retractionstatus",
        "withdrawn",
        "status",
        "loeuf",
        "loeuflower",
        "loeufupper",
        "pli",
        "oelof",
        "oelofupper",
        "lofoe",
        "lofoeciupper",
        "lofoecilower",
    }

    def visit(item, path=""):
        nonlocal omitted
        if not isinstance(item, Mapping):
            return
        header, row = item.get("header", item.get("headers")), item.get("row")
        if isinstance(header, (list, tuple)) and isinstance(row, str):
            values = next(csv.reader([row], delimiter="\t"))
            if len(header) == len(values):
                visit(dict(zip(header, values, strict=True)), path + "row.")
        for name, raw in item.items():
            normalized = re.sub(r"[^a-z0-9]", "", name.lower())
            if normalized in fields or normalized.startswith("geneconstraint"):
                facts[path + name] = _ranking_value(raw)
            elif isinstance(raw, Mapping):
                visit(raw, path + name + ".")
            else:
                omitted += 1

    visit(data)
    value.update(facts=facts, omitted_resource_field_count=omitted)
    return value


def prepared_variant_projection(variant, *, include_observation=True):
    value = _present(
        {
            name: variant.get(name)
            for name in (
                "variant_ref",
                "variant_key",
                "transcript",
                "hgvs_c",
                "hgvs_p",
                "inheritance_hypotheses",
            )
        }
    )
    reference = variant.get("reference")
    value["reference"] = (
        {
            name: reference[name]
            for name in (
                "status",
                "normalized_hgvs",
                "reference_validation",
                "normalization_method",
                "issues",
            )
            if name in reference
        }
        if isinstance(reference, Mapping)
        else reference
    )

    def resource(item):
        if isinstance(item, Mapping):
            return {
                key: resource(raw)
                for key, raw in item.items()
                if key not in {"target", "requested_key", "original_key", "normalized_key", "variant_key"}
                and raw is not None
                and raw not in ("", (), [], {})
            }
        if isinstance(item, (list, tuple)):
            return [resource(raw) for raw in item]
        return item

    for name in ("population", "annotation", "predictions", "pvs1"):
        value[name] = resource(variant.get(name))
    bindings = variant.get("bindings") or ()
    value["bindings"] = [
        {name: row[name] for name in ("source_id", "role", "match_scope") if name in row}
        for row in bindings
        if row.get("role") == "patient_observation"
    ]
    value["omitted_binding_count"] = len(bindings) - len(value["bindings"])
    observation = variant.get("observation") or {}
    value["observation"] = {
        name: observation[name]
        for name in (
            "status",
            "origin",
            "parentage_status",
            "site_passed",
            "quality_status",
            "family_conflict",
            "limitations",
        )
        if name in observation
    }
    if include_observation:
        calls = []
        for call in observation.get("calls", ()):
            original = call.get("original", {})
            projected = {
                name: call.get(name, original.get(name))
                for name in (
                    "role",
                    "gt",
                    "ad",
                    "dp",
                    "gq",
                    "pl",
                    "genotype_state",
                    "quality_state",
                    "usable",
                    "selected_alt_dosage",
                )
                if name in call or name in original
            }
            calls.append(projected)
        if calls:
            value["observation"]["calls"] = calls
    else:
        value["observation"]["call_details_in_selected_genotype"] = True
    return value


def genotype_projection(assessment, displayed_variants):
    """Scope every genotype collection to the alleles actually shown by this stage."""
    configuration = assessment.get("genotype_configuration", {})
    refs = tuple(item["variant_ref"] for item in displayed_variants)
    keys = {tuple(item["identity"]["variant_key"]) for item in displayed_variants}
    ref_set = set(refs)
    keys_by_ref = {item["variant_ref"]: list(item["identity"]["variant_key"]) for item in displayed_variants}
    value = {
        name: _ranking_value(configuration[name])
        for name in (
            "inheritance_mode",
            "normalized_inheritance_mode",
            "compatibility",
            "genotype_status",
            "phase",
            "selected_variant_refs_role",
            "hypothesis_selection_reason",
        )
        if name in configuration
    }
    omitted = {}
    selected = tuple(
        configuration.get("selected_variant_refs") or assessment.get("selected_variant_refs", ())
    )
    unresolved = configuration.get("selected_variant_refs_role") == "observed_unresolved_candidates"
    value["selected_variant_refs"] = [] if unresolved else [ref for ref in selected if ref in ref_set]
    if not unresolved and len(selected) > len(value["selected_variant_refs"]):
        omitted["selected_variant_refs"] = len(selected) - len(value["selected_variant_refs"])
    observed = tuple(
        configuration.get("observed_variant_refs") or (item["variant_ref"] for item in assessment["variants"])
    )
    value["observed_variant_refs"] = [ref for ref in observed if ref in ref_set]
    value["omitted_observed_variant_count"] = len(observed) - len(value["observed_variant_refs"])
    for name in ("family_observations", "parental_origins", "proband_states"):
        source = configuration.get(name, {})
        value[name] = {ref: _ranking_value(source[ref]) for ref in refs if ref in source}
        if name == "family_observations":
            call_fields = {
                "genotype": "GT",
                "allele_depths": "AD",
                "depth": "DP",
                "genotype_quality": "GQ",
                "pl": "PL",
            }
            value[name] = {
                ref: {
                    role: {
                        call_fields.get(key, key): raw
                        for key, raw in call.items()
                        if key != "role" and raw is not None
                    }
                    for role, call in calls.items()
                }
                for ref, calls in value[name].items()
            }
        if len(source) > len(value[name]):
            omitted[name] = len(source) - len(value[name])
    for name in ("homozygous_alternate_refs", "hemizygous_alternate_refs"):
        source = configuration.get(name, ())
        value[name] = [ref for ref in source if ref in ref_set]
        if len(source) > len(value[name]):
            omitted[name] = len(source) - len(value[name])

    def hypothesis_view(hypothesis):
        result = {
            name: _ranking_value(hypothesis[name])
            for name in (
                "inheritance_mode",
                "normalized_inheritance_mode",
                "compatibility",
                "phase",
                "source",
                "status",
                "hypothesis_index",
                "kind",
            )
            if name in hypothesis
        }
        dropped = {}
        for name in ("variant_refs", "selected_variant_refs"):
            source = hypothesis.get(name, ())
            if name in hypothesis:
                result[name] = [ref for ref in source if ref in ref_set]
                dropped[name] = len(source) - len(result[name])
        for name in ("variant_keys", "observed_variant_keys"):
            source = hypothesis.get(name, ())
            if name in hypothesis:
                result[name] = _ranking_value([key for key in source if tuple(key) in keys])
                dropped[name] = len(source) - len(result[name])
        source = hypothesis.get("variant_disease_modes", ())
        if "variant_disease_modes" in hypothesis:
            result["variant_disease_modes"] = _ranking_value(
                [row for row in source if row.get("variant_ref") in ref_set]
            )
            dropped["variant_disease_modes"] = len(source) - len(result["variant_disease_modes"])
        for name, limit in (
            ("source_ids", 3),
            ("disease_ids", 2),
            ("missing_variant_keys", 4),
            ("conflicts", 4),
            ("uncertainties", 4),
        ):
            source = hypothesis.get(name, ())
            if name in hypothesis:
                candidates = source
                if name == "source_ids":
                    source_keys = assessment.get("_source_variant_keys", {})
                    candidates = [
                        source_id
                        for source_id in source
                        if not source_keys.get(source_id)
                        or any(tuple(key) in keys for key in source_keys[source_id])
                    ]
                result[name] = _ranking_value(candidates[:limit])
                dropped[name] = len(source) - len(result[name])
        result["omitted_counts"] = {name: count for name, count in dropped.items() if count}
        if result.get("observed_variant_keys") == result.get("variant_keys"):
            result.pop("observed_variant_keys", None)
        if result.get("variant_refs") and result.get("variant_keys") == [
            keys_by_ref[ref] for ref in result["variant_refs"]
        ]:
            result.pop("variant_keys")
        result = {
            name: raw for name, raw in result.items() if raw is not None and raw not in ("", (), [], {})
        }
        return result

    value["selected_hypothesis"] = (
        hypothesis_view(configuration["selected_hypothesis"])
        if configuration.get("selected_hypothesis")
        else {}
    )
    value["phase_evidence"] = (
        hypothesis_view(configuration["phase_evidence"]) if configuration.get("phase_evidence") else {}
    )
    conflicts = configuration.get("hypothesis_conflicts", ())
    value["hypothesis_conflicts"] = [hypothesis_view(row) for row in conflicts[:4]]
    if len(conflicts) > 4:
        omitted["hypothesis_conflicts"] = len(conflicts) - 4
    notes = configuration.get("unresolved", ())
    value["unresolved"] = _ranking_value(notes[:4])
    if len(notes) > 4:
        omitted["unresolved"] = len(notes) - 4
    value["omitted_counts"] = omitted
    value["displayed_variant_refs"] = list(refs)
    return value


def layered_evidence(assessment, *, recall=False):
    if assessment is None:
        return None
    configuration = assessment.get("genotype_configuration", {})
    selected = tuple(
        configuration.get("selected_variant_refs") or assessment.get("selected_variant_refs", ())
    )
    unresolved_selection = configuration.get("selected_variant_refs_role") == "observed_unresolved_candidates"
    by_ref = {item["variant_ref"]: item for item in assessment["variants"]}
    selected_variants = (
        [by_ref[ref] for ref in selected if ref in by_ref][:2] if not unresolved_selection else []
    )
    status = configuration.get("genotype_status", "uncertain")
    variants = selected_variants if selected_variants else list(assessment["variants"][:4]) if recall else []
    decisive, supporting, counter, expert, external = [], [], [], [], []
    for variant in variants:
        for criterion in variant["criteria"]:
            if criterion["status"] != "met" or criterion.get("counted") is not True:
                continue
            item = {
                "variant_ref": variant["variant_ref"],
                **{
                    key: _ranking_value(criterion.get(key))
                    for key in (
                        "code",
                        "strength",
                        "direction",
                        "independence_group",
                    )
                },
                "source_refs": list(criterion.get("source_refs", ()))[:3],
                "source_excerpt": str(next(iter(criterion.get("excerpts", ())), ""))[:500],
                "limitations": list(criterion.get("limitations", ()))[:2],
            }
            code = criterion["code"].upper()
            facts = criterion.get("raw_facts", {})
            mode = str(configuration.get("inheritance_mode", "")).upper()
            if criterion["direction"] == "benign" and (
                code in {"BA1", "BS1", "BS2", "BS3", "BS4"}
                or criterion["strength"] in {"strong", "very_strong", "stand_alone"}
            ):
                counter.append(item)
            elif criterion["direction"] == "pathogenic" and (
                code in {"PS1", "PS3", "PS4"}
                or code == "PS2"
                and mode
                in {
                    "AD",
                    "X",
                    "XR",
                    "XD",
                    "AUTOSOMAL_DOMINANT",
                    "X_LINKED_RECESSIVE",
                    "X_LINKED_DOMINANT",
                    "X_LINKED",
                }
                and facts.get("de_novo_confirmed_count", 0) > 0
                or code == "PM3"
                and mode
                in {
                    "AR",
                    "AUTOSOMAL_RECESSIVE",
                    "AUTOSOMAL_RECESSIVE_COMP_HET",
                    "AUTOSOMAL_RECESSIVE_HOM_ALT",
                }
                and status == "complete"
            ):
                decisive.append(item)
            else:
                supporting.append(item)
        diseases = {
            item.get("disease_id") for item in variant.get("disease_contexts", ()) if item.get("disease_id")
        }
        for assertion in variant.get("external_assertions", ()):
            projected = {"variant_ref": variant["variant_ref"], **_ranking_value(assertion)}
            external.append(projected)
            classification = str(assertion.get("classification", "")).lower().replace(" ", "_")
            if (
                assertion.get("expert_panel")
                and assertion.get("disease_id") in diseases
                and classification in {"pathogenic", "likely_pathogenic"}
                and not assertion.get("withdrawn")
                and not assertion.get("retracted")
                and not assertion.get("conflicting_evidence")
            ):
                expert.append(projected)
    if recall:
        return {
            "gene": assessment["gene"],
            "purpose": "candidate_recall_not_final_causal_discrimination",
            "selected_variant_refs": [item["variant_ref"] for item in selected_variants],
            "genotype_status": status,
            "variant_summaries": [variant_projection(item) for item in variants],
            "human_or_expert_recall_signals": [*decisive, *expert],
            "hard_counterevidence": counter,
            "external_variant_assertions": external,
            "supporting_signal_codes": sorted({item["code"] for item in supporting}),
            "missing_information_is_not_negative": True,
        }
    review_variants = (
        ([by_ref[ref] for ref in selected if ref in by_ref] or list(assessment["variants"]))[:4]
        if unresolved_selection
        else selected_variants
    )
    genotype = genotype_projection(assessment, review_variants)
    return {
        "gene": assessment["gene"],
        "patient_causal_chain": {"genotype_status": status},
        "selected_genotype": genotype,
        "selected_variants": [variant_projection(item) for item in variants],
        "decisive_rank_evidence": decisive,
        "supporting_context_not_independent_votes": supporting,
        "decisive_counterevidence": counter,
        "expert_variant_assertions": expert,
        "external_variant_assertions": external,
        "unresolved": list(assessment.get("unresolved", ()))[:4],
        "omitted_unresolved_count": max(0, len(assessment.get("unresolved", ())) - 4),
        "omitted_variant_count": max(0, len(assessment["variants"]) - len(review_variants)),
        "observed_unresolved_variant_evidence": [
            {
                "variant_ref": item["variant_ref"],
                "identity": _ranking_value(item["identity"]),
                "classification": item.get("classification"),
                "patient_causality_not_established": True,
                "criteria": [
                    {
                        "code": criterion["code"],
                        "strength": criterion["strength"],
                        "status": criterion["status"],
                        "counted": criterion["counted"],
                        "source_refs": list(criterion["source_refs"])[:3],
                        "omitted_source_ref_count": max(0, len(criterion["source_refs"]) - 3),
                        "source_excerpt": str(next(iter(criterion.get("excerpts", ())), ""))[:500],
                    }
                    for criterion in item["criteria"]
                ],
            }
            for item in review_variants
        ]
        if unresolved_selection
        else [],
    }
