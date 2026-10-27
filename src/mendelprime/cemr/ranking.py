"""Evidence Brief, Global, Pairwise, and Reflection record contracts."""

import hashlib
import json
import re
from collections import Counter

from mendelprime import prompting
from mendelprime.domain import MendelPrimeError, Ranking

from . import hypotheses, views

COMPARISONS = {"already_discriminating", "partially_discriminating", "non_discriminating"}
QUALITIES = {"strong", "partial", "weak"}


def fail(code):
    raise MendelPrimeError(code, category="response")


def evidence_pool(evidence, gene):
    return evidence.get(
        gene,
        {
            "gene": gene,
            "documents": (),
            "omim_context": (),
            "retrieval_errors": (),
        },
    )


def brief_get(briefs, gene):
    return next((item for item in briefs["results"] if item["gene"] == gene), None)


def _chunks(value, limit=420):
    output = []
    for sentence in re.split(r"(?<=[.!?])\s+(?=[A-Z0-9(\[])", " ".join(value.split())):
        remaining = sentence.strip()
        while len(remaining) > limit:
            split = remaining.rfind(" ", 0, limit + 1)
            split = split if split >= limit // 2 else limit
            output.append(remaining[:split].strip())
            remaining = remaining[split:].strip()
        if remaining:
            output.append(remaining)
    return tuple(output)


def _omim(document):
    segments = [{"segment_id": f"{document['prompt_id']}#T0", "text": document["title"]}]
    segments.extend(
        {"segment_id": f"{document['prompt_id']}#S{index}", "text": content}
        for index, content in enumerate(_chunks(document["prompt_text"]), 1)
    )
    return {
        "evidence_id": document["prompt_id"],
        "title": document["title"],
        "report_units": [
            {
                "report_unit_id": f"{document['prompt_id']}#RU_CURATED",
                "report_unit_kind": "curated_source",
                "segments": segments,
            }
        ],
        "section": document.get("section") or "molecular_genetics",
        "pmids": list(document.get("cited_pmids", ())),
    }


def brief_prompt(profile, gene, pool, repair_error=""):
    if pool["gene"] != gene:
        raise MendelPrimeError("gene_mismatch", field="pool.gene")
    high = set(views.high_information_hpo_ids(profile))
    terms = sorted(
        views.present_terms(profile), key=lambda term: (-term["information_content"], term["hpo_id"])
    )
    value = prompting.template("brief_prompt", "brief")
    phenotype = value["phenotype_context"]
    phenotype["high_information_phenotype_ids"] = list(views.high_information_hpo_ids(profile))
    phenotype["all_present_phenotypes"] = [
        {
            "hpo_id": term["hpo_id"],
            "label": term["label"],
            "definition": term["definition"] or None,
            "phenotype_specificity_ic": round(term["information_content"], 6),
            "information_role": "high_information_clue" if term["hpo_id"] in high else "additional_context",
        }
        for term in terms
    ]
    views.add_clinical_context(phenotype, profile)
    if not profile["clinical_description"].strip():
        phenotype.pop("clinical_description", None)
    value["candidate_gene"] = gene
    value["evidence_contract"]["allowed_article_evidence_ids"] = [
        item["prompt_id"] for item in pool["documents"]
    ]
    value["organized_evidence"]["article_sources"] = [
        {
            "evidence_id": item["prompt_id"],
            "title": item["title"],
            "publication_types": list(item.get("publication_types", ())),
            "abstract_text": item["prompt_text"],
        }
        for item in pool["documents"]
    ]
    value["organized_evidence"]["omim_gene_context"] = [_omim(item) for item in pool["omim_context"]]
    value["required_json_schema"]["gene"] = gene
    if repair_error:
        value["repair_instruction"] = {
            "previous_error": repair_error,
            "instruction": (
                "Return a complete replacement JSON object containing only the gene and "
                "valid PMID evidence summaries."
            ),
        }
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def with_similar(prompt, evidence):
    if not evidence:
        return prompt
    value = json.loads(prompt)
    value["similar_case_evidence"] = [views.similar_payload(item) for item in evidence]
    value["similar_case_instructions"] = [
        "This is retrieved patient-level evidence, not a publication.",
        "known_causal_gene is the reported causal gene for the similar case.",
        "Use the complete case HPO terms and explicit HPO match fields as context.",
        "Do not convert a similar case into a PMID or claim it proves the current patient diagnosis.",
    ]
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def parse_brief(value, gene, allowed, limit):
    if str(value.get("gene") or "").strip().upper() != gene:
        fail("gene_contract_mismatch")
    raw = value.get("selected_evidence")
    if not isinstance(raw, list):
        fail("selected_evidence_must_be_list")
    selected, seen = [], set()
    for item in raw:
        if not isinstance(item, dict):
            continue
        prompt_id = str(item.get("evidence_id") or "").strip()
        summary = " ".join(str(item.get("evidence_summary") or "").split())[:800]
        if prompt_id in allowed and prompt_id not in seen and summary:
            selected.append({"document": allowed[prompt_id], "summary": summary})
            seen.add(prompt_id)
        if len(selected) == limit:
            break
    return tuple(selected)


def source_rankings(hypothesis):
    return {method: item for method, item in hypothesis["sources"].items() if item.genes}


def brief_payload(value, limit=None):
    summary = value["summary"] if limit is None else value["summary"][:limit].rstrip()
    return {"evidence_id": value["document"]["prompt_id"], "evidence_summary": summary}


def evidence_card(gene, briefs, limit=None):
    brief = brief_get(briefs, gene)
    value = {
        "selected_evidence": [brief_payload(item, limit) for item in brief["selected_evidence"]]
        if brief
        else []
    }
    if brief and brief["similar_case_evidence"]:
        value["similar_case_evidence"] = [
            views.similar_payload(item) for item in brief["similar_case_evidence"]
        ]
    return value


def _candidate_card(gene, position, admission, hypothesis, briefs):
    fusion = admission.details.get("fusion_genes", ())
    return {
        "input_rank": position,
        "fusion_rank": fusion.index(gene) + 1 if gene in fusion else None,
        "gene": gene,
        "source_llm_ranks": {
            f"{source}_llm": item.position(gene) for source, item in source_rankings(hypothesis).items()
        },
        **evidence_card(gene, briefs, 800),
    }


def global_prompt(profile, admission, hypothesis, briefs, repair_error=""):
    sources = source_rankings(hypothesis)
    patient = {
        "present_phenotypes": [views.term_payload(term, True) for term in views.present_terms(profile)],
        "ic_prioritized_hpo_ids": list(views.high_information_hpo_ids(profile)),
        "information_content_note": "Information content indicates phenotype specificity, not diagnostic importance.",
        "negated_phenotypes": [term["hpo_id"] for term in views.absent_terms(profile)],
    }
    views.add_clinical_context(patient, profile)
    value = prompting.template("global_prompt", "global")
    value["patient"] = patient
    value["source_llm_candidate_counts"] = {
        f"{source}_llm": len(item.genes) for source, item in sources.items()
    }
    value["baseline_candidate_order"] = list(admission.genes)
    value["candidates_in_baseline_order"] = [
        _candidate_card(gene, position, admission, hypothesis, briefs)
        for position, gene in enumerate(admission.genes, 1)
    ]
    value["candidate_contract"].update(
        allowed_genes=list(admission.genes), candidate_count=len(admission.genes)
    )
    value["required_json_schema"]["ranked_candidate_genes"] = list(admission.genes)
    if repair_error:
        instruction = (
            "Return every allowed gene exactly once. The top decision and every comparison "
            "must use allowed genes and cite only supplied PMID: or SC: evidence for the corresponding gene."
        )
        value["repair_instruction"] = {"previous_error": repair_error, "instruction": instruction}
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def _refs(value, first, second, allowed, prefixes):
    left = str(value.get(f"{prefixes[0]}_primary_source_ref") or "")
    right = str(value.get(f"{prefixes[1]}_primary_source_ref") or "")
    return (not left or left in allowed[first]) and (not right or right in allowed[second])


def parse_global(value, genes, briefs):
    raw = value.get("ranked_candidate_genes", ())
    ranked = tuple(str(gene).strip().upper() for gene in raw) if isinstance(raw, (list, tuple)) else ()
    if not ranked or len(ranked) != len(genes) or len(set(ranked)) != len(genes) or set(ranked) != set(genes):
        fail("ranked_candidate_genes_contract_mismatch")
    allowed = {}
    for gene in genes:
        brief = brief_get(briefs, gene)
        allowed[gene] = (
            {item["document"]["prompt_id"] for item in brief["selected_evidence"]} if brief else set()
        )
        if brief:
            allowed[gene].update(item["evidence_id"] for item in brief["similar_case_evidence"])
    top = value.get("top_decision")
    if not isinstance(top, dict) or str(top.get("selected_gene") or "").upper() != ranked[0]:
        fail("top_decision_invalid")
    alternative = str(top.get("strongest_alternative_gene") or "").upper()
    if alternative not in allowed or alternative == ranked[0]:
        fail("top_decision_alternative_invalid")
    if not _refs(top, ranked[0], alternative, allowed, ("selected", "alternative")):
        fail("top_decision_source_ref_invalid")
    comparisons = value.get("key_comparisons") or []
    if not isinstance(comparisons, list):
        fail("key_comparisons_invalid")
    positions = {gene: rank for rank, gene in enumerate(ranked)}
    for item in comparisons:
        if not isinstance(item, dict):
            fail("key_comparison_invalid")
        higher, lower = str(item.get("higher_gene") or "").upper(), str(item.get("lower_gene") or "").upper()
        if higher not in positions or lower not in positions or positions[higher] >= positions[lower]:
            fail("key_comparison_order_invalid")
        if not _refs(item, higher, lower, allowed, ("higher", "lower")):
            fail("key_comparison_source_ref_invalid")
    groups, used = [], set()
    raw_groups = value.get("indistinguishable_groups") or []
    if not isinstance(raw_groups, list):
        fail("indistinguishable_group_invalid")
    for raw_group in raw_groups:
        group = tuple(str(gene).upper() for gene in raw_group) if isinstance(raw_group, list) else ()
        if len(group) < 2 or set(group) - set(genes) or used.intersection(group):
            fail("indistinguishable_group_invalid")
        groups.append(group)
        used.update(group)
    summary = " ".join(str(value.get("case_pattern_summary") or "").split())
    if not summary:
        fail("case_pattern_summary_missing")
    return Ranking(
        genes=ranked,
        stage="global",
        details={
            "assessment": {
                "top_decision": top,
                "key_comparisons": tuple(comparisons),
                "indistinguishable_groups": tuple(groups),
                "case_pattern_summary": summary,
                "uncertainty": " ".join(str(value.get("uncertainty") or "").split()),
            }
        },
    )


def pair_prompt(profile, first, second, admission, hypothesis, briefs, include_clinical=False):
    sources = source_rankings(hypothesis)

    def card(gene):
        admitted = gene in admission.genes
        return {
            "gene": gene,
            "source_specific_ranks": {
                f"{source}_llm": {
                    "rank": item.position(gene) if admitted else None,
                    "candidate_count": len(item.genes) if admitted else 0,
                }
                for source, item in sources.items()
            },
            **evidence_card(gene, briefs),
        }

    patient = {
        "present_phenotypes": [views.term_payload(term) for term in views.present_terms(profile)],
        "high_information_hpo_ids": list(views.high_information_hpo_ids(profile)),
        "negated_phenotypes": [term["hpo_id"] for term in views.absent_terms(profile)],
    }
    if include_clinical:
        views.add_clinical_context(patient, profile)
    return json.dumps(
        {
            "task": "Which candidate gene better explains this patient?",
            "patient": patient,
            "candidate_a": card(first),
            "candidate_b": card(second),
            "required_output": "A or B",
        },
        ensure_ascii=False,
        separators=(",", ":"),
    )


# VCF views consume the canonical prepared facts without altering HPO projections.
def vcf_dump(value):
    return json.dumps(views._ranking_value(value), ensure_ascii=False, separators=(",", ":"))


def _bounded_vcf_dump(value, limit, preserve_strings=()):
    """Keep established Prompt bytes unless the complete request exceeds its frozen limit."""
    prompt = vcf_dump(value)
    if len(prompt) <= limit:
        return prompt
    packed = views.compact_prompt_payload(value, preserve_strings)
    compact = vcf_dump(packed)
    if len(compact) > limit:
        compact = min(compact, vcf_dump(views.compact_prompt_columns(packed)), key=len)
    return compact if len(compact) < len(prompt) else prompt


def vcf_patient(profile):
    high = set(views.high_information_hpo_ids(profile))
    terms = views.present_terms(profile)
    value = {
        "present_phenotypes": [
            {
                "hpo_id": term["hpo_id"],
                "label": term["label"],
                "definition": term["definition"] or None,
                "information_role": "high_information_clue"
                if term["hpo_id"] in high
                else "additional_context",
            }
            for term in terms
        ],
        "high_information_hpo_ids": list(views.high_information_hpo_ids(profile)),
        "information_content_note": "High-information marks ontology specificity, not diagnostic importance.",
        "phenotype_input_quality": {
            "present_phenotype_count": len(terms),
            "sparse_input": len(terms) < 3,
            "absence_of_unreported_features_is_not_negative": True,
        },
        "negated_phenotypes": [
            {"hpo_id": term["hpo_id"], "label": term["label"]} for term in views.absent_terms(profile)
        ],
    }
    if profile["clinical_description"].strip():
        views.add_clinical_context(value, profile)
    return value


def vcf_source_prompt(profile, tool, observations):
    genes = tool["ranking"].genes[:100]
    return _bounded_vcf_dump(
        {
            "task": f"Rerank all {len(genes)} {tool['tool_id']} genes from most to least likely.",
            "patient": vcf_patient(profile),
            "source": tool["tool_id"],
            "baseline_order": genes,
            "candidate_cards_alphabetical": views.compact_source_cards(tool, observations),
            "candidate_contract": {
                "allowed_genes": genes,
                "candidate_count": len(genes),
                "return_ranked_subset_without_duplicates": True,
            },
            "required_json_schema": {"ranked_genes": genes},
        },
        160000,
    )


def _vcf_science(profile, genes, prepared, acmg, evidence):
    if prepared is None or acmg is None:
        raise MendelPrimeError("prepared_ranking_evidence_required", field="ranking_context")
    return views.ranking_context_views(
        profile["case_key"], genes, {"prepared": prepared, "acmg": acmg, "evidence": evidence or {}}
    )


def vcf_repair(prompt, error):
    value = json.loads(prompt)
    value["repair_instruction"] = {
        "previous_error": error,
        "instruction": "Return a complete replacement satisfying the supplied gene, source and prefix contracts.",
    }
    return vcf_dump(value)


def vcf_parse_tail(value, allowed, required):
    raw = value.get("ranked_tail_genes")
    if not isinstance(raw, list):
        fail("ranked_tail_genes_must_be_list")
    genes = tuple(dict.fromkeys(str(gene).strip().upper() for gene in raw if str(gene).strip()))
    if len(genes) != required:
        fail(f"expected_{required}_genes_got_{len(genes)}")
    if set(genes) - set(allowed):
        fail("output_contains_disallowed_gene")
    return genes


def vcf_selected_prepared(context, *, limit=4, causal=False, include_observation=True):
    assessment = context["ranking_acmg"] or {}
    configuration = assessment.get("genotype_configuration", {})
    role = configuration.get("selected_variant_refs_role")
    selected = tuple(
        configuration.get("selected_variant_refs") or assessment.get("selected_variant_refs", ())
    )
    if causal and role == "observed_unresolved_candidates":
        return []
    by_ref = {item["variant_ref"]: item for item in context["prepared_variants"] if item.get("variant_ref")}
    ordered = [by_ref[ref] for ref in selected if ref in by_ref]
    if not causal:
        ordered.extend(item for item in context["prepared_variants"] if item not in ordered)
    return [
        views.prepared_variant_projection(row, include_observation=include_observation)
        for row in (ordered if limit is None else ordered[:limit])
    ]


def _vcf_fusion_variant(card, context):
    value = dict(card)
    assessment = context["ranking_acmg"] or {}
    selected = set(assessment.get("selected_variant_refs", ()))
    keys = {
        ":".join(map(str, item["identity"]["variant_key"]))
        for item in assessment.get("variants", ())
        if item["variant_ref"] in selected
    }
    alleles = [
        item
        for item in card["alleles"]
        if not keys
        or item.get("vcf_allele_key") in keys
        or item["allele_identity"] in keys
        or any(annotation.get("variant_key") in keys for annotation in item["tool_annotations"])
    ]
    value["alleles"] = [
        {
            key: raw
            for key, raw in row.items()
            if key != "family_observation" and not (key == "vcf_allele_key" and raw == row["allele_identity"])
        }
        for row in alleles[:2]
    ]
    value["omitted_allele_count"] = max(0, len(card["alleles"]) - len(value["alleles"]))
    value["native_tool_records"] = [
        {**item, "records": item["records"][:2], "omitted_record_count": max(0, len(item["records"]) - 2)}
        for item in card["native_tool_records"]
    ]
    shown_keys = {item["allele_identity"] for item in value["alleles"]}
    shown_keys.update(item.get("vcf_allele_key") for item in value["alleles"] if item.get("vcf_allele_key"))
    displayed = [
        item
        for item in assessment.get("variants", ())
        if ":".join(map(str, item["identity"]["variant_key"])) in shown_keys
    ]
    value["selected_variant_refs"] = (
        views.genotype_projection(assessment, displayed).get("selected_variant_refs", ())
        if assessment
        else ()
    )
    value["unparsed_variant_source_ids"] = card.get("unparsed_variant_source_ids", [])[:4]
    value["omitted_unparsed_variant_count"] = max(0, len(card.get("unparsed_variant_source_ids", ())) - 4)
    return value


def _vcf_fusion_acmg(context, family_columns):
    assessment = context["ranking_acmg"]
    if assessment is None:
        return None
    selected = tuple(
        assessment.get("genotype_configuration", {}).get("selected_variant_refs")
        or assessment.get("selected_variant_refs", ())
    )
    by_ref = {item["variant_ref"]: item for item in assessment["variants"]}
    variants, displayed = [], []
    selected_rows = [by_ref[ref] for ref in selected if ref in by_ref]
    for item in selected_rows or assessment["variants"]:
        projected = views._present(views.variant_projection(item))
        # Identity is already shown in structured_variant_evidence under the same variant_ref.
        projected.pop("identity", None)
        if projected.get("independent_classification") == projected.get("classification"):
            projected.pop("independent_classification", None)
        variants.append(
            {
                **projected,
                "criteria": [
                    {
                        "code": criterion["code"],
                        "strength": criterion["strength"],
                        "direction": criterion["direction"],
                        "status": criterion["status"],
                        "counted": criterion["counted"],
                        "source_refs": list(criterion["source_refs"])[:3],
                        "source_excerpt": str(next(iter(criterion.get("excerpts", ())), ""))[:300],
                        "method": criterion.get("method", ""),
                        "limitations": list(criterion.get("limitations", ()))[:2],
                        "independence_group": criterion.get("independence_group", ""),
                    }
                    for criterion in item["criteria"]
                ],
            }
        )
        displayed.append(item)
        if len(variants) == 2:
            break
    genotype = views.genotype_projection(assessment, displayed)
    family = genotype["family_observations"]
    if family:
        genotype["family_observations"] = {
            ref: {role: [call.get(key) for key in family_columns] for role, call in calls.items()}
            for ref, calls in family.items()
        }
    return {
        "gene": assessment["gene"],
        "phase": assessment.get("phase"),
        "selected_variant_refs": genotype.get("selected_variant_refs", []),
        "genotype_configuration": genotype,
        "variants": variants,
        "unresolved": assessment.get("unresolved", ())[:4],
        "omitted_unresolved_count": max(0, len(assessment.get("unresolved", ())) - 4),
        "omitted_variant_count": max(0, len(assessment["variants"]) - len(displayed)),
        "computation_ineligible_variant_count": assessment.get("computation_ineligible_variant_count", 0),
    }


def _vcf_stage_context(context, displayed, *, limit):
    """Project auxiliary facts at the same allele scope as the stage's variant cards."""
    keys = {tuple(item["variant_key"]) for item in displayed}
    clinical = context.get("clinical_context", {})
    omitted = {"computation_ineligible_variants": context.get("computation_ineligible_variant_count", 0)}

    def scoped(rows, name):
        relevant = []
        for row in rows:
            key = row.get("variant_key") or row.get("identity", {}).get("variant_key")
            if key is None or not keys or tuple(key) in keys:
                relevant.append(row)
        result = relevant[:limit]
        omitted[name] = len(rows) - len(result)
        return result

    clinical_view = {
        **clinical,
        "variant_assertions": scoped(clinical.get("variant_assertions", ()), "clinical_variant_assertions"),
        "lookups": scoped(clinical.get("lookups", ()), "clinical_lookups"),
        "gene_context": list(clinical.get("gene_context", ()))[:2],
    }

    def document_view(row):
        return row if len(row.get("text", "")) <= 400 else views.clinical_document_projection(row)

    clinical_view["gene_context"] = [document_view(row) for row in clinical_view["gene_context"]]
    clinical_view["variant_assertions"] = [
        {
            **{key: raw for key, raw in row.items() if key != "document"},
            "document": document_view(row["document"]),
        }
        for row in clinical_view["variant_assertions"]
    ]
    omitted["clinical_gene_context"] = max(0, len(clinical.get("gene_context", ())) - 2)
    source_ids = views.projection_source_ids(
        [row for row in context["prepared_variants"] if tuple(row["variant_key"]) in keys]
    )
    structured = {}
    for name, rows in context.get("structured_gene_evidence", {}).items():
        preferred = [row for row in rows if row.get("source_id") in source_ids]
        preferred.extend(row for row in rows if row.get("source_id") not in source_ids)
        structured[name] = [
            {
                **row,
                "facts": views.clinical_document_projection(
                    {"text": json.dumps(row["facts"]), "source_id": row["source_id"]}
                ).get("facts", {}),
            }
            if name == "gene_disease_assertions"
            else row
            for row in preferred[:2]
        ]
        omitted[name] = max(0, len(rows) - 2)
    unusable = scoped(context.get("unusable_observations", ()), "unusable_observations")
    return {
        "clinical_database_context": clinical_view,
        "structured_gene_evidence": structured,
        "unusable_observations": unusable,
        "omitted_context_counts": {name: count for name, count in omitted.items() if count},
    }


def vcf_fusion_prompt(profile, hypothesis, native, pool, prepared, acmg, evidence):
    genes = pool["genes"]
    science = _vcf_science(profile, genes, prepared, acmg, evidence)
    native_cards = {
        item["gene"]: item for item in views.candidate_cards(native, genes, include_unprioritized=True)
    }
    displayed = {
        gene: vcf_selected_prepared(science[gene], limit=2, include_observation=False) for gene in genes
    }
    auxiliary = {gene: _vcf_stage_context(science[gene], displayed[gene], limit=2) for gene in genes}
    family_columns = (
        "GT",
        "AD",
        "DP",
        "GQ",
        "PL",
        "call_status",
        "quality_status",
        "usable_for_inheritance",
        "state",
        "carries_alternate",
        "ploidy",
    )
    return _bounded_vcf_dump(
        {
            "task": "Fuse only the supplied genotype ToolLLM Top10 rankings and return an exact full ranking of their deduplicated union.",
            "patient": vcf_patient(profile),
            "source_contract": {
                "allowed_sources": list(pool["seed_genes"]),
                "source_rankings": pool["seed_genes"],
                "source_kinds": pool["source_kinds"],
                "phenotype_only_sources_present": False,
                "ranks_are_correlated_priors": True,
            },
            "structured_evidence_contract": {
                "clingen_vcep_and_matching_clinvar_are_one_assertion": True,
                "gene_disease_validity_is_not_variant_pathogenicity": True,
                "patient_variant_observation_is_not_pathogenicity_evidence": True,
                "missing_evidence_is_uncertainty_not_benign_evidence": True,
                "acmg_variant_identity_lookup": "structured_variant_evidence.variants by variant_ref",
                "family_observation_columns": family_columns,
            },
            "deterministic_rrf_baseline": genes,
            "candidate_cards_alphabetical": [
                {
                    "gene": gene,
                    **pool["source_rank_context"][gene],
                    "variant_context": _vcf_fusion_variant(native_cards[gene], science[gene]),
                    "structured_variant_evidence": {
                        "variants": displayed[gene],
                        **auxiliary[gene]["structured_gene_evidence"],
                    },
                    "clinical_database_context": auxiliary[gene]["clinical_database_context"],
                    "unusable_observations": auxiliary[gene]["unusable_observations"],
                    "omitted_context_counts": auxiliary[gene]["omitted_context_counts"],
                    "mgagent_preliminary_acmg": _vcf_fusion_acmg(science[gene], family_columns),
                }
                for gene in sorted(genes)
            ],
            "candidate_contract": {
                "allowed_genes": genes,
                "candidate_count": len(genes),
                "return_every_gene_exactly_once": True,
            },
            "required_json_schema": {"ranked_genes": genes},
        },
        220000,
    )


def vcf_anchor_document(pool):
    documents = [
        entry["document"]
        for entry in vcf_article_catalog(pool).values()
        if entry["document"]["prompt_text"].strip()
    ]

    def priority(item):
        provenance = item.get("provenance", ())
        return next(
            (
                index
                for index, name in enumerate(("variant", "must_review", "amelie", "omim"))
                if name in provenance
            ),
            4,
        ), item["prompt_id"]

    return min(documents, key=priority) if documents else None


def _vcf_recall_view(context):
    recall = views.layered_evidence(context["ranking_acmg"], recall=True) or {
        "gene": context["gene"],
        "purpose": "candidate_recall_not_final_causal_discrimination",
        "genotype_status": "not_assessed",
        "selected_variant_refs": [],
        "missing_information_is_not_negative": True,
    }
    auxiliary = _vcf_stage_context(context, vcf_selected_prepared(context), limit=4)
    return {
        **recall,
        "clinical_database_context": auxiliary["clinical_database_context"],
        "unusable_observations": auxiliary["unusable_observations"],
        "omitted_context_counts": auxiliary["omitted_context_counts"],
    }


def vcf_admission_prompt(profile, fusion, native, pool, count, prepared, acmg, evidence):
    genes = fusion.genes
    science = _vcf_science(profile, genes, prepared, acmg, evidence)
    return vcf_dump(
        {
            "task": f"Rank the complete VCF candidate pool and return an exact full ranking; the first {count} genes are admitted to full evidence review.",
            "patient": vcf_patient(profile),
            "rank_contract": {
                "source_kinds": pool["source_kinds"],
                "ranks_are_correlated_priors": True,
                "missing_rank_is_not_negative_evidence": True,
            },
            "vcf_genotype_fusion_baseline": genes,
            "candidate_cards_alphabetical": [
                {
                    "gene": gene,
                    "rank_context": pool["source_rank_context"][gene],
                    "recall_evidence": _vcf_recall_view(science[gene]),
                    "anchor_article": views.document_view(document)
                    if (document := vcf_anchor_document((evidence or {}).get(gene, {})))
                    else None,
                }
                for gene in sorted(genes)
            ],
            "candidate_contract": {
                "allowed_genes": genes,
                "candidate_count": len(genes),
                "admission_count": count,
                "admitted_genes_must_equal_full_ranking_prefix": True,
            },
            "required_json_schema": {
                "ranked_candidate_genes": genes,
                "top15_genes": genes[:count],
                "uncertainty": "main limitation or empty string",
            },
        }
    )


def vcf_source_anchors(document):
    text = document["prompt_text"]
    start, spans = 0, []
    for boundary in re.finditer(r"(?<=[.!?])\s+(?=[A-Z0-9(\[])", text):
        spans.append((start, boundary.start()))
        start = boundary.end()
    spans.append((start, len(text)))
    bounded = []
    for start, end in spans:
        while start < end and text[start].isspace():
            start += 1
        while end > start and text[end - 1].isspace():
            end -= 1
        while end - start > 1200:
            split = text.rfind(" ", start + 600, start + 1201)
            split = split if split >= start + 600 else start + 1200
            bounded.append((start, split))
            start = split
            while start < end and text[start].isspace():
                start += 1
        if start < end:
            bounded.append((start, end))
    digest = hashlib.sha256(text.encode()).hexdigest()
    return tuple(
        {
            "anchor_id": f"{document['prompt_id']}#A{index}",
            "verbatim_text": text[start:end],
            "start_char": start,
            "end_char": end,
            "source_id": _vcf_source_fields(document)["source_id"],
            "text_sha256": digest,
        }
        for index, (start, end) in enumerate(bounded, 1)
    )


def _vcf_source_fields(document):
    if "source" in document:
        source = document["source"]
        return {
            "source_id": source.source_id,
            "kind": source.kind,
            "provider": source.provider,
            "variant_keys": [":".join(map(str, key)) for key in source.variant_keys],
            "parent_source_ids": list(source.parent_source_ids),
        }
    return {
        key: document.get(key, [] if key in {"variant_keys", "parent_source_ids"} else "")
        for key in ("source_id", "kind", "provider", "variant_keys", "parent_source_ids")
    }


def vcf_article_catalog(pool):
    documents, lineage = {}, {}
    for document in pool.get("documents", ()):
        source = _vcf_source_fields(document)
        if document.get("section") == "article_lineage" or (
            source["kind"] == "derived" and source["provider"] == "pmc"
        ):
            if source["kind"] != "derived":
                fail("invalid_article_lineage_source")
            try:
                record = json.loads(document["prompt_text"])
            except ValueError:
                fail("invalid_article_lineage_record")
            if record.get("record_type") != "pmc_article_lineage_v1" or set(source["parent_source_ids"]) != {
                record.get("abstract_source_id"),
                record.get("passage_source_id"),
            }:
                fail("invalid_article_lineage_parents")
            lineage[record["passage_source_id"]] = record
            continue
        if source["kind"] != "literature":
            fail("vcf_literature_document_required")
        identifier = document["prompt_id"]
        if not re.fullmatch(r"PMID:\d+|PMC:\d+#[^\s]+", identifier):
            fail("vcf_literature_identifier_required")
        if identifier in documents and documents[identifier]["source_id"] != source["source_id"]:
            fail("ambiguous_document_identifier")
        documents[identifier] = {**source, "document": document}
    for identifier, entry in documents.items():
        metadata = entry["document"].get("metadata", {})
        link = lineage.get(entry["source_id"], {})
        if (
            link.get("pmid")
            and metadata.get("parent_pmid")
            and str(link["pmid"]).removeprefix("PMID:") != str(metadata["parent_pmid"]).removeprefix("PMID:")
        ):
            fail("pmc_parent_article_mismatch")
        parent = str(
            link.get("pmid")
            or metadata.get("parent_pmid")
            or (
                entry["document"].get("article_id", "")
                if str(entry["document"].get("article_id", "")).startswith("PMID:")
                else ""
            )
        )
        parent = "PMID:" + parent.removeprefix("PMID:") if parent else ""
        if parent and not re.fullmatch(r"PMID:\d+", parent):
            fail("invalid_parent_article_identifier")
        if identifier.startswith("PMC:"):
            pmcid = "PMC" + identifier.split("#", 1)[0].split(":", 1)[1]
            if metadata.get("pmcid", pmcid) != pmcid or link.get("pmcid", pmcid) != pmcid:
                fail("pmc_article_identity_mismatch")
            entry.update(
                article_id=parent or identifier.split("#", 1)[0],
                parent_pmid=parent or None,
                pmcid=pmcid,
                evidence_type="pmc_fulltext_passage",
                content_type=metadata.get("content_type", "fulltext_passage"),
            )
        else:
            entry.update(
                article_id=identifier,
                parent_pmid=identifier,
                pmcid=None,
                evidence_type="pubmed_abstract" if entry["provider"] == "pubmed" else "literature_excerpt",
                content_type=metadata.get("content_type", "abstract"),
            )
        entry["license"] = metadata.get("license")
        entry["dependency_group"] = entry["article_id"]
    return documents


def _vcf_omim(document):
    metadata = document.get("metadata", {})
    return {
        **views.document_view(document),
        **{
            key: views._ranking_value(metadata.get(key, ()))
            for key in ("variant_number", "linked_patient_variant_refs", "matched_identifiers")
        },
        "pmids": list(document.get("cited_pmids", metadata.get("pmids", ()))),
        "evidence_type": "omim_context_not_article",
    }


def vcf_brief_prompt(profile, gene, pool, context, error=""):
    catalog = vcf_article_catalog(pool)
    sources, allowed = [], {}
    groups = {}
    for identifier, entry in catalog.items():
        document = entry["document"]
        anchors = vcf_source_anchors(document)
        allowed[document["prompt_id"]] = [item["anchor_id"] for item in anchors]
        marked, cursor = [], 0
        for anchor in anchors:
            marked.extend(
                (
                    document["prompt_text"][cursor : anchor["start_char"]],
                    f"[{anchor['anchor_id']}] ",
                    anchor["verbatim_text"],
                )
            )
            cursor = anchor["end_char"]
        marked.append(document["prompt_text"][cursor:])
        sources.append(
            {
                "evidence_id": document["prompt_id"],
                "title": document["title"],
                "source_text_with_inline_source_anchors": "".join(marked),
                **{key: value for key, value in entry.items() if key != "document"},
                "publication_types": list(document.get("publication_types", ())),
                "provenance": list(document.get("provenance", ())),
                "linked_patient_variant_refs": list(
                    document.get("metadata", {}).get("linked_patient_variant_refs", ())
                ),
            }
        )
        groups.setdefault(entry["article_id"], []).append(identifier)
    displayed = vcf_selected_prepared(context)
    auxiliary = _vcf_stage_context(context, displayed, limit=4)
    assessment = context["ranking_acmg"]
    displayed_refs = {item["variant_ref"] for item in displayed}
    genotype = (
        views.genotype_projection(
            assessment, [item for item in assessment["variants"] if item["variant_ref"] in displayed_refs]
        )
        if assessment is not None
        else {}
    )
    value = {
        "task": "Select zero to three distinct parent articles relevant to this gene and patient. PubMed abstracts and PMC passages from the same article form one dependent source, not independent publications. Summarize only supplied source text and return its exact document and anchor IDs.",
        "patient": vcf_patient(profile),
        "candidate_gene": gene,
        "variant_context": {
            "variants": displayed,
            "selected_genotype": genotype,
            "missing_information_is_uncertainty": True,
            "clinical_database_context": auxiliary["clinical_database_context"],
            "unusable_observations": auxiliary["unusable_observations"],
            "omitted_context_counts": auxiliary["omitted_context_counts"],
            "omitted_variant_count": max(0, len(context["prepared_variants"]) - len(displayed)),
        },
        "evidence_contract": {
            "allowed_article_evidence_ids": list(allowed),
            "allowed_source_anchor_ids_by_evidence_id": allowed,
            "article_groups": groups,
            "maximum_distinct_articles": 3,
            "maximum_anchors_per_article": 4,
            "each_item_anchors_belong_to_its_exact_document_id": True,
            "lineage_records_are_not_medical_evidence": True,
        },
        "organized_evidence": {
            "article_sources": sources,
            "omim_molecular_genetics_context": [
                _vcf_omim(item)
                for item in pool["omim_context"]
                if item.get("section", "molecular_genetics") == "molecular_genetics"
            ],
            "omim_matched_allelic_variant_context": [
                _vcf_omim(item) for item in pool["omim_context"] if item.get("section") == "allelic_variant"
            ],
            "retrieval_errors": list(pool["retrieval_errors"]),
        },
        "required_json_schema": {
            "gene": gene,
            "selected_evidence": [
                {
                    "evidence_id": "PMID:123",
                    "evidence_summary": "source-grounded summary; at most 800 characters",
                    "source_anchor_ids": ["PMID:123#A1"],
                }
            ],
        },
        "selection_rules": [
            "Prefer exact patient-variant human reports, family segregation, relevant functional assays, and affected-human phenotype-genotype evidence.",
            "Do not select a source merely for gene mention, a review, general mechanism, or broad disease description.",
            "Do not merge patients or publications. Preserve variant identifiers and original source type.",
            "Choose no more than three distinct parent articles and four supplied anchors per article across all its documents.",
            "Use a separate selected_evidence item for each document. Its source_anchor_ids must all belong to its evidence_id; never attach PMC passage anchors to a PMID document, even for the same parent article.",
            "OMIM is contextual evidence, not an article selection. Article-lineage records are provenance only.",
            "Do not rank genes, make a diagnosis or classify a variant in this filtering step.",
        ],
    }
    prompt = vcf_dump(value)
    return vcf_repair(prompt, error) if error else prompt


def vcf_parse_brief(value, gene, allowed):
    if str(value.get("gene", "")).strip().upper() != gene or not isinstance(
        value.get("selected_evidence"), list
    ):
        fail("brief_contract_mismatch")
    selected, seen, articles, anchor_counts = [], set(), set(), Counter()
    catalog = vcf_article_catalog({"documents": tuple(allowed.values())})
    if len(value["selected_evidence"]) > 12:
        fail("too_many_selected_evidence")
    for item in value["selected_evidence"]:
        if not isinstance(item, dict):
            fail("selected_evidence_record_required")
        ref, ids = str(item.get("evidence_id") or "").strip(), item.get("source_anchor_ids")
        summary = " ".join(str(item.get("evidence_summary") or "").split())[:800]
        ids = (
            [str(identifier).strip() for identifier in ids if str(identifier).strip()]
            if isinstance(ids, list)
            else None
        )
        if (
            ref not in allowed
            or ref in seen
            or not summary
            or not isinstance(ids, list)
            or not ids
            or len(set(ids)) != len(ids)
        ):
            fail("brief_source_contract_mismatch")
        anchors = {anchor["anchor_id"]: anchor for anchor in vcf_source_anchors(allowed[ref])}
        if any(identifier not in anchors for identifier in ids):
            fail("brief_source_anchor_invalid")
        article = catalog[ref]["article_id"]
        ids = ids[:4]
        articles.add(article)
        anchor_counts[article] += len(ids)
        if len(articles) > 3 or anchor_counts[article] > 4:
            fail("brief_article_budget_exceeded")
        selected.append(
            {
                "document": allowed[ref],
                "summary": summary,
                "source_anchors": tuple(anchors[identifier] for identifier in ids),
                "article_id": article,
                "dependency_group": article,
            }
        )
        seen.add(ref)
    return tuple(selected)


def _vcf_brief_evidence(brief):
    articles = {}
    for item in brief["selected_evidence"] if brief else ():
        document = item["document"]
        article = (
            item.get("article_id")
            or document.get("article_id")
            or document.get("metadata", {}).get("parent_pmid")
            or document["prompt_id"].split("#", 1)[0]
        )
        row = articles.setdefault(
            article,
            {
                "article_id": article,
                "dependency_group": article,
                "evidence_id": document["prompt_id"],
                "evidence_summary": "",
                "source_excerpt": "",
                "documents": [],
            },
        )
        row["evidence_summary"] = (row["evidence_summary"] + " " + item["summary"]).strip()[:600]
        row["source_excerpt"] = (
            row["source_excerpt"]
            + " "
            + " ".join(anchor["verbatim_text"] for anchor in item.get("source_anchors", ()))
        ).strip()[:400]
        row["documents"].append({"evidence_id": document["prompt_id"], **_vcf_source_fields(document)})
    return list(articles.values())


def vcf_fallback_anchor(gene, pool, terms=(), anchor_id=""):
    catalog = vcf_article_catalog(pool)
    publication = catalog.get(anchor_id, {}).get("document") or vcf_anchor_document(pool)
    omim = next(iter(pool.get("omim_context", ())), None)
    if publication is None and omim is None:
        return {"selected_pmid": None, "publication": None, "omim_context": None, "status": "unavailable"}

    def compact(document, limit, kind):
        if document is None:
            return None
        mentioned = bool(document.get("metadata", {}).get("candidate_gene_mentioned")) or bool(
            re.search(rf"\b{re.escape(gene)}\b", document["title"] + " " + document["prompt_text"], re.I)
        )
        metadata = document.get("metadata", {})
        linked = (
            bool(metadata.get("linked_patient_variant_refs"))
            or any("variant" in route.lower() for route in document.get("provenance", ()))
            or bool(_vcf_source_fields(document)["variant_keys"])
        )
        sentences = re.split(r"(?<=[.!?])\s+", " ".join(document["prompt_text"].split()))
        needles = tuple(str(term).lower() for term in (gene, *terms) if len(str(term)) >= 3)
        ranked = sorted(
            enumerate(sentences),
            key=lambda item: (-sum(term in item[1].lower() for term in needles), item[0]),
        )
        relevant = " ".join(text for _, text in sorted(ranked[:3]))
        return {
            "evidence_id": document["prompt_id"],
            "source_id": document["source_id"] if "source_id" in document else document["source"].source_id,
            "title": document["title"],
            "text": relevant[:limit] if mentioned or linked or kind == "gene_mechanism_context" else "",
            "provenance": list(document.get("provenance", ())),
            "evidence_quality": "gene_mechanism_context"
            if kind == "gene_mechanism_context"
            else "variant_linked"
            if linked
            else "gene_specific"
            if mentioned
            else "broad_context",
            "variant_linked": linked,
            "candidate_gene_mentioned": mentioned,
            "article_id": next(
                (
                    entry["article_id"]
                    for entry in catalog.values()
                    if entry["document"]["prompt_id"] == document["prompt_id"]
                ),
                None,
            ),
            "metadata": metadata,
        }

    return {
        "selected_pmid": publication["prompt_id"] if publication else None,
        "publication": compact(
            publication,
            800,
            "variant_linked"
            if publication
            and (
                publication.get("metadata", {}).get("linked_patient_variant_refs")
                or any("variant" in route.lower() for route in publication.get("provenance", ()))
            )
            else "gene_specific",
        ),
        "omim_context": compact(omim, 700, "gene_mechanism_context"),
    }


def vcf_discrimination(context, brief):
    value = views.layered_evidence(context["ranking_acmg"])
    if value is None:
        value = {
            "gene": context["gene"],
            "patient_causal_chain": {"genotype_status": "not_assessed"},
            "selected_genotype": {},
            "selected_variants": [],
            "decisive_rank_evidence": [],
            "supporting_context_not_independent_votes": [],
            "decisive_counterevidence": [],
            "expert_variant_assertions": [],
            "unresolved": ["No eligible patient variant assessment"],
        }
    displayed = vcf_selected_prepared(context, limit=2, causal=True)
    if not displayed:
        displayed = vcf_selected_prepared(context)
    auxiliary = _vcf_stage_context(context, displayed, limit=len(displayed) or 4)
    value["clinical_database_context"] = auxiliary["clinical_database_context"]
    value["unusable_observations"] = auxiliary["unusable_observations"]
    value["omitted_context_counts"] = auxiliary["omitted_context_counts"]
    value["exact_patient_variant_literature"] = []
    value["structured_variant_facts"] = vcf_selected_prepared(
        context, limit=2, causal=True, include_observation=False
    )
    value["gene_disease_validity"] = auxiliary["structured_gene_evidence"].get("gene_disease_assertions", [])
    value["gene_constraint"] = auxiliary["structured_gene_evidence"].get("gene_constraint", [])
    value["ranking_features"] = auxiliary["structured_gene_evidence"].get("ranking_features", [])
    pool = context["ranking_literature"]
    raw_catalog = vcf_article_catalog(
        {"documents": tuple(_vcf_view_document(item) for item in pool["documents"])}
    )
    for variant in value["selected_variants"]:
        identity = variant["identity"]
        tokens = [
            token
            for raw in (identity.get("hgvs_c"), identity.get("hgvs_p"))
            if len(token := re.sub(r"[^a-z0-9.>*_]", "", str(raw or "").lower())) >= 6
        ]
        candidates = list(brief["selected_evidence"]) if brief else []
        pool = context["ranking_literature"]
        known = {item["document"]["prompt_id"] for item in candidates}
        for document in (_vcf_view_document(item) for item in pool["documents"]):
            if document["prompt_id"] not in known and document.get("kind") == "literature":
                candidates.append({"document": document, "source_anchors": vcf_source_anchors(document)})
        for item in candidates:
            document = item["document"]
            metadata = document.get("metadata", {})
            linked_refs = metadata.get("linked_patient_variant_refs", ())
            match_status = metadata.get("match_status", ())
            match_status = (match_status,) if isinstance(match_status, str) else match_status
            source_fields = _vcf_source_fields(document)
            key_text = ":".join(map(str, identity["variant_key"]))
            exact_bound = key_text in source_fields["variant_keys"] or (
                variant["variant_ref"] in linked_refs
                and any(status in {"coordinate_exact", "variation_id_exact"} for status in match_status)
            )
            if not exact_bound:
                continue
            anchors = item.get("source_anchors", ())
            matches = [
                anchor
                for anchor in anchors
                if any(
                    token in re.sub(r"[^a-z0-9.>*_]", "", anchor["verbatim_text"].lower()) for token in tokens
                )
            ]
            if matches:
                value["exact_patient_variant_literature"].append(
                    {
                        "variant_ref": variant["variant_ref"],
                        "evidence_id": item["document"]["prompt_id"],
                        "match_type": "patient_variant_identifier_in_source",
                        "verbatim_source_excerpts": [item["verbatim_text"][:800] for item in matches],
                        "source_id": source_fields["source_id"],
                        "article_id": item.get("article_id")
                        or raw_catalog.get(document["prompt_id"], {}).get("article_id")
                        or document["prompt_id"].split("#", 1)[0],
                        "dependency_group": item.get("article_id")
                        or raw_catalog.get(document["prompt_id"], {}).get("article_id")
                        or document["prompt_id"].split("#", 1)[0],
                        "affected_human_evidence_requires_source_confirmation": True,
                    }
                )
                break
    return value


def vcf_global_prompt(profile, admission, hypothesis, briefs, contexts, error=""):
    genes = admission.genes
    if len(genes) > 10:
        raise MendelPrimeError("vcf_shortlist_requires_top10", field="admission")
    cards = []
    for rank, gene in enumerate(genes, 1):
        brief = brief_get(briefs, gene)
        card = {
            "input_rank": rank,
            "gene": gene,
            "molecular_evidence": vcf_discrimination(contexts[gene], brief),
            "rank_context": {
                "genotype_source_ranks": {
                    f"{method}_llm" if result.status == "ok" else method: {
                        "rank": result.position(gene),
                        "candidate_count": len(result.genes),
                    }
                    for method, result in hypothesis["sources"].items()
                },
                "genotype_fusion": {
                    "rank": admission.details.get("fusion_genes", ()).index(gene) + 1
                    if gene in admission.details.get("fusion_genes", ())
                    else None,
                    "candidate_count": len(admission.details.get("fusion_genes", ())),
                },
            },
            "variant_context": vcf_selected_prepared(contexts[gene]),
            "selected_evidence": _vcf_brief_evidence(brief),
        }
        refs, variants = vcf_reference_catalog(card)
        card["primary_source_options"] = sorted(refs)
        card["primary_variant_options"] = sorted(variants)
        cards.append(card)
    value = {
        "task": "Rank the recalled candidates using only patient-specific causal evidence; the first five enter Pairwise refinement.",
        "patient": vcf_patient(profile),
        "rank_contract": {
            "ranks_are_correlated_priors": True,
            "missing_rank_is_not_negative_evidence": True,
            "molecular_evidence_is_compressed_not_exhaustive": True,
            "evidence_categories_are_not_independent_votes": True,
            "acmg_strengths_are_not_fixed_gene_ranking_weights": True,
        },
        "baseline_candidate_order": genes,
        "candidates_in_admission_order": cards,
        "candidate_contract": {
            "allowed_genes": genes,
            "candidate_count": len(genes),
            "pairwise_shortlist_count": min(5, len(genes)),
            "pairwise_shortlist_is_ranked_prefix": True,
        },
        "required_json_schema": {
            "ranked_candidate_genes": genes,
            "case_pattern_summary": "patient comparison",
            "top_decision": {
                "selected_gene": genes[0] if genes else "",
                "strongest_alternative_gene": genes[1] if len(genes) > 1 else "",
                "selected_primary_source_ref": "exact selected-gene primary_source_options entry or empty",
                "alternative_primary_source_ref": "exact alternative-gene primary_source_options entry or empty",
                "selected_primary_variant_ref": "supplied variant_ref or empty",
                "alternative_primary_variant_ref": "supplied variant_ref or empty",
                "comparative_reason": "bilateral evidence-grounded comparison",
            },
            "key_comparisons": [
                {
                    "higher_gene": genes[0] if genes else "",
                    "lower_gene": genes[1] if len(genes) > 1 else "",
                    "higher_primary_source_ref": "exact higher-gene primary_source_options entry or empty",
                    "lower_primary_source_ref": "exact lower-gene primary_source_options entry or empty",
                    "higher_primary_variant_ref": "supplied higher-gene variant ref or empty",
                    "lower_primary_variant_ref": "supplied lower-gene variant ref or empty",
                    "comparative_reason": "concise bilateral source-grounded comparison",
                }
            ],
            "indistinguishable_groups": [],
            "uncertainty": "main limitation or empty",
        },
    }
    output_refs = {
        ref
        for card in cards
        for name in ("primary_source_options", "primary_variant_options")
        for ref in card[name]
    }
    prompt = _bounded_vcf_dump(value, 160000, output_refs)
    return vcf_repair(prompt, error) if error else prompt


def vcf_reference_catalog(value):
    sources, variants = set(), set()

    def visit(item):
        if isinstance(item, dict):
            for key, child in item.items():
                if key in {"source_id", "source_ref", "evidence_id"} and isinstance(child, str) and child:
                    sources.add(child)
                elif key in {"source_ids", "source_refs"} and isinstance(child, (list, tuple)):
                    sources.update(ref for ref in child if isinstance(ref, str) and ref)
                elif key == "variant_ref" and isinstance(child, str) and child:
                    variants.add(child)
                elif isinstance(child, (dict, tuple, list)):
                    visit(child)
        elif isinstance(item, (tuple, list)):
            for child in item:
                visit(child)

    visit(value)
    return sources, variants


def vcf_parse_global(value, genes, catalogs):
    ranked, _ = hypotheses.decode_vcf_permutation(
        {"ranked_genes": value.get("ranked_candidate_genes")}, genes
    )
    top, comparisons = value.get("top_decision"), value.get("key_comparisons", [])
    if (
        not isinstance(top, dict)
        or str(top.get("selected_gene") or "").upper() != ranked[0]
        or not isinstance(comparisons, list)
    ):
        fail("top_decision_invalid")
    alternative = str(top.get("strongest_alternative_gene") or "").upper()
    if len(genes) > 1 and (alternative not in genes or alternative == ranked[0]):
        fail("top_decision_alternative_invalid")
    allowed = {gene: values[0] for gene, values in catalogs.items()}
    variant_refs = {gene: values[1] for gene, values in catalogs.items()}

    def validate_refs(item, pairs):
        for prefix, gene in pairs:
            source, variant = (
                item.get(prefix + "_primary_source_ref", ""),
                item.get(prefix + "_primary_variant_ref", ""),
            )
            if (
                source
                and (gene not in allowed or source not in allowed[gene])
                or variant
                and (gene not in variant_refs or variant not in variant_refs[gene])
            ):
                fail("global_source_or_variant_scope_mismatch")

    validate_refs(top, (("selected", ranked[0]), ("alternative", alternative)))
    positions = {gene: index for index, gene in enumerate(ranked)}
    for comparison in comparisons:
        if not isinstance(comparison, dict):
            fail("key_comparison_invalid")
        higher, lower = (
            str(comparison.get("higher_gene") or "").upper(),
            str(comparison.get("lower_gene") or "").upper(),
        )
        if higher not in positions or lower not in positions or positions[higher] >= positions[lower]:
            fail("key_comparison_order_invalid")
        validate_refs(comparison, (("higher", higher), ("lower", lower)))
    summary = str(value.get("case_pattern_summary") or "").strip()
    if not summary:
        fail("case_pattern_summary_missing")
    raw_groups = value.get("indistinguishable_groups") or []
    if not isinstance(raw_groups, list):
        fail("indistinguishable_group_invalid")
    groups, used = [], set()
    for raw_group in raw_groups:
        if (
            not isinstance(raw_group, list)
            or len(raw_group) < 2
            or len(set(raw_group)) != len(raw_group)
            or set(raw_group) - set(genes)
            or used.intersection(raw_group)
        ):
            fail("indistinguishable_group_invalid")
        groups.append(tuple(raw_group))
        used.update(raw_group)
    return Ranking(
        genes=ranked,
        stage="global",
        details={
            "assessment": {
                "top_decision": top,
                "key_comparisons": tuple(comparisons),
                "case_pattern_summary": summary,
                "indistinguishable_groups": tuple(groups),
                "uncertainty": str(value.get("uncertainty") or ""),
            }
        },
    )


def vcf_pair_cards(profile, genes, admission, briefs, contexts):
    """Prepare each candidate once; both comparison directions consume the same evidence."""
    cards = {}
    for gene in genes:
        brief = brief_get(briefs, gene)
        selected = _vcf_brief_evidence(brief)
        pool = contexts[gene]["ranking_literature"]
        assessment = contexts[gene]["ranking_acmg"] or {}
        terms = tuple(
            value
            for variant in assessment.get("variants", ())
            for value in (variant["identity"].get("hgvs_c"), variant["identity"].get("hgvs_p"))
            if value
        )
        terms += tuple(term["label"] for term in views.present_terms(profile))
        cards[gene] = {
            "gene": gene,
            "admission_rank_prior": admission.position(gene),
            "causal_discrimination_evidence": vcf_discrimination(contexts[gene], brief),
            "selected_patient_relevant_literature": selected,
            "fallback_literature_anchor": None
            if selected
            else vcf_fallback_anchor(
                gene,
                {
                    "documents": tuple(_vcf_view_document(item) for item in pool["documents"]),
                    "omim_context": tuple(_vcf_view_document(item) for item in pool["omim_context"]),
                },
                terms,
                admission.details.get("anchor_evidence_ids", {}).get(gene, ""),
            ),
        }

    return cards


def vcf_pair_prompt(profile, first, second, cards):
    return vcf_dump(
        {
            "task": "Which candidate gene better explains this patient?",
            "patient": vcf_patient(profile),
            "decision_contract": {
                "task_is_final_patient_level_causal_discrimination": True,
                "admission_rank_is_tiebreaker_only": True,
                "decisive_rank_evidence_can_change_the_winner": True,
                "supporting_context_cannot_independently_change_the_winner": True,
                "more_text_or_more_evidence_items_is_not_stronger_evidence": True,
                "missing_information_is_uncertainty_not_negative_evidence": True,
                "exact_identifier_occurrence_alone_does_not_prove_affected_human_evidence": True,
            },
            "candidate_a": cards[first],
            "candidate_b": cards[second],
            "required_output": "A or B",
        }
    )


def _vcf_view_document(value):
    return {
        "prompt_id": value["evidence_id"],
        "title": value["title"],
        "prompt_text": value["text"],
        "source_id": value["source_id"],
        "provenance": value.get("provenance", ()),
        "metadata": value.get("metadata", {}),
        "section": value.get("section", ""),
        "publication_types": value.get("publication_types", ()),
        **{key: value.get(key) for key in ("kind", "provider", "variant_keys", "parent_source_ids")},
    }


def score_directions(candidates, directions, clinical_consumed=False, stage="pairwise"):
    by_pair = {(item["gene_a"], item["gene_b"]): item for item in directions}
    wins, ties, losses, inconsistent = Counter(), Counter(), Counter(), 0
    for index, first in enumerate(candidates):
        for second in candidates[index + 1 :]:
            forward, reverse = by_pair[(first, second)], by_pair[(second, first)]
            if (
                forward["status"] != "ok"
                or reverse["status"] != "ok"
                or forward["selected_gene"] != reverse["selected_gene"]
            ):
                inconsistent += int(forward["status"] == reverse["status"] == "ok")
                ties.update((first, second))
            else:
                winner = forward["selected_gene"]
                wins[winner] += 1
                losses[second if winner == first else first] += 1
    scores = tuple(
        {
            "gene": gene,
            "score": wins[gene] + 0.5 * ties[gene],
            "wins": wins[gene],
            "ties": ties[gene],
            "losses": losses[gene],
        }
        for gene in candidates
    )
    ranked = tuple(item["gene"] for item in sorted(scores, key=lambda item: -item["score"]))
    invalid = sum(item["status"] != "ok" for item in directions)
    return Ranking(
        genes=ranked,
        stage=stage,
        status="degraded" if invalid else "ok",
        details={
            "scores": scores,
            "directions": directions,
            "inconsistent_pairs": inconsistent,
            "invalid_directions": invalid,
            "clinical_consumed": clinical_consumed,
        },
        issues=(f"{invalid}_invalid_direction_outputs",) if invalid else (),
    )


def pairwise_score_gap(result):
    genes = result.genes[:2]
    scores = {item["gene"]: item["score"] for item in result.details.get("scores", ())}
    if len(genes) != 2 or any(gene not in scores for gene in genes):
        return None
    return round(abs(float(scores[genes[0]]) - float(scores[genes[1]])) * 2) / 2


def reflection_patient(profile):
    value = {
        "present_phenotypes": [
            {
                "hpo_id": term["hpo_id"],
                "label": term["label"],
                "definition": term["definition"] or None,
                "information_content": round(term["information_content"], 6),
            }
            for term in views.present_terms(profile)
        ],
        "negated_phenotypes": [
            {"hpo_id": term["hpo_id"], "label": term["label"]} for term in views.absent_terms(profile)
        ],
    }
    views.add_clinical_context(value, profile)
    return value


def reflection_brief_card(gene, briefs):
    value = {"gene": gene, **evidence_card(gene, briefs)}
    value.setdefault("similar_case_evidence", [])
    return value


def gap_prompt(profile, genes, briefs):
    return json.dumps(
        {
            "patient": reflection_patient(profile),
            "candidate_genes_alphabetical": sorted(genes),
            "initial_evidence_briefs": [reflection_brief_card(gene, briefs) for gene in sorted(genes)],
        },
        ensure_ascii=False,
        separators=(",", ":"),
    )


def _text(value, limit):
    return " ".join(str(value or "").split())[:limit]


def _strings(value, limit, length):
    if not isinstance(value, list):
        return ()
    return tuple(dict.fromkeys(_text(item, length) for item in value if _text(item, length)))[:limit]


def invalid_query(value, genes):
    upper = value.upper()
    return (
        "PMID" in upper
        or any(token in upper for token in (" AND ", " OR ", "[TITLE", "[ABSTRACT"))
        or any(re.search(rf"(?<![A-Z0-9]){re.escape(gene)}(?![A-Z0-9])", upper) for gene in genes)
    )


def parse_gap(value, genes, briefs):
    comparison = str(value.get("comparison_status") or "").strip().lower()
    needed, assessments, raw_queries = (
        value.get("retrieval_needed"),
        value.get("candidate_assessments"),
        value.get("queries_by_gene"),
    )
    reason, gap = _text(value.get("reason"), 1200), _text(value.get("comparison_gap"), 1000)
    if (
        comparison not in COMPARISONS
        or not isinstance(needed, bool)
        or not isinstance(assessments, dict)
        or not isinstance(raw_queries, dict)
        or not reason
    ):
        fail("invalid_gap_shape")
    parsed = {}
    for gene in genes:
        item = assessments.get(gene)
        if not isinstance(item, dict) or str(item.get("brief_quality") or "").lower() not in QUALITIES:
            fail("invalid_candidate_assessment")
        brief = brief_get(briefs, gene)
        allowed = (
            (
                {item["document"]["prompt_id"] for item in brief["selected_evidence"]}
                | {item["evidence_id"] for item in brief["similar_case_evidence"]}
            )
            if brief
            else set()
        )
        parsed[gene] = {
            "brief_quality": str(item["brief_quality"]).lower(),
            "useful_evidence_ids": tuple(
                identifier
                for identifier in _strings(item.get("useful_evidence_ids"), 8, 100)
                if identifier in allowed
            ),
            "redundant_evidence_ids": tuple(
                identifier
                for identifier in _strings(item.get("redundant_evidence_ids"), 8, 100)
                if identifier in allowed
            ),
            "missing_discriminator": _text(item.get("missing_discriminator"), 600),
        }
    queries = {}
    for raw_gene, items in raw_queries.items():
        gene, phrases = str(raw_gene).strip().upper(), _strings(items, 2, 180)
        if gene not in genes or any(invalid_query(phrase, genes) for phrase in phrases):
            fail("invalid_query_phrase")
        if phrases:
            queries[gene] = phrases
    if needed and (not gap or not queries):
        fail("retrieval_requires_gap_and_queries")
    return {
        "status": "ok",
        "comparison_status": comparison,
        "candidate_assessments": parsed,
        "retrieval_needed": needed,
        "comparison_gap": gap if needed else "",
        "queries_by_gene": queries if needed else {},
        "reason": reason,
    }


def fill_prompt(profile, genes, briefs, gap, evidence):
    value = json.loads(gap_prompt(profile, genes, briefs))
    value["identified_gap"] = {
        "comparison_status": gap["comparison_status"],
        "comparison_gap": gap["comparison_gap"],
        "queries_by_gene": {gene: list(items) for gene, items in gap["queries_by_gene"].items()},
    }
    value["new_pubmed_abstracts"] = {
        gene: [
            {
                "evidence_id": item["prompt_id"],
                "title": item["title"],
                "abstract_text": item["prompt_text"],
                "already_in_initial_raw_pool": False,
                "matched_query_count": len(item.get("query_hits", ())),
            }
            for item in evidence_pool(evidence, gene)["documents"]
        ]
        for gene in genes
    }
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def parse_fill(value, genes, evidence):
    found, raw = value.get("comparison_discriminator_found"), value.get("gene_results")
    summary, remaining = (
        _text(value.get("discriminator_summary"), 1200),
        _text(value.get("remaining_gap"), 1000),
    )
    if not isinstance(found, bool) or not isinstance(raw, dict):
        fail("invalid_fill_shape")
    results = {}
    for gene in genes:
        item = raw.get(gene)
        if (
            not isinstance(item, dict)
            or not isinstance(item.get("gap_filled"), bool)
            or not (reason := _text(item.get("reason"), 700))
        ):
            fail("invalid_gene_fill_result")
        allowed = {document["prompt_id"] for document in evidence_pool(evidence, gene)["documents"]}
        useful = tuple(
            identifier
            for identifier in _strings(item.get("useful_new_evidence_ids"), 8, 100)
            if identifier in allowed
        )
        results[gene] = {
            "gap_filled": bool(item["gap_filled"] and useful),
            "useful_new_evidence_ids": useful,
            "reason": reason,
        }
    if found and (not summary or not any(item["useful_new_evidence_ids"] for item in results.values())):
        fail("found_discriminator_requires_new_evidence")
    if not found and not remaining:
        fail("unresolved_requires_remaining_gap")
    return {
        "status": "ok",
        "comparison_discriminator_found": found,
        "gene_results": results,
        "discriminator_summary": summary,
        "remaining_gap": remaining,
    }


def useful_evidence(genes, retrieved, results):
    return {
        gene: {
            "gene": gene,
            "documents": tuple(
                item
                for item in evidence_pool(retrieved, gene)["documents"]
                if item["prompt_id"] in set(results[gene]["useful_new_evidence_ids"])
            ),
            "omim_context": (),
            "retrieval_errors": (),
        }
        for gene in genes
    }


def reflection_pair_prompt(profile, pair, briefs, gap, fill, evidence):
    def card(gene):
        value = reflection_brief_card(gene, briefs)
        value["initial_brief_evidence"] = value.pop("selected_evidence")
        value["validated_new_pubmed_evidence"] = [
            {
                "evidence_id": item["prompt_id"],
                "title": item["title"],
                "abstract_text": " ".join(item["prompt_text"].split())[:2800],
                "already_in_initial_raw_pool": False,
            }
            for item in evidence_pool(evidence, gene)["documents"]
        ]
        return value

    return json.dumps(
        {
            "patient": reflection_patient(profile),
            "fixed_comparison_gap": gap["comparison_gap"],
            "validated_discriminator_summary": fill["discriminator_summary"],
            "remaining_gap": fill["remaining_gap"],
            "candidate_a": card(pair[0]),
            "candidate_b": card(pair[1]),
            "required_output": "A or B",
        },
        ensure_ascii=False,
        separators=(",", ":"),
    )
