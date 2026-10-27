"""Fixed-head tail completion and bounded Similar Case challenge."""

import json

from mendelprime.domain import MendelPrimeError, Ranking

from . import views

PHENOTYPE_SOURCES = (
    "phen2disease",
    "phen2gene",
    "local_overlap",
    "amelie_article",
    "direct_hpo",
)


def candidate_rows(candidates, rankings, depth):
    positions = {
        source: {gene: rank for rank, gene in enumerate(values, 1)} for source, values in rankings.items()
    }
    return tuple(
        {
            "gene": gene,
            **{f"{source}_llm_rank": ranks.get(gene) for source, ranks in positions.items()},
            "source_ranks": {f"{source}_llm": ranks.get(gene) for source, ranks in positions.items()},
            "top50_source_count": sum(ranks.get(gene, 51) <= min(depth, 50) for ranks in positions.values()),
        }
        for gene in sorted(candidates)
    )


def tail_prompt(profile, fixed, rows, required, error=""):
    patient = {
        "present_phenotypes": [
            {
                "hpo_id": term["hpo_id"],
                "name": term["label"],
                "specificity": term["information_content"],
                "role": term["role"],
            }
            for term in views.present_terms(profile)
        ]
    }
    views.add_clinical_context(patient, profile)
    payload = {
        "patient_present_phenotypes": patient.pop("present_phenotypes"),
        "fixed_top20_genes_excluded_from_task": list(fixed),
        "tail_candidate_count": len(rows),
        "required_output_count": required,
        "tail_candidates_alphabetical": list(rows),
    }
    if patient:
        payload["patient_context"] = patient
    if error:
        payload["repair_instruction"] = (
            f"The previous response was invalid: {error}. Return exactly {required} distinct allowed genes."
        )
    return json.dumps(payload, ensure_ascii=False, separators=(",", ":"))


def parse_tail(value, allowed, required):
    raw = value.get("ranked_tail_genes")
    if not isinstance(raw, list):
        raw = value.get("ranked_genes")
    if not isinstance(raw, list):
        raise MendelPrimeError("ranked_tail_genes_must_be_list", category="response")
    selected = []
    for row in raw:
        gene = row.get("gene") if isinstance(row, dict) else row
        if isinstance(gene, str):
            gene = gene.strip().upper()
            if gene in allowed and gene not in selected:
                selected.append(gene)
    genes = tuple(selected[:required])
    if required and not genes:
        raise MendelPrimeError("no_valid_allowed_tail_genes", category="response")
    return genes


def tail_rrf(candidates, rankings, k):
    positions = [{gene: rank for rank, gene in enumerate(values, 1)} for values in rankings.values()]

    def key(gene):
        ranks = [mapping[gene] for mapping in positions if gene in mapping]
        return -sum(1.0 / (k + rank) for rank in ranks), -len(ranks), min(ranks), gene

    return tuple(sorted(candidates, key=key))


def tail_inputs(fixed, hypothesis, mode):
    fixed = tuple(dict.fromkeys(gene.upper() for gene in fixed if gene))[:20]
    outcomes = {method: source.genes for method, source in hypothesis["sources"].items()}
    source_ids = PHENOTYPE_SOURCES if mode == "phenotype" else tuple(outcomes)
    rankings = {source: outcomes.get(source, ()) for source in source_ids}
    rrf_depth = 100 if mode == "phenotype" else 50
    output_count = 80 if mode == "phenotype" else max(0, 50 - len(fixed))
    candidates = {gene for values in rankings.values() for gene in values[:50] if gene not in set(fixed)}
    rrf_pool = {gene for values in rankings.values() for gene in values[:rrf_depth] if gene not in set(fixed)}
    rows = candidate_rows(candidates, rankings, 50)
    return {
        "fixed": fixed,
        "rankings": rankings,
        "candidates": candidates,
        "rrf_pool": rrf_pool,
        "rows": rows,
        "required": min(output_count, len(rows)),
        "output_count": output_count,
        "limit": 100 if mode == "phenotype" else 50,
    }


def complete_tail(inputs, llm_genes, fusion, error=""):
    selected = list(llm_genes)
    for gene in (*tail_rrf(inputs["rrf_pool"], inputs["rankings"], 60.0), *fusion):
        if len(selected) >= inputs["output_count"]:
            break
        if gene not in set(inputs["fixed"]) and gene not in selected:
            selected.append(gene)
    tail_genes = tuple(selected[: inputs["output_count"]])
    complete = tuple(dict.fromkeys((*inputs["fixed"], *tail_genes)))[: inputs["limit"]]
    return Ranking(
        genes=complete,
        stage="tail",
        status="ok" if llm_genes or inputs["required"] == 0 else "degraded",
        details={"fixed_genes": inputs["fixed"], "tail_genes": tail_genes, "llm_genes": tuple(llm_genes)},
        issues=(error,) if error else (),
    )


def challenger_card(value):
    hit, decision = value["hit"], value["decision"]
    return {
        "case_id": hit["record_id"],
        "known_causal_gene": hit["gene"],
        "match_type": decision["match_type"],
        "matched_patient_hpo_ids": list(decision["matched_patient_hpo_ids"]),
        "matched_case_hpo_ids": list(decision["matched_case_hpo_ids"]),
        "summary": decision["summary"],
    }


def challenger_prompt(profile, fixed, baseline, evidence, maximum):
    value = {
        "task": "Evaluate bounded Similar Case challengers for insertion into positions 21-100.",
        "patient_present_hpo": [
            {"hpo_id": term["hpo_id"], "label": term["label"]} for term in views.present_terms(profile)
        ],
        "frozen_top20_genes": list(fixed),
        "baseline_positions_21_100": [
            {"rank": rank, "gene": gene} for rank, gene in enumerate(baseline[20:], 21)
        ],
        "max_admitted_challengers": maximum,
        "challengers": [
            {"gene": gene, "similar_cases": [challenger_card(item) for item in cards]}
            for gene, cards in evidence.items()
        ],
        "rules": [
            "Return one decision per challenger; admit or reject.",
            "Never change Top20 or the relative order of non-challenger genes.",
            "Use only supplied case IDs and matched HPO IDs; each decisive_hpo_matches item must be an "
            "object containing case_id, patient_hpo_id, case_hpo_id, and relation; gene identity alone "
            "is insufficient.",
            "An admitted existing gene must move upward and every insertion anchor must be a rank 21-100 "
            "non-challenger.",
        ],
        "required_output": {
            "challenger_decisions": [
                {
                    "gene": "GENE",
                    "decision": "admit | reject",
                    "insert_before_gene": "GENE or empty",
                    "supporting_case_ids": ["supplied case_id"],
                    "decisive_hpo_matches": [
                        {
                            "case_id": "one supporting case_id",
                            "patient_hpo_id": "supplied matched patient HPO",
                            "case_hpo_id": "supplied matched case HPO",
                            "relation": "brief supplied-evidence relation",
                        }
                    ],
                    "reason": "brief reason",
                }
            ]
        },
    }
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def _band(decision, anchor, positions):
    if decision == "reject":
        return "reject"
    position = positions[anchor]
    return (
        "21-30" if position <= 30 else "31-50" if position <= 50 else "51-80" if position <= 80 else "81-100"
    )


def parse_challengers(value, profile, baseline, evidence, maximum):
    rows = value.get("challenger_decisions")
    if not isinstance(rows, list):
        raise MendelPrimeError("challenger_decisions_must_be_list", category="response")
    positions = {gene: rank for rank, gene in enumerate(baseline, 1)}
    anchors, patient_hpo = set(baseline[20:]) - set(evidence), set(views.present_hpo_ids(profile))
    output, seen, admitted = [], set(), 0
    for row in rows:
        if not isinstance(row, dict):
            raise MendelPrimeError("challenger_decision_must_be_object", category="response")
        gene = str(row.get("gene") or "").upper()
        decision = str(row.get("decision") or "").lower()
        anchor = str(row.get("insert_before_gene") or "").upper()
        reason = " ".join(str(row.get("reason") or "").split())[:1200]
        cards = {item["hit"]["record_id"]: item for item in evidence.get(gene, ())}
        raw_ids, matches = row.get("supporting_case_ids", ()), row.get("decisive_hpo_matches")
        case_ids = (
            tuple(dict.fromkeys(str(item) for item in raw_ids if str(item) in cards))
            if isinstance(raw_ids, (list, tuple))
            else ()
        )
        valid_matches = isinstance(matches, list) and all(
            isinstance(item, dict)
            and item.get("case_id") in cards
            and item.get("patient_hpo_id") in patient_hpo
            and item.get("patient_hpo_id")
            in cards[str(item.get("case_id"))]["decision"]["matched_patient_hpo_ids"]
            and item.get("case_hpo_id") in cards[str(item.get("case_id"))]["decision"]["matched_case_hpo_ids"]
            for item in matches
        )
        if (
            not valid_matches
            or gene not in evidence
            or gene in seen
            or decision not in {"admit", "reject"}
            or not reason
        ):
            raise MendelPrimeError("challenger_decision_invalid", category="response")
        if decision == "admit":
            admitted += 1
            if (
                admitted > maximum
                or anchor not in anchors
                or not case_ids
                or not matches
                or positions.get(anchor, 101) >= positions.get(gene, 101)
            ):
                raise MendelPrimeError("challenger_admission_invalid", category="response")
        elif anchor:
            raise MendelPrimeError("rejected_challenger_has_insertion", category="response")
        output.append(
            {
                "gene": gene,
                "decision": decision,
                "insert_before_gene": anchor,
                "target_band": _band(decision, anchor, positions),
                "supporting_case_ids": case_ids,
                "decisive_hpo_matches": tuple(
                    {
                        "case_id": str(item.get("case_id")),
                        "patient_hpo_id": str(item.get("patient_hpo_id")),
                        "case_hpo_id": str(item.get("case_hpo_id")),
                        "relation": str(item.get("relation") or ""),
                    }
                    for item in matches
                ),
                "reason": reason,
            }
        )
        seen.add(gene)
    if seen != set(evidence):
        raise MendelPrimeError("challenger_decisions_incomplete", category="response")
    return tuple(output)


def challenger_result(baseline, decisions=(), status="ok", issue=""):
    admitted = tuple(row["gene"] for row in decisions if row["decision"] == "admit")
    insertions = {}
    for row in decisions:
        if row["decision"] == "admit":
            insertions.setdefault(row["insert_before_gene"], []).append(row["gene"])
    rebuilt = list(baseline[:20])
    for gene in baseline[20:]:
        rebuilt.extend(insertions.get(gene, ()))
        if gene not in admitted:
            rebuilt.append(gene)
    final = tuple(dict.fromkeys(rebuilt))[: len(baseline)]
    return Ranking(
        genes=final,
        stage="tail_challenger",
        status=status,
        details={
            "decisions": decisions,
            "admitted_genes": admitted,
            "inserted_external_genes": tuple(gene for gene in admitted if gene not in baseline),
            "promoted_existing_genes": tuple(gene for gene in admitted if gene in baseline),
            "displaced_genes": tuple(gene for gene in baseline if gene not in final),
            "baseline": baseline,
        },
        issues=(issue,) if issue else (),
    )
