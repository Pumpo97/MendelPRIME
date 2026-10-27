"""Candidate-source, Fusion, and Admission science on plain records."""

import json
import re
from collections.abc import Mapping

from mendelprime import prompting
from mendelprime.domain import MendelPrimeError, Ranking

from . import views

PRIMARY = ("phen2disease", "phen2gene", "local_overlap")
PHENOTYPE_SOURCE_ORDER = (*PRIMARY, "amelie_article", "direct_hpo")
SOURCE_STAGES = {
    "phen2disease": "llm_phen2disease",
    "phen2gene": "llm_phen2gene",
    "local_overlap": "llm_local_overlap",
    "amelie_article": "llm_amelie_article",
    "direct_hpo": "llm_direct",
}


def fail(code):
    raise MendelPrimeError(code, category="response")


def tool_genes(record):
    return tuple(item["gene"] for item in record["entries"])


def nominations(rankings, depth):
    output = {}
    for source, genes in rankings.items():
        for rank, gene in enumerate(genes[:depth], 1):
            output.setdefault(gene, {})[source] = rank
    return output


def nomination_order(values):
    return tuple(
        sorted(
            values,
            key=lambda gene: (
                min(values[gene].values()),
                -len(values[gene]),
                -sum(1.0 / (60.0 + rank) for rank in values[gene].values()),
                gene,
            ),
        )
    )


def admission_order(fusion, rankings, source_depth=10, limit=15):
    head = tuple(fusion[:10])
    positions = {
        source: {gene: rank for rank, gene in enumerate(genes, 1)} for source, genes in rankings.items()
    }
    external = {
        gene
        for values in positions.values()
        for gene, rank in values.items()
        if rank <= source_depth and gene not in set(head)
    }

    def priority(gene):
        ranks = [
            values[gene] for values in positions.values() if gene in values and values[gene] <= source_depth
        ]
        return min(ranks), -len(ranks), -sum(1.0 / rank for rank in ranks), gene

    return (*head, *sorted(external, key=priority))[:limit]


def decode_allowed_ranking(value, allowed):
    ranked = []
    raw = value.get("ranked_genes")
    for item in raw if isinstance(raw, list) else ():
        gene = str(item.get("gene") if isinstance(item, dict) else item or "").strip().upper()
        if gene in allowed and gene not in ranked:
            ranked.append(gene)
    if not ranked:
        fail("no_valid_allowed_genes")
    ranked.extend(gene for gene in allowed if gene not in ranked)
    reasons = value.get("reasoning_by_gene")
    reasoning = (
        tuple((str(gene).upper(), str(reason)) for gene, reason in reasons.items())
        if isinstance(reasons, dict)
        else ()
    )
    return tuple(ranked), reasoning


def decode_direct(value):
    raw = value.get("ranked_genes")
    genes = tuple(
        dict.fromkeys(
            gene
            for item in (raw if isinstance(raw, list) else ())
            if re.fullmatch(r"[A-Z0-9][A-Z0-9._-]*", gene := str(item or "").strip().upper())
        )
    )[:100]
    if not genes:
        fail("no_valid_genes")
    return genes


def decode_vcf_subset(value, allowed):
    raw = value.get("ranked_genes")
    if not isinstance(raw, list):
        fail("ranked_genes_must_be_list")
    genes = tuple(
        dict.fromkeys(str(item).strip().upper() for item in raw if str(item).strip().upper() in allowed)
    )
    if not genes:
        fail("ranked_genes_contains_no_valid_candidate")
    return genes, ()


def decode_vcf_permutation(value, allowed):
    raw = value.get("ranked_genes")
    genes = (
        tuple(str(item).strip().upper() for item in raw if str(item).strip()) if isinstance(raw, list) else ()
    )
    if len(genes) != len(allowed) or len(set(genes)) != len(genes) or set(genes) != set(allowed):
        fail("ranked_genes_must_be_an_exact_candidate_permutation")
    return genes, ()


def decode_vcf_admission(value, allowed, count):
    genes, _ = decode_vcf_permutation({"ranked_genes": value.get("ranked_candidate_genes")}, allowed)
    raw = value.get("top15_genes")
    prefix = (
        tuple(str(item).strip().upper() for item in raw if str(item).strip()) if isinstance(raw, list) else ()
    )
    if prefix != genes[:count]:
        fail("top15_prefix_contract_mismatch")
    return prefix, genes


def genotype_pool(hypothesis, native):
    """Build the Top10 union with an RRF prior from full effective-source ranks."""
    rankings, origins = {}, {}
    for tool in native["tools"]:
        method = tool["tool_id"]
        result = hypothesis["sources"][method]
        successful = result.status == "ok" and bool(result.genes)
        source = f"{method}_llm" if successful else method
        rankings[source] = result.genes if successful else tool["ranking"].genes
        origins[source] = "genotype_tool_llm" if successful else "genotype"
    selected = {gene for genes in rankings.values() for gene in genes[:10]}
    positions = {
        source: {gene: rank for rank, gene in enumerate(genes, 1)} for source, genes in rankings.items()
    }

    def priority(gene):
        ranks = [values[gene] for values in positions.values() if gene in values]
        return -sum(1.0 / (60.0 + rank) for rank in ranks), -len(ranks), min(ranks), gene

    return {
        "genes": tuple(sorted(selected, key=priority))[:20],
        "seed_genes": {source: genes[:10] for source, genes in rankings.items()},
        "source_rank_context": {
            gene: {
                "genotype_source_ranks": {
                    source: {"rank": values.get(gene), "candidate_count": len(rankings[source])}
                    for source, values in positions.items()
                }
            }
            for gene in sorted(selected)
        },
        "source_kinds": origins,
    }


def decode_admission(value, allowed, count):
    rows = value.get("ranked_candidate_genes")
    genes = tuple(dict.fromkeys(str(gene).strip().upper() for gene in rows)) if isinstance(rows, list) else ()
    if len(genes) != count:
        fail(f"expected_{count}_genes")
    if set(genes) - set(allowed):
        fail("disallowed_gene")
    return genes


def _observed(profile):
    return [
        {"hpo_id": term["hpo_id"], "label": term["label"] or term["hpo_id"]}
        for term in views.present_terms(profile)
    ]


def _clinical(value, profile):
    if profile["clinical_description"].strip():
        views.add_clinical_context(value, profile)
    else:
        value.pop("clinical_description", None)


def direct_prompt(profile, count=100):
    value = prompting.template("source_prompts", "direct")
    value["patient_phenotypes"] = _observed(profile)
    _clinical(value, profile)
    value["task"] = f"Return exactly {count} distinct candidate genes, ranked from most to least likely."
    return json.dumps(value, ensure_ascii=False, indent=2)


def tool_prompt(profile, genes, raw, reason_limit=10):
    value = prompting.template("source_prompts", "tool")
    source = raw["tool_id"] if raw else ""
    entries = {item["gene"]: item for item in raw.get("entries", ())} if raw else {}
    names = {term["hpo_id"]: term["label"] for term in views.present_terms(profile)}
    rows = []
    for rank, gene in enumerate(genes, 1):
        row = {
            "consensus_rank": rank,
            "gene": gene,
            "matched_tool_count": 1,
            "mean_tool_rank": float(rank),
            "best_tool_rank": rank,
            "tool_presence": [source] if source else [],
            "has_amelie_pmids": False,
            "source_count": 1,
            "tool_ranks": {source: rank} if source else {},
            "average_tool_rank": float(rank),
        }
        context = entries.get(gene, {}).get("context")
        if isinstance(context, Mapping) and source == "phen2disease":
            row["phen2disease_context"] = {
                "top_disease": context.get("disease_id"),
                "integrated_score": context.get("integrated_score"),
            }
        elif isinstance(context, Mapping) and source == "local_overlap":
            hpo_ids = tuple(context.get("hpo_ids", ()))
            row["local_overlap_context"] = {
                "matched_hpo_ids": list(hpo_ids),
                "matched_phenotypes": [{"hpo_id": item, "name": names.get(item, item)} for item in hpo_ids],
                "matched_term_count": len(hpo_ids),
                "query_term_count": context.get("query_term_count"),
                "gene_term_count": context.get("gene_term_count"),
            }
        rows.append(row)
    patient = value["patient"]
    patient["present_phenotypes"] = [term["label"] for term in views.present_terms(profile)]
    patient["hpo_ids"] = list(views.present_hpo_ids(profile))
    patient["specificity_prioritized_phenotypes"] = [
        term["label"] for term in views.present_terms(profile) if term["role"] == "specificity_prioritized"
    ]
    _clinical(patient, profile)
    value["candidate_contract"].update(allowed_genes=list(genes), candidate_count=len(genes))
    value["candidate_rows"] = rows
    if reason_limit > 0:
        value["ranking_rules"].append(f"Provide reasoning_by_gene for at most the top {reason_limit} genes.")
        value["required_json_schema"].update(
            reasoning_by_gene={"GENE1": "brief reason"},
            case_level_reasoning="brief case-level logic",
            uncertainty="main uncertainty",
        )
    return json.dumps(value, ensure_ascii=False, indent=2)


def _rank_schema():
    return {
        "ranked_genes": ["GENE1", "GENE2"],
        "final_top10_genes": ["GENE1", "GENE2"],
        "reasoning_by_gene": {"GENE1": "brief article-grounded reason for at most Top10"},
        "case_level_reasoning": "brief comparative logic",
        "uncertainty": "main evidence limitation",
    }


def amelie_prompt(profile, genes, articles):
    value = prompting.template("source_prompts", "amelie")
    value["patient"]["observed_phenotypes"] = _observed(profile)
    _clinical(value["patient"], profile)
    value["candidate_contract"]["allowed_genes"] = list(genes)
    value["candidate_packets"] = [
        {"gene": gene, "article": views.article_payload(articles.get(gene, {"gene": gene}))} for gene in genes
    ]
    value["required_json_schema"] = _rank_schema()
    return json.dumps(value, ensure_ascii=False, indent=2)


def fusion_prompt(profile, genes, source_nominations, articles):
    value = prompting.template("source_prompts", "fusion")
    value["patient"]["observed_phenotypes"] = _observed(profile)
    _clinical(value["patient"], profile)
    value["candidate_contract"]["allowed_genes"] = list(genes)
    value["candidate_packets"] = [
        {
            "gene": gene,
            "source_nominations": [
                {"source": f"{method}_llm", "rank": rank}
                for method, rank in sorted(source_nominations[gene].items())
            ],
            "article": views.article_payload(articles.get(gene, {"gene": gene})),
        }
        for gene in sorted(genes)
    ]
    value["baseline_candidate_order"] = list(genes)
    value["required_json_schema"] = _rank_schema()
    return json.dumps(value, ensure_ascii=False, indent=2)


def source_ranking(method, genes, status="ok", error_code="", reasoning=(), kind="tool_llm"):
    stage = f"llm_{method}" if kind == "genotype_tool_llm" else SOURCE_STAGES[method]
    return Ranking(
        genes=tuple(genes),
        stage=stage,
        status=status,
        details={"kind": kind, "reasoning": tuple(reasoning)},
        issues=(error_code,) if error_code else (),
    )
