"""Public V4 projections, with one canonical scientific record graph."""

import json
import math
from collections.abc import Mapping
from copy import deepcopy
from pathlib import PurePath

from mendelprime.domain import (
    MendelPrimeError,
    Source,
    case_record,
    family_origin,
    is_evidence_limitation,
    profile_term,
    require,
    validate_case_key,
    variant_record,
    variant_text,
)
from mendelprime.harness.memory import OTEM

_PRIVATE_TOKENS = frozenset(
    "accesskey accesstoken apikey apitoken attempts auth authorization clientsecret config "
    "credential credentials credentialscope gold goldgene goldgenes idtoken password passwd "
    "prompt prompts rawoutput refreshtoken retries retry retrycount secret secretkey settings "
    "token transportretries signature privatekey sessiontoken xamzcredential xamzsignature "
    "xamzsecuritytoken xgoogsignature xgoogcredential".split()
)


def public_value(value, field="result"):
    if value is None or type(value) in {str, bool, int}:
        return value
    if type(value) is float:
        return value if math.isfinite(value) else None
    if isinstance(value, Mapping):
        output = {}
        for key, item in value.items():
            require(type(key) is str, field, "string_mapping_key_required")
            token = "".join(character for character in key.casefold() if character.isalnum())
            require(
                token not in _PRIVATE_TOKENS and not key.casefold().endswith(("_path", "_paths")),
                key,
                "private_runtime_field_forbidden",
            )
            output[key] = public_value(item, key)
        return output
    if isinstance(value, (list, tuple)):
        return [public_value(item, field) for item in value]
    raise MendelPrimeError(
        "private_path_forbidden" if isinstance(value, PurePath) else "unsupported_public_value", field=field
    )


def _compact(value):
    return {key: item for key, item in value.items() if item is not None and item not in ("", (), [], {})}


def _ranking(ranking, reasoning=False):
    value = {"genes": ranking.genes}
    if reasoning and ranking.details.get("reasoning"):
        value["reasoning"] = dict(ranking.details["reasoning"])
    return value


def _rank_rows(genes, scores=None):
    rows = []
    for rank, gene in enumerate(genes, 1):
        row = {"gene": gene, "rank": rank}
        if (score := (scores or {}).get(gene)) is not None:
            row["score"] = public_value(score)
        rows.append(row)
    return rows


def _hit(hit):
    return {
        **{key: value for key, value in hit.items() if key != "source"},
        "source": hit["source"].provider,
        "source_ref": hit["source"].source_id,
    }


def _similar(item):
    return {**item, "hit": _hit(item["hit"])}


def _similar_result(value):
    def section(prefix, hits):
        return {
            "status": value[prefix + "status"],
            "hits": [_hit(item) for item in hits],
            "decisions": value[prefix + "decisions"],
            "evidence_by_gene": {
                gene: [_similar(item) for item in items]
                for gene, items in value[prefix + "evidence_by_gene"].items()
            },
        }

    return {
        **section("", value["hits"]),
        "top15_support": section("", value["judged_hits"]),
        "tail_challengers": section("tail_", value["tail_hits"]),
    }


def _optional_rankings(stages, final):
    output = {"tail_challenger": {"status": "disabled"}, "reflection": {"status": "disabled"}}
    if (result := stages["challenger"]) is not None:
        details = dict(result.details)
        details["baseline_top100"] = details.pop("baseline")
        output["tail_challenger"] = {
            "status": result.status,
            **details,
            "final_top100": result.genes,
            "ranking": _ranking(result),
        }
    if (result := stages["reflection"]) is not None:
        output["reflection"] = {
            "status": result.status,
            **{key: item for key, item in result.details.items() if key != "evidence_updates"},
            "revised_top2": result.genes,
            "ranking": _ranking(final),
        }
    return output


_AUTO_CONTEXT = object()


def project_result(run, explanation=None, *, context=_AUTO_CONTEXT):
    """Project one completed ranking run to the public V4 contract."""
    require(hasattr(run, "patient") and hasattr(run, "final"), "result.run", "completed_ranking_required")
    value = run.explanation if explanation is None else explanation
    require(
        isinstance(value, Mapping) and value.get("status") in {"disabled", "ok", "degraded"},
        "result.explanation",
        "valid_explanation_required",
    )
    if context is _AUTO_CONTEXT:
        from mendelprime.runtime.explanation import freeze_context

        try:
            context = freeze_context(run) if run.final.genes else None
        except MendelPrimeError as error:
            context = None
            value = {
                "status": "degraded",
                "issues": (f"explanation_context_invalid:{error.field}:{error.code}",),
                "markdown": "",
            }
    if context is None and run.final.genes:
        require(
            value["status"] != "ok"
            and not any(value.get(k) for k in ("report", "assessment", "assessment_progress", "markdown")),
            "result.explanation",
            "validated_context_required",
        )
        value = {
            "status": "degraded",
            "issues": value.get("issues") or ("explanation_context_unavailable",),
            "markdown": "",
        }
    result = (
        _phenotype_result(run, context) if run.patient.genomic is None else _genotype_result(run, context)
    )
    if value["status"] != "disabled":
        if context is None:
            result["status"] = "degraded"
            result["stages"]["final_explanation"] = project_explanation(
                run.patient.key, run.final.genes[:5], value
            )
        else:
            result = with_explanation(result, run.patient.key, context, value)
    issues = tuple(
        dict.fromkeys(
            [code for row in run.stages.get("trace", ()) for code in row["error_codes"]]
            + [
                code
                for item in run.evidence
                if item.kind == "gene_documents"
                for code in item.facts.get("retrieval_errors", ())
            ]
            + (list(value.get("issues", ())) if context is None else [])
        )
    )
    limits = [issue for issue in issues if is_evidence_limitation(issue)]
    if issues:
        target = "genotype" if run.patient.genomic is not None else "evidence_brief"
        result["stages"][target]["evidence_limitations"] = limits
        result["stages"][target]["execution_issues"] = [i for i in issues if not is_evidence_limitation(i)]
    return result


def _phenotype_result(run, context):
    """Project completed HPO/clinical ranking without running any report work."""
    from mendelprime.runtime.explanation import encode_context

    require(run.patient.genomic is None, "result.ranking", "phenotype_run_required")
    require(run.explanation["status"] == "disabled", "result.explanation", "ranking_projection_required")
    patient, stages = run.patient, run.stages
    profile, hypothesis, fusion = stages["profile"], stages["hypotheses"], stages["fusion"]
    admission, pairwise, tail = (stages[key] for key in ("admission", "pairwise", "tail"))
    global_result, briefs = stages["global_result"], stages["briefs"]["results"]
    inputs = {
        "phenotypes": [
            _compact(
                {
                    "hpo_id": term["hpo_id"],
                    "label": term["label"],
                    "status": "present" if term["present"] else "absent",
                }
            )
            for term in patient.phenotypes
        ]
    }
    projected_profile = {"terms": [_compact(profile_term(term)) for term in profile["terms"]]}
    if patient.clinical_description:
        inputs["clinical_description"] = projected_profile["clinical_description"] = (
            patient.clinical_description
        )
    output = {
        "phenotype_profile": projected_profile,
        "candidate_generation": {
            "tool_llm": {
                method: _compact(
                    {
                        "status": value.status,
                        "ranking": {"genes": _rank_rows(value.genes)},
                        "reasoning": dict(value.details["reasoning"]),
                    }
                )
                for method, value in hypothesis["sources"].items()
            },
            "tools": {
                value["tool_id"]: {
                    "genes": _rank_rows(
                        tuple(row["gene"] for row in value["entries"]),
                        {row["gene"]: row["score"] for row in value["entries"]},
                    )
                }
                for value in stages["candidates"]["rankings"]
            },
        },
        "fusion": {
            "status": fusion.status,
            "ranking": _ranking(fusion, True),
            "candidate_union": fusion.details["nominated_genes"],
        },
        "evidence_admission": {
            "candidate_genes": admission.genes,
            "fusion_top10": fusion.genes[:10],
            "external_genes": tuple(gene for gene in admission.genes if gene not in fusion.genes[:10]),
            "source_kinds": admission.details["source_kinds"],
        },
    }
    if stages["similar"] is not None:
        output["similar_cases"] = _similar_result(stages["similar"])
    output.update(
        {
            "evidence_brief": {
                "briefs": {
                    row["gene"]: _compact(
                        {
                            "gene": row["gene"],
                            "selected_evidence": [
                                {"evidence_id": item["document"]["prompt_id"], "summary": item["summary"]}
                                for item in row["selected_evidence"]
                            ],
                            "similar_case_evidence": [
                                _similar(item) for item in row["similar_case_evidence"]
                            ],
                        }
                    )
                    for row in briefs
                },
                "status_by_gene": {row["gene"]: row["status"] for row in briefs},
            },
            "global": {
                "status": global_result.status,
                **(global_result.details.get("assessment") or {}),
                "ranking": _ranking(global_result),
            },
            "pairwise": {
                "status": pairwise.status,
                "ranking": _ranking(pairwise),
                "scores_by_gene": {
                    row["gene"]: {key: value for key, value in row.items() if key != "gene"}
                    for row in pairwise.details.get("scores", ())
                },
                "directions": pairwise.details.get("directions", ()),
                "inconsistent_pair_count": pairwise.details.get("inconsistent_pairs", 0),
                "clinical_description_consumed": pairwise.details.get("clinical_consumed", False),
            },
            "tail": {
                "status": tail.status,
                "llm_ranked_genes": tail.details["llm_genes"],
                "ranking": {"genes": tail.details["tail_genes"]},
            },
            "global_boundary": {"status": "disabled"},
            **_optional_rankings(stages, run.final),
            "final_explanation": {"status": "disabled"},
        }
    )
    return public_value(
        {
            "schema_version": "mgagent_case_result_v4",
            "case_id": patient.case_id,
            "status": run.status,
            "input_mode": patient.input_kind,
            "input": inputs,
            "stages": output,
            "final": {f"top{count}": run.final.genes[:count] for count in (10, 20, 50, 100)},
            "provenance": encode_context(context) if context is not None else None,
        }
    )


def _genotype_result(run, context):
    patient, stages = run.patient, run.stages
    profile, hypothesis, fusion = stages["profile"], stages["hypotheses"], stages["fusion"]
    admission, pairwise, tail = (stages[key] for key in ("admission", "pairwise", "tail"))
    global_result, briefs = stages["global_result"], stages["briefs"]["results"]
    inputs = {
        "phenotypes": [
            _compact(
                {
                    "hpo_id": term["hpo_id"],
                    "label": term["label"],
                    "status": "present" if term["present"] else "absent",
                }
            )
            for term in patient.phenotypes
        ],
        "genomic_input": _compact(
            {
                "assembly": stages["acquisition"]["family"]["assembly"],
                "family_mode": patient.genomic["family_mode"],
                "sex": patient.genomic["sex"],
                "assay_type": patient.genomic["assay_type"],
            }
        ),
    }
    projected_profile = {"terms": [_compact(profile_term(term)) for term in profile["terms"]]}
    if patient.clinical_description:
        inputs["clinical_description"] = projected_profile["clinical_description"] = (
            patient.clinical_description
        )
    output = {
        "phenotype_profile": projected_profile,
        "candidate_generation": {"tool_llm": _tool_llm(hypothesis, ranked_rows=True)},
        "fusion": {"status": fusion.status, "ranking": _ranking(fusion, True)},
        "evidence_admission": _genotype_admission(admission, fusion),
        "similar_cases": {"status": "disabled"},
        "evidence_brief": {
            "briefs": {
                row["gene"]: _compact(
                    {
                        "gene": row["gene"],
                        "selected_evidence": [
                            {"evidence_id": item["document"]["prompt_id"], "summary": item["summary"]}
                            for item in row["selected_evidence"]
                        ],
                        "similar_case_evidence": [_similar(item) for item in row["similar_case_evidence"]],
                    }
                )
                for row in briefs
            },
            "status_by_gene": {row["gene"]: row["status"] for row in briefs},
        },
        "global": {
            "status": global_result.status,
            **(dict(global_result.details.get("assessment") or {})),
            "ranking": _ranking(global_result),
        },
        "pairwise": {
            "status": pairwise.status,
            "ranking": _ranking(pairwise),
            "scores_by_gene": {
                row["gene"]: {key: value for key, value in row.items() if key != "gene"}
                for row in pairwise.details.get("scores", ())
            },
            "directions": pairwise.details.get("directions", ()),
            "inconsistent_pair_count": pairwise.details.get("inconsistent_pairs", 0),
            "clinical_description_consumed": pairwise.details.get("clinical_consumed", False),
        },
        "tail": {
            "status": tail.status,
            "llm_ranked_genes": tail.details["llm_genes"],
            "ranking": {"genes": tail.details["tail_genes"]},
        },
        "global_boundary": {"status": "disabled"},
        **_optional_rankings(stages, run.final),
        "genotype": _genotype_stage(run),
        "final_explanation": {"status": "disabled"},
    }
    from mendelprime.runtime.explanation import encode_context

    return public_value(
        {
            "schema_version": "mgagent_case_result_v4",
            "case_id": patient.case_id,
            "status": run.status,
            "input_mode": patient.input_kind,
            "input": inputs,
            "stages": output,
            "final": {f"top{count}": run.final.genes[:count] for count in (10, 20, 50, 100)},
            "provenance": encode_context(context) if context is not None else None,
        }
    )


def _tool_llm(hypothesis, *, ranked_rows=False):
    return {
        method: _compact(
            {
                "status": ranking.status,
                "ranking": {"genes": _rank_rows(ranking.genes) if ranked_rows else ranking.genes},
                "reasoning": dict(ranking.details.get("reasoning", ())),
            }
        )
        for method, ranking in hypothesis["sources"].items()
    }


def _genotype_admission(admission, fusion):
    return {
        "status": admission.status,
        "candidate_genes": admission.genes,
        "candidate_pool_genes": fusion.genes,
        "ranked_pool_genes": admission.details.get("full_ranking", admission.genes),
        "source_kinds": admission.details["source_kinds"],
    }


def _genotype_stage(run):
    stages, acquisition = run.stages, run.stages["acquisition"]
    native, hypothesis = acquisition["evidence"], stages["hypotheses"]
    failed = {row["tool_id"] for row in acquisition["failures"]}
    source_rankings = {
        tool["tool_id"]: {
            "status": "degraded" if tool["issues"] or tool["tool_id"] in failed else "ok",
            "version": tool["version"],
            "ranking": tool["ranking"].genes,
        }
        for tool in native["tools"]
    }
    admission, fusion = stages["admission"], stages["fusion"]
    ranking_acmg = stages.get("ranking_acmg")
    return {
        "status": "ok"
        if acquisition["status"] == fusion.status == admission.status == "ok"
        and all(value.status == "ok" for value in hypothesis["sources"].values())
        and (ranking_acmg is None or ranking_acmg["status"] == "ok")
        else "degraded",
        "analysis": {"status": acquisition["status"], "source_rankings": source_rankings},
        "variant_evidence": {
            "family": _family(acquisition["family"]),
            "family_observations": _family_observations(native["observations"]),
        },
        "structured_variant_evidence": {"by_gene": _native_variants(native["tools"])},
        "tool_llm": _tool_llm(hypothesis),
        "candidate_pool": {"genes": fusion.genes},
        "candidate_admission": _genotype_admission(admission, fusion),
        **({"ranking_acmg": ranking_acmg} if ranking_acmg is not None else {}),
    }


def _family(family):
    roles = tuple(name for name in ("proband", "father", "mother") if family["roles"].get(name))
    return _compact(
        {
            "assembly": family["assembly"],
            "assembly_basis": family["assembly_basis"],
            "mode": family["mode"],
            "roles": roles,
            "sample_roles": {role: family["roles"][role] for role in roles},
            "sex": family["sex"],
            "role_basis": family["role_basis"],
            "issues": family["issues"],
        }
    )


def _family_observations(observations):
    output = {}
    for observation in observations["observations"]:
        calls = {}
        for call in observation["samples"]:
            calls[call["role"]] = _compact(
                {
                    "sample_id": call["sample_id"],
                    "gt": call["gt"],
                    "alleles": call["alleles"],
                    "ad": call["ad"],
                    "dp": call["dp"],
                    "gq": call["gq"],
                    "pl": call["pl"],
                    "phase_set": call["phase_set"],
                    "filter": call["filter_text"],
                }
            )
        key = variant_text(observation["variant_key"])
        origin = observation.get(
            "origin", family_origin(observation["samples"], observation["alternate_index"])
        )
        output[key] = _compact(
            {
                "variant_key": key,
                "vcf_variant_key": variant_text(observation["vcf_variant_key"]),
                "site_quality": observation["quality"],
                "filter": observation["filter_text"],
                "observations": calls,
                "origin": origin,
                "parentage_status": observation.get("parentage_status", "unavailable"),
            }
        )
    return output


def _native_variants(tools):
    by_gene = {}
    for tool in tools:
        for native_gene in tool["genes"]:
            row = by_gene.setdefault(native_gene["gene"], {"variants": []})
            for variant in native_gene["variants"]:
                row["variants"].append(
                    _compact(
                        {
                            "variant_key": variant_text(variant["key"]) if variant["key"] else None,
                            "chromosome": variant["chromosome"],
                            "position": variant["position"],
                            "reference": variant["reference"],
                            "alternate": variant["alternate"],
                            "transcript": variant["transcript"],
                            "hgvs_c": variant["hgvs_c"],
                            "hgvs_p": variant["hgvs_p"],
                            "genotype": variant["genotype_raw"],
                            "tool_annotation": {
                                "source": tool["tool_id"],
                                "version": tool["version"],
                                **dict(variant["details"]),
                            },
                        }
                    )
                )
    return by_gene


def project_assessment(assessment, seed):
    if isinstance(assessment, Mapping) and assessment.get("workflow") == "clinical_profile_v1":
        return _project_clinical_profile_assessment(assessment, seed)
    return _project_legacy_assessment(assessment, seed)


def _project_legacy_assessment(assessment, seed):
    records = _assessment_record(assessment)
    from mendelprime.runtime.explanation import assessment_origin

    return {
        "status": "ok",
        "origin": assessment_origin(seed),
        "by_gene": _assessment_by_gene(records["final"]),
        **records,
    }


def _project_clinical_profile_assessment(assessment, seed):
    from mendelprime.runtime.explanation import assessment_origin
    from mendelprime.tools.acmg import validate_result

    validate_result(assessment)
    profiles = [_assessment_value(row, public=True, citation_index=False) for row in assessment["profiles"]]
    final = [_assessment_value(row, public=True, citation_index=False) for row in assessment["final"]]
    return {
        "status": "ok",
        "origin": assessment_origin(seed),
        "by_gene": _assessment_by_gene(final),
        "workflow": "clinical_profile_v1",
        "profiles": profiles,
        "final": final,
        "sources": [source.source_id for source in assessment["sources"]],
    }


def _assessment_by_gene(decisions):
    by_gene = {}
    for decision in decisions:
        gene = decision["identity"]["gene"]
        disease = decision["disease"]
        by_gene.setdefault(gene, {"variants": []})["variants"].append(
            {
                **decision,
                "disease_id": disease["identifier"],
                "disease_label": disease["label"],
                "mode_of_inheritance": disease["inheritance"],
            }
        )
    return by_gene


def _assessment_record(assessment):
    require(isinstance(assessment, Mapping), "assessment", "object_required")
    output = {}
    for name in ("independent", "partners", "final"):
        require(type(assessment.get(name)) is tuple, f"assessment.{name}", "tuple_required")
        output[name] = [_assessment_value(row, public=True) for row in assessment[name]]
    return output


def _assessment_value(value, *, public=False, citation_index=True):
    if isinstance(value, Source):
        return {"source_id": value.source_id}
    if type(value) is tuple:
        if len(value) == 5 and isinstance(value[0], str) and value[0] in {"GRCh37", "GRCh38"}:
            return variant_text(value) if public else variant_record(value)
        return [_assessment_value(item, public=public, citation_index=citation_index) for item in value]
    if isinstance(value, Mapping):
        result = {
            key: _assessment_value(item, public=public, citation_index=citation_index)
            for key, item in value.items()
        }
        if public and citation_index and "citations" in result:
            result["source_ids"] = list(dict.fromkeys(row["source_id"] for row in result["citations"]))
        return result
    return public_value(value)


def project_explanation(key, fixed_genes, explanation):
    if explanation["status"] == "disabled":
        return {"status": "disabled"}
    validate_case_key(key)
    genes = list(fixed_genes)
    value = {"status": explanation["status"], "fixed_ranked_genes": genes}
    if genes:
        value["top_gene"] = genes[0]
    structured = {}
    report = explanation.get("report")
    if report is not None:
        require(report["case_key"] == key, "result.explanation.report", "case_mismatch")
        require(report["fixed_genes"] == fixed_genes, "result.explanation.report", "ranking_mismatch")
        value["overall_conclusion"] = " ".join(row["text"] for row in report["summary"])
        value["gene_explanations"] = [
            {"rank": rank, "gene": section["gene"], "paragraphs": public_value(section["paragraphs"])}
            for rank, section in enumerate(report["gene_sections"], 1)
        ]
        value["top_gene_comparisons"] = [
            _compact(
                {
                    "first_gene": row["first_gene"],
                    "second_gene": row["second_gene"],
                    "comparison": row["text"],
                    "evidence_claim_ids": row["evidence_claim_ids"],
                    "favored_gene": row["favored_gene"],
                }
            )
            for row in report["comparisons"]
        ]
        structured["report"] = _report_record(report)
    if "historical_audit" in explanation:
        structured["semantic_audit"] = deepcopy(explanation["historical_audit"])
    if "review_status" in explanation:
        value["review_status"] = explanation["review_status"]
    if structured:
        value["structured_explanation"] = structured
    if explanation.get("markdown"):
        value["report_markdown_en"] = explanation["markdown"]
    if explanation.get("issues"):
        value["issues"] = explanation["issues"]
    for name in ("evidence_limitations", "execution_issues"):
        if name in explanation:
            value[name] = explanation[name]
    return public_value(value)


def _report_record(value):
    if isinstance(value, Source):
        return {"source_id": value.source_id}
    if type(value) is tuple:
        if len(value) == 5 and isinstance(value[0], str) and value[0] in {"GRCh37", "GRCh38"}:
            return variant_text(value)
        return [_report_record(item) for item in value]
    if isinstance(value, Mapping):
        output = {}
        for key, item in value.items():
            if key == "case_key" and "claim_id" in value:
                continue
            if key == "case_key":
                output[key] = case_record(item)
            else:
                output[key] = _report_record(item)
        if "citations" in output:
            output["source_ids"] = list(dict.fromkeys(row["source_id"] for row in output["citations"]))
        return output
    return public_value(value)


def _frozen_ranking_status(public, context):
    """Recover report-only failures only with a case/order-bound ranking receipt."""
    def failed(value):
        if isinstance(value, Mapping):
            if value.get("status") in ("failed", "degraded") or value.get("execution_issues"):
                return True
            statuses = value.get("status_by_gene", {})
            if isinstance(statuses, Mapping) and any(s in ("failed", "degraded") for s in statuses.values()):
                return True
            return any(failed(item) for name, item in value.items()
                       if name not in {"final_explanation", "final_acmg", "assessment_progress"})
        return isinstance(value, (list, tuple)) and any(failed(item) for item in value)

    if failed(public["stages"]):
        return "degraded"
    base = context["base"]
    genes = tuple(base["fixed_genes"])
    receipts = [source for source in base["sources"]
                if source.provider == "mendelprime" and source.kind == "tool_observation"
                and source.locator == "method:ranking-observations/frozen-order"
                and source.version == "ranking-observations-v1" and source.case_key == base["case_key"]
                and set(source.genes) == set(genes) and not source.variant_keys
                and source.source_id in base["available_source_ids"]]
    if len(receipts) == 1:
        try:
            payload = json.loads(receipts[0].content)
            rows = payload["gene_observations"]
            valid = payload["fixed_genes"] == list(genes) and isinstance(rows, list) and len(rows) == len(genes)
            valid = valid and all(isinstance(row, dict) and row.get("gene") == gene
                                  and type(row.get("final_rank")) is int and row["final_rank"] == rank
                                  for rank, (row, gene) in enumerate(zip(rows, genes), 1))
            states = {row.get("ranking_status") for row in rows} if valid else set()
            if states in ({"ok"}, {"degraded"}):
                return states.pop()
        except (ValueError, KeyError, TypeError):
            pass
    return public["status"]


def with_explanation(public, key, context, explanation):
    """Return a validated report-bearing copy; ranking data is never recomputed."""
    from mendelprime.runtime.explanation import (
        _encode,
        assessment_availability,
        encode_context,
        validate_assessment_origin,
        validate_context,
    )

    validate_public_result(public, require_provenance=True)
    validate_case_key(key)
    validate_context(context)
    require(public["case_id"] == key[1] == context["base"]["case_key"][1], "result", "case_mismatch")
    genes = tuple(public["final"]["top100"][:5])
    require(genes == context["base"]["fixed_genes"], "result", "ranking_mismatch")
    copied = deepcopy(public)
    copied["stages"]["final_explanation"] = project_explanation(key, genes, explanation)
    assessment = explanation.get("assessment")
    require(
        context["assessment_seed"] is None
        or assessment is not None
        or explanation.get("assessment_progress") is not None
        or explanation.get("report") is None,
        "result.explanation",
        "vcf_assessment_required_for_report",
    )
    if assessment is not None:
        require(context["assessment_seed"] is not None, "result.explanation", "assessment_only_for_vcf")
        validate_assessment_origin(
            context["assessment_seed"],
            assessment["sources"],
            tuple((row["identity"]["gene"], row["identity"]["variant_key"]) for row in assessment["final"]),
        )
        copied["stages"]["genotype"]["final_acmg"] = project_assessment(
            assessment, context["assessment_seed"]
        )
        copied["stages"]["genotype"].pop("assessment_progress", None)
    progress = explanation.get("assessment_progress")
    require(
        progress is not None or "assessment_availability" not in explanation,
        "result.explanation",
        "coverage_requires_partial_assessment",
    )
    if progress is not None:
        from mendelprime.tools.acmg import validate_progress

        require(
            context["assessment_seed"] is not None
            and assessment is None
            and explanation["status"] == "degraded",
            "result.assessment_progress",
            "partial_assessment_requires_degraded_vcf",
        )
        validate_progress(progress)
        coverage = public_value(assessment_availability(progress))
        require(
            "assessment_availability" not in explanation
            or public_value(explanation["assessment_availability"]) == coverage,
            "result.explanation",
            "classification_coverage_mismatch",
        )
        copied["stages"]["final_explanation"]["assessment_availability"] = coverage
        validate_assessment_origin(
            context["assessment_seed"],
            progress["sources"],
            tuple((row["gene"], row["variant_key"]) for row in progress["expected_identities"]),
        )
        copied["stages"]["genotype"]["assessment_progress"] = _encode(progress)
        copied["stages"]["genotype"].pop("final_acmg", None)
    if explanation["status"] == "degraded":
        copied["status"] = "degraded"
    elif public["stages"].get("final_explanation", {}).get("status") == "degraded":
        copied["status"] = _frozen_ranking_status(public, context)
    registry = OTEM(key)
    registry.add(context["base"]["sources"])
    registry.add(explanation.get("sources", ()))
    required = {source.source_id for source in context["base"]["sources"]}
    if progress is not None:
        registry.add(progress["sources"])
        required.update(source.source_id for source in progress["sources"])
    required.update(_source_references(copied["stages"]["final_explanation"]))
    if assessment is not None:
        required.update(_source_references(copied["stages"]["genotype"]["final_acmg"]))
        if assessment.get("workflow") == "clinical_profile_v1":
            required.update(source.source_id for source in assessment["sources"])
        else:
            required.update(
                source.source_id
                for source in assessment["sources"]
                if source.kind == "guideline"
                or (source.kind != "derived" and set(source.genes).intersection(genes))
            )
    frozen = {
        "base": {**context["base"], "sources": registry.closure(tuple(required))},
        "assessment_seed": context["assessment_seed"],
    }
    copied["provenance"] = encode_context(frozen)
    return public_value(copied)


def _source_references(value):
    references = set()
    if isinstance(value, Mapping):
        for name, item in value.items():
            if name == "source_id" and isinstance(item, str):
                references.add(item)
            elif (name == "source_ids" or name.endswith("_source_ids")) and isinstance(item, (list, tuple)):
                references.update(item)
            references.update(_source_references(item))
    elif isinstance(value, (list, tuple)):
        for item in value:
            references.update(_source_references(item))
    return references


def validate_public_shape(value, *, decoded_json):
    """Common V4 shape only; callers retain their distinct legacy/privacy policies."""
    require(
        type(value) is dict and value.get("schema_version") == "mgagent_case_result_v4",
        "result.schema_version",
        "unsupported_result_schema",
    )
    status = value.get("status")
    require(
        (not decoded_json or type(status) is str) and status in {"ok", "degraded"},
        "result.status",
        "unsupported_result_status",
    )
    require(
        value.get("input_mode") in {"hpo", "hpo_clinical", "hpo_vcf"},
        "result.input_mode",
        "unsupported_input_mode",
    )
    require(
        type(value.get("input")) is dict and type(value.get("stages")) is dict,
        "result",
        "input_and_stages_required",
    )
    final = value.get("final")
    require(type(final) is dict, "result.final", "rank_array_required")
    legacy_vcf = decoded_json and value.get("provenance") is None and value["input_mode"] == "hpo_vcf"
    genes = final.get("top100", final.get("top50") if legacy_vcf else None)
    require(type(genes) is list, "result.final", "rank_array_required")
    require(
        len(genes) <= 100
        and all(
            (type(gene) is str if decoded_json else isinstance(gene, str)) and gene.strip() for gene in genes
        ),
        "result.final.top100",
        "gene_array_required",
    )
    for count in (10, 20, 50):
        require(final.get(f"top{count}") == genes[:count], "result.final", "rank_prefix_mismatch")
    return genes


def validate_public_result(value, *, require_provenance=False):
    require(type(value) is dict, "result", "object_required")
    require(
        set(value)
        == {"schema_version", "case_id", "status", "input_mode", "input", "stages", "final", "provenance"},
        "result",
        "schema_mismatch",
    )
    validate_public_shape(value, decoded_json=False)
    if require_provenance:
        require(
            value["provenance"] is not None, "result.provenance", "standalone_explanation_provenance_missing"
        )
    public_value(value)
    return value


__all__ = [
    "project_assessment",
    "project_explanation",
    "project_result",
    "public_value",
    "validate_public_result",
    "with_explanation",
]
