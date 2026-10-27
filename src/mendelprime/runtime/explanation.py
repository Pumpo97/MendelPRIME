"""Freeze, hydrate, and execute post-ranking explanations without rerunning ranking."""

import hashlib
import json
import math
import re
from collections.abc import Mapping

from mendelprime.domain import (
    MendelPrimeError,
    Source,
    case_record,
    is_evidence_limitation,
    require,
    source_from_record,
    source_record,
    validate_case_key,
    validate_gene,
    validate_verified_paper_view,
    variant_key,
    variant_record,
    variant_text,
    verified_paper_binding_basis,
)
from mendelprime.harness import report as report_workflow
from mendelprime.harness.memory import OTEM

_ENVELOPE = {
    "schema_version",
    "case_key",
    "fixed_genes",
    "sources",
    "base_available_source_ids",
    "assessment_seed",
}
_SEED_FIELDS = {
    "profile",
    "fixed_genes",
    "family",
    "native_genes",
    "observations",
    "family_sources",
    "variant_knowledge",
    "evidence",
    "sources",
}
_VARIANT_FIELDS = {"assembly", "chromosome", "position", "reference", "alternate"}
_PROFILE_FIELDS = {
    "case_key",
    "terms",
    "clinical_description",
    "ontology_manifest",
    "statistics_manifest",
    "consulted_hpo_ids",
}
_TERM_FIELDS = {
    "hpo_id",
    "label",
    "present",
    "information_content",
    "depth",
    "priority_score",
    "role",
    "definition",
    "parent_ids",
    "parent_labels",
    "ontology_known",
}
_FAMILY_FIELDS = {
    "assembly",
    "assembly_basis",
    "mode",
    "roles",
    "sex",
    "role_basis",
    "pedigree",
    "issues",
}
_NATIVE_GENE_FIELDS = {"gene", "records", "variants", "unparsed_variants"}
_NATIVE_VARIANT_FIELDS = {
    "gene",
    "key",
    "chromosome",
    "position",
    "reference",
    "alternate",
    "transcript",
    "hgvs_c",
    "hgvs_p",
    "genotype_raw",
    "details",
    "source",
}
_EXOMISER_RECORD_FIELDS = {
    "source",
    "rank",
    "moi",
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
    "human_phenotype_evidence",
    "mouse_phenotype_evidence",
    "fish_phenotype_evidence",
    "human_ppi_evidence",
    "mouse_ppi_evidence",
    "fish_ppi_evidence",
}
_LIRICAL_RECORD_FIELDS = {
    "source",
    "json_source",
    "gene_mapping_source",
    "rank",
    "disease_id",
    "disease_name",
    "pretest_probability",
    "posttest_probability",
    "composite_lr",
    "genotype_match_type",
    "genotype_explanation",
}
_EXOMISER_DETAILS = {
    "moi",
    "consequence",
    "quality",
    "filter",
    "contributing",
    "whitelist",
    "variant_score",
    "acmg_classification",
    "acmg_evidence",
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
    "all_frequencies",
    "max_pathogenicity_source",
    "max_pathogenicity_score",
    "all_pathogenicity_scores",
}
_LIRICAL_DETAILS = {"pathogenicity", "disease_id", "raw_annotation"}
_OBSERVATION_FIELDS = {
    "variant_key",
    "alternate_index",
    "samples",
    "quality",
    "filter_text",
    "raw_record",
    "vcf_variant_key",
}
_CALL_FIELDS = {
    "role",
    "sample_id",
    "gt",
    "alleles",
    "ad",
    "dp",
    "gq",
    "pl",
    "phase_set",
    "filter_text",
    "raw_fields",
}
_DOCUMENT_FIELDS = {
    "source",
    "prompt_id",
    "title",
    "prompt_text",
    "provenance",
    "publication_types",
    "query_hits",
    "cited_pmids",
    "section",
}


def freeze_context(run):
    """Freeze only evidence needed to explain an already fixed ranking."""
    require(hasattr(run, "patient") and hasattr(run, "final"), "report.run", "completed_ranking_required")
    patient, stages, genes = run.patient, run.stages, run.final.genes[:5]
    require(bool(genes), "report.run", "ranking_required_for_explanation")
    memory, available, metadata = OTEM(patient.key), set(), set()
    memory.add(run.sources)
    seed = _assessment_seed(run) if patient.genomic is not None else None

    if patient.genomic is not None:
        receipt = Source(
            provider="patient",
            kind="patient_observation",
            locator="input:patient/genomic-settings",
            version="submitted-genomic-settings-v1",
            case_key=patient.key,
            content=json.dumps(
                _encode(
                    {
                        "assembly": patient.genomic["assembly"],
                        "family_mode": patient.genomic["family_mode"],
                        "sex": patient.genomic["sex"],
                        "assay_type": patient.genomic["assay_type"],
                    }
                ),
                ensure_ascii=False,
                sort_keys=True,
            ),
        )
        memory.add((receipt,))
        metadata.add(receipt.source_id)

    def expose(source):
        memory.add((source,))
        available.add(source.source_id)

    expose(memory.register_patient(patient))
    term_ids = {term["hpo_id"] for term in stages["profile"]["terms"]}
    for source in run.sources:
        if source.kind == "ontology" and source.locator.rsplit("#", 1)[-1] in term_ids:
            expose(source)
    expose(
        Source(
            provider="patient",
            kind="patient_observation",
            locator="input:patient/data-availability",
            version="submitted-input-types-v1",
            case_key=patient.key,
            content=json.dumps(
                {
                    "genomic_input_supplied": patient.genomic is not None,
                    "clinical_description_supplied": bool(patient.clinical_description.strip()),
                    "interpretation": "Availability of submitted input, not evidence of normal clinical tests",
                },
                sort_keys=True,
            ),
        )
    )
    evidence = {item.genes[0]: item.facts for item in run.evidence if item.kind == "gene_documents"}
    briefs = {row["gene"]: row for row in stages["briefs"]["results"]}
    scores = {row["gene"]: row for row in stages["pairwise"].details.get("scores", ())}
    observations = []
    for gene in genes:
        pool, brief = evidence.get(gene), briefs.get(gene)
        if pool:
            for document in (*pool["documents"], *pool["omim_context"]):
                if document.get("section") != "article_lineage":
                    expose(document["source"])
        if brief:
            for item in brief["selected_evidence"]:
                expose(item["document"]["source"])
            for item in brief["similar_case_evidence"]:
                expose(item["hit"]["source"])
        observations.append(
            {
                "gene": gene,
                "final_rank": run.final.position(gene),
                "global_rank": stages["global_result"].position(gene),
                "pairwise_rank": stages["pairwise"].position(gene),
                "pairwise_score": scores.get(gene),
                "ranking_status": run.status,
                "evidence_retrieval_issues": pool["retrieval_errors"] if pool else (),
            }
        )
    expose(
        Source(
            provider="mendelprime",
            kind="tool_observation",
            locator="method:ranking-observations/frozen-order",
            version="ranking-observations-v1",
            genes=genes,
            case_key=patient.key,
            content=json.dumps(
                _encode(
                    {
                        "interpretation": "Recorded software behavior, not independent medical evidence",
                        "fixed_genes": genes,
                        "gene_observations": observations,
                        "global_reasoning": stages["global_result"].details.get("assessment"),
                    }
                ),
                ensure_ascii=False,
                sort_keys=True,
            ),
        )
    )
    required = (
        available
        | metadata
        | ({source.source_id for source in seed["sources"]} if seed is not None else set())
        | set(stages.get("ranking_acmg", {}).get("source_ids", ()))
    )
    context = {
        "base": {
            "case_key": patient.key,
            "fixed_genes": genes,
            "sources": memory.closure(tuple(required)),
            "available_source_ids": tuple(
                source.source_id for source in memory.sources if source.source_id in available
            ),
            "assessment": None,
        },
        "assessment_seed": seed,
    }
    if seed is not None and stages.get("ranking_acmg"):
        context["review"] = _review_context(
            seed, stages["ranking_acmg"], patient.genomic, context["base"]["sources"]
        )
    validate_context(context)
    return context


def _assessment_seed(run):
    patient, stages, genes = run.patient, run.stages, run.final.genes[:5]
    native = stages["acquisition"]["evidence"]
    native_genes = tuple(row for tool in native["tools"] for row in tool["genes"] if row["gene"] in genes)
    keys = {
        variant["key"] for row in native_genes for variant in row["variants"] if variant["key"] is not None
    }
    observations = {
        **native["observations"],
        "observations": tuple(
            {
                name: value
                for name, value in row.items()
                if name not in {"origin", "parentage_status", "match_method"}
            }
            for row in native["observations"]["observations"]
            if row["variant_key"] in keys
        ),
        "missing_targets": tuple(key for key in native["observations"]["missing_targets"] if key in keys),
    }
    family_sources = tuple(
        source
        for source in native["family_sources"]
        if len(source.genes) == 1 and source.genes[0] in genes and set(source.variant_keys).intersection(keys)
    )
    clinical = {gene: value for gene, value in stages["variant_knowledge"].items() if gene in genes}
    evidence = {
        item.genes[0]: item.facts
        for item in run.evidence
        if item.kind == "gene_documents" and item.genes[0] in genes
    }
    # Expansion coverage is execution state. Source text and parent links retain
    # article provenance in the unchanged standalone-report document contract.
    evidence = {
        gene: {
            **{name: value for name, value in row.items() if name != "variant_expansion"},
            **{
                name: tuple(
                    {key: value for key, value in doc.items() if key != "metadata"} for doc in row[name]
                )
                for name in ("documents", "omim_context")
            },
        }
        for gene, row in evidence.items()
    }
    required = set()
    for row in native_genes:
        required.update(source.source_id for source in _sources_in(row))
    required.update(source.source_id for source in family_sources)
    for value in clinical.values():
        required.update(source.source_id for source in _sources_in(value))
        required.update(seed["discovered_from_source_id"] for seed in value["citation_seeds"])
    for value in evidence.values():
        required.update(source.source_id for source in _sources_in(value))
    registry = OTEM(patient.key)
    registry.add(run.sources)
    seed = {
        "profile": stages["profile"],
        "fixed_genes": genes,
        "family": stages["acquisition"]["family"],
        "native_genes": native_genes,
        "observations": observations,
        "family_sources": family_sources,
        "variant_knowledge": clinical,
        "evidence": evidence,
        "sources": registry.closure(tuple(required)),
    }
    _validate_seed(seed, patient.key, genes)
    return seed


def _review_context(seed, ledger, genomic, sources):
    """Transfer allele nominations and cited facts, never a rank or final grade."""
    native = {
        (row["gene"], v["key"]) for row in seed["native_genes"] for v in row["variants"] if v.get("key")
    }
    registry = {source.source_id: source for source in sources}
    nominations, facts = {}, {}
    for gene in seed["fixed_genes"]:
        entry = ledger.get("by_gene", {}).get(gene, {})
        variants = {row["variant_ref"]: row for row in entry.get("variants", ())}
        require(len(variants) == len(entry.get("variants", ())), "review", "duplicate_variant_reference")
        refs = tuple(entry.get("selected_variant_refs", ()))
        require(
            len(set(refs)) == len(refs) and all(ref in variants for ref in refs),
            "review",
            "unknown_or_duplicate_nomination",
        )
        addresses = {
            ref: _identity({"identity": {name: row["identity"][name] for name in ("gene", "variant_key")}})
            for ref, row in variants.items()
        }
        selected = tuple(addresses[ref][1] for ref in refs)
        require(all((gene, key) in native for key in selected), "review", "nominated_allele_not_observed")
        nominations[gene] = selected
        for row in variants.values():
            identity = addresses[row["variant_ref"]]
            require(identity in native and identity[0] == gene, "review", "ranking_allele_identity_mismatch")
            retained = []
            for criterion in row.get("criteria", ()):
                citations = criterion.get("citations", ())
                if not citations:
                    continue
                require(
                    all(
                        c.get("source_id") in registry
                        and c.get("quote")
                        and c["quote"] in registry[c["source_id"]].content
                        for c in citations
                    ),
                    "review",
                    "ranking_quote_not_in_source",
                )
                retained.append(
                    {
                        name: criterion.get(name)
                        for name in (
                            "code",
                            "status",
                            "strength",
                            "rationale",
                            "citations",
                            "independence_group",
                        )
                    }
                )
                retained[-1]["rationale"] = criterion.get("rationale") or criterion.get("method") or ""
                retained[-1]["citations"] = tuple(dict(c) for c in citations)
                retained[-1]["counted"] = criterion.get("counted", True)
            facts[identity] = tuple(retained)
    return {
        "case_key": seed["profile"]["case_key"],
        "nominations": nominations,
        "facts": facts,
        "vcf_path": str((genomic or {}).get("vcf_path") or ""),
        "sources": sources,
    }


def _sources_in(value):
    if isinstance(value, Source):
        yield value
    elif isinstance(value, Mapping):
        for item in value.values():
            yield from _sources_in(item)
    elif isinstance(value, (tuple, list)):
        for item in value:
            yield from _sources_in(item)


def encode_context(context):
    validate_context(context)
    base = context["base"]
    return {
        "schema_version": "mendelprime_explanation_input_v1",
        "case_key": case_record(base["case_key"]),
        "fixed_genes": list(base["fixed_genes"]),
        "sources": [source_record(source) for source in base["sources"]],
        "base_available_source_ids": list(base["available_source_ids"]),
        "assessment_seed": _encode(context["assessment_seed"]),
    }


def load_context(public, key, *, for_assessment=False):
    """Hydrate current or old public V4 evidence, rejecting branch or origin drift."""
    from mendelprime.runtime.results import validate_public_result

    try:
        json.dumps(public, allow_nan=False)
    except (TypeError, ValueError, RecursionError) as error:
        raise MendelPrimeError("strict_json_required", field="result") from error
    validate_public_result(public, require_provenance=True)
    validate_case_key(key)
    require(public["case_id"] == key[1], "result.case_id", "case_mismatch")
    genes = tuple(public["final"]["top100"][:5])
    raw = public["provenance"]
    require(type(raw) is dict and set(raw) == _ENVELOPE, "result.provenance", "schema_mismatch")
    require(
        raw["schema_version"] == "mendelprime_explanation_input_v1",
        "result.provenance",
        "unsupported_provenance_schema",
    )
    require(raw["case_key"] == case_record(key), "result.provenance", "case_mismatch")
    require(raw["fixed_genes"] == list(genes), "result.provenance", "ranking_mismatch")
    require(type(raw["sources"]) is list, "result.provenance.sources", "array_required")
    sources = tuple(source_from_record(row) for row in raw["sources"])
    require(
        len({source.source_id for source in sources}) == len(sources),
        "result.provenance.sources",
        "duplicate_source",
    )
    registry = {source.source_id: source for source in sources}
    available = raw["base_available_source_ids"]
    require(
        type(available) is list and len(set(available)) == len(available),
        "result.provenance.base_available_source_ids",
        "unique_array_required",
    )
    require(
        all(source_id in registry for source_id in available),
        "result.provenance.base_available_source_ids",
        "unknown_source",
    )
    seed = _decode(raw["assessment_seed"], registry) if raw["assessment_seed"] is not None else None
    context = {
        "base": {
            "case_key": key,
            "fixed_genes": genes,
            "sources": sources,
            "available_source_ids": tuple(available),
            "assessment": None,
        },
        "assessment_seed": seed,
    }
    ledger = public["stages"].get("genotype", {}).get("ranking_acmg")
    if (
        seed is not None
        and ledger
        and (for_assessment or not public["stages"].get("genotype", {}).get("final_acmg"))
    ):
        context["review"] = _review_context(seed, ledger, public.get("input", {}).get("genomic"), sources)
    validate_context(context)
    _validate_input_origin(public, context)
    return context


def validate_context(context):
    require(
        isinstance(context, Mapping)
        and set(context) in ({"base", "assessment_seed"}, {"base", "assessment_seed", "review"}),
        "explanation.context",
        "schema_mismatch",
    )
    base = context["base"]
    require(
        isinstance(base, Mapping)
        and set(base) == {"case_key", "fixed_genes", "sources", "available_source_ids", "assessment"},
        "explanation.context.base",
        "schema_mismatch",
    )
    validate_case_key(base["case_key"])
    genes = base["fixed_genes"]
    require(
        type(genes) is tuple and 1 <= len(genes) <= 5 and len(set(genes)) == len(genes),
        "fixed_genes",
        "top1_to_top5_tuple_required",
    )
    for gene in genes:
        validate_gene(gene)
    sources = base["sources"]
    require(
        type(sources) is tuple and all(isinstance(source, Source) for source in sources),
        "sources",
        "source_tuple_required",
    )
    source_ids = {source.source_id for source in sources}
    require(len(source_ids) == len(sources), "sources", "duplicate_source")
    memory = OTEM(base["case_key"])
    memory.add(sources)
    available = base["available_source_ids"]
    require(
        type(available) is tuple and available and len(set(available)) == len(available),
        "available_source_ids",
        "nonempty_tuple_required",
    )
    for source_id in available:
        require(source_id in source_ids, "available_source_ids", "unknown_source")
        source = memory.get_source(source_id)
        require(
            not source.genes or not set(source.genes).isdisjoint(genes),
            "available_source_ids",
            "unrelated_gene_source",
        )
    require(
        base["assessment"] is None or isinstance(base["assessment"], Mapping),
        "assessment",
        "assessment_record_invalid",
    )
    seed = context["assessment_seed"]
    if "review" in context:
        require(
            seed is not None and context["review"]["case_key"] == base["case_key"], "review", "case_mismatch"
        )
    if seed is not None:
        _validate_seed(seed, base["case_key"], genes)
        recorded = {source.source_id: source for source in sources}
        require(
            all(recorded.get(source.source_id) == source for source in seed["sources"]),
            "report.assessment",
            "current_run_source_changed_or_missing",
        )
    return context


def _validate_seed(seed, case_key, genes):
    require(isinstance(seed, Mapping) and set(seed) == _SEED_FIELDS, "assessment_seed", "schema_mismatch")
    require(seed["fixed_genes"] == genes, "assessment_seed.fixed_genes", "ranking_mismatch")
    profile, family, observations = seed["profile"], seed["family"], seed["observations"]
    _keys(profile, _PROFILE_FIELDS, "assessment_seed.profile")
    require(profile["case_key"] == case_key, "assessment_seed.profile", "case_mismatch")
    require(isinstance(profile["clinical_description"], str), "assessment_seed.profile", "text_required")
    require(
        type(profile["terms"]) is tuple and type(profile["consulted_hpo_ids"]) is tuple,
        "assessment_seed.profile",
        "tuple_required",
    )
    seen_hpo = set()
    for term in profile["terms"]:
        _keys(term, _TERM_FIELDS, "assessment_seed.profile.terms")
        require(
            isinstance(term["hpo_id"], str)
            and term["hpo_id"].startswith("HP:")
            and type(term["present"]) is bool,
            "assessment_seed.profile.terms",
            "invalid_term",
        )
        require(term["hpo_id"] not in seen_hpo, "assessment_seed.profile.terms", "duplicate_hpo_identifier")
        require(
            type(term["parent_ids"]) is tuple and type(term["parent_labels"]) is tuple,
            "assessment_seed.profile.terms",
            "tuple_required",
        )
        seen_hpo.add(term["hpo_id"])
    manifests = (profile["ontology_manifest"],)
    if profile["statistics_manifest"] is not None:
        manifests += (profile["statistics_manifest"],)
    for manifest in manifests:
        _keys(manifest, {"locator", "version", "sha256"}, "assessment_seed.profile.manifest")
        require(
            all(isinstance(manifest[name], str) and manifest[name] for name in manifest),
            "assessment_seed.profile.manifest",
            "nonempty_text_required",
        )

    _keys(family, _FAMILY_FIELDS, "assessment_seed.family")
    require(family["assembly"] in {"GRCh37", "GRCh38"}, "assessment_seed.family", "assembly_enum_required")
    require(family["mode"] in {"trio", "proband_only"}, "assessment_seed.family", "family_mode_invalid")
    require(family["sex"] in {"MALE", "FEMALE", "UNKNOWN"}, "assessment_seed.family", "sex_invalid")
    require(
        isinstance(family["roles"], Mapping) and set(family["roles"]) == {"proband", "father", "mother"},
        "assessment_seed.family.roles",
        "schema_mismatch",
    )
    role_values = tuple(sample for sample in family["roles"].values() if sample is not None)
    require(
        bool(family["roles"]["proband"]) and len(role_values) == len(set(role_values)),
        "assessment_seed.family.roles",
        "roles_must_be_distinct",
    )
    require(
        (family["mode"] == "trio") == bool(family["roles"]["father"] and family["roles"]["mother"]),
        "assessment_seed.family.roles",
        "family_mode_role_mismatch",
    )

    _keys(observations, {"case_key", "observations", "missing_targets"}, "assessment_seed.observations")
    require(observations["case_key"] == case_key, "assessment_seed.observations", "case_mismatch")
    require(
        type(observations["observations"]) is tuple and type(observations["missing_targets"]) is tuple,
        "assessment_seed.observations",
        "tuple_required",
    )
    observed_keys = set()
    for row in observations["observations"]:
        _validate_observation(row, family)
        require(row["variant_key"] not in observed_keys, "assessment_seed.observations", "duplicate_variant")
        observed_keys.add(row["variant_key"])
    for key in observations["missing_targets"]:
        _assembly(key, family["assembly"], "assessment_seed.observations.missing_targets")

    for name in ("native_genes", "family_sources", "sources"):
        require(type(seed[name]) is tuple, f"assessment_seed.{name}", "tuple_required")
    require(
        all(isinstance(source, Source) for source in (*seed["family_sources"], *seed["sources"])),
        "assessment_seed.sources",
        "source_tuple_required",
    )
    require(
        isinstance(seed["variant_knowledge"], Mapping) and isinstance(seed["evidence"], Mapping),
        "assessment_seed",
        "mapping_required",
    )
    source_registry = {source.source_id: source for source in seed["sources"]}
    require(len(source_registry) == len(seed["sources"]), "assessment_seed.sources", "duplicate_source")
    for source in seed["sources"]:
        require(
            all(parent in source_registry for parent in source.parent_source_ids),
            "assessment_seed.sources",
            "unknown_source",
        )
        require(
            all(key[0] == family["assembly"] for key in source.variant_keys),
            "assessment_seed.sources",
            "assembly_mismatch",
        )
    for source in _sources_in(seed):
        require(
            source_registry.get(source.source_id) == source,
            "assessment_seed.sources",
            "current_run_source_changed_or_missing",
        )

    native_keys = set()
    for row in seed["native_genes"]:
        _validate_native_gene(row, genes, case_key, family["assembly"], native_keys)
    require(
        observed_keys <= native_keys | set(observations["missing_targets"]),
        "assessment_seed.observations",
        "unknown_variant",
    )
    for source in seed["family_sources"]:
        require(
            source.case_key == case_key
            and len(source.genes) == 1
            and source.genes[0] in genes
            and bool(set(source.variant_keys).intersection(native_keys)),
            "assessment_seed.family_sources",
            "source_scope_mismatch",
        )
    for mapping_name in ("variant_knowledge", "evidence"):
        mapping = seed[mapping_name]
        require(set(mapping) <= set(genes), f"assessment_seed.{mapping_name}", "unrelated_gene")
        for gene, value in mapping.items():
            if mapping_name == "variant_knowledge":
                _validate_knowledge(gene, value, family["assembly"])
            else:
                _validate_evidence(gene, value)


def _keys(value, expected, field):
    require(isinstance(value, Mapping) and set(value) == expected, field, "schema_mismatch")


def _assembly(key, assembly, field):
    require(type(key) is tuple and len(key) == 5, field, "variant_key_invalid")
    variant_key(*key)
    require(key[0] == assembly, field, "assembly_mismatch")


def _validate_observation(row, family):
    field = "assessment_seed.observations.observations"
    _keys(row, _OBSERVATION_FIELDS, field)
    _assembly(row["variant_key"], family["assembly"], field)
    _assembly(row["vcf_variant_key"], family["assembly"], field)
    require(
        type(row["alternate_index"]) is int and row["alternate_index"] > 0,
        field,
        "positive_alternate_index_required",
    )
    require(row["quality"] is None or type(row["quality"]) in {int, float}, field, "quality_invalid")
    require(
        isinstance(row["filter_text"], str) and isinstance(row["raw_record"], str), field, "text_required"
    )
    require(type(row["samples"]) is tuple, field, "sample_tuple_required")
    calls = {}
    for call in row["samples"]:
        _keys(call, _CALL_FIELDS, f"{field}.samples")
        role = call["role"]
        require(
            role in {"proband", "father", "mother"} and role not in calls,
            f"{field}.samples",
            "sample_role_invalid",
        )
        require(call["sample_id"] == family["roles"].get(role), f"{field}.samples", "family_role_mismatch")
        require(call["gt"] is None or isinstance(call["gt"], str), f"{field}.samples", "gt_invalid")
        require(
            type(call["alleles"]) is tuple
            and all(allele is None or type(allele) is int and allele >= 0 for allele in call["alleles"]),
            f"{field}.samples",
            "alleles_invalid",
        )
        if call["gt"] is not None:
            tokens = re.split(r"[/|]", call["gt"])
            require(
                bool(tokens) and all(token == "." or token.isdigit() for token in tokens),
                f"{field}.samples.gt",
                "gt_invalid",
            )
            parsed = tuple(None if token == "." else int(token) for token in tokens)
            require(parsed == call["alleles"], f"{field}.samples", "gt_allele_mismatch")
        for name in ("ad", "pl", "raw_fields"):
            require(
                call[name] is None or type(call[name]) is tuple,
                f"{field}.samples.{name}",
                "tuple_or_none_required",
            )
        for name in ("dp", "gq"):
            require(
                call[name] is None or type(call[name]) in {int, float} and call[name] >= 0,
                f"{field}.samples.{name}",
                "nonnegative_number_or_none_required",
            )
        require(
            call["phase_set"] is None or isinstance(call["phase_set"], str),
            f"{field}.samples.phase_set",
            "text_or_none_required",
        )
        calls[role] = call
    require(
        set(calls) == {role for role, sample in family["roles"].items() if sample},
        f"{field}.samples",
        "family_role_mismatch",
    )


def _validate_native_gene(row, genes, case_key, assembly, native_keys):
    field = "assessment_seed.native_genes"
    _keys(row, _NATIVE_GENE_FIELDS, field)
    gene = row["gene"]
    require(gene in genes, field, "unrelated_gene")
    require(type(row["records"]) is tuple and row["records"], field, "nonempty_record_tuple_required")
    providers = set()
    for record in row["records"]:
        keys = set(record) if isinstance(record, Mapping) else set()
        require(keys == _EXOMISER_RECORD_FIELDS or keys == _LIRICAL_RECORD_FIELDS, field, "schema_mismatch")
        provider = "exomiser" if keys == _EXOMISER_RECORD_FIELDS else "lirical"
        providers.add(provider)
        for name in ("source", "json_source", "gene_mapping_source"):
            source = record.get(name)
            if source is not None:
                require(isinstance(source, Source), field, "source_required")
                if name == "gene_mapping_source":
                    require(
                        source.provider == "hgnc"
                        and source.kind == "annotation"
                        and source.genes == (gene,)
                        and source.case_key is None,
                        field,
                        "gene_mapping_scope_mismatch",
                    )
                else:
                    require(
                        source.provider == provider
                        and source.genes == (gene,)
                        and source.case_key == case_key
                        and not source.variant_keys,
                        field,
                        "source_scope_mismatch",
                    )
    require(len(providers) == 1, field, "mixed_tool_records")
    provider = next(iter(providers))
    require(
        type(row["variants"]) is tuple and type(row["unparsed_variants"]) is tuple, field, "tuple_required"
    )
    for variant in row["variants"]:
        _keys(variant, _NATIVE_VARIANT_FIELDS, f"{field}.variants")
        require(
            variant["gene"] == gene and isinstance(variant["source"], Source),
            f"{field}.variants",
            "gene_variant_mismatch",
        )
        require(
            variant["source"].provider == provider
            and variant["source"].case_key == case_key
            and variant["source"].genes == (gene,),
            f"{field}.variants",
            "source_scope_mismatch",
        )
        details = variant["details"]
        require(
            isinstance(details, Mapping)
            and (set(details) == _EXOMISER_DETAILS or set(details) == _LIRICAL_DETAILS),
            f"{field}.variants.details",
            "schema_mismatch",
        )
        key = variant["key"]
        if key is not None:
            _assembly(key, assembly, f"{field}.variants.key")
            require(variant["source"].variant_keys == (key,), f"{field}.variants", "allele_scope_mismatch")
            chromosome = variant["chromosome"].removeprefix("chr")
            chromosome = "MT" if chromosome == "M" else chromosome
            require(
                (chromosome, variant["position"], variant["reference"].upper(), variant["alternate"].upper())
                == key[1:],
                f"{field}.variants",
                "reported_allele_mismatch",
            )
            native_keys.add(key)
        else:
            require(not variant["source"].variant_keys, f"{field}.variants", "allele_scope_mismatch")
    require(
        all(
            isinstance(source, Source)
            and source.provider == provider
            and source.case_key == case_key
            and source.genes == (gene,)
            and not source.variant_keys
            for source in row["unparsed_variants"]
        ),
        field,
        "unparsed_source_invalid",
    )


def _validate_document(value, gene, field):
    _keys(value, _DOCUMENT_FIELDS, field)
    require(
        isinstance(value["source"], Source) and value["source"].genes == (gene,),
        field,
        "source_scope_mismatch",
    )
    for name in ("provenance", "publication_types", "query_hits", "cited_pmids"):
        require(type(value[name]) is tuple, f"{field}.{name}", "tuple_required")


def _validate_knowledge(gene, value, assembly):
    field = "assessment_seed.variant_knowledge"
    _keys(value, {"gene", "variant_evidence", "gene_context", "citation_seeds", "lookups"}, field)
    require(value["gene"] == gene, field, "gene_result_mismatch")
    for name in ("variant_evidence", "gene_context", "citation_seeds", "lookups"):
        require(type(value[name]) is tuple, f"{field}.{name}", "tuple_required")
    for item in value["variant_evidence"]:
        _keys(
            item,
            {
                "gene",
                "variant_key",
                "match_basis",
                "assertion_scope",
                "document",
                "cited_pmids",
                "normalized_hgvs",
                "linked_source_id",
                "retracted",
            },
            field,
        )
        require(item["gene"] == gene, field, "gene_result_mismatch")
        _assembly(item["variant_key"], assembly, field)
        _validate_document(item["document"], gene, field)
    for document in value["gene_context"]:
        _validate_document(document, gene, field)
    for item in value["citation_seeds"]:
        _keys(
            item,
            {"gene", "pmid", "variant_key", "discovered_from_source_id", "route", "source_fields"},
            field,
        )
        require(item["gene"] == gene and type(item["source_fields"]) is tuple, field, "citation_seed_invalid")
        _assembly(item["variant_key"], assembly, field)
    for item in value["lookups"]:
        _keys(item, {"provider", "gene", "variant_key", "status", "error_code"}, field)
        require(item["gene"] == gene, field, "gene_result_mismatch")
        if item["variant_key"] is not None:
            _assembly(item["variant_key"], assembly, field)


def _validate_evidence(gene, value):
    field = "assessment_seed.evidence"
    _keys(value, {"gene", "documents", "omim_context", "retrieval_errors"}, field)
    require(value["gene"] == gene, field, "gene_result_mismatch")
    for name in ("documents", "omim_context", "retrieval_errors"):
        require(type(value[name]) is tuple, f"{field}.{name}", "tuple_required")
    for document in (*value["documents"], *value["omim_context"]):
        _validate_document(document, gene, field)


def _validate_input_origin(public, context):
    raw_input = public["input"]
    rows, clinical = raw_input.get("phenotypes"), raw_input.get("clinical_description", "")
    require(type(rows) is list and type(clinical) is str, "result.provenance", "patient_input_mismatch")
    terms = []
    for row in rows:
        require(
            type(row) is dict and row.get("status", "present") in {"present", "absent"},
            "result.provenance",
            "patient_input_mismatch",
        )
        terms.append(
            {
                "hpo_id": row["hpo_id"],
                "label": row.get("label", ""),
                "present": row.get("status", "present") == "present",
            }
        )
    expected = json.dumps(
        {"phenotypes": terms, "clinical_description": clinical},
        ensure_ascii=False,
        separators=(",", ":"),
    )
    observed = tuple(
        source
        for source in context["base"]["sources"]
        if source.provider == "patient"
        and source.kind == "patient_observation"
        and source.locator == "input:patient/observations"
        and source.version == "canonical-case-v1"
        and not source.genes
        and not source.variant_keys
    )
    require(
        len(observed) == 1 and observed[0].content == expected, "result.provenance", "patient_input_mismatch"
    )
    seed = context["assessment_seed"]
    require(
        (seed is not None) == (public["input_mode"] == "hpo_vcf"),
        "result.provenance",
        "input_branch_mismatch",
    )
    if seed is None:
        require(
            "genomic_input" not in raw_input
            and public["input_mode"] == ("hpo_clinical" if clinical.strip() else "hpo"),
            "result.provenance",
            "input_branch_mismatch",
        )
        return
    genomic, family, profile = raw_input.get("genomic_input"), seed["family"], seed["profile"]
    require(
        type(genomic) is dict and {"family_mode", "sex", "assay_type"} <= set(genomic),
        "result.provenance",
        "genomic_input_mismatch",
    )
    require(
        genomic.get("assembly") in {None, family["assembly"]}
        and genomic["family_mode"] in {"auto", family["mode"]}
        and genomic["sex"] in {"UNKNOWN", family["sex"]},
        "result.provenance",
        "genomic_input_mismatch",
    )
    receipt = json.dumps(
        {field: genomic.get(field) for field in ("assembly", "family_mode", "sex", "assay_type")},
        ensure_ascii=False,
        sort_keys=True,
    )
    require(
        any(
            source.provider == "patient"
            and source.kind == "patient_observation"
            and source.locator == "input:patient/genomic-settings"
            and source.version == "submitted-genomic-settings-v1"
            and source.content == receipt
            and not source.genes
            and not source.variant_keys
            for source in context["base"]["sources"]
        ),
        "result.provenance",
        "genomic_input_mismatch",
    )
    require(
        profile["clinical_description"] == clinical
        and tuple((term["hpo_id"], term["present"]) for term in profile["terms"])
        == tuple((term["hpo_id"], term["present"]) for term in terms),
        "result.provenance",
        "patient_input_mismatch",
    )


def assessment_origin(seed):
    encoded = json.dumps(
        _encode(seed), ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return {
        "case_key": case_record(seed["profile"]["case_key"]),
        "seed_sha256": hashlib.sha256(encoded).hexdigest(),
    }


def validate_assessment_origin(seed, sources, identities):
    recorded = {source.source_id: source for source in sources}
    require(
        all(recorded.get(source.source_id) == source for source in seed["sources"]),
        "report.assessment",
        "current_run_source_changed_or_missing",
    )
    observed = {
        (row["gene"], variant["key"])
        for row in seed["native_genes"]
        for variant in row["variants"]
        if variant["key"] is not None
    }
    require(
        all(identity in observed for identity in identities),
        "report.assessment",
        "assessment_allele_not_in_current_run",
    )


def _profile_citations(profile):
    for section in ("allele_disease",):
        for citation in profile[section]["citations"]:
            yield section, citation
    for section in ("key_support", "key_limits", "retained_observations"):
        for item in profile.get(section, ()):
            for citation in item["citations"]:
                yield f"{section}:{item['id']}", citation
    anchor = profile["validated_pvs1_anchor"]
    if anchor is not None:
        for citation in anchor["citations"]:
            yield "validated_pvs1_anchor", citation


def _verified_exact_paper_bindings(seed, assessment, sources, case_key):
    """Create identity-only views when the saved provenance proves an exact paper/allele link."""
    knowledge_by_gene = seed.get("variant_knowledge", {})
    if not any(value.get("citation_seeds") for value in knowledge_by_gene.values()):
        return ()
    registry = {source.source_id: source for source in sources}
    native = {}
    for row in seed["native_genes"]:
        for variant in row["variants"]:
            key, source = variant.get("key"), variant.get("source")
            if (
                key is not None
                and isinstance(source, Source)
                and row["gene"] == variant.get("gene")
                and variant.get("transcript")
                and variant.get("hgvs_c")
                and source.genes == (row["gene"],)
                and source.variant_keys == (key,)
            ):
                native.setdefault((row["gene"], key), []).append(
                    (variant["transcript"], variant["hgvs_c"], source)
                )
    profiles = {_identity(profile): profile for profile in assessment["profiles"]}
    output = []
    for gene, knowledge in knowledge_by_gene.items():
        for seed_record in knowledge["citation_seeds"]:
            key, pmid = seed_record["variant_key"], seed_record["pmid"]
            if not (
                seed_record["gene"] == gene
                and seed_record["route"] == "clinvar_citation"
                and "citations.pmid" in seed_record["source_fields"]
                and (gene, key) in profiles
            ):
                continue
            clinvar = registry.get(seed_record["discovered_from_source_id"])
            if clinvar is None:
                continue
            for transcript, hgvs, identity_source in native.get((gene, key), ()):
                for section, citation in _profile_citations(profiles[(gene, key)]):
                    paper = registry.get(citation["source_id"])
                    if not (
                        paper is not None
                        and paper.provider == "pubmed"
                        and paper.kind == "literature"
                        and paper.genes == (gene,)
                        and not paper.variant_keys
                        and re.fullmatch(
                            rf"https://pubmed\.ncbi\.nlm\.nih\.gov/{re.escape(pmid)}/?",
                            paper.locator,
                        )
                        and citation["quote"] in paper.content
                    ):
                        continue
                    try:
                        basis = verified_paper_binding_basis(
                            clinvar,
                            paper,
                            identity_source,
                            gene,
                            key,
                            case_key,
                            pmid,
                            transcript,
                            hgvs,
                            citation["quote"],
                        )
                    except MendelPrimeError:
                        continue
                    content = {
                        "record_type": "verified_exact_paper_binding_v1",
                        "scope": "allele_identity_only",
                        "binding_basis": {**basis, "validated_profile_section": section},
                        "does_not_establish": [
                            "phase",
                            "genotype_state",
                            "independent_case_count",
                            "PS4_strength",
                            "pathogenic_classification",
                        ],
                    }
                    receipt = Source(
                        provider="mendelprime",
                        kind="derived",
                        locator=(
                            f"method:verified-exact-paper-binding/{gene}/{variant_text(key)}/PMID:{pmid}"
                        ),
                        version="verified-exact-paper-binding-v1",
                        content=json.dumps(content, ensure_ascii=False, sort_keys=True),
                        genes=(gene,),
                        variant_keys=(key,),
                        case_key=case_key,
                        parent_source_ids=(
                            clinvar.source_id,
                            paper.source_id,
                            identity_source.source_id,
                        ),
                    )
                    scoped = Source(
                        provider="pubmed",
                        kind="literature",
                        locator=(
                            paper.locator.rstrip("/")
                            + f"/#verified-allele-{receipt.source_id.removeprefix('SRC_')[:16]}"
                        ),
                        version=(
                            f"{paper.version};scope=verified-exact-paper-binding-v1;"
                            f"receipt={receipt.source_id}"
                        ),
                        content=paper.content,
                        genes=(gene,),
                        variant_keys=(key,),
                        case_key=case_key,
                        parent_source_ids=(
                            receipt.source_id,
                            clinvar.source_id,
                            paper.source_id,
                            identity_source.source_id,
                        ),
                    )
                    validate_verified_paper_view(
                        scoped,
                        {
                            **registry,
                            receipt.source_id: receipt,
                            scoped.source_id: scoped,
                        },
                    )
                    output.extend((receipt, scoped))
                    break
                else:
                    continue
                break
    return tuple(dict.fromkeys(output))


def assessment_availability(progress):
    """Project coverage from validated progress without assigning missing grades."""
    completed = {_identity(row) for row in progress["final"]}
    return {
        "expected_alleles": len(progress["expected_identities"]),
        "classified_alleles": len(completed),
        "unclassified_alleles": tuple(
            {"gene": row["gene"], "variant_key": variant_text(row["variant_key"])}
            for row in progress["expected_identities"]
            if (row["gene"], row["variant_key"]) not in completed
        ),
    }


def attach_assessment(context, assessment=None, *, progress=None):
    validate_context(context)
    seed, base = context["assessment_seed"], context["base"]
    availability = {}
    if progress is not None:
        from mendelprime.tools.acmg import completed_assessment

        require(seed is not None and assessment is None, "report.assessment", "partial_without_full_required")
        assessment = completed_assessment(progress)
        validate_assessment_origin(
            seed,
            progress["sources"],
            tuple((row["gene"], row["variant_key"]) for row in progress["expected_identities"]),
        )
        availability = {"assessment_availability": assessment_availability(progress)}
        unresolved_locators = {
            f"method:post-ranking-acmg/{row['gene']}/{row['variant_key']}"
            for row in availability["assessment_availability"]["unclassified_alleles"]
        }
        retained = tuple(
            source
            for source in base["sources"]
            if not (
                source.provider == "mendelprime"
                and source.kind == "derived"
                and source.locator in unresolved_locators
            )
        )
        retained_ids = {source.source_id for source in retained}
        base = {
            **base,
            "sources": retained,
            "available_source_ids": tuple(
                source_id for source_id in base["available_source_ids"] if source_id in retained_ids
            ),
        }
    require(
        (seed is not None) == (assessment is not None or progress is not None),
        "report.assessment",
        "assessment_required_only_for_vcf",
    )
    if assessment is None:
        if progress is None:
            return base
        index = OTEM(base["case_key"])
        index.add(base["sources"])
        index.add(progress["sources"])
        return {**base, "sources": index.sources, **availability}
    clinical_profile = assessment.get("workflow") == "clinical_profile_v1"
    if clinical_profile:
        from mendelprime.tools.acmg import validate_result

        validate_result(assessment)
    identities = tuple(_identity(row) for row in assessment["final"])
    validate_assessment_origin(seed, assessment["sources"], identities)
    identities = tuple(_identity(row) for row in assessment["final"])
    derived_locators = {f"method:post-ranking-acmg/{gene}/{variant_text(key)}" for gene, key in identities}
    index = OTEM(base["case_key"])
    index.add(
        tuple(
            source
            for source in base["sources"]
            if not (source.kind == "derived" and source.locator in derived_locators)
        )
    )
    index.add(assessment["sources"])
    available = set(base["available_source_ids"])
    paper_bindings = (
        _verified_exact_paper_bindings(seed, assessment, index.sources, base["case_key"])
        if clinical_profile
        else ()
    )

    def expose(source):
        index.add((source,))
        available.add(source.source_id)

    if clinical_profile:
        rows = zip(assessment["profiles"], assessment["final"], strict=True)
    else:
        rows = zip(assessment["independent"], assessment["partners"], assessment["final"], strict=True)
    for row in rows:
        if clinical_profile:
            profile, final = row
            require(_identity(profile) == _identity(final), "report.assessment", "profile_final_mismatch")
            criteria = ()
        else:
            review, partner, final = row
            criteria = tuple(
                item
                for item in (*review["criteria"], partner["criterion"])
                if item["status"] in {"met", "uncertain"}
            )
        gene, key = _identity(final)
        originals = tuple(
            source
            for source in assessment["sources"]
            if source.genes == (gene,) and key in source.variant_keys and source.kind != "derived"
        )
        if clinical_profile:
            citations = (
                *profile["allele_disease"]["citations"],
                *profile["patient_fit"]["citations"],
                *(citation for item in profile["key_support"] for citation in item["citations"]),
                *(citation for item in profile["key_limits"] for citation in item["citations"]),
                *(
                    citation
                    for item in profile.get("retained_observations", ())
                    for citation in item["citations"]
                ),
                *(
                    profile["validated_pvs1_anchor"]["citations"]
                    if profile["validated_pvs1_anchor"] is not None
                    else ()
                ),
            )
            parents = tuple(
                dict.fromkeys(
                    (
                        *(source.source_id for source in originals),
                        *profile["input_source_ids"],
                        *(citation["source_id"] for citation in citations),
                        *final["decisive_source_ids"],
                        *final["contextual_source_ids"],
                    )
                )
            )
        else:
            parents = tuple(
                dict.fromkeys(
                    (
                        *(source.source_id for source in originals),
                        *final["decisive_source_ids"],
                        *final["contextual_source_ids"],
                        *(
                            citation["source_id"]
                            for criterion in criteria
                            for citation in criterion["citations"]
                        ),
                    )
                )
            )
        require(bool(parents), "report.assessment", "allele_originals_required")
        for source_id in parents:
            expose(index.get_source(source_id))
        content = (
            {
                "interpretation": (
                    "Post-ranking clinical evidence profile and five-tier assessment; "
                    "did not cause the frozen gene ranking"
                ),
                "workflow": "clinical_profile_v1",
                "profile": _encode(profile, omit_source_ids=True),
                "decision": _encode(final, omit_source_ids=True),
            }
            if clinical_profile
            else {
                "interpretation": "Post-ranking model assessment; did not cause the frozen gene ranking",
                "decision": _encode(final, omit_source_ids=True),
                "reviewed_criteria": [_encode(item, omit_source_ids=True) for item in criteria],
            }
        )
        serialized = json.dumps(content, ensure_ascii=False, sort_keys=True)
        expose(
            Source(
                provider="mendelprime",
                kind="derived",
                locator=f"method:post-ranking-acmg/{gene}/{variant_text(key)}",
                version="post-ranking-clinical-profile-v1;sha256="
                + hashlib.sha256(serialized.encode()).hexdigest()
                if clinical_profile
                else "post-ranking-acmg-v1",
                genes=(gene,),
                variant_keys=(key,),
                case_key=base["case_key"],
                parent_source_ids=parents,
                content=serialized,
            )
        )
    for source in paper_bindings:
        expose(source)
    report_assessment = (
        assessment
        if clinical_profile
        else {name: assessment[name] for name in ("independent", "partners", "final")}
    )
    return {
        "case_key": base["case_key"],
        "fixed_genes": base["fixed_genes"],
        "sources": index.sources,
        "available_source_ids": tuple(
            source.source_id for source in index.sources if source.source_id in available
        ),
        "assessment": report_assessment,
        **availability,
    }


def explain(model, context, readers=None, saved_assessment=None, saved_progress=None, *, acmg_workers=1):
    """Run only post-ranking ACMG/report work over a frozen context."""
    validate_context(context)
    frozen = context["base"]["fixed_genes"]
    assessment = saved_assessment
    progress = None
    issues = ()
    try:
        seed = context["assessment_seed"]
        if assessment is not None:
            require(seed is not None, "explanation", "assessment_only_for_vcf")
            if assessment.get("workflow") == "clinical_profile_v1":
                from mendelprime.tools.acmg import validate_result

                validate_result(assessment)
            validate_assessment_origin(
                seed, assessment["sources"], tuple(_identity(row) for row in assessment["final"])
            )
        elif seed is not None:
            require(readers is not None, "explanation", "vcf_explanation_resources_required")
            from mendelprime.harness import assessment as assessment_workflow
            from mendelprime.tools import acmg

            actions = readers.get("actions", readers)
            prepared = assessment_workflow.build(
                seed, actions, **({"review": context["review"]} if "review" in context else {})
            )
            require(
                prepared["input"]["profile"] == seed["profile"]
                and prepared["input"]["family"] == seed["family"],
                "explanation",
                "assessment_patient_context_mismatch",
            )
            validate_assessment_origin(
                seed,
                prepared["input"]["sources"],
                tuple(
                    (row["primary"]["gene"], row["primary"]["key"]) for row in prepared["input"]["targets"]
                ),
            )
            issues = prepared["issues"]
            try:
                assessment = acmg.run(
                    model,
                    prepared["input"],
                    workers=acmg_workers,
                    **({"progress": saved_progress} if saved_progress is not None else {}),
                )
            except MendelPrimeError as error:
                progress = getattr(error, "partial_result", None)
                if progress is None:
                    raise
                acmg.validate_progress(progress)
                issues = (*issues, *(row["code"] for row in progress["issues"]))
        report_context = attach_assessment(context, assessment, progress=progress)
        result = report_workflow.run(model, report_context)
        combined_issues = (*issues, *result["issues"])
        result = {
            **result,
            "status": "degraded"
            if any(not is_evidence_limitation(issue) for issue in combined_issues)
            else result["status"],
            "assessment": assessment,
            "sources": report_context["sources"],
            "issues": combined_issues,
            "evidence_limitations": tuple(i for i in combined_issues if is_evidence_limitation(i)),
            "execution_issues": tuple(i for i in combined_issues if not is_evidence_limitation(i)),
            **(
                {
                    "assessment_progress": progress,
                    "status": "degraded",
                    "assessment_availability": report_context["assessment_availability"],
                }
                if progress is not None
                else {}
            ),
        }
        require(context["base"]["fixed_genes"] == frozen, "explanation", "frozen_ranking_changed")
        return result
    except (MendelPrimeError, OSError) as error:
        require(context["base"]["fixed_genes"] == frozen, "explanation", "frozen_ranking_changed")
        progress = (
            getattr(error, "partial_result", progress or saved_progress) if assessment is None else None
        )
        issues = (*issues, "explanation_resource_io_error" if isinstance(error, OSError) else error.code)
        return {
            "status": "degraded",
            "assessment": assessment,
            "report": None,
            "sources": (
                progress["sources"]
                if progress is not None
                else assessment["sources"]
                if assessment is not None
                else context["base"]["sources"]
            ),
            **({"assessment_progress": progress} if progress is not None else {}),
            "markdown": "",
            "issues": issues,
            "evidence_limitations": tuple(i for i in issues if is_evidence_limitation(i)),
            "execution_issues": tuple(i for i in issues if not is_evidence_limitation(i)),
        }


def load_saved(public, context):
    """Read saved assessments and report text without rerendering or re-auditing history."""
    from mendelprime.runtime.results import public_value, validate_public_result

    validate_public_result(public, require_provenance=True)
    validate_context(context)
    stage = public["stages"].get("final_explanation")
    require(
        type(stage) is dict
        and type(stage.get("status")) is str
        and stage["status"] in {"disabled", "ok", "degraded"},
        "result.explanation",
        "valid_status_required",
    )
    genotype = public["stages"].get("genotype")
    raw_assessment = genotype.get("final_acmg") if type(genotype) is dict else None
    raw_progress = genotype.get("assessment_progress") if type(genotype) is dict else None
    if stage["status"] == "disabled":
        require(
            raw_assessment is None and raw_progress is None and set(stage) == {"status"},
            "result.explanation",
            "disabled_must_be_empty",
        )
        return {
            "status": "disabled",
            "assessment": None,
            "report": None,
            "sources": context["base"]["sources"],
            "markdown": "",
            "issues": (),
        }
    require(
        type(stage.get("issues", [])) is list
        and all(type(issue) is str and issue for issue in stage.get("issues", [])),
        "result.explanation.issues",
        "error_code_array_required",
    )
    assessment = _load_assessment(raw_assessment, context) if raw_assessment is not None else None
    issue_views = {}
    for name in ("evidence_limitations", "execution_issues"):
        if name in stage:
            require(
                type(stage[name]) is list and all(type(i) is str and i for i in stage[name]),
                "result.explanation." + name,
                "error_code_array_required",
            )
            issue_views[name] = tuple(stage[name])
    progress = None
    if raw_progress is not None:
        from mendelprime.tools.acmg import validate_progress

        require(
            assessment is None and stage["status"] == "degraded" and context["assessment_seed"] is not None,
            "result.assessment_progress",
            "partial_assessment_requires_degraded_vcf",
        )
        registry = {source.source_id: source for source in context["base"]["sources"]}
        progress = _decode(raw_progress, registry)
        validate_progress(progress)
        validate_assessment_origin(
            context["assessment_seed"],
            progress["sources"],
            tuple((row["gene"], row["variant_key"]) for row in progress["expected_identities"]),
        )
    structured = stage.get("structured_explanation", {})
    require(type(structured) is dict, "result.explanation", "structured_report_required")
    if "assessment_availability" in stage or progress is not None and "report" in structured:
        require(
            progress is not None
            and stage.get("assessment_availability") == public_value(assessment_availability(progress)),
            "result.explanation.assessment_availability",
            "classification_coverage_mismatch",
        )
    if "report" not in structured:
        require(
            stage["status"] == "degraded"
            and not stage.get("report_markdown_en")
            and "semantic_audit" not in structured,
            "result.explanation",
            "structured_report_required",
        )
        return {
            "status": "degraded",
            "assessment": assessment,
            "report": None,
            "sources": context["base"]["sources"],
            "markdown": "",
            "issues": tuple(stage.get("issues", ())),
            **issue_views,
            **({"assessment_progress": progress} if progress is not None else {}),
        }
    require(
        context["assessment_seed"] is None or assessment is not None or progress is not None,
        "result.explanation",
        "vcf_assessment_required_for_report",
    )
    report_context = attach_assessment(context, assessment, progress=progress)
    report = _load_report(structured["report"], report_context["case_key"])
    markdown = stage.get("report_markdown_en", "")
    require(type(markdown) is str, "result.explanation", "saved_markdown_string_required")
    return {
        "status": stage["status"],
        "assessment": assessment,
        "report": report,
        **(
            {
                "assessment_progress": progress,
                "assessment_availability": report_context["assessment_availability"],
            }
            if progress is not None
            else {}
        ),
        **({"historical_audit": structured["semantic_audit"]} if "semantic_audit" in structured else {}),
        **({"review_status": stage["review_status"]} if "review_status" in stage else {}),
        "sources": report_context["sources"],
        "markdown": markdown,
        "issues": tuple(stage.get("issues", ())),
        **issue_views,
    }


def _load_assessment(value, context):
    if type(value) is dict and value.get("workflow") == "clinical_profile_v1":
        return _load_clinical_profile_assessment(value, context)
    return _load_legacy_assessment(value, context)


def _load_legacy_assessment(value, context):
    require(
        type(value) is dict
        and set(value) == {"status", "origin", "by_gene", "independent", "partners", "final"},
        "result.assessment",
        "schema_mismatch",
    )
    seed = context["assessment_seed"]
    require(
        seed is not None and value["origin"] == assessment_origin(seed),
        "result.assessment.origin",
        "assessment_origin_mismatch",
    )
    assessment = {
        name: tuple(_decode_public_record(row) for row in value[name])
        for name in ("independent", "partners", "final")
    }
    seed_ids = {source.source_id for source in seed["sources"]}
    genes = set(context["base"]["fixed_genes"])
    assessment["sources"] = tuple(
        source
        for source in context["base"]["sources"]
        if source.source_id in seed_ids
        or source.kind == "guideline"
        or (source.kind != "derived" and bool(set(source.genes).intersection(genes)))
    )
    identities = tuple(_identity(row) for row in assessment["final"])
    validate_assessment_origin(seed, assessment["sources"], identities)
    expected = {}
    for decision in value["final"]:
        gene = decision["identity"]["gene"]
        expected.setdefault(gene, {"variants": []})["variants"].append(
            {
                **decision,
                "disease_id": decision["disease"]["identifier"],
                "disease_label": decision["disease"]["label"],
                "mode_of_inheritance": decision["disease"]["inheritance"],
            }
        )
    require(value["by_gene"] == expected, "result.assessment.by_gene", "projection_mismatch")
    return assessment


def _load_clinical_profile_assessment(value, context):
    require(
        type(value) is dict
        and set(value) == {"status", "origin", "by_gene", "workflow", "profiles", "final", "sources"},
        "result.assessment",
        "schema_mismatch",
    )
    require(value["status"] == "ok", "result.assessment.status", "unsupported_status")
    require(
        type(value["profiles"]) is list and type(value["final"]) is list,
        "result.assessment",
        "profile_and_final_arrays_required",
    )
    seed = context["assessment_seed"]
    require(
        seed is not None and value["origin"] == assessment_origin(seed),
        "result.assessment.origin",
        "assessment_origin_mismatch",
    )
    source_ids = value["sources"]
    require(
        type(source_ids) is list
        and all(type(source_id) is str and source_id for source_id in source_ids)
        and len(set(source_ids)) == len(source_ids),
        "result.assessment.sources",
        "unique_source_id_array_required",
    )
    registry = {source.source_id: source for source in context["base"]["sources"]}
    require(
        all(source_id in registry for source_id in source_ids),
        "result.assessment.sources",
        "unknown_source",
    )
    assessment = {
        "workflow": "clinical_profile_v1",
        "profiles": tuple(_decode_public_record(row) for row in value["profiles"]),
        "final": tuple(_decode_public_record(row) for row in value["final"]),
        "sources": tuple(registry[source_id] for source_id in source_ids),
    }
    from mendelprime.tools.acmg import validate_result

    validate_result(assessment)
    validate_assessment_origin(
        seed,
        assessment["sources"],
        tuple(_identity(row) for row in assessment["final"]),
    )
    expected = {}
    for decision in value["final"]:
        gene = decision["identity"]["gene"]
        disease = decision["disease"]
        expected.setdefault(gene, {"variants": []})["variants"].append(
            {
                **decision,
                "disease_id": disease["identifier"],
                "disease_label": disease["label"],
                "mode_of_inheritance": disease["inheritance"],
            }
        )
    require(value["by_gene"] == expected, "result.assessment.by_gene", "projection_mismatch")
    return assessment


def _load_report(value, case_key):
    require(
        type(value) is dict
        and set(value) == {"case_key", "fixed_genes", "claims", "summary", "gene_sections", "comparisons"},
        "result.explanation.report",
        "schema_mismatch",
    )
    report = _decode_public_record(value)
    require(report["case_key"] == case_key, "result.explanation.report", "case_mismatch")
    for claim in report["claims"]:
        claim["case_key"] = case_key
    return report


def _decode_public_record(value):
    if type(value) is list:
        return tuple(_decode_public_record(item) for item in value)
    if type(value) is dict:
        if set(value) == {"cohort", "case_id"}:
            return value["cohort"], value["case_id"]
        return {key: _decode_public_record(item) for key, item in value.items() if key != "source_ids"}
    if isinstance(value, str):
        parts = value.split(":")
        if len(parts) == 5 and parts[0] in {"GRCh37", "GRCh38"}:
            try:
                return variant_key(parts[0], parts[1], int(parts[2]), parts[3], parts[4])
            except (TypeError, ValueError):
                pass
    return value


def _identity(record):
    row = record.get("identity") if type(record) is dict else None
    require(
        type(row) is dict and set(row) == {"gene", "variant_key"}, "assessment", "assessment_identity_invalid"
    )
    validate_gene(row["gene"], "assessment.identity.gene")
    key = row["variant_key"]
    if isinstance(key, str):
        parts = key.split(":")
        require(len(parts) == 5, "assessment", "variant_key_invalid")
        key = variant_key(parts[0], parts[1], int(parts[2]), parts[3], parts[4])
    else:
        key = variant_key(*key)
    return row["gene"], key


def _encode(value, *, field="", omit_source_ids=False):
    if value is None or type(value) in {str, bool, int}:
        return value
    if type(value) is float:
        return value if math.isfinite(value) else None
    if isinstance(value, Source):
        return {"source_id": value.source_id}
    if isinstance(value, (tuple, list)):
        if len(value) == 5 and isinstance(value[0], str) and value[0] in {"GRCh37", "GRCh38"}:
            return variant_record(value)
        if field == "case_key" and len(value) == 2:
            return case_record(value)
        return [_encode(item, omit_source_ids=omit_source_ids) for item in value]
    if isinstance(value, Mapping):
        return {
            key: _encode(item, field=key, omit_source_ids=omit_source_ids)
            for key, item in value.items()
            if not (omit_source_ids and key == "source_ids")
        }
    raise MendelPrimeError("unsupported_public_value", field=field or "provenance")


def _decode(value, registry, *, field=""):
    if type(value) is list:
        return tuple(_decode(item, registry) for item in value)
    if type(value) is dict:
        if set(value) == {"source_id"}:
            require(value["source_id"] in registry, field or "provenance", "unknown_source")
            return registry[value["source_id"]]
        if set(value) == _VARIANT_FIELDS:
            return variant_key(**value)
        if set(value) == {"cohort", "case_id"} and field == "case_key":
            return value["cohort"], value["case_id"]
        return {key: _decode(item, registry, field=key) for key, item in value.items()}
    return value


__all__ = [
    "assessment_origin",
    "attach_assessment",
    "encode_context",
    "explain",
    "freeze_context",
    "load_context",
    "load_saved",
    "validate_assessment_origin",
    "validate_context",
]
