"""Case-local evidence registry; retrieval and stage views keep source identity."""

import hashlib
import json
import re
from collections import defaultdict

from mendelprime import prompting
from mendelprime.domain import (
    Evidence,
    MendelPrimeError,
    Source,
    require,
    validate_case_key,
    validate_gene,
    validate_verified_paper_view,
    variant_key,
    variant_text,
)
from mendelprime.knowledge.hpo import ancestors, consulted_terms

_CAUSAL_ASSERTION = re.compile(
    r"\b(pathogenic|benign|causal|causative|causes?|disease[- ]causing|diagnos(?:is|tic)|etiolog(?:y|ic)|"
    r"explains? the (?:patient|phenotype)|responsible for|proves?|confirms? (?:the )?(?:cause|diagnosis)|"
    r"establish(?:es|ed) (?:pathogenicity|benignity|diagnosis|causality))\b",
    re.IGNORECASE,
)


def _patient_observation_context(claim, cited):
    """Qualify observations by their typed scope and records, not prose keywords."""
    if (
        claim["scope"] != "allele"
        or claim["polarity"] not in {"background", "limitation"}
        or _CAUSAL_ASSERTION.search(claim["text"])
    ):
        return False
    genes, variants = set(claim["genes"]), set(claim["variant_keys"])
    observed = set()
    for source in cited:
        if source.kind != "patient_observation" or source.case_key != claim["case_key"]:
            return False
        if not source.genes and not source.variant_keys:
            continue
        if set(source.genes) != genes:
            return False
        observed.update(source.variant_keys)
        try:
            record = json.loads(source.content)
        except (json.JSONDecodeError, TypeError):
            return False
        if not isinstance(record, dict) or not isinstance(record.get("samples"), list):
            return False
    return bool(cited) and variants <= observed


def validate_claim_shape(claim, case_key):
    """Validate an immutable claim's intrinsic shape before binding it to Sources."""
    require(claim["case_key"] == case_key, "claim.case_key", "case_mismatch")
    genes, variants, scope = claim["genes"], claim["variant_keys"], claim["scope"]
    require(
        isinstance(claim["claim_id"], str) and claim["claim_id"].strip(),
        "claim_id",
        "nonempty_text_required",
    )
    require(
        isinstance(claim["text"], str) and claim["text"].strip(),
        "text",
        "nonempty_text_required",
    )
    require(len(set(genes)) == len(genes), "genes", "duplicate_claim_gene")
    require(len(set(variants)) == len(variants), "variant_keys", "duplicate_variant")
    for gene in genes:
        validate_gene(gene)
    for key in variants:
        variant_key(*key)
    require(scope in {"case", "gene", "comparison", "allele"}, "scope", "claim_scope_enum_required")
    require(
        claim["polarity"] in {"support", "limitation", "background"},
        "polarity",
        "claim_polarity_enum_required",
    )
    require(bool(claim["citations"]), "citations", "nonempty_tuple_required")
    require(
        len({(citation["source_id"], citation["quote"]) for citation in claim["citations"]})
        == len(claim["citations"]),
        "citations",
        "duplicate_citation",
    )
    if scope == "case":
        require(not genes and not variants, "scope", "case_scope_must_be_global")
    elif scope == "gene":
        require(len(genes) == 1 and not variants, "scope", "gene_scope_requires_one_gene")
    elif scope == "comparison":
        require(len(genes) >= 2, "scope", "comparison_scope_requires_two_genes")
        require(not variants, "scope", "comparison_scope_forbids_variants")
    else:
        require(len(genes) == 1 and variants, "scope", "allele_scope_requires_gene_and_variant")


def bound_alleles(sources, genes):
    """Original, single-gene allele records; derived summaries are not new identity proof."""
    return {
        key
        for source in sources
        if source.kind != "derived" and len(source.genes) == 1 and source.genes[0] in genes
        for key in source.variant_keys
    }


class OTEM:
    def __init__(self, case_key):
        validate_case_key(case_key)
        self.case_key = case_key
        self._sources, self._logical, self._evidence = {}, {}, {}

    @property
    def sources(self):
        return tuple(self._sources.values())

    @property
    def evidence(self):
        return tuple(self._evidence.values())

    def add(self, sources=(), evidence=()):
        for source in sources:
            require(isinstance(source, Source), "source", "source_required")
            require(source.case_key in {None, self.case_key}, "source.case_key", "case_mismatch")
            for parent_id in source.parent_source_ids:
                require(parent_id in self._sources, "parent_source_ids", "parent_not_registered")
                parent = self._sources[parent_id]
                require(
                    parent.case_key in {None, source.case_key}, "parent_source_ids", "parent_case_mismatch"
                )
            logical = (
                source.provider,
                source.kind,
                source.locator,
                source.version,
                source.genes,
                source.variant_keys,
                source.case_key,
                source.parent_source_ids,
            )
            require(
                self._logical.get(logical, source.source_id) == source.source_id,
                "source",
                "logical_record_content_changed",
            )
            require(self._sources.get(source.source_id, source) == source, "source_id", "content_collision")
            if source.kind == "literature" and source.parent_source_ids:
                validate_verified_paper_view(source, {**self._sources, source.source_id: source})
            self._sources[source.source_id] = source
            self._logical[logical] = source.source_id
        for value in evidence:
            require(isinstance(value, Evidence), "evidence", "evidence_required")
            require(value.case_key == self.case_key, "evidence.case_key", "case_mismatch")
            require(
                all(key in self._sources for key in value.source_ids), "evidence.source_ids", "unknown_source"
            )
            require(
                all(key in self._evidence for key in value.derived_from),
                "evidence.derived_from",
                "parent_not_registered",
            )
            require(self._evidence.get(value.evidence_id, value) == value, "evidence_id", "content_collision")
            self._evidence[value.evidence_id] = value

    def get_source(self, source_id):
        return self._sources[source_id]

    def closure(self, source_ids):
        selected, pending = set(), list(source_ids)
        while pending:
            source_id = pending.pop()
            if source_id not in selected:
                require(source_id in self._sources, "source_ids", "unknown_source")
                selected.add(source_id)
                pending.extend(self._sources[source_id].parent_source_ids)
        return tuple(source for key, source in self._sources.items() if key in selected)

    def register_molecular(self, native, genes, family, knowledge):
        from mendelprime.cemr.views import candidate_cards, document_view, projection_source_ids
        from mendelprime.knowledge.literature import evidence_document
        from mendelprime.knowledge.variants import validate_result

        require(self.case_key == native["case_key"], "molecular_context", "case_mismatch")
        require(
            all(tool["assembly"] == family["assembly"] for tool in native["tools"]),
            "molecular_context",
            "assembly_mismatch",
        )
        roles = {role: sample for role, sample in family["roles"].items() if sample}
        observations = native["observations"]["observations"]
        require(
            all(
                {call["role"]: call["sample_id"] for call in row["samples"]} == roles for row in observations
            ),
            "molecular_context",
            "family_role_mismatch",
        )
        family_source = Source(
            provider="patient",
            kind="patient_observation",
            locator="input:patient/genomic-context",
            version="family-resolution-v1",
            case_key=self.case_key,
            content=json.dumps(family, ensure_ascii=False, sort_keys=True),
        )
        self.add((family_source,))
        observed = {variant_text(row["vcf_variant_key"]): row for row in observations}
        annotations = {}
        for gene in genes:
            result = knowledge.get(gene)
            require(
                isinstance(result, dict) and result.get("gene") == gene,
                "molecular_context",
                "variant_knowledge_gene_mismatch",
            )
            validate_result(result)
            annotations[gene] = {
                "gene": gene,
                "documents": tuple(
                    row["document"] for row in result["variant_evidence"] if not row["retracted"]
                )
                + result["gene_context"],
                "omim_context": (),
                "retrieval_errors": tuple(
                    f"{row['provider']}:{row['status']}:{row['error_code']}"
                    for row in result["lookups"]
                    if row["status"] not in {"found", "not_found", "not_applicable"}
                ),
            }
        output = {}
        for card in candidate_cards(native, genes, annotations, include_unprioritized=True):
            gene = card["gene"]
            card["case_context"] = {
                "source_id": family_source.source_id,
                "assembly": family["assembly"],
                "family_mode": family["mode"],
                "proband_sex": family["sex"],
            }
            for number, allele in enumerate(card["alleles"], 1):
                allele["variant_ref"] = f"VAR:{gene}:{number}"
                observation = observed.get(allele["vcf_allele_key"])
                if observation is None or allele["family_observation"] is None:
                    continue
                for call in observation["samples"]:
                    raw, gt = dict(call["raw_fields"]), call["gt"]
                    allele["family_observation"]["observations"][call["role"]].update(
                        genotype_ploidy=len(call["alleles"]) if gt is not None else None,
                        genotype_phased=bool(gt and "|" in gt and "/" not in gt),
                        phase_set=call["phase_set"],
                        sample_filter=call["filter_text"],
                        physical_genotype=None if raw.get("PGT") in {None, "."} else raw["PGT"],
                        physical_phase_id=None if raw.get("PID") in {None, "."} else raw["PID"],
                    )
            card["variant_assertions"], card["withdrawn_assertions"] = [], []
            card["assertion_semantics"] = {
                "variant_assertions": "Quoted database assertions, not a new MendelPRIME ACMG classification.",
                "withdrawn_assertions": "Historical records only; not usable as supporting assertions.",
                "linked_curations": "Records sharing a dependency group are not independent patient or experiment observations.",
            }
            for row in knowledge[gene]["variant_evidence"]:
                source_id, linked = row["document"]["source"].source_id, row["linked_source_id"]
                assertion = {
                    "variant_key": variant_text(row["variant_key"]),
                    "match_basis": row["match_basis"],
                    "assertion_scope": row["assertion_scope"],
                    "source_id": source_id,
                    "dependency_group_id": linked or source_id,
                }
                if linked:
                    assertion["linked_source_id"] = linked
                assertion["usable_as_assertion"] = not row["retracted"]
                if row["retracted"]:
                    assertion["record"] = document_view(row["document"])
                card["withdrawn_assertions" if row["retracted"] else "variant_assertions"].append(assertion)
            card["knowledge_lookups"] = [
                {
                    "provider": row["provider"],
                    "status": row["status"],
                    "error_code": row["error_code"],
                    "variant_key": variant_text(row["variant_key"]) if row["variant_key"] else None,
                }
                for row in knowledge[gene]["lookups"]
            ]
            parents = projection_source_ids(card)
            for identifier in parents:
                require(identifier in self._sources, "molecular_context", "parent_not_registered")
                parent = self.get_source(identifier)
                require(not parent.genes or gene in parent.genes, "molecular_context", "parent_gene_mismatch")
            content = json.dumps(card, ensure_ascii=False, separators=(",", ":"))
            source = Source(
                provider="mendelprime",
                kind="derived",
                locator=f"method:molecular-context/{gene}",
                version="molecular-context-v2",
                content=content,
                genes=(gene,),
                case_key=self.case_key,
                parent_source_ids=tuple(sorted(parents)),
            )
            self.add((source,))
            output[gene] = evidence_document(
                source,
                f"VCF:{gene}",
                f"Observed molecular context for {gene}",
                content,
                provenance=("native_tools", "family_observations", "resource_annotations"),
            )
        return output

    def retrieve_similar(self, profile, candidate_genes, retrieve, hpo, model, max_cases=20):
        from mendelprime.knowledge.retrieval import unique_case_hits, unique_upper_gene_groups

        try:
            hits = retrieve(profile)
        except Exception as error:
            return {
                "status": "degraded",
                "error_code": f"case_retrieval_failed:{type(error).__name__}",
                "hits": (),
                "judged_hits": (),
                "decisions": (),
                "evidence_by_gene": {},
                "tail_status": "ok",
                "tail_error_code": "",
                "tail_hits": (),
                "tail_decisions": (),
                "tail_evidence_by_gene": {},
            }
        allowed = set(candidate_genes)
        ordered = sorted(
            (hit for hit in hits if hit["gene"].upper() in allowed),
            key=lambda hit: (
                hit["rerank_rank"],
                hit["gene_case_rank"],
                hit["source"].provider,
                hit["retrieval_rank"],
            ),
        )
        selected = tuple(unique_case_hits(ordered))
        if max_cases > 0:
            selected = selected[:max_cases]
        ordered = sorted(
            hits,
            key=lambda hit: (
                hit["rerank_rank"] or 10_000,
                hit["gene_case_rank"],
                hit["retrieval_rank"],
                hit["source"].provider,
                hit["record_id"],
            ),
        )
        genes = tuple(
            dict.fromkeys(hit["gene"].upper() for hit in ordered if hit["gene"].upper() not in allowed)
        )[:10]
        counts, tail = defaultdict(int), []
        for hit in unique_upper_gene_groups(ordered):
            gene = hit["gene"].upper()
            if gene in genes and counts[gene] < 2:
                tail.append(hit)
                counts[gene] += 1
        top_result = _judge_cases(profile, selected, hpo, model)
        tail_result = _judge_cases(profile, tuple(tail), hpo, model)
        self.add(tuple(hit["source"] for hit in hits))
        return {
            "status": top_result["status"],
            "hits": tuple(hits),
            "judged_hits": selected,
            "decisions": top_result["decisions"],
            "evidence_by_gene": top_result["evidence"],
            "error_code": top_result["error_code"],
            "tail_status": tail_result["status"],
            "tail_hits": tuple(tail),
            "tail_decisions": tail_result["decisions"],
            "tail_evidence_by_gene": tail_result["evidence"],
            "tail_error_code": tail_result["error_code"],
        }

    def register_patient(self, patient):
        require(patient.key == self.case_key, "patient.case_key", "case_mismatch")
        content = json.dumps(
            {
                "phenotypes": [dict(term) for term in patient.phenotypes],
                "clinical_description": patient.clinical_description,
            },
            ensure_ascii=False,
            separators=(",", ":"),
        )
        source = Source(
            provider="patient",
            kind="patient_observation",
            locator="input:patient/observations",
            version="canonical-case-v1",
            content=content,
            case_key=patient.key,
        )
        self.add((source,))
        return source

    def register_profile(self, hpo, profile):
        require(profile["case_key"] == self.case_key, "profile.case_key", "case_mismatch")
        require(
            profile["ontology_manifest"] == hpo["ontology_manifest"]
            and profile["statistics_manifest"] == hpo["statistics_manifest"],
            "profile",
            "ontology_manifest_mismatch",
        )
        require(
            profile["consulted_hpo_ids"]
            == consulted_terms(hpo, tuple(term["hpo_id"] for term in profile["terms"])),
            "profile.consulted_hpo_ids",
            "ontology_terms_mismatch",
        )
        manifest = hpo["ontology_manifest"]
        for identifier in profile["consulted_hpo_ids"]:
            original = hpo["terms"].get(identifier, {}).get("original_record")
            require(original, "profile.consulted_hpo_ids", "original_term_record_missing")
            source = Source(
                provider="hpo",
                kind="ontology",
                locator=f"resource:hpo/terms/{manifest['sha256']}#{identifier}",
                version=manifest["version"],
                content=original,
            )
            self.add((source,))

    def validate_claim(self, claim):
        """Prove exact citation and scope binding, not clinical entailment."""
        validate_claim_shape(claim, self.case_key)
        genes, variants, scope = claim["genes"], claim["variant_keys"], claim["scope"]
        cited = []
        for citation in claim["citations"]:
            require(citation["source_id"] in self._sources, "citation.source_id", "unknown_source")
            source = self._sources[citation["source_id"]]
            require(source.case_key in {None, self.case_key}, "citation.source_id", "case_mismatch")
            require(
                isinstance(citation["quote"], str) and citation["quote"].strip(),
                "quote",
                "nonempty_text_required",
            )
            require(citation["quote"] in source.content, "citation.quote", "quote_not_in_source")
            cited.append(source)
        if scope == "case":
            require(
                not any(s.genes or s.variant_keys for s in cited), "citations", "case_scope_source_forbidden"
            )
            return
        for source in cited:
            require(
                not source.genes or set(genes).intersection(source.genes),
                "citations",
                "source_gene_scope_mismatch",
            )
            require(
                not source.variant_keys or set(variants).intersection(source.variant_keys),
                "citations",
                "source_variant_scope_mismatch",
            )
        observation_context = _patient_observation_context(claim, cited)
        support = cited if claim["polarity"] == "support" else self.closure(tuple(s.source_id for s in cited))
        for gene in genes:
            require(
                observation_context
                or any(gene in s.genes and s.kind not in {"patient_observation", "derived"} for s in support),
                f"genes[{gene}]",
                "non_patient_source_required",
            )
        anchored = bound_alleles(support, genes)
        for variant in variants:
            require(
                variant in anchored,
                f"variant_keys[{':'.join(map(str, variant))}]",
                "allele_source_required",
            )


def _hpo_name(hpo, identifier):
    return hpo["terms"].get(identifier, {}).get("name") or identifier


def _relations(patient, case, hpo):
    output = []
    for left in patient:
        for right in case:
            if left == right:
                relation = "exact"
            elif right in hpo["terms"].get(left, {}).get("parent_ids", ()) or left in hpo["terms"].get(
                right, {}
            ).get("parent_ids", ()):
                relation = "direct_ontology_relation"
            elif left in ancestors(hpo, right):
                relation = "patient_term_is_ancestor_of_case_term"
            elif right in ancestors(hpo, left):
                relation = "case_term_is_ancestor_of_patient_term"
            else:
                continue
            output.append(
                {
                    "patient_hpo_id": left,
                    "patient_label": _hpo_name(hpo, left),
                    "case_hpo_id": right,
                    "case_label": _hpo_name(hpo, right),
                    "relation": relation,
                }
            )
    return output


def judge_prompt(profile, hits, hpo):
    patient = [
        {"hpo_id": term["hpo_id"], "label": term["label"]}
        for term in profile["terms"]
        if term["present"] and term["hpo_id"]
    ]
    value = prompting.template("memory_prompt", "judge")
    value["patient_hpo"] = patient
    value["retrieved_cases"] = [
        {
            "case_index": index,
            "hpo_terms": [{"hpo_id": key, "label": _hpo_name(hpo, key)} for key in hit["hpo_ids"]],
            "explicit_ontology_relations": _relations(
                tuple(x["hpo_id"] for x in patient), hit["hpo_ids"], hpo
            ),
        }
        for index, hit in enumerate(hits, 1)
    ]
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def _matched_ids(values, allowed):
    return (
        tuple(dict.fromkeys(str(value).strip() for value in values if str(value).strip() in allowed))
        if isinstance(values, list)
        else ()
    )


def parse_decisions(value, profile, hits):
    def check(condition, code):
        if not condition:
            raise MendelPrimeError(code, category="response")

    rows = value.get("case_decisions")
    check(isinstance(rows, list), "similar_case_hpo_decisions_must_be_list")
    indexed = dict(enumerate(hits, 1))
    patient_ids = {term["hpo_id"] for term in profile["terms"] if term["present"]}
    output, seen = [], set()
    for row in rows:
        check(
            isinstance(row, dict)
            and isinstance(row.get("case_index"), int)
            and row["case_index"] in indexed
            and row["case_index"] not in seen,
            "similar_case_hpo_case_index_invalid",
        )
        index = row["case_index"]
        hit = indexed[index]
        match = row.get("hpo_match")
        kind = str(row.get("match_type") or "").strip()
        patient_match = _matched_ids(row.get("matched_patient_hpo_ids"), patient_ids)
        case_match = _matched_ids(row.get("matched_case_hpo_ids"), set(hit["hpo_ids"]))
        summary = " ".join(str(row.get("summary") or "").split())[:1200]
        check(
            isinstance(match, bool)
            and kind
            in {
                "exact_or_equivalent",
                "ancestor_descendant",
                "direct_ontology_relation",
                "semantic_related",
                "none",
            }
            and summary,
            "similar_case_hpo_decision_content_invalid",
        )
        check(
            not match or (kind != "none" and patient_match and case_match),
            "similar_case_hpo_positive_match_invalid",
        )
        check(match or kind == "none", "similar_case_hpo_negative_match_invalid")
        output.append(
            {
                "record_id": hit["record_id"],
                "gene": hit["gene"],
                "hpo_match": match,
                "match_type": kind,
                "matched_patient_hpo_ids": patient_match,
                "matched_case_hpo_ids": case_match,
                "summary": summary,
            }
        )
        seen.add(index)
    check(seen == set(indexed), "similar_case_hpo_decisions_incomplete")
    return tuple(output)


def _judge_cases(profile, hits, hpo, model):
    if not hits:
        return {"status": "ok", "error_code": "", "decisions": (), "evidence": {}}
    request = {
        "case_key": profile["case_key"],
        "stage_id": "similar_case_hpo_semantic_judge_v2",
        "system_prompt": prompting.text("similar_case_judge"),
        "user_prompt": judge_prompt(profile, hits, hpo),
        "max_tokens": 3000,
        "timeout_seconds": 180.0,
        "cache_choice": False,
    }
    try:
        decisions = prompting.call_json(
            model, request, lambda value: parse_decisions(value, profile, hits), 2
        )
    except MendelPrimeError as error:
        return {"status": "degraded", "error_code": error.code, "decisions": (), "evidence": {}}
    paired = {(hit["record_id"], hit["gene"]): hit for hit in hits}
    grouped = defaultdict(list)
    ordered = sorted(
        (row for row in decisions if row["hpo_match"]),
        key=lambda row: (
            paired[(row["record_id"], row["gene"])]["rerank_rank"],
            paired[(row["record_id"], row["gene"])]["retrieval_rank"],
            row["record_id"],
        ),
    )
    for decision in ordered:
        hit = paired[(decision["record_id"], decision["gene"])]
        identifier = hashlib.sha256(hit["record_id"].encode("utf-8")).hexdigest()[:20]
        grouped[hit["gene"]].append(
            {
                "evidence_id": f"SC:{hit['source'].provider}:{identifier}",
                "hit": hit,
                "decision": decision,
                "hpo_terms": tuple((key, _hpo_name(hpo, key)) for key in hit["hpo_ids"]),
            }
        )
    return {
        "status": "ok",
        "error_code": "",
        "decisions": decisions,
        "evidence": {gene: tuple(values) for gene, values in grouped.items()},
    }
