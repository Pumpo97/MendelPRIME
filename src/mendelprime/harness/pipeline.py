"""The scientific stage sequence; runtime owns batching and publication."""

import logging
import time

from mendelprime.cemr import views
from mendelprime.domain import (
    CaseResult,
    Evidence,
    MendelPrimeError,
    Ranking,
    is_evidence_limitation,
    require,
)
from mendelprime.tools.phenotype import require_positive_hpo

from .memory import OTEM

LOG = logging.getLogger("mendelprime.pipeline")


class Harness:
    def __init__(self, tools, knowledge, cemr, *, reflection=False, similar=True, ranking_readers=None):
        self.tools, self.knowledge, self.cemr = tools, knowledge, cemr
        self.reflection, self.similar = reflection, similar
        self.ranking_readers = ranking_readers

    def run(self, patient, excluded_pmids=(), *, progress=None):
        require_positive_hpo(patient)
        genotype = patient.genomic is not None
        memory, trace = OTEM(patient.key), []
        observer_failures = 0

        def emit(name, event, record=None):
            nonlocal observer_failures
            LOG.info("[%s/%s] %s %s", *patient.key, event, name)
            if progress is not None:
                try:
                    progress({"case_key": patient.key, "stage_id": name, "event": event, "record": record})
                except Exception:
                    observer_failures += 1
                    LOG.warning("Stage observer failed; scientific execution continues")

        def step(name, action, quality=None):
            emit(name, "started")
            started = time.perf_counter()
            status, issues, fallback = "failed", ("stage_execution_failed",), False
            try:
                value = action()
                status, issues = quality(value) if quality else _quality(value)
                fallback = status != "ok"
                return value
            finally:
                record = {
                    "stage_id": name,
                    "status": status,
                    "elapsed_seconds": max(0.0, time.perf_counter() - started),
                    "error_codes": issues,
                    "fallback_used": fallback,
                }
                trace.append(record)
                emit(name, "finished", record)

        memory.register_patient(patient)

        def prepare_profile():
            profile = self.tools.profile(patient)
            memory.register_profile(self.knowledge.hpo(), profile)
            return profile

        profile = step("phenotype_profile", prepare_profile)
        if genotype:
            stages = self._genotype(patient, profile, memory, step, excluded_pmids)
        else:
            stages = self._phenotype(profile, memory, step, excluded_pmids)
        hypothesis, fusion, admission, evidence = (
            stages[name] for name in ("hypotheses", "fusion", "admission", "evidence")
        )
        similar, molecular = stages.get("similar"), stages.get("molecular_context")
        ranking_options = {}
        review = admission.genes[:10] if genotype else admission.genes

        def review_evidence():
            if genotype:
                ranking_options["vcf_views"] = views.ranking_context_views(
                    profile["case_key"], review, stages["ranking_context"]
                )
                views.validate_molecular_context(molecular, profile["case_key"], review)
            return self.cemr.brief(
                profile,
                evidence,
                review,
                similar=similar["evidence_by_gene"] if similar else None,
                **ranking_options,
            )

        briefs = step("evidence_brief", review_evidence)
        global_result = step(
            "global_evidence_ranking",
            lambda: self.cemr.global_rank(
                profile,
                admission,
                hypothesis,
                briefs,
                **ranking_options,
            ),
        )
        pairwise_count = min(5, len(global_result.genes)) if genotype else 10
        pairwise = step(
            "pairwise_evidence_ranking",
            lambda: self.cemr.pairwise(
                profile,
                global_result.genes[:pairwise_count],
                admission,
                hypothesis,
                briefs,
                include_clinical=bool(profile["clinical_description"].strip()),
                required_count=pairwise_count,
                timeout_seconds=120,
                **ranking_options,
            ),
        )
        head = pairwise.genes if pairwise.status == "ok" else global_result.genes[:pairwise_count]
        fixed = tuple(dict.fromkeys((*head, *global_result.genes, *admission.genes, *fusion.genes)))[:20]
        tail = step(
            "tail_candidate_ranking",
            lambda: self.cemr.tail(
                profile,
                fixed,
                hypothesis,
                fusion.genes,
                mode="genotype" if genotype else "phenotype",
                **({"native": stages["acquisition"]["evidence"]} if genotype else {}),
            ),
        )
        challenger, reflection, final = None, None, tail.genes
        if similar is not None:
            challenger = step(
                "tail_challenger_ranking", lambda: self.cemr.challenge_tail(profile, tail.genes, similar)
            )
            final = challenger.genes
        if self.reflection and not genotype:

            def refine():
                nonlocal final
                if pairwise.genes[:2] != final[:2]:
                    return self.cemr._reflection(
                        pairwise.genes[:2], "pairwise_top2_not_current_final_top2", -1.0
                    )
                result = self.cemr.reflect(
                    profile,
                    pairwise,
                    evidence,
                    briefs,
                    self.knowledge.retrieve_reflection,
                    excluded_pmids=excluded_pmids,
                )
                if result.status == "ok":
                    updates = result.details["evidence_updates"]
                    _register_documents(memory, updates)
                    for gene, update in updates.items():
                        old = evidence.get(
                            gene, {"documents": (), "omim_context": (), "retrieval_errors": ()}
                        )
                        documents = {
                            doc["prompt_id"]: doc for doc in (*old["documents"], *update["documents"])
                        }
                        evidence[gene] = {
                            "gene": gene,
                            "documents": tuple(documents.values()),
                            "omim_context": old["omim_context"],
                            "retrieval_errors": (*old["retrieval_errors"], *update["retrieval_errors"]),
                        }
                    if result.genes != result.details["original_top2"]:
                        final = (*result.genes, *(gene for gene in final if gene not in result.genes))
                return result

            reflection = step("top2_reflection", refine)
        stages.update(
            profile=profile,
            similar=similar,
            briefs=briefs,
            global_result=global_result,
            pairwise=pairwise,
            tail=tail,
            challenger=challenger,
            reflection=reflection,
            trace=tuple(trace),
            observer_failures=observer_failures,
        )
        observations = tuple(
            Evidence(
                evidence_id=f"gene-documents:{gene}",
                kind="gene_documents",
                case_key=patient.key,
                genes=(gene,),
                facts=value,
                source_ids=tuple(
                    doc["source"].source_id for doc in (*value["documents"], *value["omim_context"])
                ),
            )
            for gene, value in evidence.items()
        )
        memory.add(evidence=observations)
        stages.pop("evidence")
        return CaseResult(
            patient=patient,
            stages=stages,
            final=Ranking(genes=final),
            sources=memory.sources,
            evidence=memory.evidence,
            status="ok" if all(row["status"] == "ok" for row in trace) else "degraded",
        )

    def _phenotype(self, profile, memory, step, excluded_pmids):
        def acquire():
            value = self.tools.phenotype_candidates(profile)
            memory.add(tuple(source for row in value["rankings"] for source in row["sources"]))
            return value

        pool = step(
            "candidate_acquisition",
            acquire,
            lambda value: (
                "degraded" if value["failures"] else "ok",
                tuple(row["code"] for row in value["failures"]),
            ),
        )

        def hypothesize():
            hypothesis = self.cemr.hypothesize(
                profile,
                pool["rankings"],
                mode="phenotype",
                article_select=self.knowledge.select_articles,
                excluded_pmids=excluded_pmids,
            )
            fusion = self.cemr.aggregate(
                profile,
                hypothesis,
                article_select=self.knowledge.select_articles,
                excluded_pmids=excluded_pmids,
            )
            articles = (*hypothesis["amelie_articles"].values(), *fusion.details.get("articles", {}).values())
            memory.add(tuple(row["source"] for row in articles if row["source"] is not None))
            return hypothesis, fusion

        hypothesis, fusion = step("multi_source_ranking", hypothesize, _hypothesis_quality)
        admission = step(
            "evidence_admission", lambda: self.cemr.admit(profile, hypothesis, fusion, mode="phenotype")
        )
        similar = None
        if self.similar:
            similar = step(
                "similar_case_evidence",
                lambda: memory.retrieve_similar(
                    profile,
                    admission.genes,
                    self.knowledge.similar_cases,
                    self.knowledge.hpo(),
                    self.cemr.model,
                ),
            )
        must_review = {
            row["gene"]: (row["pmid"],) for row in fusion.details.get("articles", {}).values() if row["pmid"]
        }
        evidence = step(
            "evidence_retrieval:0",
            lambda: self._ground(
                memory, profile, admission.genes, hypothesis["paper_scores"], must_review, excluded_pmids
            ),
            _retrieval_quality,
        )
        return {
            "candidates": pool,
            "hypotheses": hypothesis,
            "fusion": fusion,
            "admission": admission,
            "similar": similar,
            "evidence": evidence,
        }

    def _genotype(self, patient, profile, memory, step, excluded_pmids):
        from mendelprime.tools.acmg import ranking_literature, ranking_preliminary

        from .assessment import prepare

        def acquire():
            value = self.tools.genotype(patient)
            require(isinstance(value, dict), "genotype_acquisition", "typed_acquisition_required")
            require(
                value["evidence"]["case_key"] == patient.key, "genotype_acquisition.case_key", "case_mismatch"
            )
            memory.add(value["evidence"]["sources"])
            return value

        acquisition = step("genotype_acquisition", acquire, _acquisition_quality)
        native = acquisition["evidence"]
        require(native["tools"], "genotype_acquisition", "native_results_required")
        hypothesis = step(
            "genotype_source_ranking",
            lambda: self.cemr.hypothesize(profile, native["tools"], mode="genotype", native=native),
            lambda value: _rankings_quality(value["sources"].values()),
        )
        pool = tuple(dict.fromkeys(gene for row in hypothesis["sources"].values() for gene in row.genes[:10]))
        require(bool(pool), "genotype_candidate_pool", "candidate_genes_required")
        require(self.ranking_readers is not None, "ranking_resources", "required_for_vcf_evidence")
        evidence, knowledge = {}, {}

        def prepare_scope(genes, previous=None):
            seed = {
                "profile": profile,
                "fixed_genes": genes,
                "family": acquisition["family"],
                "native_genes": tuple(row for tool in native["tools"] for row in tool["genes"]),
                "observations": native["observations"],
                "family_sources": native["family_sources"],
                "variant_knowledge": knowledge,
                "evidence": evidence,
                "sources": memory.sources,
            }
            prepared = prepare(seed, self.ranking_readers, genes=genes, previous=previous)
            memory.add(prepared["input"]["sources"])
            knowledge.update(prepared["clinical"])
            _validate_seeds(memory, knowledge)
            return prepared

        prepared = step("ranking_variant_preparation", lambda: prepare_scope(pool), _preparation_quality)
        ledger = step("preliminary_acmg", lambda: ranking_preliminary(prepared))
        context = {"prepared": prepared, "acmg": ledger, "evidence": evidence}
        fusion = step(
            "genotype_fusion",
            lambda: self.cemr.aggregate(profile, hypothesis, native=native, ranking_context=context),
        )
        initial_review = fusion.genes[:10]

        def retrieve(genes):
            seeds = tuple(seed for gene in genes for seed in knowledge[gene]["citation_seeds"])
            return self._ground(
                memory,
                profile,
                genes,
                {},
                {},
                excluded_pmids,
                seeds,
                variant_targets=prepared["input"]["targets"],
            )

        evidence = step(
            "evidence_retrieval:0",
            lambda: retrieve(initial_review),
            _retrieval_quality,
        )

        def expand(genes):
            updates = self.knowledge.expand_variants(
                evidence,
                prepared["input"]["targets"],
                genes,
                excluded_pmids=excluded_pmids,
                ranking_acmg=ledger,
            )
            _register_documents(memory, updates)
            evidence.update(updates)
            return updates

        def assess_literature(genes):
            nonlocal ledger
            update = ranking_literature(
                self.cemr.model,
                prepared,
                ledger,
                evidence,
                targets=genes,
                workers=self.cemr.workers["literature"],
            )
            memory.add(update["sources"])
            _register_documents(memory, update["evidence_updates"])
            evidence.update(update["evidence_updates"])
            ledger = update["ledger"]
            context.update(prepared=prepared, acmg=ledger, evidence=evidence)
            return update

        step("variant_literature_expansion:0", lambda: expand(initial_review), _retrieval_quality)
        step("literature_acmg:0", lambda: assess_literature(initial_review), _preparation_quality)
        admission = step(
            "genotype_admission",
            lambda: self.cemr.admit(
                profile, hypothesis, fusion, mode="genotype", native=native, ranking_context=context
            ),
        )
        admitted = admission.genes[:10]
        # Existing gene names do not establish allele/resource coverage. Rebuild
        # bindings while the preparation cache reuses completed reader facts.
        prepared = step(
            "admitted_variant_preparation", lambda: prepare_scope(admitted, prepared), _preparation_quality
        )
        missing = tuple(
            gene
            for gene in admitted
            if gene not in evidence or _retrieval_quality({gene: evidence[gene]})[0] != "ok"
        )
        if missing:
            updates = step("evidence_retrieval:1", lambda: retrieve(missing), _retrieval_quality)
            for gene, update in updates.items():
                old = evidence.get(gene)
                if old:
                    documents = {
                        doc["source"].source_id: doc for doc in (*old["documents"], *update["documents"])
                    }
                    update = {**update, "documents": tuple(documents.values())}
                evidence[gene] = update
        step("variant_literature_expansion:1", lambda: expand(admitted), _retrieval_quality)
        step("literature_acmg:1", lambda: assess_literature(admitted), _preparation_quality)
        molecular = step(
            "molecular_context",
            lambda: memory.register_molecular(native, admitted, acquisition["family"], knowledge),
        )
        return {
            "acquisition": acquisition,
            "hypotheses": hypothesis,
            "fusion": fusion,
            "admission": admission,
            "variant_knowledge": knowledge,
            "evidence": evidence,
            "molecular_context": molecular,
            "ranking_context": context,
            "ranking_acmg": ledger,
        }

    def _ground(self, memory, profile, genes, papers, must_review, excluded_pmids, seeds=(), **options):
        evidence = self.knowledge.ground(
            profile, genes, papers, must_review, excluded_pmids, citation_seeds=seeds, **options
        )
        _register_documents(memory, evidence)
        return evidence


def _register_documents(memory, evidence):
    memory.add(
        tuple(
            doc["source"]
            for value in evidence.values()
            for doc in (*value["documents"], *value["omim_context"])
        )
    )


def _quality(value):
    if isinstance(value, Ranking):
        return value.status, value.issues
    return value.get("status", "ok"), value.get("issues", ()) or (
        (value["error_code"],) if value.get("error_code") else ()
    )


def _hypothesis_quality(value):
    hypothesis, fusion = value
    return _rankings_quality((fusion, *hypothesis["sources"].values()))


def _rankings_quality(results):
    results = tuple(results)
    return (
        "ok" if all(row.status == "ok" for row in results) else "degraded",
        tuple(issue for row in results for issue in row.issues),
    )


def _acquisition_quality(value):
    errors = tuple(
        dict.fromkeys(
            (
                *(row["code"] for row in value["failures"]),
                *(
                    row["error_code"]
                    for row in value["executions"]
                    if row["status"] != "ok" and row["error_code"]
                ),
                *(issue["code"] for tool in value["evidence"]["tools"] for issue in tool["issues"]),
            )
        )
    )
    return (value["status"] if value["evidence"]["tools"] else "failed"), errors


def _preparation_quality(value):
    issues = value.get("issues", ())
    return ("degraded" if any(not is_evidence_limitation(issue) for issue in issues) else "ok"), issues


def _retrieval_quality(evidence):
    errors = tuple(
        dict.fromkeys(
            error
            for row in evidence.values()
            for error in row["retrieval_errors"]
            if error != "article_or_abstract_missing" and not error.endswith(":article_or_abstract_missing")
        )
    )
    return ("degraded" if errors else "ok"), errors


def _validate_seeds(memory, knowledge):
    cited = {
        document["source"].source_id: set(document["cited_pmids"])
        for value in knowledge.values()
        for document in (*[row["document"] for row in value["variant_evidence"]], *value["gene_context"])
    }
    for value in knowledge.values():
        for seed in value["citation_seeds"]:
            try:
                source = memory.get_source(seed["discovered_from_source_id"])
            except KeyError:
                raise MendelPrimeError("source_not_registered", field="citation_seeds") from None
            require(
                source.genes == (seed["gene"],) and seed["variant_key"] in source.variant_keys,
                "citation_seeds",
                "source_scope_mismatch",
            )
            require(seed["pmid"] in cited.get(source.source_id, ()), "citation_seeds", "pmid_not_in_source")
