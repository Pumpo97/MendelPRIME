"""One per-run owner for shared knowledge readers and query caches."""

import threading

from mendelprime.domain import freeze_record
from mendelprime.runtime.storage import private_directory

from . import hpo, literature, retrieval


class KnowledgeBase:
    def __init__(self, resources, cache_root, device="cuda"):
        self.resources = freeze_record(resources)
        self.cache_root, self.device = cache_root, device
        self._readers = {}
        self._lock = threading.RLock()

    def _once(self, name, loader, *args, **kwargs):
        with self._lock:
            if name not in self._readers:
                self._readers[name] = loader(*args, **kwargs)
            return self._readers[name]

    def hpo(self):
        resources = self.resources["resources"]
        return self._once("hpo", hpo.load_hpo, resources["hpo_obo"], resources["hpo_stats"])

    def omim(self):
        resources = self.resources["resources"]
        return self._once("omim", literature.load_omim, resources["omim_jsonl"], resources["mim2gene"])

    def pubmed(self):
        private_directory(self.cache_root)
        return self._once("pubmed", literature.new_pubmed, self.cache_root / "pubmed")

    def select_articles(self, genes, paper_scores, excluded_pmids=()):
        return literature.select_articles(self.pubmed(), genes, paper_scores, excluded_pmids)

    def ground(
        self,
        profile,
        genes,
        paper_scores,
        must_review,
        excluded_pmids=(),
        *,
        citation_seeds=(),
        variant_targets=None,
    ):
        variant_options = {}
        if variant_targets is not None:
            variant_options = {
                "variant_aware": True,
                "variant_search_terms_by_gene": {
                    gene: tuple(
                        dict.fromkeys(
                            term
                            for target in variant_targets
                            if target["primary"]["gene"] == gene
                            for term in literature._variant_search_terms(target)
                        )
                    )
                    for gene in genes
                },
            }
        return literature.retrieve_evidence(
            self.pubmed(),
            self.omim(),
            profile,
            genes,
            paper_scores,
            must_review,
            excluded_pmids,
            citation_seeds=citation_seeds,
            **variant_options,
        )

    def retrieve_reflection(self, genes, queries, excluded_pmids=()):
        return literature.retrieve_reflection(self.pubmed(), genes, queries, excluded_pmids)

    def expand_variants(self, evidence, targets, genes, *, excluded_pmids=(), ranking_acmg=None):
        selected = {}
        for gene, assessment in (ranking_acmg or {}).get("by_gene", {}).items():
            genotype = assessment["genotype_configuration"]
            if genotype.get("selected_variant_refs_role") == "selected_hypothesis":
                refs = set(assessment["selected_variant_refs"])
                selected[gene] = tuple(
                    variant["identity"]["variant_key"]
                    for variant in assessment["variants"]
                    if variant["variant_ref"] in refs
                )
        return {
            gene: literature.expand_variant_evidence(
                self.pubmed(),
                gene,
                evidence.get(
                    gene, {"gene": gene, "documents": (), "omim_context": (), "retrieval_errors": ()}
                ),
                tuple(target for target in targets if target["primary"]["gene"] == gene),
                excluded_pmids=excluded_pmids,
                selected_variant_keys=selected.get(gene),
            )
            for gene in genes
        }

    def variant_lookup(self, gene, keys, reference=None):
        from . import variants

        state = self._once("variants", variants.new_variant_store, self.resources["variant_knowledge"])
        return variants.lookup(state, gene, keys, reference)

    def scope_native_transcripts(self, tools):
        from . import variants

        state = self._once("variants", variants.new_variant_store, self.resources["variant_knowledge"])
        return variants.scope_native_transcripts(state, tools)

    def similar_cases(self, profile):
        resources = self.resources["resources"]
        index = self._once("cases", retrieval.load_index, resources["similar_case_root"])
        dense = retrieval.dense_hits(
            index, profile, resources["embedding_model"], resources["cross_encoder"], self.device
        )
        semantic = retrieval.semantic_hits(index, self.hpo(), profile, True, 10, 50)
        return retrieval.hybrid_hits(dense, semantic)

    def close(self):
        with self._lock:
            self._readers.clear()
