"""The sole stateful CEMR entry point."""

import json
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations

from mendelprime import prompting
from mendelprime.domain import MendelPrimeError, Ranking, source_diverse_rrf

from . import hypotheses, ranking, tail, views

LOG = logging.getLogger("mendelprime.cemr")

DEFAULT_WORKERS = {"sources": 4, "briefs": 8, "pairwise": 8, "literature": 4}
PAIRWISE_RETRY_WORKERS = 2
# Frozen request policy: asset, output tokens, timeout, structural attempts.
REQUEST_POLICY = {
    **{method: ("tool_ranking", 4096, 180, 1) for method in hypotheses.PRIMARY},
    "amelie_article": ("amelie_ranking", 4096, 180, 2),
    "direct_hpo": ("direct_hpo_ranking", 4096, 180, 2),
    "genotype_source": ("genotype_source_ranking", 5000, 300, 2),
    "genotype_fusion": ("vcf_genotype_fusion", 3000, 180, 2),
    "candidate_fusion": ("fusion_ranking", 4096, 180, 2),
    "vcf_candidate_admission": ("vcf_candidate_admission", 3000, 180, 2),
    "vcf_evidence_brief": ("vcf_evidence_brief", 1800, 150, 2),
    "vcf_causal_shortlist": ("vcf_global_ranking", 3200, 180, 2),
    "vcf_structured_pairwise": ("vcf_pairwise", 200, 120, 1),
    "vcf_tail_ranking": ("vcf_tail_ranking", 2400, 180, 2),
    "evidence_brief": ("evidence_brief", 2600, 150, 2),
    "global_evidence_ranking": ("global_ranking", 4500, 180, 2),
    "pairwise_ranking": ("pairwise", 800, 120, 1),
    "tail_ranking": ("tail_ranking", 2400, 180, 2),
    "tail_challenger_ranking_v1": ("tail_challenger", 3000, 180, 2),
    "top2_evidence_gap_analysis": ("reflection_gap", 3200, 300, 1),
    "top2_evidence_gap_fill_assessment": ("reflection_fill", 3200, 300, 1),
    "top2_reflection_pairwise": ("reflection_pairwise", 800, 300, 1),
}


class CEMR:
    """Execute individual scientific stages; Harness owns their order."""

    def __init__(self, model, workers=None):
        self.model = model
        self.workers = {**DEFAULT_WORKERS, **(workers or {})}
        if set(self.workers) != set(DEFAULT_WORKERS):
            raise ValueError("unknown CEMR worker setting")
        if any(type(value) is not int or value < 1 for value in self.workers.values()):
            raise ValueError("CEMR worker counts must be positive integers")

    @staticmethod
    def _request(profile, stage, prompt, *, system=None, timeout=None, cache=False):
        asset, tokens, default_timeout, _ = REQUEST_POLICY[stage.split(":", 1)[0]]
        return {
            "case_key": profile["case_key"],
            "stage_id": stage,
            "system_prompt": prompting.text(asset) if system is None else system,
            "user_prompt": prompt,
            "max_tokens": tokens,
            "timeout_seconds": default_timeout if timeout is None else timeout,
            "cache_choice": cache,
        }

    def _call(self, profile, stage, prompt, decode, *, repair=None, system=None):
        key = stage.split(":", 1)[0]
        limit = 220000 if key == "genotype_fusion" else 160000

        def within_budget(value):
            if not key.startswith(("genotype_", "vcf_")) or len(value) <= limit:
                return value
            try:
                payload = json.loads(value)
                sections = sorted(
                    (
                        (name, len(json.dumps(value, ensure_ascii=False, separators=(",", ":"))))
                        for name, value in payload.items()
                    ),
                    key=lambda item: item[1],
                    reverse=True,
                )[:4]
                summary = ",".join(f"{name}:{size}" for name, size in sections)
            except (TypeError, ValueError):
                summary = "unparseable"
            diagnostics = f"characters={len(value)};limit={limit};largest_sections={summary}"
            case_key = profile.get("case_key", ("unknown", "unknown"))
            LOG.warning("[%s/%s] %s prompt rejected: %s", *case_key, key, diagnostics)
            raise MendelPrimeError(
                "vcf_prompt_exceeds_limit",
                field="prompt",
                source=diagnostics,
            )

        request = self._request(profile, stage, within_budget(prompt), system=system)
        attempts = REQUEST_POLICY[stage.split(":", 1)[0]][3]
        bounded_repair = (lambda code: within_budget(repair(code))) if repair is not None else None
        return prompting.call_json(self.model, request, decode, attempts, bounded_repair)

    @staticmethod
    def _parallel(calls, workers, failure):
        output = {}
        with ThreadPoolExecutor(max_workers=min(workers, len(calls)) or 1) as executor:
            futures = {executor.submit(call): key for key, call in calls.items()}
            for future in as_completed(futures):
                key = futures[future]
                try:
                    output[key] = future.result()
                except Exception:
                    output[key] = failure(key)
        return output

    @staticmethod
    def _source(method, genes, kind, error="", reasoning=()):
        result = hypotheses.source_ranking(
            method,
            genes,
            status="degraded" if error else "ok",
            error_code=error,
            reasoning=reasoning,
            kind=kind,
        )
        if kind == "genotype_tool_llm":
            return Ranking(
                genes=result.genes,
                stage=result.stage,
                status=result.status,
                issues=result.issues,
                details={**result.details, "origin": "raw_tool_fallback" if error else "tool_llm"},
            )
        return result

    def _rank_source(self, profile, method, prompt, allowed, kind="tool_llm", stage=None):
        if not allowed:
            return self._source(method, (), kind, "no_valid_allowed_genes")
        try:
            genes, reasons = self._call(
                profile,
                stage or method,
                prompt,
                lambda value: (
                    hypotheses.decode_vcf_subset(value, allowed)
                    if kind == "genotype_tool_llm"
                    else hypotheses.decode_allowed_ranking(value, allowed)
                ),
            )
            return self._source(method, genes, kind, reasoning=reasons)
        except MendelPrimeError as error:
            return self._source(method, allowed, kind, error.code)

    def _phenotype_source(self, profile, method, raw):
        genes = hypotheses.tool_genes(raw)[:100] if raw else ()
        if not genes:
            return self._source(method, (), "tool_llm", "no_valid_allowed_genes")
        for reason_limit in (10, 10, 10, 0):
            result = self._rank_source(
                profile,
                method,
                hypotheses.tool_prompt(profile, genes, raw, reason_limit),
                genes,
            )
            if result.status == "ok":
                break
        return result

    def _direct(self, profile):
        try:
            genes = self._call(
                profile,
                "direct_hpo",
                hypotheses.direct_prompt(profile),
                hypotheses.decode_direct,
            )
            return self._source("direct_hpo", genes, "direct_llm")
        except MendelPrimeError as error:
            return self._source("direct_hpo", (), "direct_llm", error.code)

    def hypothesize(
        self, profile, tool_records, *, mode, article_select=None, native=None, excluded_pmids=()
    ):
        if mode == "genotype":
            return self._genotype_hypotheses(profile, tool_records, native)
        if mode != "phenotype":
            raise MendelPrimeError("mode_required", field="mode")
        raw = {record["tool_id"]: record for record in tool_records}
        papers = {
            item["gene"]: tuple(item.get("context") or ())
            for item in raw.get("amelie", {}).get("entries", ())
            if isinstance(item.get("context"), (tuple, list))
        }
        amelie_genes = hypotheses.tool_genes(raw["amelie"])[:20] if "amelie" in raw else ()
        articles = (
            article_select(amelie_genes, papers, excluded_pmids)
            if article_select is not None and amelie_genes
            else {}
        )
        calls = {
            method: (lambda method=method: self._phenotype_source(profile, method, raw.get(method)))
            for method in hypotheses.PRIMARY
        }
        calls["amelie_article"] = lambda: self._rank_source(
            profile,
            "amelie_article",
            hypotheses.amelie_prompt(profile, amelie_genes, articles),
            amelie_genes,
        )
        calls["direct_hpo"] = lambda: self._direct(profile)

        def failed(method):
            fallback = (
                amelie_genes
                if method == "amelie_article"
                else hypotheses.tool_genes(raw[method])
                if method in raw
                else ()
            )
            kind = "direct_llm" if method == "direct_hpo" else "tool_llm"
            return self._source(method, fallback, kind, "source_execution_failed")

        output = self._parallel(calls, self.workers["sources"], failed)
        raw_rankings = {record["tool_id"]: hypotheses.tool_genes(record) for record in tool_records}
        consensus = source_diverse_rrf(
            raw_rankings, 100, 30, 10, ("phen2disease", "phen2gene", "local_overlap", "amelie")
        )
        return {
            "mode": mode,
            "sources": {method: output[method] for method in hypotheses.PHENOTYPE_SOURCE_ORDER},
            "consensus": consensus,
            "paper_scores": papers,
            "amelie_articles": articles,
            "issues": tuple(issue for record in tool_records for issue in record.get("issues", ())),
        }

    def _genotype_hypotheses(self, profile, tools, native):
        if native is None or native.get("case_key") != profile["case_key"]:
            raise MendelPrimeError("case_or_evidence_mismatch", field="native")

        by_method = {tool["tool_id"]: tool for tool in tools}
        calls = {
            method: lambda tool=tool: self._rank_source(
                profile,
                tool["tool_id"],
                ranking.vcf_source_prompt(profile, tool, native["observations"]),
                tool["ranking"].genes[:100],
                kind="genotype_tool_llm",
                stage=f"genotype_source:{tool['tool_id']}",
            )
            for method, tool in by_method.items()
        }
        output = self._parallel(
            calls,
            self.workers["sources"],
            lambda method: self._source(
                method,
                by_method[method]["ranking"].genes[:100],
                "genotype_tool_llm",
                "source_execution_failed",
            ),
        )
        return {
            "mode": "genotype",
            "sources": {method: output[method] for method in calls},
            "consensus": (),
            "paper_scores": {},
            "amelie_articles": {},
            "issues": (),
        }

    @staticmethod
    def _source_genes(hypothesis):
        return {method: source.genes for method, source in hypothesis["sources"].items()}

    def _validate_native(self, hypothesis, native):
        if native is None or native.get("case_key") is None:
            raise MendelPrimeError("genotype_evidence_required", field="native")
        expected = tuple(tool["tool_id"] for tool in native["tools"])
        actual = tuple(hypothesis["sources"])
        if actual != expected:
            raise MendelPrimeError("source_order_or_coverage_mismatch", field="sources")
        for source, tool in zip(hypothesis["sources"].values(), native["tools"], strict=True):
            if len(set(source.genes)) != len(source.genes) or not set(source.genes) <= set(
                tool["ranking"].genes[:100]
            ):
                raise MendelPrimeError("source_candidate_mismatch", field="sources")

    def aggregate(
        self,
        profile,
        hypothesis,
        *,
        article_select=None,
        excluded_pmids=(),
        native=None,
        ranking_context=None,
    ):
        if hypothesis["mode"] == "genotype":
            if native is None or native.get("case_key") != profile["case_key"]:
                raise MendelPrimeError("case_mismatch", field="native")
            prepared, ranking_acmg, evidence = (
                ranking_context[key] for key in ("prepared", "acmg", "evidence")
            )
            self._validate_native(hypothesis, native)
            pool = hypotheses.genotype_pool(hypothesis, native)
            nominated = pool["genes"]
            if not nominated:
                return Ranking(
                    genes=(),
                    stage="fusion",
                    status="degraded",
                    details={"nominated_genes": (), "reasoning": ()},
                    issues=("no_source_candidates",),
                )
            if not any(source.status == "ok" for source in hypothesis["sources"].values()):
                return Ranking(
                    genes=nominated,
                    stage="fusion",
                    status="degraded",
                    details={"nominated_genes": nominated, "reasoning": ()},
                    issues=("all_genotype_toolllm_rankings_degraded",),
                )
            prompt = ranking.vcf_fusion_prompt(
                profile, hypothesis, native, pool, prepared, ranking_acmg, evidence
            )
            stage = "genotype_fusion"
        else:
            eligible = {
                method: source.genes
                for method, source in hypothesis["sources"].items()
                if method != "direct_hpo" or source.status == "ok"
            }
            votes = hypotheses.nominations(eligible, 5)
            nominated = hypotheses.nomination_order(votes)[:100]
            if not nominated:
                return Ranking(
                    genes=hypothesis["consensus"],
                    stage="fusion",
                    status="degraded",
                    details={
                        "nominated_genes": (),
                        "reasoning": (),
                        "articles": {},
                        "amelie_articles": hypothesis["amelie_articles"],
                    },
                    issues=("no_source_candidates",),
                )
            articles = (
                article_select(nominated, hypothesis["paper_scores"], excluded_pmids)
                if article_select is not None
                else {}
            )
            prompt = hypotheses.fusion_prompt(profile, nominated, votes, articles)
            stage = "candidate_fusion"
        try:
            genes, reasons = self._call(
                profile,
                stage,
                prompt,
                lambda value: (
                    hypotheses.decode_vcf_permutation(value, nominated)
                    if hypothesis["mode"] == "genotype"
                    else hypotheses.decode_allowed_ranking(value, nominated)
                ),
                repair=(lambda error: ranking.vcf_repair(prompt, error))
                if hypothesis["mode"] == "genotype"
                else None,
            )
            status, issues = "ok", ()
        except MendelPrimeError as error:
            genes, reasons, status, issues = nominated, (), "degraded", (error.code,)
        details = {"nominated_genes": nominated, "reasoning": reasons}
        if hypothesis["mode"] == "phenotype":
            details.update(articles=articles, amelie_articles=hypothesis["amelie_articles"])
            genes = tuple(dict.fromkeys((*genes, *hypothesis["consensus"], *eligible.get("direct_hpo", ()))))
        return Ranking(genes=genes, stage="fusion", status=status, details=details, issues=issues)

    def admit(
        self,
        profile,
        hypothesis,
        fusion,
        *,
        mode,
        native=None,
        ranking_context=None,
    ):
        source_genes = self._source_genes(hypothesis)
        if mode == "phenotype":
            usable = {
                method: genes
                for method, genes in source_genes.items()
                if genes and (method != "direct_hpo" or hypothesis["sources"][method].status == "ok")
            }
            genes, status, issues = hypotheses.admission_order(fusion.genes, usable), "ok", ()
        else:
            prepared, ranking_acmg, evidence = (
                ranking_context[key] for key in ("prepared", "acmg", "evidence")
            )
            if native is None or native.get("case_key") != profile["case_key"]:
                raise MendelPrimeError("case_mismatch", field="native")
            self._validate_native(hypothesis, native)
            candidate_pool = hypotheses.genotype_pool(hypothesis, native)
            expected = candidate_pool["genes"]
            if fusion.details.get("nominated_genes") != expected or set(fusion.genes) != set(expected):
                raise MendelPrimeError("fusion_candidate_mismatch", field="fusion")
            pool, count = fusion.genes[:20], min(15, len(fusion.genes))
            if not count:
                return Ranking(
                    genes=(),
                    stage="admission",
                    status="degraded",
                    details={"fusion_genes": fusion.genes, "source_kinds": {}},
                    issues=("no_fusion_candidates",),
                )
            try:
                prompt = ranking.vcf_admission_prompt(
                    profile, fusion, native, candidate_pool, count, prepared, ranking_acmg, evidence
                )
                genes, full_ranking = self._call(
                    profile,
                    "vcf_candidate_admission",
                    prompt,
                    lambda value: hypotheses.decode_vcf_admission(value, pool, count),
                    repair=lambda error: ranking.vcf_repair(prompt, error),
                )
                status, issues = "ok", ()
            except MendelPrimeError as error:
                genes, status, issues = pool[:count], "degraded", (error.code,)
                full_ranking = pool
        kinds = {
            method: source.details["kind"] for method, source in hypothesis["sources"].items() if source.genes
        }
        if mode == "genotype":
            kinds = candidate_pool["source_kinds"]
        return Ranking(
            genes=genes,
            stage="admission",
            status=status,
            details={
                "fusion_genes": fusion.genes,
                "source_kinds": kinds,
                **({"full_ranking": full_ranking} if mode == "genotype" else {}),
                **(
                    {
                        "anchor_evidence_ids": {
                            gene: document["prompt_id"]
                            for gene in full_ranking
                            if (document := ranking.vcf_anchor_document((evidence or {}).get(gene, {})))
                            is not None
                        }
                    }
                    if mode == "genotype"
                    else {}
                ),
            },
            issues=issues,
        )

    def brief(
        self,
        profile,
        evidence,
        genes,
        *,
        similar=None,
        system_prompt=None,
        vcf_views=None,
    ):
        if vcf_views is not None and len(genes) > 10:
            raise MendelPrimeError("vcf_brief_requires_top10", field="genes")
        if vcf_views is not None:
            for gene in genes:
                ranking.vcf_article_catalog(ranking.evidence_pool(evidence, gene))
        system = system_prompt or prompting.text(
            "vcf_evidence_brief" if vcf_views is not None else "evidence_brief"
        )

        def outcome(gene, selected=(), status="ok", code=""):
            return {
                "gene": gene,
                "selected_evidence": selected,
                "similar_case_evidence": tuple((similar or {}).get(gene, ())),
                "status": status,
                "error_code": code,
            }

        def one(gene):
            pool = ranking.evidence_pool(evidence, gene)
            allowed = (
                {
                    identifier: {**entry["document"], "article_id": entry["article_id"]}
                    for identifier, entry in ranking.vcf_article_catalog(pool).items()
                }
                if vcf_views is not None
                else {item["prompt_id"]: item for item in pool["documents"]}
            )
            code, transport = "", False
            stage = "vcf_evidence_brief" if vcf_views is not None else "evidence_brief"
            for attempt in range(REQUEST_POLICY[stage][3]):
                prompt = (
                    ranking.vcf_brief_prompt(profile, gene, pool, vcf_views[gene], code if attempt else "")
                    if vcf_views is not None
                    else ranking.with_similar(
                        ranking.brief_prompt(
                            profile,
                            gene,
                            pool,
                            code if attempt else "",
                        ),
                        tuple((similar or {}).get(gene, ())),
                    )
                )
                request = self._request(profile, f"{stage}:{gene}", prompt, system=system)
                try:
                    reply = self.model.complete(request)
                except MendelPrimeError as error:
                    code, transport = error.code, True
                    continue
                except Exception:
                    code, transport = "brief_model_exception", True
                    continue
                try:
                    parsed = prompting.parse_object(reply["content"])
                    selected = (
                        ranking.vcf_parse_brief(parsed, gene, allowed)
                        if vcf_views is not None
                        else ranking.parse_brief(parsed, gene, allowed, 3)
                    )
                except MendelPrimeError as error:
                    code, transport = error.code, False
                    continue
                return outcome(gene, selected)
            return outcome(gene, status="failed" if transport else "degraded", code=code)

        results = self._parallel(
            {gene: lambda gene=gene: one(gene) for gene in genes},
            self.workers["briefs"],
            lambda gene: outcome(gene, status="failed", code="brief_worker_exception"),
        )
        ordered = tuple(results[gene] for gene in genes)
        issues = tuple(item["error_code"] for item in ordered if item["error_code"])
        return {
            "stage": "evidence_brief",
            "status": "degraded" if issues else "ok",
            "issues": issues,
            "results": ordered,
        }

    def global_rank(
        self,
        profile,
        admission,
        hypothesis,
        briefs,
        *,
        system_prompt=None,
        vcf_views=None,
    ):
        if vcf_views is not None:
            admission = Ranking(
                genes=admission.genes[:10],
                stage=admission.stage,
                status=admission.status,
                details=admission.details,
                issues=admission.issues,
            )
        system = system_prompt or prompting.text(
            "vcf_global_ranking" if vcf_views is not None else "global_ranking"
        )
        if any(item["similar_case_evidence"] for item in briefs["results"]):
            system += prompting.text("global_similar_overlay")

        def build(error=""):
            if vcf_views is not None:
                return ranking.vcf_global_prompt(profile, admission, hypothesis, briefs, vcf_views, error)
            return ranking.global_prompt(profile, admission, hypothesis, briefs, error)

        try:
            prompt = build()
            catalogs = (
                {
                    card["gene"]: (set(card["primary_source_options"]), set(card["primary_variant_options"]))
                    for card in views.expand_prompt_payload(json.loads(prompt))[
                        "candidates_in_admission_order"
                    ]
                }
                if vcf_views is not None
                else None
            )
            return self._call(
                profile,
                "vcf_causal_shortlist" if vcf_views is not None else "global_evidence_ranking",
                prompt,
                lambda value: (
                    ranking.vcf_parse_global(value, admission.genes, catalogs)
                    if vcf_views is not None
                    else ranking.parse_global(value, admission.genes, briefs)
                ),
                repair=build,
                system=system,
            )
        except MendelPrimeError as error:
            return Ranking(
                genes=admission.genes,
                stage="global",
                status="degraded",
                issues=(error.code,),
            )

    def _direction(self, request, pair):
        try:
            choice = self.model.complete(request)["content"].strip().upper()
        except Exception:
            choice = ""
        return {
            "gene_a": pair[0],
            "gene_b": pair[1],
            "selected_gene": pair[0] if choice == "A" else pair[1] if choice == "B" else "",
            "status": "ok" if choice in {"A", "B"} else "invalid",
        }

    def _directions(self, requests, retries=True):
        jobs, output = tuple(requests), {}
        limits = (self.workers["pairwise"], PAIRWISE_RETRY_WORKERS) if retries else (len(jobs),)
        pending = jobs
        for limit in limits:
            if not pending:
                break
            output.update(
                self._parallel(
                    {pair: lambda pair=pair: self._direction(requests[pair], pair) for pair in pending},
                    limit,
                    lambda pair: {
                        "gene_a": pair[0],
                        "gene_b": pair[1],
                        "selected_gene": "",
                        "status": "invalid",
                    },
                )
            )
            pending = tuple(pair for pair in pending if output[pair]["status"] != "ok")
        return tuple(output[pair] for pair in jobs)

    def pairwise(
        self,
        profile,
        genes,
        admission,
        hypothesis,
        briefs,
        *,
        required_count=10,
        include_clinical=False,
        system_prompt=None,
        timeout_seconds=120,
        vcf_views=None,
    ):
        candidates = tuple(dict.fromkeys(gene.upper() for gene in genes if gene))
        if len(candidates) != required_count:
            return Ranking(
                genes=candidates,
                stage="pairwise",
                status="degraded",
                issues=(f"pairwise_rank_requires_exactly_{required_count}_candidates",),
            )
        requests = {}
        base = system_prompt or prompting.text("vcf_pairwise" if vcf_views is not None else "pairwise")
        cards = (
            ranking.vcf_pair_cards(profile, candidates, admission, briefs, vcf_views)
            if vcf_views is not None
            else None
        )
        for pair in (direction for item in combinations(candidates, 2) for direction in (item, item[::-1])):
            similar = any(
                (brief := ranking.brief_get(briefs, gene)) and brief["similar_case_evidence"] for gene in pair
            )
            system = base + (prompting.text("pairwise_similar_overlay") if similar else "")
            prompt = (
                ranking.vcf_pair_prompt(profile, *pair, cards)
                if vcf_views is not None
                else ranking.pair_prompt(profile, *pair, admission, hypothesis, briefs, include_clinical)
            )
            requests[pair] = self._request(
                profile,
                "vcf_structured_pairwise" if vcf_views is not None else "pairwise_ranking",
                prompt,
                system=system,
                timeout=120 if vcf_views is not None else timeout_seconds,
                cache=True,
            )
        return ranking.score_directions(
            candidates,
            self._directions(requests),
            bool(include_clinical and profile["clinical_description"].strip()),
        )

    def tail(
        self,
        profile,
        fixed,
        hypothesis,
        fusion,
        *,
        mode,
        system_prompt=None,
        native=None,
    ):
        inputs = tail.tail_inputs(fixed, hypothesis, mode)
        if mode == "genotype" and native is not None:
            if native.get("case_key") != profile["case_key"]:
                raise MendelPrimeError("case_mismatch", field="native")
            raw = {tool["tool_id"]: tool["ranking"].genes for tool in native["tools"]}
            llm = {
                f"{method}_llm": result.genes
                for method, result in hypothesis["sources"].items()
                if result.status == "ok"
            }
            rankings = {**raw, **llm}
            fixed_set = set(inputs["fixed"])
            candidates = {
                gene for values in rankings.values() for gene in values[:50] if gene not in fixed_set
            }
            reservoir = {
                gene for values in rankings.values() for gene in values[:100] if gene not in fixed_set
            }
            positions = {
                source: {gene: rank for rank, gene in enumerate(values, 1)}
                for source, values in rankings.items()
            }
            rows = tuple(
                {
                    "gene": gene,
                    **{f"{source}_llm_rank": None for source in tail.PHENOTYPE_SOURCES},
                    **{f"{source}_rank": values.get(gene) for source, values in positions.items()},
                    "source_ranks": {
                        **{f"{source}_llm": None for source in tail.PHENOTYPE_SOURCES},
                        **{source: values.get(gene) for source, values in positions.items()},
                    },
                    "top50_source_count": sum(gene in values[:50] for values in rankings.values()),
                    "variant_context": None,
                }
                for gene in sorted(candidates)
            )
            inputs.update(
                rankings=rankings,
                candidates=candidates,
                rrf_pool=reservoir,
                rows=rows,
                required=min(inputs["output_count"], len(rows)),
            )
            if not inputs["required"]:
                return tail.complete_tail(inputs, (), fusion)

        def build(error=""):
            if mode == "genotype" and native is not None:
                payload = {
                    "task": "Rank the supplied VCF-aware tail genes without changing the fixed prefix.",
                    "patient": ranking.vcf_patient(profile),
                    "fixed_prefix_excluded_from_task": inputs["fixed"],
                    "tail_candidates_alphabetical": inputs["rows"],
                    "required_output_count": inputs["required"],
                    "required_json_schema": {"ranked_tail_genes": ["GENE1", "GENE2"]},
                }
                prompt = ranking.vcf_dump(payload)
                return ranking.vcf_repair(prompt, error) if error else prompt
            return tail.tail_prompt(profile, inputs["fixed"], inputs["rows"], inputs["required"], error)

        system = system_prompt or prompting.text(
            "tail_ranking" if mode == "phenotype" else "vcf_tail_ranking"
        )
        try:
            genes = self._call(
                profile,
                "vcf_tail_ranking" if mode == "genotype" and native is not None else "tail_ranking",
                build(),
                lambda value: (
                    ranking.vcf_parse_tail(value, inputs["candidates"], inputs["required"])
                    if mode == "genotype" and native is not None
                    else tail.parse_tail(value, inputs["candidates"], inputs["required"])
                ),
                repair=build,
                system=system,
            )
            error = ""
        except MendelPrimeError as failure:
            genes, error = (), failure.code
        return tail.complete_tail(inputs, genes, fusion, error)

    def challenge_tail(self, profile, baseline, memory, *, maximum=5, fixed_top20=None):
        baseline, fixed = tuple(baseline[:100]), tuple(baseline[:20])
        if fixed_top20 is not None and fixed != tuple(fixed_top20[:20]):
            return tail.challenger_result(
                baseline, status="degraded", issue="fixed_top20_does_not_match_baseline"
            )
        evidence = {
            gene: cards
            for gene, cards in memory.get("tail_evidence_by_gene", {}).items()
            if gene not in set(fixed) and cards
        }
        if memory.get("tail_status") == "degraded":
            return tail.challenger_result(baseline, status="degraded", issue="tail_challenger_judge_degraded")
        if not evidence:
            return tail.challenger_result(baseline)
        prompt = tail.challenger_prompt(profile, fixed, baseline, evidence, maximum)

        try:
            decisions = self._call(
                profile,
                "tail_challenger_ranking_v1",
                prompt,
                lambda value: tail.parse_challengers(value, profile, baseline, evidence, maximum),
                repair=lambda error: self._repair(prompt, error),
            )
        except MendelPrimeError as error:
            return tail.challenger_result(baseline, status="degraded", issue=error.code)
        return tail.challenger_result(baseline, decisions)

    @staticmethod
    def _repair(prompt, error):
        value = json.loads(prompt)
        value["repair_instruction"] = {
            "previous_error": error,
            "instruction": "Return a complete replacement object satisfying every supplied constraint.",
        }
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))

    @staticmethod
    def _reflection(
        genes,
        reason,
        score,
        status="ok",
        triggered=False,
        gap=None,
        fill=None,
        updates=None,
        issue="",
        original=None,
    ):
        return Ranking(
            genes=genes,
            stage="reflection",
            status=status,
            details={
                "triggered": triggered,
                "trigger_reason": reason,
                "score_gap": score,
                "original_top2": original or genes,
                "gap_analysis": gap or {},
                "fill_assessment": fill or {},
                "evidence_updates": updates or {},
            },
            issues=(issue,) if issue else (),
        )

    def reflect(self, profile, pairwise, evidence, briefs, retrieve, *, excluded_pmids=()):
        genes, score = pairwise.genes[:2], ranking.pairwise_score_gap(pairwise)
        if pairwise.status != "ok":
            return self._reflection(genes, "pairwise_not_ok", score if score is not None else -1.0)
        if len(genes) != 2 or score is None:
            return self._reflection(genes, "pairwise_top2_or_scores_unavailable", -1.0)
        if round(score * 2) not in {0, 1}:
            return self._reflection(genes, "pairwise_score_gap_above_0_5", score)
        try:
            gap = self._call(
                profile,
                "top2_evidence_gap_analysis",
                ranking.gap_prompt(profile, genes, briefs),
                lambda value: ranking.parse_gap(value, genes, briefs),
            )
        except MendelPrimeError as error:
            return self._reflection(genes, "gap_analysis_failed", score, "degraded", True, issue=error.code)
        if not gap["retrieval_needed"]:
            return self._reflection(genes, "existing_evidence_sufficient", score, triggered=True, gap=gap)
        seen = tuple(
            document["prompt_id"].removeprefix("PMID:")
            for gene in genes
            for document in ranking.evidence_pool(evidence, gene)["documents"]
            if document["prompt_id"].startswith("PMID:")
        )
        try:
            retrieved = retrieve(
                genes, gap["queries_by_gene"], tuple(dict.fromkeys((*excluded_pmids, *seen)))
            )
        except Exception:
            return self._reflection(
                genes,
                "evidence_retrieval_failed",
                score,
                "degraded",
                True,
                gap=gap,
                issue="reflection_evidence_exception",
            )
        if any(ranking.evidence_pool(retrieved, gene)["retrieval_errors"] for gene in genes):
            return self._reflection(
                genes,
                "evidence_retrieval_failed",
                score,
                "degraded",
                True,
                gap=gap,
                issue="evidence_retrieval_failed",
            )
        try:
            fill = self._call(
                profile,
                "top2_evidence_gap_fill_assessment",
                ranking.fill_prompt(profile, genes, briefs, gap, retrieved),
                lambda value: ranking.parse_fill(value, genes, retrieved),
            )
        except MendelPrimeError as error:
            return self._reflection(
                genes, "evidence_admission_failed", score, "degraded", True, gap=gap, issue=error.code
            )
        updates = ranking.useful_evidence(genes, retrieved, fill["gene_results"])
        if not fill["comparison_discriminator_found"] or not any(
            item["documents"] for item in updates.values()
        ):
            return self._reflection(
                genes, "no_discriminating_evidence", score, triggered=True, gap=gap, fill=fill
            )
        requests = {
            pair: self._request(
                profile,
                "top2_reflection_pairwise",
                ranking.reflection_pair_prompt(profile, pair, briefs, gap, fill, updates),
            )
            for pair in (genes, genes[::-1])
        }
        compared = ranking.score_directions(
            genes, self._directions(requests, retries=False), stage="reflection_pairwise"
        )
        if compared.status != "ok":
            return self._reflection(
                genes,
                "pairwise_failed",
                score,
                "degraded",
                True,
                gap=gap,
                fill=fill,
                issue=compared.issues[0],
            )
        return self._reflection(
            compared.genes[:2],
            "reranked",
            score,
            triggered=True,
            gap=gap,
            fill=fill,
            updates=updates,
            original=genes,
        )
