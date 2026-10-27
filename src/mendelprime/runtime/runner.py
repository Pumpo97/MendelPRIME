"""Own shared services and schedule inference or saved-result explanation once."""

import csv
import json
import logging
import os
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path
from uuid import uuid4

from mendelprime.cemr.ranker import CEMR
from mendelprime.domain import MendelPrimeError, Patient, require, validate_case_key, validate_gene
from mendelprime.harness.pipeline import Harness
from mendelprime.knowledge.base import KnowledgeBase
from mendelprime.tools.phenotype import require_positive_hpo
from mendelprime.tools.service import MedicalTools

from .io import prepare_output_directory, publish_case, read_public
from .llm import LLMClient, RequestBudget
from .storage import atomic_replace_bytes, private_directory, stat_fingerprint

LOG = logging.getLogger("mendelprime.runner")
PHENOTYPE_RESOURCES = (
    "genes_to_phenotype",
    "phen2disease_root",
    "similar_case_root",
    "embedding_model",
    "cross_encoder",
)
PHEN2DISEASE_FILES = (
    "data/matrix/lin_similarity_matrix.json",
    "data/association/Disease-Hpo/disease2hpo20210413.json",
    "data/association/Gene-Disease/disease2genecard2021.json",
    "data/obo/hp20210413.obo",
)


def result_paths(root):
    root = Path(root)
    if root.is_file():
        return (root,)
    require(root.is_dir(), "results", "directory_or_file_required")
    if (root / "result.json").is_file():
        return (root / "result.json",)
    found = []
    for parent, directories, files in os.walk(root, followlinks=False):
        directories[:] = sorted(
            name
            for name in directories
            if not name.startswith(".") and not (Path(parent) / name).is_symlink()
        )
        if "result.json" in files:
            found.append(Path(parent) / "result.json")
        require(len(found) <= 10000, "results", "too_many_case_results")
    require(bool(found), "results", "no_case_results")
    return tuple(sorted(found))


def check_targets(root, keys):
    require(isinstance(root, Path) and root.is_absolute(), "output_root", "absolute_path_required")
    require(bool(keys) and len(set(keys)) == len(keys), "cases", "nonempty_unique_case_keys_required")
    prepare_output_directory(root)
    for key in keys:
        validate_case_key(key)
        parent = root / key[0]
        prepare_output_directory(parent)
        require(not os.path.lexists(parent / key[1]), "output", "case_already_exists")
        require(
            not os.path.lexists(parent / f".{key[1]}.publishing"), "output", "case_publication_in_progress"
        )


def _resource(path, field, *, directory=False, code=None):
    try:
        valid = isinstance(path, Path) and (path.is_dir() if directory else path.is_file())
        valid = valid and os.access(path, os.R_OK | (os.X_OK if directory else 0))
    except OSError:
        valid = False
    require(valid, field, code or ("directory_required" if directory else "regular_file_required"))


def _preflight(config, cases):
    require(
        type(cases) is tuple and bool(cases) and all(isinstance(case, Patient) for case in cases),
        "cases",
        "nonempty_case_tuple_required",
    )
    require(len({case.key for case in cases}) == len(cases), "cases", "duplicate_case_key")
    for case in cases:
        require_positive_hpo(case)
    phenotype = any(case.genomic is None for case in cases)
    genotype = any(case.genomic is not None for case in cases)
    assemblies = tuple(
        dict.fromkeys(
            case.genomic["assembly"]
            for case in cases
            if case.genomic and case.genomic["assembly"] is not None
        )
    )
    resources = config["resources"]
    require(resources is not None, "resources", "required_for_application")
    similar = config["ranking"].get("similar_cases", True)
    phenotype_resources = PHENOTYPE_RESOURCES if similar else PHENOTYPE_RESOURCES[:2]
    if phenotype:
        for name in phenotype_resources:
            require(resources[name] is not None, f"resources.{name}", "required_for_phenotype_application")
    if genotype:
        for name in ("genotype", "variant_knowledge"):
            require(config[name] is not None, name, "required_for_vcf_application")
        available = {value["assembly"] for value in config["genotype"]["assemblies"]}
        require(
            all(value in available for value in assemblies),
            "genotype.assemblies",
            "assembly_resources_not_configured",
        )
    if genotype:
        _assessment_settings(config, assemblies)
    for name in ("hpo_obo", "hpo_stats", "omim_jsonl", "mim2gene"):
        _resource(resources[name], f"resources.{name}")
    if phenotype:
        _resource(resources["genes_to_phenotype"], "resources.genes_to_phenotype")
        for name in phenotype_resources[1:]:
            _resource(resources[name], f"resources.{name}", directory=True)
        for relative in PHEN2DISEASE_FILES:
            _resource(
                resources["phen2disease_root"] / relative,
                "resources.phen2disease_root",
                code="required_layout_missing",
            )
        if similar:
            root = resources["similar_case_root"] / "rarea"
            _resource(root, "resources.similar_case_root", directory=True, code="required_layout_missing")
            for name in ("records.jsonl", "record_ids.json", "embeddings.npy"):
                _resource(root / name, "resources.similar_case_root", code="required_layout_missing")
        if config["external"]["ca_bundle"] is not None:
            _resource(config["external"]["ca_bundle"], "external.ca_bundle")
    if genotype:
        _variant_paths(config)
    return phenotype, genotype


def _variant_paths(config):
    require(config["variant_knowledge"] is not None, "variant_knowledge", "required_for_vcf_application")
    for name, path in config["variant_knowledge"].items():
        _resource(path, f"variant_knowledge.{name}")


def _assessment_settings(config, assemblies):
    settings = config["explanation_resources"]
    require(settings is not None, "explanation_resources", "required_for_vcf_evidence")
    require(
        all(assembly == settings["assembly"] for assembly in assemblies),
        "explanation_resources",
        "assembly_mismatch",
    )


class Runner:
    """Keep configuration, global budgeting and resource lifetime outside science."""

    def __init__(self, config, cache_root, *, on_call=None, responses=None, on_stage=None, audit_root=None):
        self.config, self.cache_root = config, Path(cache_root).absolute()
        self.on_call, self.responses, self.on_stage = on_call, responses, on_stage
        self.audit_root = Path(audit_root).absolute() if audit_root is not None else None
        self.model = self.knowledge = self.resources = None

    def _model(self):
        if self.model is None:
            execution = self.config["execution"]
            private_directory(self.cache_root)
            choices = private_directory(self.cache_root / "choices") if execution["cache_choices"] else None
            if self.audit_root is not None:
                private_directory(self.audit_root)

            def observe(event):
                self._audit("timing", event)
                (self.on_call or _log_call)(event)

            def respond(request, reply):
                self._audit("response", {"request": request, "reply": reply})
                if self.responses is not None:
                    self.responses(request, reply)

            self.model = LLMClient(
                self.config["llm"],
                RequestBudget(execution["max_llm_concurrency"]),
                execution["queue_timeout_seconds"],
                cache_root=choices,
                observer=observe,
                responses=respond,
            )
        return self.model

    def _audit(self, kind, value):
        if self.audit_root is not None:
            # Only validated messages/replies and sanitized timing events enter
            # this optional private sink, never HTTP headers or credentials.
            data = json.dumps(value, ensure_ascii=False, allow_nan=False).encode("utf-8")
            atomic_replace_bytes(self.audit_root / f"{kind}-{uuid4().hex}.json", data, mode=0o600)

    def _knowledge(self):
        if self.knowledge is None:
            self.knowledge = KnowledgeBase(
                {name: self.config[name] for name in ("resources", "variant_knowledge")},
                self.cache_root,
                self.config["execution"]["device"],
            )
        return self.knowledge

    def _assessment_resources(self, assemblies, *, purpose="report"):
        if self.resources is None and assemblies:
            from mendelprime.knowledge.assessment_resources import open_resources

            _assessment_settings(self.config, assemblies)
            _variant_paths(self.config)
            root = private_directory(self.cache_root / "explanation")
            self.resources = open_resources(
                self.config["explanation_resources"],
                self._knowledge().variant_lookup,
                root,
                purpose=purpose,
                workers=self.config["workers"]["literature"],
            )
        return self.resources

    def run(self, cases, output_root):
        from .explanation import explain, freeze_context
        from .results import project_result

        started = time.monotonic()
        output_root = Path(output_root).absolute()
        require(
            self.audit_root is None or not self.audit_root.is_relative_to(output_root),
            "audit_root",
            "audit_must_be_outside_formal_results",
        )
        check_targets(output_root, tuple(case.key for case in cases))
        phenotype, genotype = _preflight(self.config, cases)
        private_directory(self.cache_root)
        for relative in ("pubmed", "pubmed/search", "pubmed/articles"):
            private_directory(self.cache_root / relative)
        if phenotype:
            private_directory(self.cache_root / "amelie")
            private_directory(self.cache_root / "amelie/by_request")
        if genotype:
            private_directory(self.cache_root / "native")
        try:
            knowledge, model = self._knowledge(), self._model()
            LOG.info("Loading shared HPO and evidence resources")
            knowledge.hpo()
            knowledge.omim()
            tools = MedicalTools(
                knowledge,
                self.cache_root,
                self.config["external"],
                phenotype=phenotype,
                genotype=self.config["genotype"],
                work_root=self.cache_root / "native",
            )
            reporting = self.config["output"]["explanation"]
            if genotype:
                assemblies = tuple(
                    case.genomic["assembly"]
                    for case in cases
                    if case.genomic and case.genomic["assembly"] is not None
                )
                self._assessment_resources(
                    assemblies or (self.config["explanation_resources"]["assembly"],),
                    purpose="report" if reporting else "ranking",
                )
            harness = Harness(
                tools,
                knowledge,
                CEMR(model, workers={k: v for k, v in self.config["workers"].items() if k != "post_acmg"}),
                reflection=self.config["ranking"]["reflection"],
                similar=self.config["ranking"].get("similar_cases", True),
                ranking_readers=self.resources["ranking_actions"] if genotype else None,
            )
            LOG.info("Shared application ready")

            def progress(event):
                self._audit("stage", event)
                if self.on_stage is not None:
                    self.on_stage(event)

            def execute(case):
                run = harness.run(
                    case,
                    progress=progress if self.audit_root is not None or self.on_stage is not None else None,
                )
                explanation = None
                if reporting:
                    try:
                        context = freeze_context(run)
                    except MendelPrimeError as error:
                        explanation = {"status": "degraded", "issues": (error.code,), "markdown": ""}
                    else:
                        explanation = explain(
                            model, context, self.resources, acmg_workers=self.config["workers"]["post_acmg"]
                        )
                payload = project_result(run, explanation)
                markdown = explanation.get("markdown", "") if explanation else ""
                return payload["status"], publish_case(output_root, case.key, payload, markdown)

            return self._batch(tuple((case.key, case) for case in cases), execute, started)
        finally:
            self.close()

    def explain(
        self,
        paths,
        output_root,
        cohort=None,
        *,
        rerender_report=False,
        recompute_assessment=False,
        case_inputs=(),
    ):
        from .explanation import explain, load_context, load_saved
        from .results import with_explanation

        started, jobs, assemblies = time.monotonic(), [], []
        original_cases = {case.key: case for case in case_inputs}
        require(len(original_cases) == len(case_inputs), "case_inputs", "duplicate_case")
        vcf_fingerprints = {}
        output_root = Path(output_root).absolute()
        require(
            self.audit_root is None or not self.audit_root.is_relative_to(output_root),
            "audit_root",
            "audit_must_be_outside_formal_results",
        )
        for raw in paths:
            path = Path(raw).absolute()
            path = path / "result.json" if path.is_dir() else path
            fingerprint = stat_fingerprint(path.stat())
            key, value = read_public(path, cohort)
            context = load_context(value, key, **({"for_assessment": True} if recompute_assessment else {}))
            saved = load_saved(value, context)
            if original_cases:
                require(
                    key in original_cases and context.get("review") is not None,
                    "case_inputs",
                    "matching_vcf_assessment_required",
                )
                case, seed = original_cases[key], context["assessment_seed"]
                require(
                    case.genomic is not None and case.genomic["assembly"] == seed["family"]["assembly"],
                    "case_inputs",
                    "assembly_mismatch",
                )
                require(
                    {(t["hpo_id"], t["present"]) for t in case.phenotypes}
                    == {(t["hpo_id"], t["present"]) for t in seed["profile"]["terms"]},
                    "case_inputs",
                    "phenotype_mismatch",
                )
                _resource(case.genomic["vcf_path"], "case_inputs.vcf_path")
                vcf_fingerprints[key] = stat_fingerprint(case.genomic["vcf_path"].stat())
            require(stat_fingerprint(path.stat()) == fingerprint, "result", "file_changed_during_preflight")
            needs_model = (
                recompute_assessment
                or rerender_report
                or saved.get("report") is None
                or saved.get("assessment_progress") is not None
            )
            if (
                needs_model
                and context["assessment_seed"]
                and (recompute_assessment or saved.get("assessment") is None)
            ):
                assemblies.append(context["assessment_seed"]["family"]["assembly"])
            jobs.append((key, (path, fingerprint, needs_model, key)))
        require(
            not original_cases or set(original_cases) == {key for key, _ in jobs},
            "case_inputs",
            "case_set_mismatch",
        )
        check_targets(output_root, tuple(key for key, _ in jobs))
        try:
            if any(job[2] for _, job in jobs):
                if assemblies:
                    _assessment_settings(self.config, assemblies)
                    _variant_paths(self.config)
                self._model()
                if assemblies:
                    self._assessment_resources(tuple(dict.fromkeys(assemblies)))

            def execute(job):
                path, fingerprint, _, expected_key = job
                require(
                    stat_fingerprint(path.stat()) == fingerprint, "result", "file_changed_after_preflight"
                )
                key, value = read_public(path, cohort)
                context = load_context(
                    value, key, **({"for_assessment": True} if recompute_assessment else {})
                )
                if key in original_cases and context.get("review") is not None:
                    case = original_cases[key]
                    require(
                        stat_fingerprint(case.genomic["vcf_path"].stat()) == vcf_fingerprints[key],
                        "case_inputs",
                        "vcf_changed_after_preflight",
                    )
                    context["review"]["vcf_path"] = str(case.genomic["vcf_path"])
                saved = load_saved(value, context)
                if recompute_assessment:
                    saved = {}
                    for field in ("final_acmg", "assessment_progress"):
                        value["stages"].get("genotype", {}).pop(field, None)
                require(
                    key == expected_key and stat_fingerprint(path.stat()) == fingerprint,
                    "result",
                    "file_changed_after_preflight",
                )
                if (
                    recompute_assessment
                    or rerender_report
                    or saved.get("report") is None
                    or saved.get("assessment_progress") is not None
                ):
                    saved = explain(
                        self.model,
                        context,
                        self.resources,
                        saved_assessment=saved.get("assessment"),
                        acmg_workers=self.config["workers"]["post_acmg"],
                        **(
                            {"saved_progress": saved["assessment_progress"]}
                            if saved.get("assessment_progress") is not None
                            else {}
                        ),
                    )
                payload = with_explanation(value, key, context, saved)
                require(
                    stat_fingerprint(path.stat()) == fingerprint, "result", "file_changed_after_preflight"
                )
                return payload["status"], publish_case(output_root, key, payload, saved["markdown"])

            return self._batch(tuple(jobs), execute, started)
        finally:
            self.close()

    def _batch(self, jobs, action, started):
        results, iterator = {}, iter(enumerate(jobs))
        with ThreadPoolExecutor(
            max_workers=self.config["execution"]["lanes"], thread_name_prefix="mendelprime-case"
        ) as pool:
            pending = {}

            def submit():
                item = next(iterator, None)
                if item is not None:
                    index, (key, job) = item
                    pending[pool.submit(_case, key, action, job)] = index

            for _ in range(min(self.config["execution"]["lanes"], len(jobs))):
                submit()
            try:
                while pending:
                    done, _ = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done:
                        results[pending.pop(future)] = future.result()
                    for _ in done:
                        submit()
            except BaseException:
                LOG.warning("Dispatch stopped; waiting for active cases before releasing shared resources")
                for future in pending:
                    future.cancel()
                raise
        outcomes = tuple(results[index] for index in range(len(jobs)))
        states = {row["status"] for row in outcomes}
        return {
            "outcomes": outcomes,
            "elapsed_seconds": time.monotonic() - started,
            "exit_code": 1 if "failed" in states else 2 if "degraded" in states else 0,
        }

    def close(self):
        self.model = None
        if self.resources is not None:
            from mendelprime.knowledge.assessment_resources import close_resources

            close_resources(self.resources)
            self.resources = None
        if self.knowledge is not None:
            self.knowledge.close()
            self.knowledge = None


def _case(key, action, job):
    started = time.monotonic()
    LOG.info("[%s/%s] case started", *key)
    status, output, code = "failed", None, ""
    try:
        status, output = action(job)
    except Exception as error:
        code = (
            error.code
            if isinstance(error, MendelPrimeError)
            else ("case_io_failed" if isinstance(error, OSError) else "case_execution_failed")
        )
        LOG.error("[%s/%s] case failed code=%s type=%s", *key, code, type(error).__name__)
    elapsed = time.monotonic() - started
    LOG.info("[%s/%s] case finished status=%s elapsed=%.2fs", *key, status, elapsed)
    return {
        "key": dict(zip(("cohort", "case_id"), key, strict=True)),
        "status": status,
        "output": str(output) if output else None,
        "elapsed_seconds": elapsed,
        "error_code": code,
    }


def _log_call(event):
    LOG.info(
        "[%s/%s] model stage=%s status=%s elapsed=%.2fs",
        *event["case_key"],
        event["stage_id"],
        event["status"],
        event["elapsed_seconds"],
    )


def evaluate(paths, truth_path, cohort=None):
    """Gold stays outside inference and public packages; missing ranks count as misses."""
    with Path(truth_path).open(encoding="utf-8-sig", newline="") as stream:
        truth = {}
        for row in csv.DictReader(stream):
            key = (row.get("cohort") or cohort or "", row.get("case_id", ""))
            require(bool(key[1]) and key not in truth, "truth", "unique_case_key_required")
            raw = row.get("gold_genes", row.get("gold_gene", row.get("gene", "")))
            genes = json.loads(raw) if raw.startswith("[") else raw.replace(";", ",").split(",")
            require(isinstance(genes, list) and bool(genes), "truth", "gold_genes_required")
            for gene in genes:
                validate_gene(gene.strip(), "truth.gold_genes")
            truth[key] = tuple(dict.fromkeys(gene.strip() for gene in genes))
    require(bool(truth), "truth", "gold_genes_required")
    records, seen = [], set()
    for path in paths:
        key, value = read_public(path, cohort)
        require(key not in seen, "results", "duplicate_case_key")
        seen.add(key)
        label_key = key if key in truth else ("", key[1])
        require(label_key in truth, "truth", "case_truth_missing")
        require(not any(row["truth_key"] == label_key for row in records), "truth", "ambiguous_case_id")
        gold, rankings = truth[label_key], _stage_rankings(value["stages"])
        rankings["final"] = value["final"].get("top100", value["final"].get("top50", ()))
        records.append(
            {
                "key": {"cohort": key[0], "case_id": key[1]},
                "truth_key": label_key,
                "status": value["status"],
                "gold_genes": gold,
                "ranks": {
                    name: min((genes.index(gene) + 1 for gene in gold if gene in genes), default=None)
                    for name, genes in rankings.items()
                },
            }
        )
    methods = sorted({method for row in records for method in row["ranks"]})

    def metrics(rows):
        result = {}
        for method in methods:
            ranks = [row["ranks"].get(method) for row in rows]
            count = len(ranks)
            result[method] = {
                "n": count,
                "available": sum(rank is not None for rank in ranks),
                "MRR": sum(1 / rank for rank in ranks if rank is not None) / count if count else 0,
                **{
                    f"R@{k}": sum(rank is not None and rank <= k for rank in ranks) / count if count else 0
                    for k in (1, 3, 5, 10, 20, 50, 100)
                },
            }
        return result

    return {
        "cases": [{key: value for key, value in row.items() if key != "truth_key"} for row in records],
        "overall": metrics(records),
        "cohorts": {
            name: metrics([row for row in records if row["key"]["cohort"] == name])
            for name in sorted({row["key"]["cohort"] for row in records})
        },
        "degraded_cases": sum(row["status"] == "degraded" for row in records),
    }


def _stage_rankings(stages):
    output = {}

    def collect(value, path):
        if not isinstance(value, dict):
            return
        genes = value.get("genes")
        if genes is None and isinstance(value.get("ranking"), list):
            genes = value["ranking"]
        if genes is None and path.endswith("candidate_admission"):
            genes = value.get("ranked_pool_genes", value.get("candidate_genes"))
        if isinstance(genes, list) and all(isinstance(item, (str, dict)) for item in genes):
            resolved = [item if isinstance(item, str) else item.get("gene") for item in genes]
            if all(isinstance(gene, str) for gene in resolved):
                output[path.removesuffix(".ranking")] = resolved
        for name, item in value.items():
            if isinstance(item, dict):
                collect(item, f"{path}.{name}" if path else name)

    collect(stages, "")
    return output
