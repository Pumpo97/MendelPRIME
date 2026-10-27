"""Read-only deployment checks; never install resources or call a model."""

import importlib
import json
import os
import re
import sys
import tempfile
import zipfile
from pathlib import Path
from urllib.parse import unquote, urlsplit

from mendelprime.domain import MendelPrimeError, require
from .process import run_process


def check_java_distribution(jar, expected_version=None):
    """Resolve the launcher's manifest class path exactly as Java does."""
    try:
        with zipfile.ZipFile(jar) as archive:
            raw = archive.read("META-INF/MANIFEST.MF").decode("utf-8")
            versions = {
                line.partition("=")[2].strip()
                for name in archive.namelist() if name.startswith("META-INF/maven/") and name.endswith("/pom.properties")
                for line in archive.read(name).decode("utf-8").splitlines() if line.startswith("version=")
            }
        lines = raw.replace("\r\n", "\n").replace("\n ", "").split("\n\n", 1)[0]
        manifest = dict(line.split(": ", 1) for line in lines.splitlines() if ": " in line)
        require(bool(manifest.get("Main-Class")), str(jar), "java_main_class_missing")
        versions.add(manifest.get("Implementation-Version", "unknown"))
        if expected_version is not None:
            require(expected_version in versions, str(jar), "java_tool_version_mismatch")
        for entry in manifest.get("Class-Path", "").split():
            require(not urlsplit(entry).scheme, str(jar), "remote_java_dependency_not_supported")
            path = jar.parent / unquote(entry)
            require(path.is_file() and path.stat().st_size and os.access(path, os.R_OK),
                    str(path), "java_dependency_missing")
        return expected_version or manifest.get("Implementation-Version", "unknown")
    except (OSError, KeyError, UnicodeError, zipfile.BadZipFile):
        raise MendelPrimeError("java_manifest_unreadable", field=str(jar)) from None


def _probe(command, expected=None, timeout=60, work_root=None):
    with tempfile.TemporaryDirectory(prefix="mendelprime-doctor-", dir=work_root) as work:
        root = Path(work)
        result = run_process(
            command=[str(part).replace("{workdir}", work) for part in command],
            root=root, timeout_seconds=timeout,
        )
        require(not result["error_code"], "tool", "tool_startup_failed")
        output = "\n".join((root / name).read_text(errors="replace") for name in ("stdout.log", "stderr.log"))
    if expected:
        require(re.search(r"(?<![\d.])" + re.escape(expected) + r"(?!\d)", output),
                "tool", "tool_version_mismatch")


def doctor(config, profile, *, check_tools=False, work_root=None, portable_root=None):
    checks = []
    similar = config.get("ranking", {}).get("similar_cases", True)

    def check(name, action):
        try:
            action()
            checks.append({"check": name, "status": "ok"})
        except Exception as error:
            checks.append({"check": name, "status": "failed",
                           "code": error.code if isinstance(error, MendelPrimeError) else type(error).__name__,
                           "field": error.field if isinstance(error, MendelPrimeError) else name})

    def exists(path, field):
        require(isinstance(path, Path) and path.exists() and os.access(path, os.R_OK),
                field, "resource_missing_or_unreadable")
        if portable_root is not None:
            require(path.resolve().is_relative_to(portable_root.resolve()), field, "resource_outside_portable_root")

    def paths(value, prefix):
        if isinstance(value, Path):
            check(prefix, lambda: exists(value, prefix))
        elif isinstance(value, dict):
            for name, item in value.items():
                paths(item, prefix + "." + name)
        elif isinstance(value, (tuple, list)):
            for index, item in enumerate(value):
                paths(item, f"{prefix}[{index}]")

    check("python", lambda: require(sys.version_info[:2] == (3, 12), "python", "python_3_12_required"))
    if config.get("llm"):
        check("llm.endpoint", lambda: require(
            not (urlsplit(config["llm"]["base_url"]).hostname or "").endswith(".invalid"),
            "llm.base_url", "replace_example_endpoint",
        ))
    modules = (("yaml", "numpy", "torch", "transformers") if similar else ("yaml",)) if profile == "hpo" else (
        "yaml", "pysam", "pyfaidx", "hgvs", "bioutils", "bs4", "psycopg2"
    )
    for name in modules:
        check("python." + name, lambda name=name: importlib.import_module(name))
    resources = config.get("resources") or {}
    paths({key: value for key, value in resources.items()
           if (profile != "vcf" or key in {"hpo_obo", "hpo_stats", "omim_jsonl", "mim2gene"})
           and (similar or key not in {"similar_case_root", "embedding_model", "cross_encoder"})}, "resources")
    if any(row["status"] == "failed" and row["check"].startswith("python.") for row in checks):
        return {"profile": profile, "checks": checks, "exit_code": 1,
                "warnings": ["Install the selected Python dependency group before checking native resources."]}
    if profile == "hpo":
        from .runner import PHENOTYPE_RESOURCES, PHEN2DISEASE_FILES

        settings = config.get("resources") or {}
        for name in PHENOTYPE_RESOURCES if similar else PHENOTYPE_RESOURCES[:2]:
            check("resources." + name, lambda name=name: exists(settings.get(name), name))
        if settings.get("phen2disease_root"):
            for relative in PHEN2DISEASE_FILES:
                path = settings["phen2disease_root"] / relative
                check(str(path), lambda path=path: exists(path, str(path)))
        if similar and settings.get("similar_case_root"):
            from mendelprime.knowledge.retrieval import load_index

            check("similar_case_index", lambda: load_index(settings["similar_case_root"]))
        for name in ("embedding_model", "cross_encoder") if similar else ():
            if settings.get(name):
                check(name, lambda name=name: _model_files(settings[name]))
        if similar and config["execution"]["device"] == "cuda":
            check("cuda", lambda: require(importlib.import_module("torch").cuda.is_available(),
                                          "execution.device", "cuda_unavailable"))
            check("cuda.bfloat16", lambda: require(importlib.import_module("torch").cuda.is_bf16_supported(),
                                                   "execution.device", "cuda_bfloat16_unavailable"))
    else:
        from mendelprime.tools.native import _validate_resources

        native = config.get("genotype")
        check("genotype", lambda: require(native is not None, "genotype", "required"))
        if native:
            paths(native, "genotype")
            for entry in native["assemblies"]:
                for tool in ("exomiser", "lirical"):
                    check(tool + "." + entry["assembly"],
                          lambda tool=tool, entry=entry: _validate_resources(config["genotype"], entry, tool))
            if check_tools:
                check("java.startup", lambda: _probe([native["java"], "-version"], "17"))
                for tool in ("exomiser", "lirical"):
                    check(tool + ".startup", lambda tool=tool: _probe(
                        [native["java"], "-jar", native[tool + "_jar"], "--version"],
                        native[tool + "_version"] if tool == "lirical" else None,
                    ))
        paths(config.get("variant_knowledge"), "variant_knowledge")
        settings = config.get("explanation_resources")
        check("explanation_resources", lambda: require(settings is not None, "explanation_resources", "required"))
        if settings:
            from mendelprime.knowledge.assessment_resources import (
                _check_resources, _load_clingen_guidance, _load_protocol, resource_settings,
            )
            reporting = config["output"]["explanation"]
            settings = resource_settings(settings, report=reporting)
            vep = reporting or settings["ranking_vep_enabled"]
            paths(settings, "explanation_resources")
            if vep and settings.get("vep_command"):
                check("vep.executable", lambda: exists(Path(settings["vep_command"][0]), "vep_command[0]"))
            check("assessment_resources", lambda: _check_resources(settings, vep_enabled=vep, report=reporting))
            if portable_root is not None:
                check("autopvs1.portable_paths", lambda: _portable_pvs1(settings["autopvs1_root"], exists))
            check("ai_cura", lambda: _load_protocol(settings, report=reporting))
            if reporting:
                check("clingen_guidance", lambda: _load_clingen_guidance(settings))
                for name in ("clingen_dosage_table", "genereviews_index", "post_evidence_catalog", "kinship_executable"):
                    if settings.get(name):
                        check(name, lambda name=name: exists(settings[name], name))
            if check_tools:
                from mendelprime.tools.pvs1_driver import _CHILD_BOOTSTRAP

                check("autopvs1.startup", lambda: _probe([
                    sys.executable, "-I", "-B", "-c", _CHILD_BOOTSTRAP,
                    str(Path(__file__).resolve().parents[2]), "--check-source", settings["autopvs1_root"],
                ], timeout=300))
                if vep:
                    check("vep.startup", lambda: _probe([*settings["vep_command"], "--help"],
                                                       settings["vep_tool_version"], timeout=300,
                                                       work_root=work_root))
                if config["output"]["explanation"] and settings.get("kinship_executable"):
                    check("plink2.startup", lambda: _probe([settings["kinship_executable"], "--version"]))
    warnings = ["Network services and API credentials were not tested; no model or patient request was sent."]
    if profile == "hpo":
        for name in ("phen2gene", "amelie"):
            if not config["external"][name + "_acknowledged"]:
                warnings.append(name + " external upload is not acknowledged; full cold ranking will be degraded.")
    return {"profile": profile, "checks": checks, "warnings": warnings,
            "exit_code": int(any(item["status"] == "failed" for item in checks))}


def _model_files(root):
    from transformers import AutoConfig, AutoTokenizer

    AutoConfig.from_pretrained(str(root), local_files_only=True, trust_remote_code=False)
    AutoTokenizer.from_pretrained(str(root), local_files_only=True, trust_remote_code=False)
    indexes = list(root.glob("*.index.json"))
    weights = list(root.glob("*.safetensors")) + list(root.glob("pytorch_model*.bin"))
    require(bool(weights), str(root), "model_weights_missing")
    for index in indexes:
        for relative in set(json.loads(index.read_text())["weight_map"].values()):
            require((root / relative).is_file(), str(root / relative), "model_shard_missing")


def _portable_pvs1(root, exists):
    import configparser

    from mendelprime.tools.pvs1_driver import _configured_path, _required_resource_paths

    config = configparser.ConfigParser()
    config.read(root / "config.ini")
    for path in _required_resource_paths(root, config):
        exists(path, str(path))
    for build in ("HG19", "HG38"):
        path = _configured_path(root, config[build]["genome"])
        if build == "HG19" or path.exists():
            exists(path, str(path))
            exists(Path(str(path) + ".fai"), str(path) + ".fai")
