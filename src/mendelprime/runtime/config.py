"""Strict loading for the run configuration and optional resource manifest."""

from __future__ import annotations

import copy
import math
import re
import tomllib
from pathlib import Path

from mendelprime.domain import MendelPrimeError, require
from mendelprime.runtime.llm import LLM_DEFAULTS, _validated_settings
from mendelprime.runtime.storage import read_regular_bounded

_MAX_CONFIG_BYTES = 64 * 1024
_ENVIRONMENT_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*", re.ASCII)
_RELEASE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*", re.ASCII)
_RESOURCE_SECTIONS = {"resources", "genotype", "variant_knowledge", "explanation_resources"}

# Schema order defines first-error field priority. Fields in
# resource sections belong to resources.toml unless explicitly marked run=True.
_SCHEMA = {
    "llm": {
        "base_url": {"required": True},
        "model": {"default": LLM_DEFAULTS["model"]},
        "api_key_env": {"default": LLM_DEFAULTS["api_key_env"]},
        "wire_api": {"default": LLM_DEFAULTS["wire_api"]},
        "transport_retries": {"default": LLM_DEFAULTS["transport_retries"]},
    },
    "execution": {
        "lanes": {"default": 1},
        "max_llm_concurrency": {"default": 8},
        "queue_timeout_seconds": {"default": 1800},
        "cache_choices": {"default": True},
        "device": {"default": "cuda"},
    },
    "workers": {
        "sources": {"default": 4},
        "briefs": {"default": 8},
        "pairwise": {"default": 8},
        "literature": {"default": 4},
        "post_acmg": {"default": 4},
    },
    "resources": {
        "hpo_obo": {"required": True, "path": True},
        "hpo_stats": {"required": True, "path": True},
        "omim_jsonl": {"required": True, "path": True},
        "mim2gene": {"required": True, "path": True},
        "genes_to_phenotype": {"default": None, "path": True},
        "phen2disease_root": {"default": None, "path": True},
        "similar_case_root": {"default": None, "path": True},
        "embedding_model": {"default": None, "path": True},
        "cross_encoder": {"default": None, "path": True},
    },
    "genotype": {
        "java": {"required": True, "path": True},
        "exomiser_jar": {"required": True, "path": True},
        "exomiser_properties": {"required": True, "path": True},
        "lirical_jar": {"required": True, "path": True},
        "lirical_data": {"required": True, "path": True},
        "gene_map": {"required": True, "path": True},
        "assemblies": {"required": True},
        "exomiser_version": {"default": "14.0.0"},
        "lirical_version": {"default": "2.4.0"},
        "lirical_transcript_db": {"default": "REFSEQ", "run": True},
        "max_heap": {"default": "12g", "run": True},
        "timeout_seconds": {"default": 1800.0, "run": True},
    },
    "variant_knowledge": {
        "clinvar_db": {"required": True, "path": True},
        "clingen_validity": {"required": True, "path": True},
        "clingen_vcep": {"required": True, "path": True},
        "clingen_manifest": {"required": True, "path": True},
    },
    "explanation_resources": {
        "assembly": {"required": True},
        "ranking_vep_enabled": {"default": False, "run": True},
        "reference_fasta": {"required": True, "path": True},
        "vep_command": {},
        "vep_cache": {"path": True},
        "vep_fasta": {"path": True},
        "vep_cache_version": {"version": True},
        "vep_tool_version": {"version": True},
        "clinvar_protein_disease_db": {"required": True, "path": True},
        "clinvar_protein_release": {"required": True, "version": True},
        "gnomad_sites": {"required": True, "path": True},
        "gnomad_coverage": {"required": True, "path": True},
        "population_release": {"required": True, "version": True},
        "gnomad_constraint": {"required": True, "path": True},
        "constraint_release": {"required": True, "version": True},
        "spliceai_snv": {"required": True, "path": True},
        "spliceai_indel": {"required": True, "path": True},
        "spliceai_release": {"required": True, "version": True},
        "revel": {"required": True, "path": True},
        "revel_release": {"required": True, "version": True},
        "mane_summary": {"path": True},
        "mane_release": {"version": True},
        "autopvs1_root": {"required": True, "path": True},
        "ai_cura_protocol": {"path": True},
        "ai_cura_workbook": {"path": True},
        "ai_cura_knowledgebase": {"required": True, "path": True},
        "ai_cura_manifest": {"required": True, "path": True},
        "clingen_guideline_manifest": {"path": True},
        "clingen_dosage_table": {"path": True},
        "genereviews_index": {"path": True},
        "post_evidence_catalog": {"path": True},
        "kinship_executable": {"path": True},
        "vep_timeout_seconds": {"default": 300, "run": True},
        "pvs1_timeout_seconds": {"default": 300, "run": True},
    },
    "external": {
        "phen2gene_acknowledged": {"default": False},
        "amelie_acknowledged": {"default": False},
        "amelie_verify_tls": {"default": True},
        "amelie_pseudonym_secret_env": {"default": "MENDELPRIME_AMELIE_PSEUDONYM_SECRET"},
        "ca_bundle": {"default": None, "path": True},
    },
    "ranking": {"reflection": {"default": False}, "similar_cases": {"default": True}},
    "output": {"explanation": {"default": False}},
}
_ASSEMBLY_SCHEMA = {
    "assembly": {},
    "exomiser_release": {"release": True},
    "phenotype_release": {"release": True},
    "exomiser_genome_data": {"path": True},
    "exomiser_phenotype_data": {"path": True},
    "lirical_exomiser_data": {"path": True},
    "lirical_variants_db": {"path": True},
    "lirical_clinvar_db": {"path": True},
}


def _read(path: Path, field: str) -> dict:
    data = None
    try:
        raw, _before, _after = read_regular_bounded(path, _MAX_CONFIG_BYTES)
        data = tomllib.loads(raw.decode("utf-8"))
    except (OSError, OverflowError, UnicodeError, ValueError):
        pass
    if not isinstance(data, dict):
        raise MendelPrimeError("unable_to_read_toml", field=field) from None
    return data


def _table(value: object, name: str, allowed: set[str]) -> dict:
    require(isinstance(value, dict), name, "table_required")
    require(not set(value) - allowed, name, "unknown_key")
    return dict(value)


def _path(value: object, field: str, base: Path) -> Path:
    require(
        isinstance(value, str) and bool(value.strip()) and "\x00" not in value, field, "path_text_required"
    )
    path = Path(value)
    return path if path.is_absolute() else (base / path).resolve()


def _positive(value: object, field: str, *, integer: bool = False) -> None:
    valid_type = type(value) is int if integer else type(value) in {int, float}
    code = "positive_integer_required" if integer else "positive_finite_number_required"
    require(valid_type and math.isfinite(value) and value > 0, field, code)


def _version(value: object, field: str) -> None:
    valid = (
        isinstance(value, str)
        and bool(value.strip())
        and value == value.strip()
        and len(value) <= 200
        and "/" not in value
        and "\\" not in value
        and not any(ord(char) < 32 or ord(char) == 127 for char in value)
    )
    require(valid, field, "declared_version_required")


def _owner(section: str, metadata: dict) -> str:
    return "run" if metadata.get("run", section not in _RESOURCE_SECTIONS) else "resource"


def _partition(data: dict, owner: str) -> None:
    allowed_sections = {
        section
        for section, fields in _SCHEMA.items()
        if any(_owner(section, metadata) == owner for metadata in fields.values())
    }
    for section in data:
        if section not in allowed_sections:
            field = section if section in _SCHEMA else "config"
            code = "misplaced_section" if section in _SCHEMA else "unknown_key"
            raise MendelPrimeError(code, field=field)
        raw = _table(data[section], section, set(_SCHEMA[section]))
        for field in raw:
            if _owner(section, _SCHEMA[section][field]) != owner:
                raise MendelPrimeError("misplaced_key", field=f"{section}.{field}")


def _normalize(name: str, raw: object, run_base: Path, resource_base: Path | None) -> dict:
    schema = _SCHEMA[name]
    values = {field: metadata["default"] for field, metadata in schema.items() if "default" in metadata}
    values.update(_table(raw, name, set(schema)))
    for field, metadata in schema.items():
        qualified = f"{name}.{field}"
        if metadata.get("required"):
            require(field in values, qualified, "required")
        if metadata.get("path") and values.get(field) is not None:
            base = resource_base if resource_base and _owner(name, metadata) == "resource" else run_base
            values[field] = _path(values[field], qualified, base)
        if metadata.get("version") and field in values:
            _version(values[field], qualified)
    if name == "llm":
        values = _validated_settings(values, wire_code="unsupported_protocol", config_order=True)
    elif name == "execution":
        _positive(values["lanes"], "execution.lanes", integer=True)
        _positive(values["max_llm_concurrency"], "execution.max_llm_concurrency", integer=True)
        _positive(values["queue_timeout_seconds"], "execution.queue_timeout_seconds")
        require(type(values["cache_choices"]) is bool, "execution.cache_choices", "boolean_required")
        require(values["device"] in {"cpu", "cuda"}, "execution.device", "cpu_or_cuda_required")
    elif name == "workers":
        for field in schema:
            _positive(values[field], f"workers.{field}", integer=True)
    elif name == "genotype":
        _normalize_genotype(values, resource_base or run_base)
    elif name == "explanation_resources":
        if isinstance(values.get("vep_command"), list) and values["vep_command"]:
            values["vep_command"][0] = str(_path(
                values["vep_command"][0], "explanation_resources.vep_command", resource_base or run_base
            ))
        _normalize_explanation(values)
    elif name == "external":
        for field in ("phen2gene_acknowledged", "amelie_acknowledged", "amelie_verify_tls"):
            require(type(values[field]) is bool, f"external.{field}", "boolean_required")
        require(
            isinstance(values["amelie_pseudonym_secret_env"], str)
            and bool(_ENVIRONMENT_NAME.fullmatch(values["amelie_pseudonym_secret_env"])),
            "external.amelie_pseudonym_secret_env",
            "environment_name_required",
        )
    elif name in {"ranking", "output"}:
        for field in schema:
            require(type(values[field]) is bool, f"{name}.{field}", "boolean_required")
    return values


def _normalize_genotype(values: dict, base: Path) -> None:
    require(
        isinstance(values["assemblies"], list) and bool(values["assemblies"]),
        "genotype.assemblies",
        "nonempty_array_of_tables_required",
    )
    assemblies = []
    for index, item in enumerate(values["assemblies"]):
        prefix = f"genotype.assemblies[{index}]"
        entry = _table(item, prefix, set(_ASSEMBLY_SCHEMA))
        for field, metadata in _ASSEMBLY_SCHEMA.items():
            qualified = f"{prefix}.{field}"
            require(field in entry, qualified, "required")
            if metadata.get("path"):
                entry[field] = _path(entry[field], qualified, base)
            if metadata.get("release"):
                require(
                    isinstance(entry[field], str) and bool(_RELEASE.fullmatch(entry[field])),
                    qualified,
                    "resource_release_required",
                )
        require(
            entry["assembly"] in {"GRCh37", "GRCh38"}, f"{prefix}.assembly", "supported_assembly_required"
        )
        assemblies.append(entry)
    require(
        len({item["assembly"] for item in assemblies}) == len(assemblies),
        "genotype.assemblies",
        "duplicate_assembly",
    )
    values["assemblies"] = tuple(assemblies)
    require(
        isinstance(values["max_heap"], str) and bool(re.fullmatch(r"[1-9]\d*[mMgG]", values["max_heap"])),
        "genotype.max_heap",
        "positive_java_heap_required",
    )
    _positive(values["timeout_seconds"], "genotype.timeout_seconds")
    require(
        values["lirical_transcript_db"] in {"REFSEQ", "REFSEQ_CURATED", "UCSC", "ENSEMBL"},
        "genotype.lirical_transcript_db",
        "unsupported_transcript_database",
    )
    for field in ("exomiser_version", "lirical_version"):
        require(
            isinstance(values[field], str) and bool(values[field].strip()),
            f"genotype.{field}",
            "tool_version_required",
        )


def _normalize_explanation(values: dict) -> None:
    require(values["assembly"] == "GRCh37", "explanation_resources.assembly", "grch37_required")
    require(
        type(values["ranking_vep_enabled"]) is bool,
        "explanation_resources.ranking_vep_enabled",
        "boolean_required",
    )
    if "vep_command" in values:
        command = values["vep_command"]
        require(
            isinstance(command, list)
            and bool(command)
            and all(isinstance(part, str) and bool(part) and "\x00" not in part for part in command),
            "explanation_resources.vep_command",
            "nonempty_argv_required",
        )
        values["vep_command"] = tuple(command)
    _positive(values["vep_timeout_seconds"], "explanation_resources.vep_timeout_seconds")
    _positive(values["pvs1_timeout_seconds"], "explanation_resources.pvs1_timeout_seconds")


def load_config(path: str | Path) -> dict:
    run_path = Path(path).resolve()
    run = _read(run_path, "config")
    resource_name = run.pop("resource_file", None)
    resource_base = None
    if resource_name is None:
        root = _table(run, "config", set(_SCHEMA))
    else:
        resource_path = _path(resource_name, "resource_file", run_path.parent)
        resources = _read(resource_path, "resource_file")
        resource_root = resources.pop("resource_root", None)
        _partition(run, "run")
        _partition(resources, "resource")
        root = dict(run)
        for name, value in resources.items():
            root[name] = dict(value) | root.get(name, {})
        resource_base = (
            _path(resource_root, "resource_root", resource_path.parent)
            if resource_root is not None else resource_path.parent
        )
    require("llm" in root, "llm.base_url", "required")
    result = {}
    for name in _SCHEMA:
        if name in root:
            result[name] = _normalize(name, root[name], run_path.parent, resource_base)
        elif name in _RESOURCE_SECTIONS:
            result[name] = None
        else:
            result[name] = _normalize(name, {}, run_path.parent, resource_base)
    return result


def apply_overrides(config: dict, **values: object) -> dict:
    require(not set(values) - {"lanes", "max_llm_concurrency", "explanation"}, "overrides", "unknown_key")
    updated = copy.deepcopy(config)
    for name, value in values.items():
        if value is None:
            continue
        section = "output" if name == "explanation" else "execution"
        updated[section][name] = value
        updated[section] = _normalize(section, updated[section], Path.cwd(), None)
    return updated
