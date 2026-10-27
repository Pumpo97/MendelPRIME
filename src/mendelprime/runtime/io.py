"""Bounded patient parsing and atomic public package I/O."""

from __future__ import annotations

import ctypes
import json
import logging
import os
import stat
from collections.abc import Mapping
from contextlib import suppress
from pathlib import Path

import yaml
from yaml.events import AliasEvent, CollectionEndEvent, CollectionStartEvent, ScalarEvent
from yaml.nodes import MappingNode, ScalarNode

from mendelprime.domain import (
    MendelPrimeError,
    Patient,
    require,
    validate_gene,
    validate_identifier,
)
from mendelprime.runtime.storage import read_regular_bounded

MAX_INPUT_BYTES = 32 * 1024 * 1024
MAX_RESULT_BYTES = 128 * 1024 * 1024
MAX_NESTING_DEPTH = 64
MAX_CASES = 10_000
MAX_DISCOVERED_FILES = 10_000
MAX_BATCH_BYTES = 256 * 1024 * 1024
MAX_YAML_NODES = 1_000_000

_LOG = logging.getLogger(__name__)
_MISSING = object()
_SCHEMAS = {"mendelprime_case_v1", "mgagent_vcf_offline_case_v1"}
_FIELDS = {
    "patient": {
        "schema_version",
        "case_id",
        "cohort",
        "phenotypes",
        "phenotype_terms",
        "clinical_description",
        "genomic",
        "genomic_input",
    },
    "phenotype": {"hpo_id", "id", "label", "name", "present"},
    "polarity": {"present", "absent"},
    "genomic": {
        "vcf_path",
        "vcf",
        "assembly",
        "family_mode",
        "sample_roles",
        "proband_sample_id",
        "sample_id",
        "pedigree_path",
        "sex",
        "assay_type",
    },
    "roles": {"proband", "father", "mother"},
}
_ASSEMBLIES = {"grch37": "GRCh37", "hg19": "GRCh37", "grch38": "GRCh38", "hg38": "GRCh38"}


def _error(code, field="$", *, source=None, index=None, category="input"):
    return MendelPrimeError(
        code,
        category=category,
        field=field,
        source=source,
        record_index=index,
    )


def _located(error, source, index=None):
    return _error(error.code, error.field, source=source, index=index)


def _mapping(value, field, schema):
    require(isinstance(value, Mapping), field, "mapping_required")
    require(all(isinstance(key, str) for key in value), f"{field}.*", "string_key_required")
    require(all(key in _FIELDS[schema] for key in value), f"{field}.*", "unknown_field")
    return value


def _text(value, field, *, strip=True, empty=False):
    require(value is not _MISSING, field, "required_field")
    require(isinstance(value, str), field, "text_required")
    normalized = value.strip() if strip else value
    require(bool(normalized) or empty, field, "nonempty_text_required")
    return normalized


def _take(data, names, field, *, required=False):
    present = [name for name in names if name in data]
    require(len(present) < 2, field, "alias_conflict")
    if present:
        return data[present[0]]
    require(not required, field, "required_field")
    return _MISSING


def read_patients(path, cohort=None):
    """Read one supported file or a deterministic directory collection."""
    try:
        source = Path(os.path.abspath(os.fspath(path)))
    except TypeError:
        raise _error("path_required", source=Path(".")) from None
    if not source.exists():
        raise _error("not_found", source=source)
    if source.is_dir():
        files = _directory_inputs(source)
        if not files:
            fallback = source / "cohort.jsonl"
            if not fallback.exists() and not fallback.is_symlink():
                raise _error("no_case_inputs", source=source)
            files = (fallback,)
    elif source.is_file():
        files = (source,)
    else:
        raise _error("unsupported_input_path", source=source)
    patients, seen, portable, total = [], set(), set(), 0
    for input_path in files:
        loaded, size = _read_patient_file(input_path, cohort, seen, portable)
        patients.extend(loaded)
        total += size
        if total > MAX_BATCH_BYTES or len(patients) > MAX_CASES:
            code = "batch_too_large" if total > MAX_BATCH_BYTES else "too_many_cases"
            raise _error(code, source=input_path)
    if not patients:
        raise _error("no_cases", source=source)
    return tuple(patients)


def parse_patient(data, base, cohort=None):
    """Normalize one decoded record into a canonical Patient."""
    if not isinstance(base, Path):
        raise _error("base_path_required", source=Path("."))
    source = Path(os.path.abspath(os.fspath(base)))
    try:
        _check_depth(data)
        return _parse_patient(data, source, cohort)
    except MendelPrimeError as error:
        raise _located(error, source) from None


def _read_patient_file(path, cohort, seen, portable):
    raw = _read_bytes(path, MAX_INPUT_BYTES, input_boundary=True)
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        raise _error("invalid_utf8", source=path) from None
    index = None
    try:
        suffix = path.suffix.lower()
        if suffix == ".jsonl":
            records = []
            for line_number, line in enumerate(text.splitlines(), 1):
                if line.strip():
                    index = line_number
                    value = _json(line, "duplicate_key", "invalid_json")
                    require(isinstance(value, Mapping), "$", "case_object_required")
                    records.append((line_number, value))
                    require(len(records) <= MAX_CASES, "$", "too_many_cases")
        else:
            require(suffix in {".json", ".yaml", ".yml"}, "$", "unsupported_format")
            value = _json(text, "duplicate_key", "invalid_json") if suffix == ".json" else _yaml(text)
            require(isinstance(value, (Mapping, list)), "$", "case_object_or_list_required")
            records = [(1, value)] if isinstance(value, Mapping) else list(enumerate(value, 1))
            require(len(records) <= MAX_CASES, "$", "too_many_cases")
        loaded = []
        for index, record in records:
            _check_depth(record)
            require(isinstance(record, Mapping), "$", "case_object_required")
            patient = _parse_patient(record, path.parent, cohort)
            require(patient.key not in seen, "case_id", "duplicate_case_key")
            folded = tuple(part.casefold() for part in patient.key)
            require(folded not in portable, "case_id", "portable_case_key_collision")
            seen.add(patient.key)
            portable.add(folded)
            loaded.append(patient)
        return loaded, len(raw)
    except MendelPrimeError as error:
        raise _located(error, path, index) from None


def _parse_patient(data, base, cohort):
    record = _mapping(data, "$", "patient")
    schema = record.get("schema_version", _MISSING)
    if schema is not _MISSING:
        require(isinstance(schema, str), "schema_version", "text_required")
        require(schema in _SCHEMAS, "schema_version", "unsupported_schema")
    case_id = _text(record.get("case_id", _MISSING), "case_id")
    selected = record.get("cohort", _MISSING)
    selected = None if selected is _MISSING else _text(selected, "cohort")
    if cohort is not None:
        requested = _text(cohort, "cohort")
        require(selected in {None, requested}, "cohort", "cohort_conflict")
        selected = requested
    raw_terms = _take(record, ("phenotypes", "phenotype_terms"), "phenotypes", required=True)
    clinical = record.get("clinical_description", "")
    require(isinstance(clinical, str), "clinical_description", "text_required")
    raw_genomic = _take(record, ("genomic", "genomic_input"), "genomic")
    genomic = None if raw_genomic is _MISSING else _genomic(raw_genomic, base)
    return Patient(
        case_id=case_id,
        cohort=selected or "cases",
        phenotypes=_phenotypes(raw_terms),
        clinical_description=clinical,
        genomic=genomic,
    )


def _phenotypes(value):
    items = []
    if isinstance(value, Mapping):
        groups = _mapping(value, "phenotypes", "polarity")
        for name, present in (("present", True), ("absent", False)):
            group = groups.get(name, [])
            require(isinstance(group, (list, tuple)), f"phenotypes.{name}", "sequence_required")
            items.extend(
                _phenotype(item, f"phenotypes.{name}[{index}]", present) for index, item in enumerate(group)
            )
    else:
        require(isinstance(value, (list, tuple)), "phenotypes", "sequence_or_polarity_mapping_required")
        items.extend(_phenotype(item, f"phenotypes[{index}]", None) for index, item in enumerate(value))
    merged, positions = [], {}
    for item in items:
        position = positions.get(item["hpo_id"])
        if position is None:
            positions[item["hpo_id"]] = len(merged)
            merged.append(item)
            continue
        previous = merged[position]
        require(previous["present"] is item["present"], "phenotypes", "conflicting_hpo_polarity")
        require(
            not previous["label"] or not item["label"] or previous["label"] == item["label"],
            "phenotypes",
            "conflicting_hpo_label",
        )
        if not previous["label"] and item["label"]:
            merged[position] = item
    return tuple(merged)


def _phenotype(value, field, forced):
    if isinstance(value, str):
        raw_id, label, present = value, "", True if forced is None else forced
    else:
        item = _mapping(value, field, "phenotype")
        raw_id = _take(item, ("hpo_id", "id"), f"{field}.hpo_id", required=True)
        raw_label = _take(item, ("label", "name"), f"{field}.label")
        label = "" if raw_label is _MISSING else _text(raw_label, f"{field}.label", empty=True)
        supplied = item.get("present", _MISSING)
        require(supplied is _MISSING or type(supplied) is bool, f"{field}.present", "boolean_required")
        present = (
            True if supplied is _MISSING and forced is None else forced if supplied is _MISSING else supplied
        )
        require(forced is None or present is forced, f"{field}.present", "polarity_conflict")
    require(isinstance(raw_id, str), f"{field}.hpo_id", "text_required")
    return {"hpo_id": raw_id.strip().upper(), "label": label, "present": present}


def _genomic(value, base):
    record = _mapping(value, "genomic", "genomic")
    raw_vcf = _take(record, ("vcf_path", "vcf"), "genomic.vcf_path", required=True)
    assembly = record.get("assembly", _MISSING)
    if assembly is not _MISSING and assembly is not None:
        assembly = _text(assembly, "genomic.assembly").lower()
        assembly = None if assembly == "unknown" else _ASSEMBLIES.get(assembly, _MISSING)
        require(assembly is not _MISSING, "genomic.assembly", "unsupported_assembly")
    else:
        assembly = None
    family = _text(record.get("family_mode", "auto"), "genomic.family_mode").lower()
    require(family in {"auto", "trio", "proband_only"}, "genomic.family_mode", "unsupported_family_mode")
    role_data = record.get("sample_roles", _MISSING)
    proband = _take(record, ("proband_sample_id", "sample_id"), "genomic.proband_sample_id")
    require(role_data is _MISSING or proband is _MISSING, "genomic.sample_roles", "role_declaration_conflict")
    roles = (
        _roles(role_data)
        if role_data is not _MISSING
        else (None if proband is _MISSING or proband is None else _roles({"proband": proband}))
    )
    pedigree = record.get("pedigree_path", _MISSING)
    pedigree = (
        None
        if pedigree is _MISSING or pedigree is None
        else _absolute_path(pedigree, base, "genomic.pedigree_path")
    )
    sex = _text(record.get("sex", "UNKNOWN"), "genomic.sex").upper()
    require(sex in {"MALE", "FEMALE", "UNKNOWN"}, "genomic.sex", "unsupported_value")
    raw_assay = _text(record.get("assay_type", "UNKNOWN"), "genomic.assay_type", strip=False)
    normalized_assay = raw_assay.strip()
    require(bool(normalized_assay), "genomic.assay_type", "nonempty_text_required")
    assay = normalized_assay.upper()
    if assay not in {"WES", "WGS", "PANEL", "UNKNOWN"}:
        assay = raw_assay
    return {
        "vcf_path": _absolute_path(raw_vcf, base, "genomic.vcf_path"),
        "assembly": assembly,
        "family_mode": family,
        "sample_roles": roles,
        "pedigree_path": pedigree,
        "sex": sex,
        "assay_type": assay,
    }


def _roles(value):
    if value is None:
        return None
    record = _mapping(value, "genomic.sample_roles", "roles")
    if not any(raw is not None for raw in record.values()):
        return None
    require(
        record.get("proband") is not None, "genomic.sample_roles.proband", "proband_required_with_parents"
    )
    return {
        name: None if record.get(name) is None else _text(record[name], f"genomic.sample_roles.{name}")
        for name in ("proband", "father", "mother")
    }


def _absolute_path(value, base, field):
    raw = _text(value, field, strip=False)
    require("\x00" not in raw, field, "path_required")
    path = Path(raw)
    return Path(os.path.abspath(os.fspath(path if path.is_absolute() else base / path)))


def _directory_inputs(root):
    found = []

    def unreadable(_):
        raise _error("directory_read_error", source=root)

    try:
        for current, directories, files in os.walk(root, topdown=True, onerror=unreadable, followlinks=False):
            parent = Path(current)
            if "case.yaml" in directories:
                raise _error("input_file_must_be_regular", source=parent / "case.yaml")
            directories[:] = sorted(name for name in directories if not (parent / name).is_symlink())
            if "case.yaml" in files:
                found.append(parent / "case.yaml")
                if len(found) > MAX_DISCOVERED_FILES:
                    raise _error("too_many_input_files", source=root)
    except MendelPrimeError:
        raise
    except OSError:
        raise _error("directory_read_error", source=root) from None
    return tuple(sorted(found, key=lambda item: item.relative_to(root).as_posix()))


def _read_bytes(path, limit, *, input_boundary=False):
    field = "$" if input_boundary else "result"
    category = "input" if input_boundary else "contract"
    try:
        if input_boundary:
            status = os.lstat(path)
            if stat.S_ISLNK(status.st_mode):
                raise _error("input_symlink_forbidden", source=path)
            if not stat.S_ISREG(status.st_mode):
                raise _error("input_file_must_be_regular", source=path)
        try:
            raw, before, after = read_regular_bounded(
                path, limit, follow_symlinks=not input_boundary, check_size=True
            )
        except OverflowError:
            raise _error("file_too_large", field, source=path, category=category) from None
        if not input_boundary and (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise _error("file_changed_during_read", "result", category="contract")
        return raw
    except MendelPrimeError:
        raise
    except OSError:
        code = "io_error" if input_boundary else "unreadable_json_result"
        raise _error(code, field, source=path, category=category) from None


def _json(text, duplicate_code, invalid_code, category="input"):
    def unique(pairs):
        value = {}
        for key, item in pairs:
            require(key not in value, "$", duplicate_code)
            value[key] = item
        return value

    try:
        return json.loads(text, object_pairs_hook=unique)
    except MendelPrimeError as error:
        if category == "input":
            raise
        raise _error(error.code, error.field, category=category) from None
    except (ValueError, OverflowError, RecursionError):
        raise _error(invalid_code, category=category) from None


def _yaml(text):
    try:
        depth = nodes = 0
        for event in yaml.parse(text, Loader=yaml.SafeLoader):
            require(not isinstance(event, AliasEvent), "$", "yaml_alias_forbidden")
            if isinstance(event, CollectionStartEvent):
                depth += 1
                nodes += 1
            elif isinstance(event, CollectionEndEvent):
                depth -= 1
            elif isinstance(event, ScalarEvent):
                nodes += 1
            require(depth <= MAX_NESTING_DEPTH, "$", "nesting_too_deep")
            require(nodes <= MAX_YAML_NODES, "$", "too_many_yaml_nodes")
        root = yaml.compose(text, Loader=yaml.SafeLoader)
        stack = [] if root is None else [root]
        while stack:
            node = stack.pop()
            if isinstance(node, MappingNode):
                keys = []
                for key, item in node.value:
                    signature = (key.tag, key.value) if isinstance(key, ScalarNode) else None
                    require(signature is None or signature not in keys, "$", "duplicate_key")
                    keys.append(signature)
                    stack.extend((key, item))
            elif hasattr(node, "value") and isinstance(node.value, list):
                stack.extend(node.value)
        return yaml.safe_load(text)
    except MendelPrimeError:
        raise
    except (yaml.YAMLError, ValueError, OverflowError, RecursionError):
        raise _error("invalid_yaml") from None


def _check_depth(value):
    stack = [(value, 1, frozenset())]
    while stack:
        current, depth, ancestors = stack.pop()
        require(depth <= MAX_NESTING_DEPTH, "$", "nesting_too_deep")
        if not isinstance(current, (Mapping, list, tuple)):
            continue
        identity = id(current)
        require(identity not in ancestors, "$", "cyclic_structure")
        children = current.values() if isinstance(current, Mapping) else current
        stack.extend((child, depth + 1, ancestors | {identity}) for child in children)


def read_public(path, cohort=None):
    """Shape-check a public V4 result without hydrating scientific records."""
    from .results import validate_public_shape

    source = Path(path)
    if source.is_dir():
        source /= "result.json"
    raw = _read_bytes(source, MAX_RESULT_BYTES)
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        raise _error("unreadable_json_result", "result", category="contract") from None
    value = _json(text, "duplicate_json_key", "unreadable_json_result", "contract")
    genes = validate_public_shape(value, decoded_json=True)
    provenance = value.get("provenance")
    recorded = provenance.get("case_key") if isinstance(provenance, dict) else None
    inferred = recorded.get("cohort") if isinstance(recorded, dict) else source.parent.parent.name
    key = (cohort if cohort is not None else inferred, value.get("case_id"))
    validate_identifier(key[0], "cohort")
    if provenance is None and isinstance(key[1], str) and key[1].count(":") == 1:
        # Legacy results used cohort:case identifiers. Reading is not permission
        # to use that identifier as a new portable output directory.
        for part in key[1].split(":"):
            validate_identifier(part, "case_id")
    else:
        validate_identifier(key[1], "case_id")
    if provenance is not None:
        require(len(set(genes)) == len(genes), "result.final.top100", "duplicate_gene")
        for gene in genes:
            validate_gene(gene, "result.final.top100")
        require(recorded == {"cohort": key[0], "case_id": key[1]}, "result.provenance", "case_mismatch")
        try:
            json.dumps(value, allow_nan=False)
        except (TypeError, ValueError, RecursionError):
            raise _error("strict_json_required", "result", category="contract") from None
    return key, value


def _rename_new_directory(source, target):
    if os.name == "nt":
        source.rename(target)
        return
    rename = getattr(ctypes.CDLL(None, use_errno=True), "renameat2", None)
    require(rename is not None, "output", "atomic_no_replace_unavailable")
    rename.argtypes = (ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint)
    rename.restype = ctypes.c_int
    if rename(-100, os.fsencode(source), -100, os.fsencode(target), 1):
        raise OSError(ctypes.get_errno(), "no-replace publication failed")


def _sync_directory(path):
    if os.name == "posix":
        descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)


def prepare_output_directory(path):
    """Create ordinary output ancestors without following symlinks or changing modes."""
    for ancestor in reversed((path, *path.parents)):
        try:
            ancestor.mkdir(mode=0o700)
        except FileExistsError:
            pass
        require(stat.S_ISDIR(ancestor.lstat().st_mode), "output_root", "ordinary_directory_required")
    return path


def publish_case(output_root, key, result_dict, markdown):
    """Atomically publish exactly result.json and traceable_explanation.md."""
    require(
        isinstance(output_root, Path) and output_root.is_absolute(), "output_root", "absolute_path_required"
    )
    require(type(key) is tuple and len(key) == 2, "case_key", "case_key_required")
    validate_identifier(key[0], "cohort")
    validate_identifier(key[1], "case_id")
    require(
        type(result_dict) is dict and result_dict.get("case_id") == key[1], "result.case_id", "case_mismatch"
    )
    require(isinstance(markdown, str), "explanation", "text_required")
    try:
        encoded = json.dumps(result_dict, ensure_ascii=False, allow_nan=False, indent=2) + "\n"
        contents = (("result.json", encoded.encode()), ("traceable_explanation.md", markdown.encode()))
    except (TypeError, ValueError, UnicodeError):
        raise _error("strict_json_and_utf8_required", "result", category="contract") from None
    try:
        parent = prepare_output_directory(output_root / key[0])
    except OSError:
        raise _error("ordinary_directory_required", "output_root", category="contract") from None
    target = parent / key[1]
    reservation = parent / f".{key[1]}.publishing"
    try:
        reservation.mkdir(mode=0o700)
    except FileExistsError:
        raise _error("case_publication_in_progress", "output", category="contract") from None
    try:
        require(not os.path.lexists(target), "output", "case_already_exists")
        for name, content in contents:
            with (reservation / name).open("xb") as stream:
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
        _sync_directory(reservation)
        try:
            _rename_new_directory(reservation, target)
        except FileExistsError:
            raise _error("case_already_exists", "output", category="contract") from None
    except BaseException:
        for name, _ in contents:
            with suppress(OSError):
                (reservation / name).unlink(missing_ok=True)
        with suppress(OSError):
            reservation.rmdir()
        raise
    try:
        _sync_directory(parent)
    except OSError:
        _LOG.warning("Case files were published; directory durability could not be confirmed")
    return target


__all__ = ["parse_patient", "publish_case", "read_patients", "read_public"]
