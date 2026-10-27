"""HPO records and cached graph queries without a parallel object model."""

import hashlib
import json
import math
import re
import threading
from collections import deque
from pathlib import Path, PurePosixPath, PureWindowsPath

from mendelprime.domain import MendelPrimeError, freeze_record, require
from mendelprime.runtime.storage import open_regular

_ID = re.compile(r"HP:\d{7}", re.ASCII)
_SECTION = re.compile(r"(?m)^\[([^\]\r\n]+)\][ \t]*(?:\r?\n|$)")
_QUOTED = re.compile(r'^"((?:\\.|[^"\\])*)"')


def _identifier(value, field="hpo_id"):
    require(isinstance(value, str) and _ID.fullmatch(value), field, "hpo_identifier_required")


def _read(path, field):
    require(isinstance(path, Path), field, "path_required")
    try:
        with open_regular(path) as (stream, _):
            return stream.read()
    except OSError:
        raise MendelPrimeError("unable_to_read_regular_file", field=field) from None


def _manifest(locator, version, content):
    require(
        isinstance(locator, str)
        and re.fullmatch(r"resource:[A-Za-z0-9][A-Za-z0-9._~:/-]*", locator)
        and ".." not in PurePosixPath(locator.removeprefix("resource:")).parts
        and not PureWindowsPath(locator.removeprefix("resource:")).drive,
        "manifest.locator",
        "logical_resource_locator_required",
    )
    require(
        isinstance(version, str)
        and version.strip() == version
        and version
        and not any(ord(c) < 32 or ord(c) == 127 for c in version),
        "manifest.version",
        "canonical_version_required",
    )
    return {"locator": locator, "version": version, "sha256": hashlib.sha256(content).hexdigest()}


def load_hpo(
    obo_path,
    statistics_path=None,
    *,
    ontology_locator="resource:hpo/hp.obo",
    statistics_locator="resource:hpo/statistics.json",
):
    raw = _read(obo_path, "hpo_obo")
    try:
        text = raw.decode("utf-8")
    except UnicodeError:
        raise MendelPrimeError("utf8_required", field="hpo_obo") from None
    sections = tuple(_SECTION.finditer(text))
    header = text[: sections[0].start()] if sections else text
    terms = {}
    for index, section in enumerate(sections):
        if section.group(1) != "Term":
            continue
        end = sections[index + 1].start() if index + 1 < len(sections) else len(text)
        term = _parse_term(text[section.start() : end])
        if term is not None:
            require(term["hpo_id"] not in terms, "hpo_ontology", "duplicate_term")
            terms[term["hpo_id"]] = term
    require(terms, "hpo_obo", "nonempty_term_records_required")
    _reject_cycles(terms)
    version = next(
        (
            line.split(":", 1)[1].strip()
            for line in header.splitlines()
            if line.startswith("data-version:") and line.split(":", 1)[1].strip()
        ),
        "sha256:" + hashlib.sha256(raw).hexdigest(),
    )
    manifest = _manifest(ontology_locator, version, raw)
    statistics, gene_ic, disease_ic = None, {}, {}
    if statistics_path is not None:
        data = _read(statistics_path, "hpo_statistics")
        try:
            payload = json.loads(data.decode("utf-8"))
        except (UnicodeError, json.JSONDecodeError):
            raise MendelPrimeError("valid_utf8_json_required", field="hpo_statistics") from None
        require(isinstance(payload, dict), "hpo_statistics", "object_required")
        gene_ic = _ic(payload.get("ic_gene", {}), "gene_information_content")
        disease_ic = _ic(payload.get("ic_disease", {}), "disease_information_content")
        version = next(
            (
                payload[key].strip()
                for key in ("data_version", "version")
                if isinstance(payload.get(key), str) and payload[key].strip()
            ),
            "sha256:" + hashlib.sha256(data).hexdigest(),
        )
        statistics = _manifest(statistics_locator, version, data)
    return {
        "terms": freeze_record(terms),
        "gene_ic": freeze_record(gene_ic),
        "disease_ic": freeze_record(disease_ic),
        "ontology_manifest": freeze_record(manifest),
        "statistics_manifest": freeze_record(statistics),
        "depths": {},
        "ancestors": {},
        "lock": threading.RLock(),
    }


def _parse_term(original):
    row = {
        "hpo_id": "",
        "name": "",
        "definition": "",
        "synonyms": [],
        "parent_ids": [],
        "original_record": original,
    }
    obsolete = False
    for raw in original.splitlines()[1:]:
        line = raw.strip()
        if line.startswith("id: "):
            row["hpo_id"] = line[4:].strip()
        elif line.startswith("name: "):
            row["name"] = line[6:]
        elif line.startswith(("def: ", "synonym: ")):
            name, value = line.split(": ", 1)
            match = _QUOTED.match(value)
            value = match.group(1).replace(r"\"", '"') if match else ""
            if name == "def":
                row["definition"] = value
            elif value and value not in row["synonyms"]:
                row["synonyms"].append(value)
        elif line.startswith("is_a: "):
            parent = line[6:].split("!", 1)[0].strip()
            if parent and parent not in row["parent_ids"]:
                row["parent_ids"].append(parent)
        elif line == "is_obsolete: true":
            obsolete = True
    if obsolete:
        return None
    require(row["hpo_id"], "hpo_obo", "term_id_required")
    _identifier(row["hpo_id"], "hpo_term.hpo_id")
    require(
        re.findall(r"(?m)^id:[ \t]*(HP:\d{7})[ \t]*\r?$", original, re.ASCII) == [row["hpo_id"]],
        "hpo_term.original_record",
        "term_identity_mismatch",
    )
    for parent in row["parent_ids"]:
        _identifier(parent, "hpo_term.parent_ids")
    return freeze_record(row)


def _ic(values, field):
    require(isinstance(values, dict), field, "mapping_required")
    result = {}
    for key, value in values.items():
        _identifier(key, field)
        require(
            type(value) in (int, float) and math.isfinite(value) and value >= 0,
            field,
            "nonnegative_finite_number_required",
        )
        result[key] = float(value)
    return result


def _reject_cycles(terms):
    indegree = dict.fromkeys(terms, 0)
    children = {key: [] for key in terms}
    for child, term in terms.items():
        for parent in term["parent_ids"]:
            if parent in terms:
                children[parent].append(child)
                indegree[child] += 1
    queue = deque(key for key, value in indegree.items() if value == 0)
    visited = 0
    while queue:
        visited += 1
        for child in children[queue.popleft()]:
            indegree[child] -= 1
            if indegree[child] == 0:
                queue.append(child)
    require(visited == len(terms), "hpo_ontology", "cycle_detected")


def ancestors(hpo, identifier):
    _identifier(identifier)
    if identifier not in hpo["terms"]:
        return frozenset()
    with hpo["lock"]:
        if identifier not in hpo["ancestors"]:
            found = set()
            queue = deque(hpo["terms"][identifier]["parent_ids"])
            while queue:
                current = queue.popleft()
                if current not in found:
                    found.add(current)
                    queue.extend(hpo["terms"].get(current, {}).get("parent_ids", ()))
            hpo["ancestors"][identifier] = frozenset(found)
        return hpo["ancestors"][identifier]


def depth(hpo, identifier):
    """Compute depth using breadth-first traversal with a shared visited set."""
    _identifier(identifier)
    if identifier not in hpo["terms"]:
        return 0
    with hpo["lock"]:
        if identifier not in hpo["depths"]:
            queue, seen, maximum = deque([(identifier, 0)]), {identifier}, 0
            while queue:
                current, level = queue.popleft()
                maximum = max(maximum, level)
                for parent in hpo["terms"].get(current, {}).get("parent_ids", ()):
                    if parent not in seen:
                        seen.add(parent)
                        queue.append((parent, level + 1))
            hpo["depths"][identifier] = maximum
        return hpo["depths"][identifier]


def consulted_terms(hpo, identifiers):
    result = set()
    for identifier in identifiers:
        if identifier in hpo["terms"]:
            result.add(identifier)
        result.update(parent for parent in ancestors(hpo, identifier) if parent in hpo["terms"])
    return tuple(sorted(result))
