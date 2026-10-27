"""Read-only VCF/PED inspection and exact family observations."""

from __future__ import annotations

import gzip
import math
import os
import re
import shutil
import tempfile
from collections.abc import Mapping
from contextlib import contextmanager, suppress
from pathlib import Path

import pysam

from mendelprime.domain import (
    MendelPrimeError,
    Patient,
    family_origin,
    is_sample_id,
    validate_case_key,
    variant_key,
)
from mendelprime.runtime.storage import open_regular, stat_fingerprint

_ALLELE = re.compile(r"[ACGTN]+", re.ASCII)
_GT = re.compile(r"(?:\d+|\.)(?:[/|](?:\d+|\.))*", re.ASCII)
_COMPRESSED = {".gz", ".bgz", ".bgzf"}
_BUILD_ALIASES = {
    "GRCh37": r"GRCh37(?:\.p\d+)?|hg19(?:un38)?|b37|hs37d5",
    "GRCh38": r"GRCh38(?:\.p\d+)?|hg38|b38",
}
_CHR1_LENGTHS = {249250621: "GRCh37", 248956422: "GRCh38"}


def _error(field, code):
    return MendelPrimeError(code, category="input", field=field)


def _check(condition, field, code):
    if not condition:
        raise _error(field, code)


def _path(value, field):
    try:
        return Path(value)
    except TypeError:
        raise _error(field, "path_required") from None


def _version(path):
    with open_regular(path) as (_, status):
        return stat_fingerprint(status)


def _unchanged(path, version):
    try:
        current = _version(path)
    except OSError:
        raise _error("vcf", "source_changed_during_read") from None
    if current != version:
        raise _error("vcf", "source_changed_during_read")


def _text(raw, field):
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        raise _error(field, "utf8_text_required") from None


@contextmanager
def _decompressed(path):
    with open_regular(path) as (raw, _):
        stream = raw
        try:
            magic = raw.read(2)
            raw.seek(0)
            if magic == b"\x1f\x8b":
                stream = gzip.GzipFile(fileobj=raw, mode="rb")
            _check(stream.read(3) != b"BCF", "vcf", "bcf_unsupported")
            stream.seek(0)
            yield stream
        finally:
            if stream is not raw:
                stream.close()


def _ordinary_gzip(path):
    with open_regular(path) as (source, _):
        fixed = source.read(12)
        if fixed[:2] != b"\x1f\x8b":
            return False
        if len(fixed) < 12 or not fixed[3] & 0x04:
            return True
        extra = source.read(int.from_bytes(fixed[10:12], "little"))
    offset = 0
    while offset + 4 <= len(extra):
        identity = extra[offset : offset + 2]
        length = int.from_bytes(extra[offset + 2 : offset + 4], "little")
        offset += 4
        end = offset + length
        if end > len(extra):
            return True
        if identity == b"BC" and length == 2:
            return False
        offset = end
    return True


@contextmanager
def _parser_source(path):
    version = _version(path)
    try:
        if not _ordinary_gzip(path):
            yield path
        else:
            with tempfile.TemporaryDirectory(prefix="mendelprime-vcf-") as directory:
                uncompressed = Path(directory) / "source.vcf"
                with _decompressed(path) as source, uncompressed.open("xb") as output:
                    shutil.copyfileobj(source, output, length=1024 * 1024)
                yield uncompressed
    finally:
        _unchanged(path, version)


@contextmanager
def _variant_file(path):
    with open_regular(path) as (source, _):
        with pysam.VariantFile(source, mode="r") as variant_file:
            yield variant_file


def _parent(value):
    return None if value in {None, "", ".", "0"} else value


def _validate_roles(roles):
    supplied = []
    for role in ("proband", "father", "mother"):
        sample = roles.get(role)
        if sample is None and role != "proband":
            continue
        _check(is_sample_id(sample), f"sample_roles.{role}", "sample_identifier_required")
        supplied.append(sample)
    _check(len(supplied) == len(set(supplied)), "sample_roles", "roles_must_be_distinct")


def _relationships(links, field):
    by_child = {}
    for link in links:
        _validate_roles(
            {"proband": link["child"], "father": link.get("father"), "mother": link.get("mother")}
        )
        child, parents = link["child"], (link.get("father"), link.get("mother"))
        _check(child not in by_child, field, "duplicate_child")
        _check(child not in parents, field, "self_parent")
        by_child[child] = parents
    while by_child:
        roots = {
            child for child, parents in by_child.items() if not any(parent in by_child for parent in parents)
        }
        _check(bool(roots), field, "pedigree_cycle")
        for child in roots:
            del by_child[child]


def _raw_header(path):
    lines, column_header = [], False
    with _decompressed(path) as stream:
        for raw_line in stream:
            line = _text(raw_line, "vcf")
            if line.startswith("#"):
                lines.append(line)
                column_header |= line.startswith("#CHROM\t")
            elif line.strip():
                break
    _check(column_header, "vcf.header", "column_header_required")
    return "".join(lines)


def _attributes(record):
    try:
        return {str(name): str(value) for name, value in record.items()}
    except (AttributeError, TypeError, ValueError):
        raise _error("vcf.header", "invalid_structured_record") from None


def _pedigree_link(attributes):
    values = {name.lower(): value for name, value in attributes.items()}
    child = values.get("child", values.get("id"))
    _check(bool(child), "vcf.pedigree", "child_required")
    return {
        "child": child,
        "father": _parent(values.get("father")),
        "mother": _parent(values.get("mother")),
    }


def _inspect(path):
    raw_text = _raw_header(path)
    with _variant_file(path) as variant_file:
        samples = tuple(variant_file.header.samples)
        references, contigs, links = [], [], []
        for record in variant_file.header.records:
            key = str(record.key)
            if key.lower() == "reference":
                value = getattr(record, "value", None)
                references.append("" if value is None else str(value))
            elif key.lower() == "contig":
                values = _attributes(record)
                name = values.get("ID")
                _check(bool(name), "vcf.contigs", "contig_id_required")
                try:
                    length = None if values.get("length") in {None, ""} else int(values["length"])
                except ValueError:
                    raise _error("vcf.contigs", "invalid_contig_length") from None
                _check(is_sample_id(name), "contig.name", "sample_identifier_required")
                _check(length is None or length >= 0, "contig.length", "nonnegative_integer_or_none_required")
                contigs.append({"name": name, "length": length, "assembly": values.get("assembly", "")})
            elif key.upper() == "PEDIGREE":
                links.append(_pedigree_link(_attributes(record)))
    valid_samples = (
        samples and len(samples) == len(set(samples)) and all(is_sample_id(sample) for sample in samples)
    )
    _check(valid_samples, "vcf.samples", "unique_nonempty_sample_tuple_required")
    _check(len({item["name"] for item in contigs}) == len(contigs), "vcf.contigs", "duplicate_contig")
    _relationships(links, "vcf.pedigree")
    return {
        "samples": samples,
        "references": tuple(references),
        "contigs": tuple(contigs),
        "pedigree_links": tuple(links),
        "raw_text": raw_text,
    }


def inspect_vcf(path):
    """Validate a text VCF and return a plain header record with exact text."""
    source = _path(path, "vcf.path")
    try:
        with _parser_source(source) as parser_source:
            return _inspect(parser_source)
    except MendelPrimeError:
        raise
    except (OSError, EOFError, NotImplementedError, ValueError):
        raise _error("vcf", "invalid_or_unreadable_vcf") from None


def read_pedigree(path):
    """Read the six leading PLINK PED fields into an exact plain record."""
    source = _path(path, "pedigree.path")
    try:
        with open_regular(source) as (stream, _):
            raw_text = _text(stream.read(), "pedigree")
    except MendelPrimeError:
        raise
    except OSError:
        raise _error("pedigree", "unreadable_pedigree") from None
    sexes = {"0": "UNKNOWN", "1": "MALE", "2": "FEMALE"}
    statuses = {"-9": 0, "0": 0, "1": 1, "2": 2}
    members = []
    for number, raw_line in enumerate(raw_text.splitlines(keepends=True), 1):
        if not raw_line.strip():
            continue
        values, field = raw_line.split(), f"pedigree.line[{number}]"
        _check(len(values) >= 6, field, "six_fields_required")
        family_id, sample_id, father, mother, sex, status = values[:6]
        _check(
            father not in {".", "-9"} and mother not in {".", "-9"},
            f"{field}.parents",
            "zero_unknown_parent_required",
        )
        _check(sample_id not in {father, mother}, field, "self_parent")
        if sex not in sexes or status not in statuses:
            suffix = "sex" if sex not in sexes else "affected"
            raise _error(f"{field}.{suffix}", "unsupported_code")
        _check(is_sample_id(family_id), "pedigree.family_id", "sample_identifier_required")
        members.append(
            {
                "family_id": family_id,
                "sample_id": sample_id,
                "father": _parent(father),
                "mother": _parent(mother),
                "sex": sexes[sex],
                "affected": statuses[status],
                "raw_line": raw_line,
            }
        )
    _check(bool(members), "pedigree.members", "nonempty_member_tuple_required")
    _check(
        len({item["sample_id"] for item in members}) == len(members),
        "pedigree.members",
        "ambiguous_sample_identity",
    )
    _relationships(
        (
            {"child": item["sample_id"], "father": item["father"], "mother": item["mother"]}
            for item in members
        ),
        "pedigree.members",
    )
    return {"members": tuple(members), "raw_text": raw_text}


def _assembly(requested, header):
    evidence = []
    if requested is not None:
        _check(requested in _BUILD_ALIASES, "assembly", "assembly_enum_required")
        evidence.append((requested, "explicit"))
    declarations = [("reference", text) for text in header["references"]]
    declarations.extend(("contig_assembly", item["assembly"]) for item in header["contigs"])
    for source, text in declarations:
        for assembly, aliases in _BUILD_ALIASES.items():
            if re.search(rf"(?<![a-z0-9])(?:{aliases})(?![a-z0-9])", text, re.IGNORECASE):
                evidence.append((assembly, source))
    for contig in header["contigs"]:
        if contig["name"] in {"1", "chr1"} and contig["length"] in _CHR1_LENGTHS:
            evidence.append((_CHR1_LENGTHS[contig["length"]], "chromosome_1_length"))
    builds = {assembly for assembly, _ in evidence}
    _check(bool(builds), "assembly", "explicit_assembly_required")
    _check(len(builds) == 1, "assembly", "conflicting_assembly_evidence")
    return evidence[0][0], tuple(dict.fromkeys(source for _, source in evidence))


def _proband(genomic, header, pedigree):
    if genomic.get("sample_roles") is not None:
        return genomic["sample_roles"]["proband"], ["explicit_roles"]
    samples = set(header["samples"])
    links, members = header["pedigree_links"], pedigree["members"] if pedigree else ()
    children = {item["child"] for item in links if item["child"] in samples}
    ped_children = {
        item["sample_id"]
        for item in members
        if (item["father"] or item["mother"]) and item["sample_id"] in samples
    }
    affected = {
        item["sample_id"] for item in members if item["affected"] == 2 and item["sample_id"] in samples
    }
    active = [
        (group, source)
        for group, source in ((children, "vcf_pedigree"), (ped_children, "pedigree_child"))
        if group
    ]
    if active:
        candidates = set.intersection(*(group for group, _ in active))
        _check(bool(candidates), "sample_roles.proband", "conflicting_proband_evidence")
        if len(candidates) > 1 and affected:
            candidates &= affected
        basis = [source for _, source in active]
    elif affected:
        candidates, basis = affected, ["pedigree_affected_subject"]
    elif len(samples) == 1:
        parents = {
            parent for item in (*links, *members) for parent in (item["father"], item["mother"]) if parent
        }
        candidates, basis = (set(), []) if samples & parents else (samples, ["single_vcf_sample"])
    else:
        candidates, basis = set(), []
    _check(len(candidates) == 1, "sample_roles.proband", "explicit_proband_required")
    return next(iter(candidates)), basis


def resolve_family(patient, header, pedigree=None):
    """Resolve assembly, proband, parents, mode, and sex without sample-order inference."""
    _check(
        isinstance(patient, Patient) and patient.genomic is not None, "patient", "genomic_patient_required"
    )
    _check(isinstance(header, Mapping), "header", "vcf_header_record_required")
    genomic = patient.genomic
    assembly, assembly_basis = _assembly(genomic.get("assembly"), header)
    _check(
        genomic.get("pedigree_path") is None or pedigree is not None,
        "pedigree",
        "declared_pedigree_not_loaded",
    )
    proband, basis = _proband(genomic, header, pedigree)
    _check(proband in header["samples"], "sample_roles.proband", "proband_missing_from_vcf")
    assertions = []
    if genomic.get("sample_roles") is not None:
        assertions.append(("explicit_roles", genomic["sample_roles"]))
    assertions.extend(("vcf_pedigree", item) for item in header["pedigree_links"] if item["child"] == proband)
    member = (
        next((item for item in pedigree["members"] if item["sample_id"] == proband), None)
        if pedigree
        else None
    )
    if member is not None:
        assertions.append(("pedigree", member))
    parents = {}
    for role in ("father", "mother"):
        values = {item.get(role) for _, item in assertions if item.get(role) is not None}
        _check(len(values) <= 1, f"sample_roles.{role}", "conflicting_parent_evidence")
        parents[role] = next(iter(values), None)
    roles = {"proband": proband, **parents}
    _validate_roles(roles)
    basis.extend(source for source, _ in assertions)
    sexes = {
        value for value in (genomic.get("sex"), member["sex"] if member else "UNKNOWN") if value != "UNKNOWN"
    }
    _check(len(sexes) <= 1, "sex", "conflicting_proband_sex")
    samples = set(header["samples"])
    complete = bool(roles["father"] and roles["mother"] and {roles["father"], roles["mother"]} <= samples)
    requested = genomic["family_mode"]
    mode = ("trio" if complete else "proband_only") if requested == "auto" else requested
    _check(mode != "trio" or complete, "sample_roles", "complete_sequenced_trio_required")
    issues = ["pedigree_proband_marked_unaffected"] if member is not None and member["affected"] == 1 else []
    if mode == "proband_only":
        if requested == "auto" and (roles["father"] or roles["mother"]):
            issues.append("incomplete_sequenced_trio")
        roles = {"proband": proband, "father": None, "mother": None}
    selected = {sample for sample in roles.values() if sample}
    families = (
        {item["family_id"] for item in pedigree["members"] if item["sample_id"] in selected}
        if pedigree
        else set()
    )
    _check(len(families) <= 1, "pedigree", "conflicting_selected_family_ids")
    if samples - selected:
        issues.append("unselected_vcf_samples")
    return {
        "assembly": assembly,
        "assembly_basis": assembly_basis,
        "mode": mode,
        "roles": roles,
        "sex": next(iter(sexes), "UNKNOWN"),
        "role_basis": tuple(dict.fromkeys(basis)),
        "pedigree": pedigree,
        "issues": tuple(issues),
    }


def _trim(position, reference, alternate):
    while len(reference) > 1 and len(alternate) > 1 and reference[-1] == alternate[-1]:
        reference, alternate = reference[:-1], alternate[:-1]
    while len(reference) > 1 and len(alternate) > 1 and reference[0] == alternate[0]:
        position, reference, alternate = position + 1, reference[1:], alternate[1:]
    return position, reference, alternate


def _vector(value, field):
    if value is None:
        return None
    values = value if isinstance(value, tuple) else (value,)
    values = tuple(None if item is None else int(item) for item in values)
    _check(
        not any(item is not None and item < 0 for item in values),
        field,
        "nonnegative_integer_or_none_required",
    )
    return values


def _call(role, sample_id, parsed, names, text):
    values = [] if not names and text == "." else text.split(":")
    _check(len(values) <= len(names), "vcf.record", "sample_has_extra_format_values")
    raw_fields = tuple(zip(names[: len(values)], values, strict=True))
    raw, gt = dict(raw_fields), dict(raw_fields).get("GT")
    _check(gt is None or _GT.fullmatch(gt), "call.gt", "vcf_genotype_required")
    gt_value = parsed.get("GT") if "GT" in raw else None
    alleles = () if gt_value is None else tuple(None if item is None else int(item) for item in gt_value)
    if gt is not None:
        decoded = tuple(None if item == "." else int(item) for item in re.split(r"[/|]", gt))
        _check(decoded == alleles, "call.gt", "gt_allele_mismatch")
    ad = None if raw.get("AD") in {None, "."} else _vector(parsed.get("AD"), "call.ad")
    pl = None if raw.get("PL") in {None, "."} else _vector(parsed.get("PL"), "call.pl")
    dp_value = None if raw.get("DP") in {None, "."} else parsed.get("DP")
    dp = None if dp_value is None else int(dp_value)
    gq = None if raw.get("GQ") in {None, "."} else float(raw["GQ"])
    _check(dp is None or dp >= 0, "call.dp", "nonnegative_integer_or_none_required")
    _check(gq is None or math.isfinite(gq) and gq >= 0, "call.gq", "nonnegative_finite_number_required")
    return {
        "role": role,
        "sample_id": sample_id,
        "gt": gt,
        "alleles": alleles,
        "ad": ad,
        "dp": dp,
        "gq": gq,
        "pl": pl,
        "phase_set": None if raw.get("PS") in {None, "."} else raw["PS"],
        "filter_text": raw.get("FT"),
        "raw_fields": raw_fields,
    }


def _observation(target, observed, alternate, calls, columns, raw_record):
    try:
        chromosome = columns[0].removeprefix("chr")
        identity = (
            "MT" if chromosome == "M" else chromosome,
            int(columns[1]),
            columns[3].upper(),
            columns[4].split(",")[alternate - 1].upper(),
        )
    except (IndexError, ValueError):
        raise _error("observation.raw_record", "raw_allele_identity_mismatch") from None
    _check(identity == observed[1:], "observation.raw_record", "raw_allele_identity_mismatch")
    return {
        "variant_key": target,
        "alternate_index": alternate,
        "samples": calls,
        "quality": None if columns[5] == "." else float(columns[5]),
        "filter_text": columns[6],
        "raw_record": raw_record,
        "vcf_variant_key": observed,
        "match_method": "exact" if target == observed else "trimmed_padding",
        "origin": family_origin(calls, alternate),
        "parentage_status": (
            "unconfirmed" if {"father", "mother"} <= {call["role"] for call in calls} else "unavailable"
        ),
    }


def _targets(values, family):
    try:
        targets = tuple(values)
    except TypeError:
        raise _error("targets", "variant_key_iterable_required") from None
    normalized = []
    for target in targets:
        if type(target) is not tuple or len(target) != 5:
            if not _ALLELE.fullmatch(str(getattr(target, "alternate", ""))):
                raise _error("targets", "unsupported_symbolic_allele")
            raise _error("targets", "variant_key_required")
        if not isinstance(target[4], str) or not _ALLELE.fullmatch(target[4]):
            raise _error("targets", "unsupported_symbolic_allele")
        key = variant_key(*target)
        if key[0] != family["assembly"]:
            raise _error("targets", "assembly_mismatch")
        normalized.append(key)
    _check(len(normalized) == len(set(normalized)), "targets", "duplicate_target_allele")
    return tuple(normalized)


def _raw_records(stream):
    for raw_line in stream:
        if not raw_line.startswith(b"#") and raw_line.strip():
            yield _text(raw_line, "vcf.record")


def read_observations(path, family, targets, *, case_key):
    """Read exact sample calls using exact or unambiguous padding-only identity."""
    _check(
        type(family) is dict and family.get("mode") in {"trio", "proband_only"},
        "family",
        "family_resolution_record_required",
    )
    validate_case_key(case_key)
    requested = _targets(targets, family)
    if not requested:
        return {"case_key": case_key, "observations": (), "missing_targets": ()}
    source = _path(path, "vcf.path")
    try:
        with _parser_source(source) as parser_source:
            header = _inspect(parser_source)
            _assembly(family["assembly"], header)
            roles = [("proband", family["roles"]["proband"])]
            roles += [
                (role, family["roles"][role]) for role in ("mother", "father") if family["roles"].get(role)
            ]
            _check(
                not any(sample not in header["samples"] for _, sample in roles),
                "family.roles",
                "selected_sample_absent_from_vcf",
            )
            identities = {}
            for target in requested:
                identity = (target[1], *_trim(target[2], target[3], target[4]))
                identities.setdefault(identity, []).append(target)
            found = {}
            sample_columns = {sample: 9 + header["samples"].index(sample) for _, sample in roles}
            with _decompressed(parser_source) as raw_stream, _variant_file(parser_source) as variant_file:
                raw_iterator = _raw_records(raw_stream)
                for record in variant_file:
                    try:
                        raw_record = next(raw_iterator)
                    except StopIteration:
                        raise _error("vcf.records", "raw_and_parsed_record_mismatch") from None
                    columns = raw_record.rstrip("\r\n").split("\t")
                    _check(
                        len(columns) == 9 + len(header["samples"]),
                        "vcf.record",
                        "sample_column_count_mismatch",
                    )
                    names = () if columns[8] in {"", "."} else tuple(columns[8].split(":"))
                    _check(len(names) == len(set(names)), "vcf.record", "duplicate_format_field")
                    _check("GT" not in names or names[0] == "GT", "vcf.record", "gt_must_be_first_format")
                    contig = str(record.contig).removeprefix("chr")
                    contig = "MT" if contig == "M" else contig
                    for alternate_index, alternate in enumerate(record.alts or (), 1):
                        reference, alternate = record.ref.upper(), alternate.upper()
                        if (
                            record.pos < 1
                            or not _ALLELE.fullmatch(reference)
                            or not _ALLELE.fullmatch(alternate)
                        ):
                            continue
                        matches = identities.get((contig, *_trim(record.pos, reference, alternate)), ())
                        if not matches:
                            continue
                        observed = variant_key(family["assembly"], contig, record.pos, reference, alternate)
                        calls = tuple(
                            _call(
                                role,
                                sample,
                                record.samples[sample],
                                names,
                                columns[sample_columns[sample]],
                            )
                            for role, sample in roles
                        )
                        for target in matches:
                            if target in found:
                                code = (
                                    "duplicate_target_allele"
                                    if found[target]["vcf_variant_key"] == observed
                                    else "ambiguous_padding_equivalence"
                                )
                                raise _error("vcf.records", code)
                            found[target] = _observation(
                                target, observed, alternate_index, calls, columns, raw_record
                            )
                try:
                    next(raw_iterator)
                except StopIteration:
                    pass
                else:
                    raise _error("vcf.records", "raw_and_parsed_record_mismatch")
    except MendelPrimeError:
        raise
    except (OSError, EOFError, NotImplementedError, ValueError, TypeError):
        raise _error("vcf", "invalid_or_unreadable_vcf") from None
    return {
        "case_key": case_key,
        "observations": tuple(found[target] for target in requested if target in found),
        "missing_targets": tuple(target for target in requested if target not in found),
    }


def _write_selected(source, output, header, samples):
    indexes = tuple(header["samples"].index(sample) for sample in samples)
    column_header = False
    with _decompressed(source) as stream:
        for raw_line in stream:
            if raw_line.startswith(b"##"):
                output.write(raw_line)
                continue
            ending = b"\r\n" if raw_line.endswith(b"\r\n") else b"\n" if raw_line.endswith(b"\n") else b""
            content = raw_line[: -len(ending)] if ending else raw_line
            columns = content.split(b"\t")
            if raw_line.startswith(b"#CHROM\t"):
                column_header = True
            elif not content:
                output.write(raw_line)
                continue
            else:
                _check(column_header, "vcf.header", "column_header_required")
            _check(len(columns) == 9 + len(header["samples"]), "vcf.record", "sample_column_count_mismatch")
            output.write(b"\t".join(columns[:9] + [columns[9 + index] for index in indexes]) + ending)
    _check(column_header, "vcf.header", "column_header_required")


def materialize_samples(source, destination, samples):
    """Atomically publish a no-overwrite VCF containing only requested samples."""
    source_path = _path(source, "vcf.path")
    destination_path = _path(destination, "destination")
    try:
        selected = tuple(samples)
    except TypeError:
        raise _error("samples", "sample_iterable_required") from None
    _check(
        bool(selected) and not any(not isinstance(sample, str) or not sample for sample in selected),
        "samples",
        "nonempty_sample_tuple_required",
    )
    _check(len(selected) == len(set(selected)), "samples", "duplicate_sample")
    _check(not os.path.lexists(destination_path), "destination", "already_exists")
    try:
        source_version, header = _version(source_path), inspect_vcf(source_path)
    except MendelPrimeError:
        raise
    except OSError:
        raise _error("vcf", "invalid_or_unreadable_vcf") from None
    _check(
        not any(sample not in header["samples"] for sample in selected), "samples", "sample_absent_from_vcf"
    )
    temporary = None
    try:
        descriptor, name = tempfile.mkstemp(
            prefix=f".{destination_path.name}.", suffix=".tmp", dir=destination_path.parent
        )
        os.close(descriptor)
        temporary = Path(name)
        if destination_path.suffix.lower() in _COMPRESSED:
            with pysam.BGZFile(str(temporary), mode="wb") as output:
                _write_selected(source_path, output, header, selected)
            with temporary.open("rb") as output:
                os.fsync(output.fileno())
        else:
            with temporary.open("wb") as output:
                _write_selected(source_path, output, header, selected)
                output.flush()
                os.fsync(output.fileno())
        _check(inspect_vcf(temporary)["samples"] == selected, "destination", "sample_projection_failed")
        _unchanged(source_path, source_version)
        try:
            os.link(temporary, destination_path)
        except FileExistsError:
            raise _error("destination", "already_exists") from None
        except OSError:
            raise _error("destination", "atomic_publication_failed") from None
        temporary.unlink()
        temporary = None
        return destination_path
    except MendelPrimeError:
        raise
    except (OSError, EOFError, ValueError):
        raise _error("destination", "materialization_failed") from None
    finally:
        if temporary is not None:
            with suppress(FileNotFoundError):
                temporary.unlink()
