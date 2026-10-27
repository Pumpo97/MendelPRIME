"""Assemble source-bound ACMG inputs from a frozen plain-dictionary seed."""

import hashlib
import json
import re
from collections.abc import Mapping

from mendelprime.domain import (
    MendelPrimeError,
    Source,
    case_tuple,
    iter_sources,
    validate_gene,
    variant_text,
    variant_tuple,
)
from mendelprime.tools.native import variant_is_candidate

_VARIANTS_PER_GENE = 8
_PASS_FILTERS = {None, "", ".", "PASS"}
_SEED_FIELDS = {
    "profile",
    "family",
    "native_genes",
    "observations",
    "family_sources",
    "variant_knowledge",
    "evidence",
    "sources",
}
_REPORT_SEED_FIELDS = _SEED_FIELDS | {"fixed_genes"}


def _key(value):
    return variant_tuple(value, MendelPrimeError("variant_key_required", field="assessment.variant_key"))


def _identity_key(value):
    return value["gene"], _key(value["variant_key"])


def _issue(resource, identity, code):
    return f"{resource}:{identity['gene']}:{variant_text(identity['variant_key'])}:{code}"


def _reader(readers, name, *args, **kwargs):
    action = readers.get(name)
    if not callable(action):
        raise MendelPrimeError("reader_required", field=f"assessment.readers.{name}")
    return action(*args, **kwargs)


def _source(value):
    if isinstance(value, Source):
        return value
    if isinstance(value, Mapping) and isinstance(value.get("source"), Source):
        return value["source"]
    raise MendelPrimeError("source_required", field="assessment.resource")


def _add_source(registry, source):
    current = registry.get(source.source_id)
    if current is not None and current != source:
        raise MendelPrimeError("source_id_collision", field="assessment.sources")
    registry.setdefault(source.source_id, source)


def _bind(registry, bindings, source, role, scope, group):
    _add_source(registry, source)
    if not any(item["source_id"] == source.source_id for item in bindings):
        bindings.append(
            {
                "source_id": source.source_id,
                "role": role,
                "match_scope": scope,
                "dependency_group": group,
            }
        )


def _selection_key(variant):
    details = variant.get("details", {})
    exomiser = variant["source"].provider == "exomiser"
    contributing = exomiser and details.get("contributing") is True
    score = details.get("variant_score") if exomiser else None
    return not exomiser, not contributing, -(score if isinstance(score, (int, float)) else 0.0)


def _select(seed, issues, limit=_VARIANTS_PER_GENE, *, ranking=False, statistics=None, nominations=None):
    output = {}
    for gene in seed["fixed_genes"]:
        variants = [
            variant
            for native_gene in seed["native_genes"]
            if native_gene["gene"] == gene
            for variant in native_gene["variants"]
        ]
        protected = tuple((nominations or {}).get(gene, ()))
        if len(protected) != len(set(protected)) or (limit is not None and len(protected) > limit):
            raise MendelPrimeError("nominated_alleles_exceed_review_budget", field="assessment.selection")
        candidates = (
            [
                variant
                for variant in variants
                if variant_is_candidate(variant) or variant.get("key") in protected
            ]
            if ranking or nominations is not None
            else variants
        )
        for variant in candidates:
            if variant.get("key") is None:
                issues.append(f"selection:{gene}:native_variant_key_missing")
        positions = {key: index for index, key in enumerate(protected)}
        ordered = sorted(
            (item for item in candidates if item.get("key") is not None),
            key=lambda item: (positions.get(_key(item["key"]), len(positions)), _selection_key(item)),
        )
        selected, seen = [], set()
        for primary in ordered:
            key = _key(primary["key"])
            if (gene, key) in seen:
                continue
            seen.add((gene, key))
            aliases = tuple(
                item
                for item in variants
                if item is not primary and item.get("key") is not None and _key(item["key"]) == key
            )
            identity = {"gene": gene, "variant_key": key}
            selected.append({"primary": primary, "aliases": aliases, "identity": identity})
            if limit is not None and len(selected) == limit:
                break
        if not selected and (not ranking or candidates):
            issues.append(f"selection:{gene}:no_keyed_native_variants")
        if statistics is not None:
            statistics[gene] = {
                "raw_records": len(variants),
                "candidate_records": len(candidates),
                "candidate_alleles": len(selected),
                "background_lirical_records": len(variants) - len(candidates),
                "background_zero_score_records": sum(
                    variant["source"].provider == "lirical"
                    and type(variant.get("details", {}).get("pathogenicity")) in (int, float)
                    and variant["details"]["pathogenicity"] == 0
                    for variant in variants
                ),
            }
        output[gene] = tuple(selected)
        if not set(protected) <= {item["identity"]["variant_key"] for item in selected}:
            raise MendelPrimeError("nominated_allele_not_observed", field="assessment.selection")
        if nominations is not None:
            omitted = {(_key(v["key"])) for v in candidates if v.get("key")} - {
                item["identity"]["variant_key"] for item in selected
            }
            issues.extend(
                f"selection:{gene}:{variant_text(key)}:deferred_by_review_budget" for key in sorted(omitted)
            )
    return output


def _assess_call(call, alternate_index, alternate_count, site_passed):
    alleles = tuple(call.get("alleles", ()))
    if any(item is not None and item > alternate_count for item in alleles):
        genotype = "invalid_genotype"
    elif alleles and all(item is not None for item in alleles):
        genotype = "called"
    elif any(item is not None for item in alleles):
        genotype = "partial_call"
    else:
        genotype = "no_call"
    if call.get("filter_text") not in _PASS_FILTERS:
        quality = "sample_filtered"
    elif call.get("dp") is None or call.get("gq") is None:
        quality = "missing_quality"
    elif call["dp"] < 10 or call["gq"] < 20:
        quality = "low_quality"
    else:
        quality = "quality_available"
    warnings = (
        ()
        if call.get("ad") is None or len(call["ad"]) == alternate_count + 1
        else ("AD_cardinality_invalid",)
    )
    return {
        "original": call,
        "genotype_state": genotype,
        "quality_state": quality,
        "usable": genotype == "called" and quality == "quality_available" and not warnings and site_passed,
        "carried": alternate_index in alleles,
        "selected_alt_dosage": alleles.count(alternate_index),
        "warnings": warnings,
    }


def _parent_sex_conflict(family):
    pedigree = family.get("pedigree")
    if not pedigree:
        return False
    members = {item["sample_id"]: item for item in pedigree["members"]}
    for role, expected in (("father", "MALE"), ("mother", "FEMALE")):
        sample = family["roles"].get(role)
        member = members.get(sample)
        if member and member["sex"] not in {"UNKNOWN", expected}:
            return True
    return False


def _assess_observation(observation, family):
    if _key(observation["variant_key"])[0] != family["assembly"]:
        raise MendelPrimeError("family_assembly_mismatch", field="observation.variant_key")
    if type(observation.get("alternate_index")) is not int or observation["alternate_index"] < 1:
        raise MendelPrimeError("positive_index_required", field="observation.alternate_index")
    expected = {"proband": family["roles"]["proband"]}
    if family["mode"] == "trio":
        expected.update(father=family["roles"]["father"], mother=family["roles"]["mother"])
    roles = [call["role"] for call in observation["samples"]]
    samples = [call["sample_id"] for call in observation["samples"]]
    if "proband" not in roles or len(roles) != len(set(roles)) or len(samples) != len(set(samples)):
        raise MendelPrimeError("unique_roles_and_samples_required", field="observation.samples")
    for call in observation["samples"]:
        if call["role"] not in expected or call["sample_id"] != expected[call["role"]]:
            raise MendelPrimeError("family_role_sample_mismatch", field="observation.samples")
    try:
        fields = observation["raw_record"].rstrip("\r\n").split("\t", 8)
        alternate_field = fields[4].split(",")
        alternate_count = len(alternate_field)
        chromosome = fields[0].removeprefix("chr")
        chromosome = "MT" if chromosome == "M" else chromosome
        observed = _key(observation.get("vcf_variant_key") or observation["variant_key"])
        raw_identity = (
            observed[0],
            chromosome,
            int(fields[1]),
            fields[3].upper(),
            alternate_field[observation["alternate_index"] - 1].upper(),
        )
    except (IndexError, ValueError, AttributeError):
        raise MendelPrimeError("raw_allele_identity_mismatch", field="observation.raw_record") from None
    if raw_identity != observed:
        raise MendelPrimeError("raw_allele_identity_mismatch", field="observation.raw_record")
    target = _key(observation["variant_key"])

    def trimmed(key):
        _, chromosome, position, reference, alternate = key
        while len(reference) > 1 and len(alternate) > 1 and reference[-1] == alternate[-1]:
            reference, alternate = reference[:-1], alternate[:-1]
        while len(reference) > 1 and len(alternate) > 1 and reference[0] == alternate[0]:
            reference, alternate, position = reference[1:], alternate[1:], position + 1
        return chromosome, position, reference, alternate

    if observed[0] != target[0] or trimmed(observed) != trimmed(target):
        raise MendelPrimeError("padding_equivalence_required", field="observation.vcf_variant_key")
    calls = tuple(
        _assess_call(
            call,
            observation["alternate_index"],
            alternate_count,
            observation.get("filter_text") in {".", "PASS"},
        )
        for call in observation["samples"]
    )
    by_role = {item["original"]["role"]: item for item in calls}
    proband = by_role["proband"]
    status = (
        "carried"
        if proband["carried"] and proband["usable"]
        else "carried_unverified"
        if proband["carried"]
        else "not_carried"
        if proband["genotype_state"] == "called"
        else proband["genotype_state"]
    )
    origin = "unresolved"
    if (
        family["mode"] == "trio"
        and set(by_role) == {"proband", "father", "mother"}
        and not _parent_sex_conflict(family)
        and all(item["usable"] for item in by_role.values())
        and proband["carried"]
    ):
        father, mother = by_role["father"]["carried"], by_role["mother"]["carried"]
        origin = (
            "paternal"
            if father and not mother
            else "maternal"
            if mother and not father
            else "candidate_de_novo"
            if not father and not mother
            else "both_parents_carry"
        )
    return {
        "original": observation,
        "calls": calls,
        "status": status,
        "origin": origin,
        "parentage_status": "unconfirmed"
        if family["mode"] == "trio" and {"father", "mother"} <= set(by_role)
        else "unavailable",
    }


def _family_source(sources, gene, key):
    matches = tuple(
        source
        for source in sources
        if source.kind == "patient_observation" and source.genes == (gene,) and key in source.variant_keys
    )
    if len(matches) > 1:
        raise MendelPrimeError("ambiguous_family_source", field="assessment.patient_observation")
    return matches[0] if matches else None


def _excluded_facts(seed, identity, assessed, reference, readers, issues, *, reference_attempted):
    if assessed is not None and not reference_attempted:
        try:
            reference = _reader(readers, "resolve", identity["variant_key"])
        except MendelPrimeError:
            raise
        except Exception:
            issues.append(_issue("reference", identity, "reference_lookup_failed"))
        else:
            if _key(reference.get("raw_key")) != identity["variant_key"]:
                raise MendelPrimeError("target_result_mismatch", field="assessment.reference")
            if reference["status"] != "sequence_verified":
                issues.append(
                    _issue("reference", identity, (reference.get("issues") or (reference["status"],))[0])
                )
    case_key = case_tuple(seed["profile"]["case_key"])
    source_ids = tuple(
        dict.fromkeys(
            source.source_id
            for source in seed["family_sources"]
            if source.kind == "patient_observation"
            and source.case_key == case_key
            and source.genes == (identity["gene"],)
            and identity["variant_key"] in source.variant_keys
        )
    )
    return {
        "observation": assessed,
        "reference": reference,
        "source_ids": source_ids,
        "identity_status": "verified"
        if assessed is not None
        and reference is not None
        and reference["status"] == "sequence_verified"
        and len(source_ids) == 1
        else "unknown",
        "absence_is_not_benign_evidence": True,
    }


def _eligible(seed, selected, readers, issues, excluded, *, retain_excluded=False):
    observations = {_key(item["variant_key"]): item for item in seed["observations"]["observations"]}
    output = []
    for gene in seed["fixed_genes"]:
        for item in selected[gene]:
            identity = item["identity"]
            assessed, reference, reference_attempted = None, None, False
            observation = observations.get(identity["variant_key"])
            if observation is None:
                reason = "observation_missing"
            else:
                assessed = _assess_observation(observation, seed["family"])
                if assessed["status"] not in {"carried", "carried_unverified"}:
                    reason = f"observation_{assessed['status']}"
                else:
                    reference_attempted = True
                    try:
                        reference = _reader(readers, "resolve", identity["variant_key"])
                    except MendelPrimeError:
                        raise
                    except Exception:
                        reason = "reference_lookup_failed"
                    else:
                        if _key(reference.get("raw_key")) != identity["variant_key"]:
                            raise MendelPrimeError("target_result_mismatch", field="assessment.reference")
                        if reference["status"] == "reference_mismatch":
                            reason = "reference_mismatch"
                        else:
                            if reference["status"] != "sequence_verified":
                                issues.append(
                                    _issue(
                                        "reference",
                                        identity,
                                        (reference.get("issues") or ("reference_unresolved",))[0],
                                    )
                                )
                            source = _family_source(seed["family_sources"], gene, identity["variant_key"])
                            if source is not None:
                                output.append(
                                    {
                                        **item,
                                        "observation": assessed,
                                        "reference": reference,
                                        "family_source": source,
                                    }
                                )
                                continue
                            reason = "patient_observation_source_missing"
            if retain_excluded:
                facts = _excluded_facts(
                    seed,
                    identity,
                    assessed,
                    reference,
                    readers,
                    issues,
                    reference_attempted=reference_attempted,
                )
                excluded.append((identity, reason, facts))
            else:
                excluded.append((identity, reason))
            issues.append(_issue("excluded", identity, reason))
    return tuple(output)


def _validated_batch(name, requested, output, row_type):
    if type(output) is not tuple or len(output) != len(requested):
        raise MendelPrimeError("complete_outcome_tuple_required", field=f"assessment.{name}")
    for target, value in zip(requested, output, strict=True):
        if not isinstance(value, row_type) or value.get("target") != target:
            raise MendelPrimeError("target_result_mismatch", field=f"assessment.{name}")
        yield value


def _batch(readers, name, values, expected, issues, accepted):
    try:
        output = _reader(readers, name, values)
    except MendelPrimeError:
        raise
    except Exception:
        for item in expected:
            issues.append(_issue(name, item, f"{name}_lookup_failed"))
        return {(_identity_key(item)): None for item in expected}
    result = {}
    for identity, value in zip(expected, _validated_batch(name, values, output, dict), strict=True):
        result[_identity_key(identity)] = value
        if value["status"] not in accepted:
            issues.append(_issue(name, identity, value.get("error_code") or value["status"]))
    return result


def _native_annotation_source(variant):
    key = _key(variant["key"])
    details = variant.get("details", {})
    consequence = details.get("consequence", details.get("raw_annotation"))
    content = json.dumps(
        {
            "record_type": "native_variant_annotation",
            "provider": variant["source"].provider,
            "variant_key": variant_text(key),
            "transcript": variant.get("transcript"),
            "hgvs_c": variant.get("hgvs_c"),
            "hgvs_p": variant.get("hgvs_p"),
            "consequence": consequence,
        },
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return Source(
        provider=variant["source"].provider,
        kind="derived",
        locator=f"method:assessment/native-annotation/{variant['source'].provider}/{variant['source'].content_sha256[:16]}",
        version=variant["source"].version,
        content=content,
        genes=(variant["gene"],),
        variant_keys=(key,),
        case_key=variant["source"].case_key,
        parent_source_ids=(variant["source"].source_id,),
    )


def _mentions(target, source):
    values = [target["reference"].get("normalized_hgvs")]
    for variant in (target["primary"], *target.get("aliases", ())):
        values.extend((variant.get("hgvs_c"), variant.get("hgvs_p")))
    text = source.content.casefold()
    return any(
        str(value).split(":")[-1].strip("() ").casefold() in text
        for value in values
        if value and len(str(value)) >= 5
    )


def _covered_absence(source):
    return {
        "code": "PM2",
        "status": "met",
        "strength": "supporting",
        "rationale": (
            "The allele was not detected in the population resource and every affected "
            "reference position met the preserved coverage threshold."
        ),
        "citations": ({"source_id": source.source_id, "quote": source.content},),
        "independence_group": "population_frequency",
        "origin": "verified_population_observation",
    }


def _knowledge_issues(resource, result, issues, *, refreshed=()):
    current = {(row["provider"], row.get("variant_key")) for row in refreshed}
    for lookup in result.get("lookups", ()):
        if (lookup["provider"], lookup.get("variant_key")) not in current and lookup["status"] not in {
            "found", "not_found", "not_applicable"
        }:
            target = (
                variant_text(_key(lookup["variant_key"])) if lookup.get("variant_key") else result["gene"]
            )
            issues.append(
                f"{resource}:{lookup['provider']}:{target}:{lookup.get('error_code') or lookup['status']}"
            )


def _clinical_binding(evidence, identity, registry, case_key):
    original = _source(evidence["document"])
    if identity["variant_key"] in original.variant_keys:
        return original
    if evidence.get("match_basis") != "sequence_equivalent" or not evidence.get("normalized_hgvs"):
        raise MendelPrimeError("allele_source_scope_mismatch", field="assessment.clinical")
    content = json.dumps(
        {
            "record_type": "verified_sequence_equivalent_clinical_evidence",
            "binding_target_key": variant_text(identity["variant_key"]),
            "source_variant_keys": [variant_text(item) for item in original.variant_keys],
            "normalized_hgvs": evidence["normalized_hgvs"],
            "original_source_id": original.source_id,
            "original_content": original.content,
        },
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    projected = Source(
        provider="mendelprime",
        kind="derived",
        locator=f"method:assessment/sequence-equivalent/{original.source_id.removeprefix('SRC_')[:16]}",
        version="reference-equivalence-v1",
        content=content,
        genes=(identity["gene"],),
        variant_keys=(identity["variant_key"],),
        case_key=case_key,
        parent_source_ids=(original.source_id,),
    )
    _add_source(registry, projected)
    return projected


def _diseases(seed, gene):
    labels, identifiers, inheritance = set(), set(), set()
    for native in seed["native_genes"]:
        if native["gene"] != gene:
            continue
        for record in native["records"]:
            for name, values in (("moi", inheritance), ("disease_name", labels), ("disease_id", identifiers)):
                value = record.get(name)
                if isinstance(value, str) and value.strip():
                    values.add(value.strip())
        for variant in native["variants"]:
            details = variant.get("details", {})
            for name, values in (
                ("acmg_disease_name", labels),
                ("acmg_disease_id", identifiers),
                ("moi", inheritance),
                ("disease_id", identifiers),
            ):
                value = details.get(name)
                if isinstance(value, str) and value.strip():
                    values.add(value.strip())
    return tuple(sorted(labels)), tuple(sorted(identifiers)), tuple(sorted(inheritance))


def _genereviews_context(fact, seed, gene):
    """A case-local, disease-specific transcript choice, never shared mutable state."""
    from mendelprime.tools.acmg import _ranking_mode

    _, identifiers, inheritance = _diseases(seed, gene)
    diseases = {"OMIM:" + value for value in fact.get("disease_ids", ())}.intersection(identifiers)
    modes = {_ranking_mode(value) for value in inheritance}
    modes = {"XL" if mode in {"XD", "XR"} else mode for mode in modes} - {"unknown"}
    if (
        fact.get("record_type") == "genereviews_mechanism_v1"
        and fact.get("gene") == gene
        and fact.get("mechanism_original") == "Mechanism of disease causation. Loss of function"
        and len(diseases) == 1
        and fact.get("inheritance") in modes
        and len(fact.get("relevant_transcripts", ())) == 1
    ):
        return {**fact, "disease_id": next(iter(diseases))}
    return None


def _disease_contexts(seed, item, registry, issues):
    gene = item["identity"]["gene"]
    contexts = []

    def add(record, source, *, variant=False):
        names = ("acmg_disease_id", "acmg_disease_name") if variant else ("disease_id", "disease_name")
        if variant and not any(record.get(name) for name in names):
            names = ("disease_id", "disease_name")
        values = {
            "disease_id": record.get(names[0]),
            "disease_name": record.get(names[1]),
            "moi": record.get("moi"),
        }
        values = {
            name: value.strip() if isinstance(value, str) and value.strip() else None
            for name, value in values.items()
        }
        if not any(values.values()):
            return
        if not isinstance(source, Source):
            issues.append(f"disease_context:{gene}:source_missing")
            return
        if gene not in source.genes:
            raise MendelPrimeError("gene_source_scope_mismatch", field="assessment.disease_contexts")
        _add_source(registry, source)
        row = {**values, "source_ids": (source.source_id,)}
        if row not in contexts:
            contexts.append(row)

    for native in seed["native_genes"]:
        if native["gene"] == gene:
            for record in native["records"]:
                add(record, record.get("source"))
    for variant in (item["primary"], *item["aliases"]):
        add(variant.get("details", {}), variant["source"], variant=True)
    return tuple(contexts)


def _target(
    seed, item, annotation, pvs1, clinical, gene_context, guidance, readers, registry, issues, *, post=False
):
    primary, identity = item["primary"], item["identity"]
    bindings = []

    def lookup(name, code, *args, **kwargs):
        try:
            return _reader(readers, name, *args, **kwargs)
        except MendelPrimeError:
            raise
        except Exception:
            issues.append(_issue(name, identity, code))
            return None

    _bind(
        registry,
        bindings,
        item["family_source"],
        "patient_observation",
        "exact_allele",
        f"patient_observation:{identity['gene']}:{variant_text(identity['variant_key'])}",
    )
    for variant in (primary, *item["aliases"]):
        _add_source(registry, variant["source"])
        source = _native_annotation_source(variant)
        _bind(
            registry,
            bindings,
            source,
            "native_annotation",
            "exact_allele",
            f"native_annotation:{variant['source'].provider}:{variant_text(identity['variant_key'])}",
        )
    for result in clinical:
        for evidence in result.get("variant_evidence", ()):
            if _key(evidence["variant_key"]) != identity["variant_key"]:
                continue
            original = _source(evidence["document"])
            _add_source(registry, original)
            if evidence.get("retracted", False):
                issues.append(_issue("clinical", identity, "retracted_assertion"))
                continue
            bound = _clinical_binding(evidence, identity, registry, case_tuple(seed["profile"]["case_key"]))
            verified = identity["variant_key"][1] != "MT" or (
                item["reference"]["status"] == "sequence_verified"
                and item["reference"]["reference_accession"] == "NC_012920.1"
            )
            _bind(
                registry,
                bindings,
                bound,
                "allele_assertion" if verified else "gene_context",
                "exact_allele" if verified else "context",
                evidence.get("linked_source_id") or original.source_id,
            )
        for document in result.get("gene_context", ()):
            source = _source(document)
            _bind(registry, bindings, source, "gene_disease_context", "context", source.source_id)
    if gene_context:
        for document in gene_context.get("gene_context", ()):
            source = _source(document)
            _bind(registry, bindings, source, "gene_context", "context", source.source_id)
    run_context = seed.get("evidence", {}).get(identity["gene"])
    mention_target = {
        "primary": primary,
        "aliases": item["aliases"],
        "reference": item["reference"],
    }
    if run_context:
        if post:
            from mendelprime.knowledge.literature import _variant_article_key, _variant_search_terms

            terms = _variant_search_terms(mention_target, report=True)
        for error in run_context.get("retrieval_errors", ()):
            issues.append(f"literature:{identity['gene']}:{error}")
        for document in run_context.get("documents", ()):
            source = _source(document)
            if post and source.variant_keys and identity["variant_key"] not in source.variant_keys:
                continue
            if post and source.kind == "literature" and not source.variant_keys:
                mentions = any(
                    re.search(
                        r"(?<![A-Za-z0-9])" + re.escape(term) + r"(?![A-Za-z0-9])", source.content, re.I
                    )
                    for term in terms
                )
                if not mentions and re.search(
                    r"\b(?:p\.\(?)?[A-Z](?:[a-z]{2})?\d+[A-Z](?:[a-z]{2})?\b", source.content
                ):
                    continue
            if source.kind == "literature":
                article = _variant_article_key(document) if post else ""
                _bind(
                    registry,
                    bindings,
                    source,
                    "literature",
                    "exact_allele"
                    if post and source.variant_keys
                    else "candidate_allele"
                    if (mentions if post else _mentions(mention_target, source))
                    else "context",
                    f"publication:{article}" if article else source.source_id,
                )
            elif post and source.kind == "annotation" and source.variant_keys:
                if source.provider == "publication-family-ledger":
                    _bind(
                        registry,
                        bindings,
                        source,
                        "family_series_ledger",
                        "exact_allele",
                        json.loads(source.content)["article_source_id"],
                    )
                else:
                    _bind(registry, bindings, source, "allele_identity", "exact_allele", source.source_id)
        for document in run_context.get("omim_context", ()):
            source = _source(document)
            _bind(registry, bindings, source, "gene_disease_context", "context", source.source_id)
    population = lookup("population", "population_lookup_failed", identity["gene"], identity["variant_key"])
    preserved = ()
    if population:
        if _key(population["original_key"]) != identity["variant_key"]:
            raise MendelPrimeError("target_result_mismatch", field="assessment.population")
        if population["status"] not in {"found", "not_found_covered"}:
            issues.append(
                _issue("population", identity, population.get("error_code") or population["status"])
            )
        if population.get("source"):
            _bind(
                registry, bindings, population["source"], "population", "exact_allele", "population_frequency"
            )
            if population["status"] == "not_found_covered":
                preserved = (_covered_absence(population["source"]),)
    if annotation and annotation.get("annotation"):
        source = annotation["annotation"]["source"]
        _bind(registry, bindings, source, "transcript_annotation", "exact_allele", source.source_id)
    predictions = lookup("predictions", "prediction_lookup_failed", identity["gene"], identity["variant_key"])
    if predictions:
        if (
            predictions["gene"] != identity["gene"]
            or _key(predictions["variant_key"]) != identity["variant_key"]
        ):
            raise MendelPrimeError("target_result_mismatch", field="assessment.predictions")
        group = f"computational_prediction:{identity['gene']}:{variant_text(identity['variant_key'])}"
        for name in ("revel", "spliceai"):
            prediction = predictions[name]
            if prediction["status"] not in {"found", "not_found"}:
                issues.append(
                    _issue(
                        "predictions",
                        identity,
                        prediction.get("error_code") or f"{name}_{prediction['status']}",
                    )
                )
            if prediction.get("source"):
                _bind(registry, bindings, prediction["source"], "prediction", "exact_allele", group)
    resolved = (pvs1 or {}).get("resolved_annotation") or ((pvs1 or {}).get("annotation") or {}).get(
        "resolved_annotation"
    )
    if pvs1 and pvs1.get("annotation"):
        source = pvs1["annotation"]["source"]
        _bind(registry, bindings, source, "pvs1", "exact_allele", f"pvs1:{source.source_id}")
    if resolved:
        group = (
            pvs1["annotation"]["source"].source_id if pvs1.get("annotation") else resolved["source"].source_id
        )
        _bind(registry, bindings, resolved["source"], "transcript_resolution", "exact_allele", group)
    if pvs1 and pvs1.get("annotation"):
        annotation_value = pvs1["annotation"]
        if (
            resolved
            and (
                (
                    pvs1["status"] == "unknown"
                    and pvs1.get("error_code") == "pvs1_strength_unset"
                    and annotation_value.get("strength") == "Unset"
                )
                or (pvs1["status"] == "assessed" and annotation_value.get("met") is True)
            )
            and annotation_value.get("decision_tree_supported") is True
            and annotation_value.get("gene_match") in {"exact", "resource_alias"}
            and annotation_value.get("strength_raw") in {"VeryStrong", "Strong", "Moderate", "Supporting"}
        ):
            from mendelprime.tools.acmg import _ranking_mode

            modes = {
                "XL" if mode in {"XD", "XR"} else mode
                for value in _diseases(seed, identity["gene"])[2]
                if (mode := _ranking_mode(value)) != "unknown"
            }
            mechanisms, dosages = [], []
            for binding in bindings:
                if binding["role"] != "gene_disease_context":
                    continue
                mechanism_source = registry[binding["source_id"]]
                if mechanism_source.provider not in (
                    {"clingen", "genereviews"} if post else {"clingen"}
                ) or mechanism_source.genes != (identity["gene"],):
                    continue
                try:
                    fact = json.loads(mechanism_source.content)
                except ValueError:
                    continue
                summary = fact.get("evidence_summary", "")
                statement = re.search(
                    r"\b(?:mechanism (?:of (?:disease|pathogenicity)|for disease)|disease mechanism)\b[^.]{0,120}\b(?:is|involves)\b[^.]{0,60}\b(?:haploinsufficiency|loss[- ]of[- ]function)\b",
                    summary,
                    re.I,
                )
                explicit = bool(
                    statement
                    and not re.search(
                        r"\b(?:not|unknown|unclear|may|might|suspected|possible)\b", statement[0], re.I
                    )
                    and fact.get("classification") in {"Strong", "Definitive"}
                    and fact.get("disease_id")
                    and fact.get("inheritance") in {"AD", "AR", "XL"}
                    and (not modes or fact["inheritance"] in modes)
                )
                if mechanism_source.provider == "genereviews":
                    context = _genereviews_context(fact, seed, identity["gene"])
                    explicit = bool(
                        context
                        and resolved["transcript_id"].split(".")[0] == context["relevant_transcripts"][0]
                    )
                    if explicit:
                        fact, statement = context, (context["mechanism_original"],)
                dosage = fact.get("row", {}) if fact.get("record_type") == "clingen_dosage" else {}
                dosage_confirmed = (
                    dosage.get("Haploinsufficiency Score") == "3"
                    and bool(dosage.get("Haploinsufficiency Disease ID"))
                    and len(modes) == 1
                    and modes <= {"AD", "XL"}
                )
                if not (explicit or dosage_confirmed):
                    continue
                disease = (
                    {
                        "identifier": fact["disease_id"],
                        "label": fact.get("disease_label", ""),
                        "inheritance": fact["inheritance"],
                    }
                    if explicit
                    else {
                        "identifier": dosage["Haploinsufficiency Disease ID"],
                        "label": "ClinGen curated dosage-sensitive disorder",
                        "inheritance": next(iter(modes)),
                    }
                )
                (mechanisms if explicit else dosages).append(
                    (
                        disease,
                        mechanism_source,
                        statement[0] if explicit else '"Haploinsufficiency Score": "3"',
                    )
                )
            choices = mechanisms or dosages
            contexts = {(d["identifier"], d["inheritance"]) for d, _, _ in choices}
            if len(contexts) == 1:
                disease, mechanism_source, quote = choices[0]
                strength_name = (
                    annotation_value.get("strength")
                    if annotation_value.get("met") is True
                    else annotation_value["strength_raw"]
                )
                strength = {
                    "VeryStrong": "very_strong",
                    "Strong": "strong",
                    "Moderate": "moderate",
                    "Supporting": "supporting",
                }[strength_name]
                preserved = (
                    *preserved,
                    {
                        "code": "PVS1",
                        "status": "met",
                        "strength": strength,
                        "disease_binding": disease,
                        "rationale": "Independent disease-specific LoF mechanism plus exact-transcript AutoPVS1 variant decision tree; an uninterpreted package gene-level suffix is not rewritten."
                        if mechanism_source.provider == "genereviews"
                        else "Independent ClinGen disease-specific LoF mechanism plus exact-transcript AutoPVS1 variant decision tree; an uninterpreted package gene-level suffix is not rewritten.",
                        "citations": (
                            {"source_id": source.source_id, "quote": source.content},
                            {
                                "source_id": mechanism_source.source_id,
                                "quote": quote,
                            },
                        ),
                        "independence_group": "validated_pvs1_mechanism_and_allele_effect",
                        "origin": "verified_mechanism_and_variant_tree",
                    },
                )
            elif choices:
                issues.append(_issue("pvs1", identity, "multiple_curated_disease_contexts"))
    for source in seed.get("sources", ()):
        if source.provider == "plink2-king" and source.case_key == case_tuple(seed["profile"]["case_key"]):
            _bind(registry, bindings, source, "observation_qc", "context", "family_relationship_qc")
    chosen = annotation.get("annotation") if annotation else None
    protein = lookup(
        "protein_context",
        "protein_context_lookup_failed",
        identity["gene"],
        variant_key=identity["variant_key"],
        transcript=chosen.get("transcript_id") if chosen else primary.get("transcript"),
        protein_hgvs=chosen.get("hgvsp") if chosen else primary.get("hgvs_p"),
    )
    if protein:
        if protein["gene"] != identity["gene"]:
            raise MendelPrimeError("gene_result_mismatch", field="assessment.protein_context")
        _knowledge_issues("protein_context", protein, issues)
        for document in protein.get("gene_context", ()):
            source = _source(document)
            _bind(registry, bindings, source, "protein_comparator", "context", source.source_id)
    for source in guidance:
        if not source.genes or identity["gene"] in source.genes:
            _bind(registry, bindings, source, "guideline", "context", source.source_id)
    labels, identifiers, inheritance = _diseases(seed, identity["gene"])
    return {
        "primary": primary,
        "observation": item["observation"],
        "reference": item["reference"],
        "population": population,
        "annotation": annotation,
        "predictions": predictions,
        "pvs1": pvs1,
        "bindings": tuple(bindings),
        "preserved_criteria": preserved,
        "aliases": item["aliases"],
        "disease_labels": labels,
        "disease_ids": identifiers,
        "inheritance_hypotheses": inheritance,
        "family_conflict": False,
    }


def validate_seed(seed):
    if not isinstance(seed, Mapping) or set(seed) != _REPORT_SEED_FIELDS:
        raise MendelPrimeError("schema_mismatch", field="assessment.seed")
    if type(seed["fixed_genes"]) is not tuple or not 1 <= len(seed["fixed_genes"]) <= 5:
        raise MendelPrimeError("top1_to_top5_tuple_required", field="assessment.seed.fixed_genes")
    if type(seed["sources"]) is not tuple or any(not isinstance(item, Source) for item in seed["sources"]):
        raise MendelPrimeError("source_tuple_required", field="assessment.seed.sources")
    return seed


def _gene_context(readers, gene, issues):
    try:
        context = _reader(readers, "gene_context", gene)
    except MendelPrimeError:
        raise
    except Exception:
        issues.append(f"gene_context:{gene}:gene_context_lookup_failed")
        return None
    if context:
        if context["gene"] != gene:
            raise MendelPrimeError("gene_result_mismatch", field="assessment.gene_context")
        _knowledge_issues("gene_context", context, issues)
    return context


def _assemble(seed, readers, *, report, review=None):
    if not isinstance(readers, Mapping):
        raise MendelPrimeError("mapping_required", field="assessment.readers")
    frozen_genes = seed["fixed_genes"]
    issues, excluded = [], []
    nomination_counts = {} if not report else None
    selected = _select(
        seed,
        issues,
        _VARIANTS_PER_GENE if report else None,
        ranking=not report,
        statistics=nomination_counts,
        nominations=review["nominations"] if report and review else None,
    )
    eligible = _eligible(seed, selected, readers, issues, excluded, retain_excluded=not report)
    if report and review and callable(readers.get("complete_report_evidence")):
        seed, extra_issues = readers["complete_report_evidence"](seed, review, eligible)
        issues.extend(extra_issues)
    if report and not any(item["primary"]["gene"] == seed["fixed_genes"][0] for item in eligible):
        raise MendelPrimeError("top1_verified_observed_allele_required", field="assessment.targets")
    genes = tuple(dict.fromkeys(item["identity"]["gene"] for item in eligible)) if report else frozen_genes
    contexts, context_issues = {}, {}
    for gene in genes if report and review else ():
        context_issues[gene] = []
        contexts[gene] = _gene_context(readers, gene, context_issues[gene])
    identities = tuple(item["identity"] for item in eligible)
    requests = tuple(
        {
            "gene": item["identity"]["gene"],
            "variant_key": item["identity"]["variant_key"],
            "transcript_id": item["primary"].get("transcript"),
        }
        for item in eligible
    )
    accepted = {"annotated"} if report else {"annotated", "not_applicable"}
    annotations = _batch(readers, "annotate", requests, identities, issues, accepted)
    preferred, mane = {}, {}
    if report and review:
        for gene, knowledge in seed["variant_knowledge"].items():
            choices = set()
            for document in knowledge.get("gene_context", ()):
                source = document["source"]
                if source.provider == "genereviews" and source.genes == (gene,):
                    context = _genereviews_context(json.loads(source.content), seed, gene)
                    if context:
                        choices.add(
                            (
                                context["disease_id"],
                                context["inheritance"],
                                context["relevant_transcripts"][0],
                            )
                        )
            if len(choices) == 1:
                preferred[gene] = (next(iter(choices))[2],)
            if not choices:
                from mendelprime.knowledge.annotation import mane_transcripts

                candidates = mane_transcripts(
                    (document["source"] for document in (contexts.get(gene) or {}).get("gene_context", ())),
                    gene,
                )
                if candidates:
                    mane[gene] = candidates
    pvs1 = _batch(
        readers,
        "pvs1",
        tuple(
            {
                **item["primary"],
                "aliases": item["aliases"],
                **(
                    {"preferred_transcripts": preferred[item["primary"]["gene"]]}
                    if item["primary"]["gene"] in preferred
                    else {}
                ),
                **(
                    {"mane_candidates": mane[item["primary"]["gene"]]}
                    if item["primary"]["gene"] in mane
                    else {}
                ),
            }
            if report and review
            else item["primary"]
            for item in eligible
        ),
        identities,
        issues,
        {"assessed", "not_applicable"},
    )
    if report and review:
        for identity, value in pvs1.items():
            resolved = (value or {}).get("resolved_annotation") or (
                (value or {}).get("annotation") or {}
            ).get("resolved_annotation")
            if resolved and (annotations[identity] or {}).get("status") != "annotated":
                annotations[identity] = {"status": "annotated", "error_code": "", "annotation": resolved}
                prefix = f"annotate:{identity[0]}:{variant_text(identity[1])}:"
                issues[:] = [code for code in issues if not code.startswith(prefix)]
    registry = {}
    for source in seed["sources"]:
        _add_source(registry, source)
    if not report:
        excluded_ids = {source_id for _, _, facts in excluded for source_id in facts["source_ids"]}
        for source in seed["family_sources"]:
            if source.source_id in excluded_ids:
                _add_source(registry, source)
    guidance_value = readers.get("guidance")
    guidance = guidance_value(tuple(seed["fixed_genes"])) if callable(guidance_value) else guidance_value
    if type(guidance) is not tuple or any(
        not isinstance(item, Source) or item.kind != "guideline" for item in guidance
    ):
        raise MendelPrimeError("guideline_source_tuple_required", field="assessment.readers.guidance")
    guidance = tuple(
        source for source in guidance if not source.genes or any(gene in source.genes for gene in genes)
    )
    for source in guidance:
        _add_source(registry, source)
    clinical = {}
    for gene in genes:
        keys = (
            tuple(item["identity"]["variant_key"] for item in eligible if item["identity"]["gene"] == gene)
            if report
            else tuple(item["identity"]["variant_key"] for item in selected[gene])
        )
        existing = seed["variant_knowledge"].get(gene)
        if existing is not None and existing.get("gene") != gene:
            raise MendelPrimeError("existing_gene_result_mismatch", field="assessment.clinical")
        values = (existing,) if existing else ()
        fresh = {}
        try:
            fresh = _reader(readers, "clinical", gene, keys)
        except MendelPrimeError:
            raise
        except Exception:
            issues.append(f"clinical:{gene}:clinical_lookup_failed")
        else:
            if fresh["gene"] != gene:
                raise MendelPrimeError("gene_result_mismatch", field="assessment.clinical")
            _knowledge_issues("clinical", fresh, issues)
            values = (*values, fresh)
        if existing:
            _knowledge_issues("clinical_existing", existing, issues, refreshed=fresh.get("lookups", ()))
        clinical[gene] = values
        if gene in contexts:
            issues.extend(context_issues[gene])
        else:
            contexts[gene] = _gene_context(readers, gene, issues)
    targets = tuple(
        _target(
            seed,
            item,
            annotations[_identity_key(item["identity"])],
            pvs1[_identity_key(item["identity"])],
            clinical[item["identity"]["gene"]],
            contexts[item["identity"]["gene"]],
            guidance,
            readers,
            registry,
            issues,
            post=bool(report and review),
        )
        for item in eligible
    )
    if not report:
        targets = tuple(
            {**target, "disease_contexts": _disease_contexts(seed, item, registry, issues)}
            for target, item in zip(targets, eligible, strict=True)
        )
        for values in clinical.values():
            for value in values:
                for source in iter_sources(value):
                    _add_source(registry, source)
    packet = {
        "profile": seed["profile"],
        "family": seed["family"],
        "targets": targets,
        "sources": tuple(registry.values()),
    }
    if report or targets:
        from mendelprime.tools.acmg import validate_input

        validate_input(packet)
    else:
        case_key = case_tuple(seed["profile"]["case_key"])
        if any(source.case_key is not None and source.case_key != case_key for source in registry.values()):
            raise MendelPrimeError("case_scope_mismatch", field="assessment.sources")
        if any(parent not in registry for source in registry.values() for parent in source.parent_source_ids):
            raise MendelPrimeError("parent_source_missing", field="assessment.sources")
    if seed["fixed_genes"] != frozen_genes:
        raise MendelPrimeError("fixed_genes_changed_during_assembly", field="assessment.seed")
    result = {
        "input": packet,
        "issues": tuple(dict.fromkeys(issues)),
        "excluded": tuple(excluded),
    }
    if not report:
        result["clinical"] = {gene: _merge_clinical(gene, values) for gene, values in clinical.items()}
        result["nomination_counts"] = nomination_counts
    return result


def build(seed, readers, *, review=None):
    """Prepare the unchanged Top5/8-variants report input without invoking a model."""
    seed = validate_seed(seed)
    result = _assemble(seed, readers, report=True, review=review)
    if review:
        registry = {source.source_id: source for source in result["input"]["sources"]}
        review_sources = {source.source_id: source for source in review["sources"]}
        for target in result["input"]["targets"]:
            target["post_evidence_recovered"] = True
            target["bindings"] = list(target["bindings"])
            identity = (target["primary"]["gene"], _key(target["primary"]["key"]))
            facts = []
            for proposal in review["facts"].get(identity, ()):
                citations = []
                for citation in proposal["citations"]:
                    source = review_sources[citation["source_id"]]
                    if source.kind == "tool_observation" and source.provider in {"exomiser", "lirical"}:
                        matches = [
                            s
                            for s in registry.values()
                            if s.locator.startswith("method:assessment/native-annotation/")
                            and source.source_id in s.parent_source_ids
                            and s.genes == (identity[0],)
                            and identity[1] in s.variant_keys
                        ]
                        if len(matches) != 1:
                            raise MendelPrimeError(
                                "native_prior_fact_projection_ambiguous", field="assessment.review"
                            )
                        source = matches[0]
                        citation = {"source_id": source.source_id, "quote": source.content}
                    citations.append(citation)
                facts.append({**proposal, "citations": tuple(citations)})
            facts = tuple(facts)
            cited = {c["source_id"] for fact in facts for c in fact["citations"]}
            for source_id in sorted(cited):
                source = registry.get(source_id) or review_sources[source_id]
                if source.genes and identity[0] not in source.genes:
                    raise MendelPrimeError("prior_fact_gene_mismatch", field="assessment.review")
                _bind(
                    registry,
                    target["bindings"],
                    source,
                    "prior_criterion_source",
                    "exact_allele" if identity[1] in source.variant_keys else "context",
                    source.source_id,
                )
            target["prior_criteria"] = facts
            target["disease_contexts"] = _disease_contexts(
                seed,
                {
                    "identity": {"gene": identity[0], "variant_key": identity[1]},
                    "primary": target["primary"],
                    "aliases": target["aliases"],
                },
                registry,
                [],
            )
            for context in target["disease_contexts"]:
                parents = tuple(registry[source_id] for source_id in context["source_ids"])
                source = Source(
                    provider=parents[0].provider,
                    kind="derived",
                    locator="method:assessment/disease-context/" + parents[0].source_id,
                    version="native-disease-context-v1",
                    genes=(identity[0],),
                    variant_keys=(identity[1],) if any(s.variant_keys for s in parents) else (),
                    case_key=case_tuple(seed["profile"]["case_key"]),
                    parent_source_ids=context["source_ids"],
                    content=json.dumps(
                        {
                            "gene": identity[0],
                            **{k: context[k] for k in ("disease_id", "disease_name", "moi")},
                            "interpretation": "Native tool disease association, not independent medical evidence or a classification.",
                        },
                        sort_keys=True,
                    ),
                )
                _bind(
                    registry,
                    target["bindings"],
                    source,
                    "gene_disease_context",
                    "exact_allele" if source.variant_keys else "context",
                    source.source_id,
                )
                context["source_ids"] = (source.source_id,)
            target["bindings"] = tuple(target["bindings"])
        for source in review["sources"]:
            if source.locator == "input:patient/genomic-settings":
                _add_source(registry, source)
        # Close only required parents; unrelated prior reports are not assessment input.
        pending = [
            parent
            for source in registry.values()
            for parent in source.parent_source_ids
            if parent not in registry
        ]
        while pending:
            source_id = pending.pop()
            if source_id not in registry:
                source = review_sources[source_id]
                _add_source(registry, source)
                pending.extend(source.parent_source_ids)
        result["input"]["sources"] = tuple(registry.values())
        from mendelprime.tools.acmg import _phase, _ranking_classify, ranking_preliminary

        packet = result["input"]
        targets = {(t["primary"]["gene"], _key(t["primary"]["key"])): t for t in packet["targets"]}
        # Opposing, unphased candidate indels cannot inherit an isolated-allele PVS1 anchor.
        for identity, target in targets.items():
            change = len(identity[1][4]) - len(identity[1][3])
            conflicts = [
                other
                for other, partner in targets.items()
                if other != identity
                and other[0] == identity[0]
                and change
                and len(other[1][4]) - len(other[1][3]) == -change
                and _phase(target, partner) != "confirmed_in_trans"
            ]
            if conflicts:
                target["local_haplotype_unresolved"] = tuple(other[1] for other in conflicts)
                target["preserved_criteria"] = tuple(
                    c for c in target["preserved_criteria"] if c["code"] != "PVS1"
                )
        ledger = ranking_preliminary(result)
        independent = {
            (g, _key(v["identity"]["variant_key"])): {
                **v,
                "classification": _ranking_classify(tuple(c for c in v["criteria"] if c["code"] != "PM3")),
            }
            for g, row in ledger["by_gene"].items()
            for v in row["variants"]
        }
        for identity, target in targets.items():
            own = independent[identity]
            target["partner_facts"] = []
            for other, partner in targets.items():
                if (
                    other == identity
                    or other[0] != identity[0]
                    or _phase(target, partner) != "confirmed_in_trans"
                ):
                    continue
                value = independent[other]
                shared = [
                    d
                    for d in own["disease_contexts"]
                    if d["disease_id"]
                    and d["mode_of_inheritance"] == "AR"
                    and any(
                        p["disease_id"] == d["disease_id"] and p["mode_of_inheritance"] == "AR"
                        for p in value["disease_contexts"]
                    )
                ]
                own_ids = {
                    (d["disease_id"], d["mode_of_inheritance"])
                    for d in own["disease_contexts"]
                    if d["disease_id"]
                }
                other_ids = {
                    (d["disease_id"], d["mode_of_inheritance"])
                    for d in value["disease_contexts"]
                    if d["disease_id"]
                }
                disease = shared[0] if shared and len(own_ids) == 1 and own_ids == other_ids else None
                if (
                    value["classification"] not in {"pathogenic", "likely_pathogenic"}
                    or disease is None
                    or partner.get("local_haplotype_unresolved")
                ):
                    continue
                criteria = tuple(
                    c
                    for c in value["criteria"]
                    if c["code"] != "PM3" and c["status"] == "met" and c.get("counted", True)
                )
                parents = tuple(
                    dict.fromkeys(c["source_id"] for criterion in criteria for c in criterion["citations"])
                )
                parents += tuple(
                    b["source_id"]
                    for t in (target, partner)
                    for b in t["bindings"]
                    if b["role"] == "patient_observation"
                )
                fact = {
                    "target_key": variant_text(identity[1]),
                    "partner_key": variant_text(other[1]),
                    "phase": "confirmed_in_trans",
                    "disease": disease,
                    "independent_partner_pathogenicity": value["classification"],
                    "PM3_excluded": True,
                    "basis": [{k: c[k] for k in ("code", "strength", "citations")} for c in criteria],
                }
                source = Source(
                    provider="mendelprime",
                    kind="derived",
                    locator=f"method:post-assessment/independent-partner/{variant_text(identity[1])}/{variant_text(other[1])}",
                    version="PM3-free-deterministic-partner-v1",
                    content=json.dumps(fact, sort_keys=True),
                    genes=(identity[0],),
                    variant_keys=(identity[1], other[1]),
                    case_key=case_tuple(seed["profile"]["case_key"]),
                    parent_source_ids=tuple(dict.fromkeys(parents)),
                )
                bindings = list(target["bindings"])
                _bind(
                    registry,
                    bindings,
                    source,
                    "partner_evidence",
                    "complex_allele",
                    "independent_partner_and_phase",
                )
                target["bindings"] = tuple(bindings)
                target["partner_facts"].append(
                    {
                        "partner_key": other[1],
                        "classification": value["classification"],
                        "source_id": source.source_id,
                        "disease": disease,
                    }
                )
            target["partner_facts"] = tuple(target["partner_facts"])
        packet["sources"] = tuple(registry.values())
        from mendelprime.tools.acmg import validate_input

        validate_input(packet)
    return result


def _convert_record(value, source_value, preserve_tuples):
    if isinstance(value, Source):
        return source_value(value)
    if isinstance(value, Mapping):
        return {name: _convert_record(item, source_value, preserve_tuples) for name, item in value.items()}
    if isinstance(value, (tuple, list)):
        items = [_convert_record(item, source_value, preserve_tuples) for item in value]
        return tuple(items) if preserve_tuples and isinstance(value, tuple) else items
    return value


def _token(value):
    return _convert_record(value, lambda source: source.source_id, False)


def _digest(value):
    return hashlib.sha256(
        json.dumps(_token(value), sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _copy_record(value):
    return _convert_record(value, lambda source: source, True)


def _merge_clinical(gene, values):
    result = {"gene": gene}
    for name in ("variant_evidence", "gene_context", "citation_seeds", "lookups"):
        rows = {}
        for value in values:
            for row in value.get(name, ()):
                rows.setdefault(_digest(row), row)
        result[name] = tuple(rows.values())
    return result


def _reusable(name, value):
    stable = {"found", "not_found", "not_applicable", "gene_scope_mismatch", "missing_prediction"}
    if name == "guidance":
        return type(value) is tuple and all(isinstance(item, Source) for item in value)
    if not isinstance(value, Mapping):
        return False
    if name == "resolve":
        return value.get("status") in {"sequence_verified", "reference_mismatch"}
    if name in {"annotate", "pvs1"}:
        return value.get("status") in {"annotated", "assessed", "not_applicable"} or (
            value.get("status") in {"unknown", "ambiguous"}
            and value.get("error_code")
            in {
                "vep_exact_transcript_not_found",
                "vep_exact_transcript_ambiguous",
                "vep_canonical_transcript_not_found",
                "vep_canonical_transcript_ambiguous",
                "pvs1_exact_transcript_unavailable",
                "pvs1_hgvsc_unavailable",
                "pvs1_strength_unset",
                "pvs1_protein_termination_unavailable",
            }
        )
    if name == "population":
        return value.get("status") in {
            "found",
            "found_filtered",
            "frequency_unavailable",
            "not_found_covered",
            "coverage_insufficient",
        }
    if name == "predictions":
        return all(value.get(key, {}).get("status") in stable for key in ("revel", "spliceai"))
    return "lookups" in value and all(row.get("status") in stable for row in value["lookups"])


def _reuse_readers(seed, readers, previous):
    case_identity = _digest({name: seed[name] for name in ("profile", "family", "observations")})
    resource = readers.get("resource_identity")
    resource_identity = _digest(
        resource()
        if callable(resource)
        else (
            resource
            if resource is not None
            else tuple((name, id(action)) for name, action in readers.items())
        )
    )
    state = {
        "version": "ranking_preparation_v1",
        "case_identity": case_identity,
        "resource_identity": resource_identity,
        "results": {},
    }
    if previous is not None:
        if not isinstance(previous, Mapping) or not isinstance(previous.get("reuse"), Mapping):
            raise MendelPrimeError("preparation_reuse_required", field="assessment.previous")
        prior = previous["reuse"]
        if prior.get("version") != state["version"] or not isinstance(prior.get("results"), Mapping):
            raise MendelPrimeError("preparation_reuse_required", field="assessment.previous")
        if prior.get("case_identity") != case_identity:
            raise MendelPrimeError("preparation_case_identity_mismatch", field="assessment.previous")
        if prior.get("resource_identity") == resource_identity:
            state["results"] = _copy_record(prior["results"])
    coverage = []

    def record(name, args, kwargs, value, reused):
        request_id = _digest((args, kwargs))
        cacheable = _reusable(name, value)
        if cacheable:
            state["results"][(name, request_id)] = _copy_record(value)
        request = args[0] if name in {"annotate", "pvs1"} else {}
        gene = request.get("gene") if isinstance(request, Mapping) else None
        if name in {"clinical", "gene_context", "protein_context", "population", "predictions"}:
            gene = args[0]
        key = request.get("variant_key", request.get("key")) if isinstance(request, Mapping) else None
        if name == "resolve":
            key = args[0]
        elif name in {"population", "predictions"}:
            key = args[1]
        elif name == "protein_context":
            key = kwargs.get("variant_key")
        outcomes = ()
        if isinstance(value, Mapping):
            outcomes = value.get("lookups", ())
            if name == "predictions":
                outcomes = tuple(
                    {"provider": provider, **value[provider]} for provider in ("revel", "spliceai")
                )
        status = "complete" if cacheable else "incomplete"
        error = ""
        if isinstance(value, Mapping):
            status = value.get("status", status)
            error = value.get("error_code", "")
        coverage.append(
            {
                "resource": name,
                "request_id": request_id,
                "gene": gene,
                "variant_key": key,
                "transcript_id": request.get("transcript_id", request.get("transcript"))
                if isinstance(request, Mapping)
                else None,
                "variant_keys": args[1] if name == "clinical" else (),
                "status": status,
                "error_code": error,
                "outcomes": tuple(
                    {
                        field: outcome.get(field)
                        for field in ("provider", "variant_key", "status", "error_code")
                    }
                    for outcome in outcomes
                ),
                "source_ids": tuple(dict.fromkeys(source.source_id for source in iter_sources(value))),
                "reused": reused,
                "cacheable": cacheable,
            }
        )
        return _copy_record(value)

    def cached(name, args, kwargs):
        value = state["results"].get((name, _digest((args, kwargs))))
        return record(name, args, kwargs, value, True) if _reusable(name, value) else None

    def invoke(name, *args, **kwargs):
        value = cached(name, args, kwargs)
        return (
            value
            if value is not None
            else record(name, args, kwargs, _reader(readers, name, *args, **kwargs), False)
        )

    def batch(name, values):
        if not values:
            return ()
        outcomes, missing, indices = [None] * len(values), [], []
        for index, value in enumerate(values):
            result = cached(name, (value,), {})
            if result is None:
                missing.append(value)
                indices.append(index)
            else:
                outcomes[index] = result
        if missing:
            result = _validated_batch(name, tuple(missing), _reader(readers, name, tuple(missing)), Mapping)
            for index, requested, value in zip(indices, missing, result, strict=True):
                outcomes[index] = record(name, (requested,), {}, value, False)
        return tuple(outcomes)

    def clinical(gene, keys):
        requested = keys or (None,)
        values, missing = [None] * len(requested), []
        for index, allele in enumerate(requested):
            arguments = (gene, (allele,) if allele is not None else ())
            value = cached("clinical", arguments, {})
            if value is None:
                missing.append((index, allele))
            else:
                values[index] = value
        if missing:
            result = _reader(readers, "clinical", gene, tuple(key for _, key in missing if key is not None))
            if not isinstance(result, Mapping) or result.get("gene") != gene:
                raise MendelPrimeError("gene_result_mismatch", field="assessment.clinical")
            for index, allele in missing:
                projected = {"gene": gene, "gene_context": result.get("gene_context", ())}
                for name in ("variant_evidence", "citation_seeds", "lookups"):
                    projected[name] = tuple(
                        row for row in result.get(name, ()) if row.get("variant_key") in {None, allele}
                    )
                arguments = (gene, (allele,) if allele is not None else ())
                values[index] = record("clinical", arguments, {}, projected, False)
        return _merge_clinical(gene, values)

    actions = {}
    for name in ("resolve", "population", "predictions", "gene_context", "protein_context"):
        actions[name] = lambda *args, _name=name, **kwargs: invoke(_name, *args, **kwargs)
    for name in ("annotate", "pvs1"):
        actions[name] = lambda values, _name=name: batch(_name, values)
    guidance = readers.get("guidance")
    actions["clinical"] = clinical
    actions["guidance"] = (lambda genes: invoke("guidance", genes)) if callable(guidance) else guidance
    return actions, state, coverage


def prepare(seed, readers, *, genes, previous=None):
    """Prepare every retained candidate allele, reusing only case/resource-bound reader facts."""
    if not isinstance(seed, Mapping) or not _SEED_FIELDS <= set(seed):
        raise MendelPrimeError("schema_mismatch", field="assessment.seed")
    if type(genes) is not tuple or len(genes) != len(set(genes)):
        raise MendelPrimeError("unique_gene_tuple_required", field="assessment.genes")
    for gene in genes:
        validate_gene(gene, field="assessment.genes")
    if not isinstance(readers, Mapping):
        raise MendelPrimeError("mapping_required", field="assessment.readers")
    if type(seed["sources"]) is not tuple or any(not isinstance(item, Source) for item in seed["sources"]):
        raise MendelPrimeError("source_tuple_required", field="assessment.seed.sources")
    case_key = case_tuple(seed["profile"]["case_key"])
    if case_tuple(seed["observations"]["case_key"]) != case_key:
        raise MendelPrimeError("case_identity_mismatch", field="assessment.observations")
    if any(source.case_key is not None and source.case_key != case_key for source in seed["sources"]):
        raise MendelPrimeError("case_identity_mismatch", field="assessment.sources")
    actions, state, coverage = _reuse_readers(seed, readers, previous)
    result = _assemble({**seed, "fixed_genes": genes}, actions, report=False)
    result.update(coverage=tuple(coverage), reuse=state)
    return result


__all__ = ["build", "prepare", "validate_seed"]
