"""The five shared values; scientific procedures consume explicit records."""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping
from dataclasses import dataclass, field, fields, replace
from pathlib import Path, PurePosixPath, PureWindowsPath
from types import MappingProxyType
from urllib.parse import parse_qsl, unquote, urlsplit, urlunsplit


class MendelPrimeError(ValueError):
    """Stable error metadata without the original private/provider payload."""

    def __init__(
        self,
        code,
        *,
        category="contract",
        field="",
        attempts=0,
        http_status=None,
        retryable=False,
        source=None,
        record_index=None,
    ):
        self.code, self.category, self.field = code, category, field
        self.attempts, self.http_status, self.retryable = attempts, http_status, retryable
        self.source, self.record_index = source, record_index
        super().__init__(f"{field}: {code}" if field else code)


def require(condition, field, code):
    if not condition:
        raise MendelPrimeError(code, field=field)


_EVIDENCE_LIMITATIONS = frozenset(
    {
        "coverage_insufficient",
        "found_filtered",
        "population_frequency_unavailable",
        "not_found",
        "not_found_covered",
        "not_applicable",
        "article_or_abstract_missing",
        "source_citation_identifier_unresolvable",
        "missing_prediction",
        "contig_unavailable",
        "population_sites_contig_unavailable",
        "population_coverage_contig_unavailable",
        "revel_contig_unavailable",
        "spliceai_snv_contig_unavailable",
        "spliceai_indel_contig_unavailable",
        "mitochondrial_reference_unsupported",
        "mitochondrial_reference_unverified",
        "observation_not_carried",
        "observation_no_call",
        "observation_partial_call",
        "no_keyed_native_variants",
        "deferred_by_review_budget",
        "pvs1_native_source_not_applicable",
        "pvs1_consequence_not_applicable",
        "pvs1_exact_transcript_unavailable",
        "pvs1_hgvsc_unavailable",
        "pvs1_protein_termination_unavailable",
        "pvs1_strength_unset",
        "vep_exact_transcript_not_found",
        "vep_canonical_transcript_not_found",
        "post_refseq_missing_or_ambiguous",
        "post_vep_exact_refseq_unavailable",
        "post_vep_coding_effect_unavailable",
        "relationship_verification_unavailable",
    }
)


def is_evidence_limitation(issue):
    """Classify known scientific uncertainty; unknown failures stay failures."""
    return str(issue).rsplit(":", 1)[-1] in _EVIDENCE_LIMITATIONS


def validate_identifier(value, field):
    require(
        isinstance(value, str) and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,119}", value),
        field,
        "safe_identifier_required",
    )
    require(
        not value.endswith(".")
        and not re.fullmatch(r"(?:CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])", value.split(".")[0].upper()),
        field,
        "portable_identifier_required",
    )


def validate_gene(gene, field="gene"):
    require(
        isinstance(gene, str) and re.fullmatch(r"[A-Z0-9][A-Z0-9._-]*", gene),
        field,
        "canonical_gene_required",
    )


def validate_case_key(key):
    require(type(key) is tuple and len(key) == 2, "case_key", "case_key_required")
    for value, name in zip(key, ("cohort", "case_id"), strict=True):
        validate_identifier(value, name)


def is_sample_id(value):
    return (
        isinstance(value, str)
        and bool(value)
        and not any(char.isspace() or ord(char) < 32 or ord(char) == 127 for char in value)
    )


def variant_key(assembly, chromosome, position, reference, alternate):
    require(assembly in {"GRCh37", "GRCh38"}, "variant.assembly", "assembly_enum_required")
    require(
        isinstance(chromosome, str)
        and re.fullmatch(
            r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT|(?:GL|KI)\d+\.\d+|(?:NC|NT|NW)_\d+\.\d+)", chromosome
        ),
        "variant.chromosome",
        "canonical_contig_required",
    )
    require(type(position) is int and position > 0, "variant.position", "positive_integer_required")
    for name, sequence in (("reference", reference), ("alternate", alternate)):
        require(
            isinstance(sequence, str) and re.fullmatch(r"[ACGTN]+", sequence),
            f"variant.{name}",
            "sequence_allele_required",
        )
    require(reference != alternate, "variant.alternate", "alternate_must_differ")
    return assembly, chromosome, position, reference, alternate


def variant_text(key):
    return ":".join(map(str, key))


def variant_tuple(value, invalid):
    """Adapt assessment representations; allele validation remains at its boundary."""
    if type(value) is tuple:
        return value
    if type(value) is list:
        return tuple(value)
    if isinstance(value, Mapping):
        return tuple(value[name] for name in ("assembly", "chromosome", "position", "reference", "alternate"))
    raise invalid


def case_tuple(value):
    if type(value) is tuple and len(value) == 2:
        return value
    if isinstance(value, Mapping) and set(value) == {"cohort", "case_id"}:
        return value["cohort"], value["case_id"]
    raise MendelPrimeError("case_key_required", field="assessment.profile.case_key")


def variant_record(key):
    return dict(zip(("assembly", "chromosome", "position", "reference", "alternate"), key, strict=True))


def case_record(key):
    return dict(zip(("cohort", "case_id"), key, strict=True))


def profile_term(term):
    """The shared scientific phenotype record in VCF messages and public results."""
    return {
        "hpo_id": term["hpo_id"],
        "label": term["label"],
        "status": "present" if term["present"] else "absent",
        "information_content": term["information_content"],
        "priority_score": term["priority_score"],
        "role": term["role"],
        "definition": term["definition"],
        "parent_terms": term["parent_ids"],
        "parent_term_labels": term["parent_labels"],
    }


def call_state(call, alternate):
    alleles = call.get("alleles", ())
    called = bool(alleles) and all(item is not None for item in alleles)
    count = alleles.count(alternate) if called else None
    state = (
        "no_call"
        if count is None
        else "homozygous_alternate"
        if count == len(alleles)
        else "heterozygous"
        if count
        else "reference"
    )
    return called, count > 0 if count is not None else None, state


def family_origin(calls, alternate):
    """Describe raw calls, without asserting parentage or quality-qualified de novo."""
    states = {call["role"]: call_state(call, alternate) for call in calls}
    called, carried, _ = states["proband"]
    mother, father = (states.get(role, (False, None, "no_call"))[1] for role in ("mother", "father"))
    if carried is not True:
        return "not_observed_in_proband" if called else "unknown"
    if mother is True and father is not True:
        return "maternal"
    if father is True and mother is not True:
        return "paternal"
    if mother is False and father is False:
        return "de_novo_candidate_unconfirmed"
    return "unresolved"


def source_diverse_rrf(rankings, limit, rrf_k, reserve, source_order=()):
    """Merge ordered gene sequences; each caller owns its fixed scientific policy."""
    for name, value in (("limit", limit), ("reserve_per_source", reserve)):
        if type(value) is not int or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer")
    if type(rrf_k) not in {int, float} or not math.isfinite(rrf_k) or rrf_k <= 0:
        raise ValueError("rrf_k must be a positive finite number")
    sources = tuple(source for source in source_order if rankings.get(source))
    sources += tuple(sorted(set(rankings) - set(sources)))
    normalized = {
        source: tuple(dict.fromkeys(gene.strip().upper() for gene in rankings[source] if gene.strip()))
        for source in sources
    }
    scores, best = {}, {}
    for source in sources:
        for rank, gene in enumerate(normalized[source], 1):
            scores[gene] = scores.get(gene, 0.0) + 1.0 / (rrf_k + rank)
            best[gene] = min(best.get(gene, rank), rank)
    ordered = sorted(scores, key=lambda gene: (-scores[gene], best[gene], gene))
    protected = {gene for source in sources for gene in normalized[source][:reserve]}
    selected = set(ordered[:limit]) | protected
    if len(selected) > limit:
        removable = [gene for gene in ordered if gene in selected and gene not in protected]
        selected.difference_update(removable[-(len(selected) - limit) :])
    return tuple(gene for gene in ordered if gene in selected)[:limit]


def freeze_record(value):
    """Freeze a record once; already-frozen mappings are shared between stages."""
    if type(value) is dict:
        return MappingProxyType({key: freeze_record(item) for key, item in value.items()})
    if type(value) is list:
        return tuple(freeze_record(item) for item in value)
    if type(value) is tuple:
        return tuple(freeze_record(item) for item in value)
    return value


_SECRET_QUERY_KEYS = frozenset(
    (
        "apikey",
        "apitoken",
        "accesskey",
        "accesstoken",
        "authtoken",
        "clientsecret",
        "idtoken",
        "refreshtoken",
        "secretkey",
        "token",
        "secret",
        "password",
        "passwd",
        "credential",
        "credentials",
        "authorization",
        "auth",
        "signature",
        "xamzcredential",
        "xamzsecuritytoken",
        "xamzsignature",
        "xgoogcredential",
        "xgoogsignature",
    )
)
SOURCE_KINDS = frozenset(
    (
        "patient_observation",
        "ontology",
        "guideline",
        "gene_disease",
        "literature",
        "tool_observation",
        "population",
        "prediction",
        "annotation",
        "similar_case",
        "derived",
    )
)


def _has_secret_query(value):
    names = [name for name, _ in parse_qsl(value, keep_blank_values=True)]
    names.extend(re.findall(r"(?:^|[&;])([^=&;:]+):", value))
    return any(re.sub(r"[^a-z0-9]", "", name.casefold()) in _SECRET_QUERY_KEYS for name in names)


def _fragment(value, logical=False):
    if not value:
        return ""
    require(
        not any(ord(c) < 32 or ord(c) == 127 for c in unquote(value)),
        "locator",
        "control_character_forbidden",
    )
    require(not _has_secret_query(unquote(value)), "locator", "credential_fragment_forbidden")
    if logical:
        require(
            "%" not in value
            and re.fullmatch(r"[A-Za-z0-9._~:/-]+", value)
            and not value.startswith(("/", "\\"))
            and not PureWindowsPath(value).is_absolute()
            and ".." not in PurePosixPath(value).parts,
            "locator",
            "logical_fragment_invalid",
        )
    return value


def source_locator(value):
    require(isinstance(value, str) and value.strip(), "locator", "nonempty_text_required")
    require(not any(ord(c) < 32 or ord(c) == 127 for c in value), "locator", "control_character_forbidden")
    try:
        parsed = urlsplit(value)
        scheme = parsed.scheme.lower()
        host, port = parsed.hostname, parsed.port
    except ValueError:
        raise MendelPrimeError("locator_parse_failed", field="locator") from None
    require(
        scheme in {"http", "https", "input", "resource", "method"},
        "locator",
        "public_locator_scheme_required",
    )
    if scheme in {"http", "https"}:
        require(
            host and parsed.username is None and parsed.password is None,
            "locator",
            "url_credentials_forbidden",
        )
        require(not any(c.isspace() for c in value), "locator", "raw_whitespace_forbidden")
        try:
            host = host.encode("idna").decode("ascii").lower().rstrip(".")
        except UnicodeError:
            raise MendelPrimeError("canonical_web_host_required", field="locator") from None
        require(
            host
            and ":" not in host
            and all(
                re.fullmatch(r"[a-z0-9](?:[a-z0-9-]{0,61}[a-z0-9])?", label) for label in host.split(".")
            ),
            "locator",
            "canonical_web_host_required",
        )
        require(not _has_secret_query(parsed.query), "locator", "credential_parameter_forbidden")
        authority = host if port in {None, 80 if scheme == "http" else 443} else f"{host}:{port}"
        return urlunsplit((scheme, authority, parsed.path, parsed.query, _fragment(parsed.fragment)))
    require(not parsed.netloc and not parsed.query, "locator", "logical_authority_or_query_forbidden")
    path = parsed.path
    require(
        path
        and not path.startswith(("/", "\\"))
        and not PurePosixPath(path).is_absolute()
        and not PureWindowsPath(path).is_absolute(),
        "locator",
        "absolute_local_path_forbidden",
    )
    require(
        "%" not in path and "\\" not in path and ".." not in PurePosixPath(path).parts,
        "locator",
        "public_logical_locator_required",
    )
    require(not any(c.isspace() for c in path), "locator", "raw_whitespace_forbidden")
    fragment = _fragment(parsed.fragment, logical=True)
    return f"{scheme}:{path}" + (f"#{fragment}" if fragment else "")


@dataclass(frozen=True, slots=True, kw_only=True)
class Patient:
    case_id: str
    phenotypes: tuple
    cohort: str = "cases"
    clinical_description: str = ""
    genomic: dict | None = None

    def __post_init__(self):
        validate_case_key(self.key)
        require(
            type(self.phenotypes) is tuple and self.phenotypes,
            "phenotypes",
            "nonempty_phenotype_tuple_required",
        )
        seen = set()
        for term in self.phenotypes:
            identifier = term.get("hpo_id")
            require(
                isinstance(identifier, str) and re.fullmatch(r"HP:\d{7}", identifier, re.ASCII),
                "hpo_id",
                "hpo_identifier_required",
            )
            require(identifier not in seen, "phenotypes", "duplicate_hpo_identifier")
            require(isinstance(term.get("label"), str), "label", "text_required")
            require(type(term.get("present")) is bool, "present", "boolean_required")
            seen.add(identifier)
        require(isinstance(self.clinical_description, str), "clinical_description", "text_required")
        if self.genomic is not None:
            for name in ("vcf_path", "pedigree_path"):
                path = self.genomic.get(name)
                require(
                    (path is None and name == "pedigree_path")
                    or (isinstance(path, Path) and path.is_absolute()),
                    name,
                    "absolute_path_required",
                )
            require(
                self.genomic.get("assembly") in {None, "GRCh37", "GRCh38"},
                "assembly",
                "assembly_enum_required",
            )
            require(
                self.genomic.get("family_mode") in {"auto", "trio", "proband_only"},
                "family_mode",
                "family_mode_enum_required",
            )
            require(self.genomic.get("sex") in {"MALE", "FEMALE", "UNKNOWN"}, "sex", "unsupported_sex")
            require(
                isinstance(self.genomic.get("assay_type"), str) and self.genomic["assay_type"].strip(),
                "assay_type",
                "nonempty_assay_description_required",
            )
            roles = self.genomic.get("sample_roles")
            if roles is not None:
                values = []
                for role in ("proband", "father", "mother"):
                    sample = roles.get(role)
                    if sample is None and role != "proband":
                        continue
                    require(
                        is_sample_id(sample),
                        f"sample_roles.{role}",
                        "sample_identifier_required",
                    )
                    values.append(sample)
                require(len(values) == len(set(values)), "sample_roles", "roles_must_be_distinct")
        object.__setattr__(self, "phenotypes", freeze_record(self.phenotypes))
        object.__setattr__(self, "genomic", freeze_record(self.genomic))

    @property
    def key(self):
        return self.cohort, self.case_id

    @property
    def input_kind(self):
        return (
            "hpo_vcf"
            if self.genomic is not None
            else ("hpo_clinical" if self.clinical_description.strip() else "hpo")
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class Source:
    provider: str
    kind: str
    locator: str
    version: str
    content: str
    genes: tuple = ()
    variant_keys: tuple = ()
    case_key: tuple | None = None
    parent_source_ids: tuple = ()
    source_id: str = field(init=False)
    content_sha256: str = field(init=False)

    def __post_init__(self):
        require(isinstance(self.provider, str), "provider", "text_required")
        provider = self.provider.strip().lower()
        require(re.fullmatch(r"[a-z0-9][a-z0-9._-]*", provider), "provider", "canonical_slug_required")
        require(self.kind in SOURCE_KINDS, "kind", "source_kind_enum_required")
        locator = source_locator(self.locator)
        require(isinstance(self.version, str) and self.version.strip(), "version", "nonempty_text_required")
        require(self.version == self.version.strip(), "version", "surrounding_whitespace_forbidden")
        require(isinstance(self.content, str) and self.content.strip(), "content", "nonempty_text_required")
        for name in ("genes", "variant_keys", "parent_source_ids"):
            require(type(getattr(self, name)) is tuple, name, "tuple_required")
        for gene in self.genes:
            validate_gene(gene)
        genes = tuple(sorted(set(self.genes)))
        variants = tuple(sorted({variant_key(*key) for key in self.variant_keys}, key=variant_text))
        parents = tuple(sorted(set(self.parent_source_ids)))
        require(
            all(isinstance(x, str) and x and not any(c.isspace() for c in x) for x in parents),
            "parent_source_ids",
            "canonical_source_id_required",
        )
        if self.case_key is not None:
            validate_case_key(self.case_key)
        verified_literature_view = (
            provider == "pubmed"
            and self.kind == "literature"
            and bool(variants)
            and self.case_key is not None
            and len(parents) >= 2
            and "#verified-allele-" in locator
            and ";scope=verified-exact-paper-binding-v1;" in self.version
        )
        require(
            genes or self.kind in {"patient_observation", "ontology", "guideline"},
            "genes",
            "explicit_gene_scope_required",
        )
        require(not variants or len(genes) == 1, "genes", "allele_source_requires_one_gene")
        require(
            self.kind not in {"patient_observation", "tool_observation", "derived"}
            or self.case_key is not None,
            "case_key",
            "case_scope_required",
        )
        require(
            (self.kind == "derived" and bool(parents))
            or (self.kind != "derived" and (not parents or verified_literature_view)),
            "parent_source_ids",
            "derived_parents_required" if self.kind == "derived" else "parents_allowed_only_for_derived",
        )
        digest = hashlib.sha256(self.content.encode("utf-8")).hexdigest()
        identity = {
            "domain": "mendelprime.source.v1",
            "provider": provider,
            "kind": self.kind,
            "locator": locator,
            "version": self.version,
            "content_sha256": digest,
            "genes": list(genes),
            "variant_keys": [variant_text(v) for v in variants],
            "case_key": list(self.case_key) if self.case_key is not None else None,
            "parent_source_ids": list(parents),
        }
        source_id = (
            "SRC_"
            + hashlib.sha256(
                json.dumps(identity, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest()
        )
        for name, value in (
            ("provider", provider),
            ("locator", locator),
            ("genes", genes),
            ("variant_keys", variants),
            ("parent_source_ids", parents),
            ("content_sha256", digest),
            ("source_id", source_id),
        ):
            object.__setattr__(self, name, value)


def _binding_hgvs(value):
    return re.sub(r"\s+", "", value).casefold() if isinstance(value, str) else ""


def _binding_quote_matches(text, hgvs, *, nucleotide_context=""):
    hgvs = _binding_hgvs(hgvs)
    pattern = re.escape(hgvs).replace(">", r"\s*>\s*")
    exact = bool(
        hgvs
        and isinstance(text, str)
        and re.search(
            rf"(?<![a-z0-9_]){pattern}(?![a-z0-9_*+.-])",
            text,
            re.IGNORECASE,
        )
    )
    if exact or not isinstance(text, str):
        return exact
    snv = re.fullmatch(r"c\.(\d+)([acgt])>([acgt])", hgvs)
    if snv is None or not nucleotide_context:
        return False
    shorthand = rf"(?<![a-z0-9_.]){snv[2]}{snv[1]}{snv[3]}(?![a-z0-9_*+.-])"
    dna = r"\b(?:genomic DNA|cDNA|nucleotide)\b"
    context = re.search(
        rf"{dna}[^.!?]{{0,180}}{shorthand}|{shorthand}[^.!?]{{0,180}}{dna}", nucleotide_context, re.I
    )
    return bool(
        re.search(shorthand, text, re.I)
        and context
        and not re.search(r"\b(?:protein|amino acid|residue)\b", context[0], re.I)
    )


def _native_binding_identity(source, gene, key, case_key, transcript, hgvs):
    if not (
        source.provider in {"exomiser", "lirical"}
        and source.kind == "tool_observation"
        and source.genes == (gene,)
        and source.variant_keys == (key,)
        and source.case_key == case_key
    ):
        return False
    try:
        record = json.loads(source.content)
    except (json.JSONDecodeError, TypeError):
        if source.provider != "lirical":
            return False
        match = re.fullmatch(
            r"([^:\s]+):(\d+)([^>\s]+)>(\S+)\s+(.+?)\s+"
            r"pathogenicity:[-+0-9.eE]+\s+\[[^]]+\]",
            source.content,
            re.IGNORECASE,
        )
        if match is None:
            return False
        chromosome = "MT" if match[1].removeprefix("chr") == "M" else match[1].removeprefix("chr")
        observed = (key[0], chromosome, int(match[2]), match[3].upper(), match[4].upper())
        identity = re.search(
            rf"(?<![A-Za-z0-9_]){re.escape(transcript)}:{re.escape(hgvs)}"
            r"(?::p\.[^:\s]+)?(?![A-Za-z0-9_])",
            match[5],
        )
        return observed == key and identity is not None
    if not isinstance(record, dict):
        return False
    if record.get("variant_key") == variant_text(key):
        return record.get("transcript") == transcript and _binding_hgvs(
            record.get("hgvs_c")
        ) == _binding_hgvs(hgvs)
    raw = record.get("HGVS")
    if isinstance(raw, str):
        parts = raw.split(":")
        return (
            len(parts) >= 3
            and parts[0] == gene
            and parts[1] == transcript
            and _binding_hgvs(parts[2]) == _binding_hgvs(hgvs)
        )
    selected = record.get("selected_transcript")
    if isinstance(selected, dict):
        return selected.get("transcript_id") == transcript and _binding_hgvs(
            selected.get("hgvsc", "").removeprefix(f"{transcript}:")
        ) == _binding_hgvs(hgvs)
    return False


def verified_paper_binding_basis(clinvar, paper, native, gene, key, case_key, pmid, transcript, hgvs, quote):
    """Return the exact identity proof shared by the builder and registry validator."""
    require(
        isinstance(hgvs, str)
        and "?" not in hgvs
        and re.search(r"\d", hgvs)
        and any(change in hgvs.casefold() for change in (">", "del", "dup", "ins", "inv")),
        "verified_paper_binding",
        "unambiguous_allele_notation_required",
    )
    require(
        clinvar.provider == "clinvar"
        and clinvar.kind == "annotation"
        and clinvar.genes == (gene,)
        and clinvar.variant_keys == (key,),
        "verified_paper_binding",
        "clinvar_exact_coordinate_source_required",
    )
    require(
        paper.provider == "pubmed"
        and paper.kind == "literature"
        and paper.genes == (gene,)
        and not paper.variant_keys
        and not paper.parent_source_ids
        and re.fullmatch(rf"https://pubmed\.ncbi\.nlm\.nih\.gov/{re.escape(pmid)}/?", paper.locator),
        "verified_paper_binding",
        "original_pubmed_source_required",
    )
    require(
        isinstance(quote, str)
        and quote in paper.content
        and _binding_quote_matches(quote, hgvs, nucleotide_context=paper.content),
        "verified_paper_binding",
        "exact_allele_quote_required",
    )
    require(
        _native_binding_identity(native, gene, key, case_key, transcript, hgvs),
        "verified_paper_binding",
        "native_identity_mismatch",
    )
    try:
        content = json.loads(clinvar.content)
        record = content["record"]
        coordinates = (
            record["assembly"],
            str(record["chromosome"]).removeprefix("chr"),
            int(record["position"]),
            str(record["reference"]).upper(),
            str(record["alternate"]).upper(),
        )
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        raise MendelPrimeError("clinvar_exact_record_required", field="verified_paper_binding") from error
    require(
        content.get("record_type") == "clinvar_sql_normalized_coordinate_match"
        and coordinates == key
        and record.get("gene_symbol") == gene
        and record.get("allele_id") not in {None, ""}
        and record.get("variation_id") not in {None, ""}
        and _binding_quote_matches(record.get("name", ""), hgvs),
        "verified_paper_binding",
        "clinvar_exact_record_required",
    )
    require(
        any(
            isinstance(row, dict)
            and str(row.get("pmid")) == pmid
            and str(row.get("allele_id")) == str(record["allele_id"])
            and str(row.get("variation_id")) == str(record["variation_id"])
            for row in content.get("citation_rows", ())
        ),
        "verified_paper_binding",
        "clinvar_exact_citation_required",
    )
    return {
        "route": "clinvar_exact_coordinate_citation",
        "gene": gene,
        "variant_key": variant_text(key),
        "transcript": transcript,
        "hgvs_c": hgvs,
        "pmid": pmid,
        "clinvar_variation_id": str(record["variation_id"]),
        "clinvar_allele_id": str(record["allele_id"]),
        "paper_quote": quote,
        "identity_source_id": native.source_id,
    }


def validate_verified_paper_view(source, index):
    """Validate one scoped PubMed view against its receipt and immutable original parents."""
    require(isinstance(source, Source) and isinstance(index, Mapping), "source", "source_registry_required")
    marker = ";scope=verified-exact-paper-binding-v1;receipt="
    require(
        source.provider == "pubmed"
        and source.kind == "literature"
        and len(source.genes) == len(source.variant_keys) == 1
        and source.case_key is not None
        and marker in source.version,
        "source",
        "verified_paper_view_required",
    )
    paper_version, receipt_id = source.version.split(marker, 1)
    require(receipt_id and ";" not in receipt_id, "source.version", "binding_receipt_required")
    require(all(parent in index for parent in source.parent_source_ids), "source", "binding_parent_missing")
    receipt = index.get(receipt_id)
    require(
        isinstance(receipt, Source)
        and receipt.provider == "mendelprime"
        and receipt.kind == "derived"
        and receipt.version == "verified-exact-paper-binding-v1"
        and receipt.genes == source.genes
        and receipt.variant_keys == source.variant_keys
        and receipt.case_key == source.case_key,
        "source",
        "binding_receipt_invalid",
    )
    try:
        proof = json.loads(receipt.content)
        basis = proof["binding_basis"]
    except (KeyError, TypeError, json.JSONDecodeError) as error:
        raise MendelPrimeError("binding_receipt_invalid", field="source") from error
    require(
        set(proof) == {"record_type", "scope", "binding_basis", "does_not_establish"}
        and proof["record_type"] == "verified_exact_paper_binding_v1"
        and proof["scope"] == "allele_identity_only"
        and set(proof["does_not_establish"])
        == {"phase", "genotype_state", "independent_case_count", "PS4_strength", "pathogenic_classification"}
        and set(basis)
        == {
            "route",
            "gene",
            "variant_key",
            "transcript",
            "hgvs_c",
            "pmid",
            "clinvar_variation_id",
            "clinvar_allele_id",
            "validated_profile_section",
            "paper_quote",
            "identity_source_id",
        },
        "source",
        "binding_receipt_invalid",
    )
    parents = [index[parent] for parent in receipt.parent_source_ids]
    papers = [parent for parent in parents if parent.provider == "pubmed" and parent.kind == "literature"]
    clinvars = [parent for parent in parents if parent.provider == "clinvar" and parent.kind == "annotation"]
    native = index.get(basis["identity_source_id"])
    require(
        len(papers) == len(clinvars) == 1
        and isinstance(native, Source)
        and native.source_id in receipt.parent_source_ids,
        "source",
        "binding_parent_roles_invalid",
    )
    paper, clinvar = papers[0], clinvars[0]
    require(
        set(receipt.parent_source_ids) == {paper.source_id, clinvar.source_id, native.source_id},
        "source",
        "binding_receipt_parent_mismatch",
    )
    expected = verified_paper_binding_basis(
        clinvar,
        paper,
        native,
        source.genes[0],
        source.variant_keys[0],
        source.case_key,
        basis["pmid"],
        basis["transcript"],
        basis["hgvs_c"],
        basis["paper_quote"],
    )
    require(
        basis == {**expected, "validated_profile_section": basis["validated_profile_section"]}
        and isinstance(basis["validated_profile_section"], str)
        and basis["validated_profile_section"].startswith(
            ("allele_disease", "key_support:", "key_limits:", "validated_pvs1_anchor")
        ),
        "source",
        "binding_basis_mismatch",
    )
    expected_parents = {receipt.source_id, paper.source_id, clinvar.source_id, native.source_id}
    expected_locator = (
        paper.locator.rstrip("/") + f"/#verified-allele-{receipt.source_id.removeprefix('SRC_')[:16]}"
    )
    require(
        set(source.parent_source_ids) == expected_parents
        and source.content == paper.content
        and source.locator == expected_locator
        and paper.version == paper_version,
        "source",
        "verified_paper_view_provenance_mismatch",
    )


def iter_sources(value):
    """Walk nested evidence records without copying, reordering or deduplicating sources."""
    if isinstance(value, Source):
        yield value
    elif isinstance(value, Mapping):
        for item in value.values():
            yield from iter_sources(item)
    elif isinstance(value, (tuple, list)):
        for item in value:
            yield from iter_sources(item)


def source_record(source):
    return {
        "provider": source.provider,
        "kind": source.kind,
        "locator": source.locator,
        "version": source.version,
        "content": source.content,
        "genes": list(source.genes),
        "variant_keys": [variant_record(key) for key in source.variant_keys],
        "case_key": case_record(source.case_key) if source.case_key else None,
        "parent_source_ids": list(source.parent_source_ids),
        "source_id": source.source_id,
        "content_sha256": source.content_sha256,
    }


def source_from_record(record):
    require(
        type(record) is dict and set(record) == {item.name for item in fields(Source)},
        "source",
        "schema_mismatch",
    )
    for name in ("genes", "variant_keys", "parent_source_ids"):
        require(type(record[name]) is list, f"source.{name}", "array_required")
    require(
        record["case_key"] is None
        or (type(record["case_key"]) is dict and set(record["case_key"]) == {"cohort", "case_id"}),
        "source.case_key",
        "schema_mismatch",
    )
    values = {key: value for key, value in record.items() if key not in {"source_id", "content_sha256"}}
    values["genes"] = tuple(record["genes"])
    values["parent_source_ids"] = tuple(record["parent_source_ids"])
    values["variant_keys"] = tuple(variant_key(**value) for value in record["variant_keys"])
    values["case_key"] = (
        tuple(record["case_key"][key] for key in ("cohort", "case_id")) if record["case_key"] else None
    )
    source = Source(**values)
    require(
        record["content_sha256"] == source.content_sha256, "source.content_sha256", "content_hash_mismatch"
    )
    require(record["source_id"] == source.source_id, "source.source_id", "source_identity_mismatch")
    return source


@dataclass(frozen=True, slots=True, kw_only=True)
class Evidence:
    evidence_id: str
    kind: str
    case_key: tuple
    facts: dict
    genes: tuple = ()
    variant_keys: tuple = ()
    source_ids: tuple = ()
    derived_from: tuple = ()
    status: str = "ok"

    def __post_init__(self):
        validate_case_key(self.case_key)
        require(
            isinstance(self.evidence_id, str) and self.evidence_id, "evidence_id", "nonempty_text_required"
        )
        for gene in self.genes:
            validate_gene(gene)
        for key in self.variant_keys:
            variant_key(*key)
        object.__setattr__(self, "facts", freeze_record(self.facts))


@dataclass(frozen=True, slots=True, kw_only=True)
class Ranking:
    genes: tuple
    stage: str = ""
    status: str = "ok"
    details: dict = field(default_factory=dict)
    issues: tuple = ()

    def __post_init__(self):
        require(type(self.genes) is tuple, "genes", "gene_tuple_required")
        for gene in self.genes:
            validate_gene(gene)
        require(len(set(self.genes)) == len(self.genes), "genes", "duplicate_gene")
        object.__setattr__(self, "details", freeze_record(self.details))

    def position(self, gene):
        validate_gene(gene)
        return self.genes.index(gene) + 1 if gene in self.genes else None


@dataclass(frozen=True, slots=True, kw_only=True)
class CaseResult:
    patient: Patient
    stages: dict
    final: Ranking
    sources: tuple = ()
    evidence: tuple = ()
    explanation: dict = field(default_factory=lambda: {"status": "disabled", "markdown": ""})
    provenance: dict | None = None
    status: str = "ok"

    def __post_init__(self):
        for name in ("stages", "explanation", "provenance"):
            object.__setattr__(self, name, freeze_record(getattr(self, name)))

    def with_explanation(self, explanation, *, sources=None, provenance=None):
        return replace(
            self,
            explanation=explanation,
            sources=self.sources if sources is None else sources,
            provenance=self.provenance if provenance is None else provenance,
            status="degraded" if explanation["status"] == "degraded" else self.status,
        )
