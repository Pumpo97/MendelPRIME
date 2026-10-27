"""Source-bound GRCh37 reference verification and HGVS normalization."""

from collections import OrderedDict
from functools import cached_property
import hashlib
from importlib.metadata import version as package_version
from pathlib import Path
from threading import RLock

import hgvs.edit
import hgvs.location
import hgvs.normalizer
import hgvs.posedit
import hgvs.sequencevariant
import pysam

from mendelprime.domain import MendelPrimeError, variant_key

HGVS_NOMENCLATURE_VERSION = "21.1.4"
HGVS_PACKAGE_VERSION = package_version("hgvs")
NORMALIZATION_METHOD = "hgvs_genomic_normalizer_with_local_fasta"
_CACHE_SIZE = 8192
_VERSIONS_37 = (10, 11, 11, 11, 9, 11, 13, 10, 11, 10, 9, 11, 10, 8, 9, 9, 10, 9, 9, 10, 8, 10)
_ACCESSIONS_37 = {str(index): f"NC_{index:06d}.{release}" for index, release in enumerate(_VERSIONS_37, 1)}
_ACCESSIONS_37.update(X="NC_000023.10", Y="NC_000024.9")
_BY_ACCESSION_37 = {accession: chromosome for chromosome, accession in _ACCESSIONS_37.items()}
# Whole-sequence identities, not contig names or length-based accession guesses.
_MT_REFERENCES = {
    "f156ff3f65bbcc80c7ebb9936dceb96b1477b4f8f535c4e1dbe7baea225cbc66": "NC_012920.1",
    "8264337d37c28c396e70204867ede43ba223090ff5d27b46625c1752d6d8bd96": "NC_001807.4",
}


class ReferenceGenome:
    """Own one local FASTA and return plain, audited allele records."""

    def __init__(self, fasta: Path, *, assembly: str):
        if not isinstance(fasta, Path) or not fasta.is_absolute():
            raise MendelPrimeError("absolute_path_required", field="reference_fasta")
        if assembly not in {"GRCh37", "GRCh38"}:
            raise MendelPrimeError("assembly_enum_required", field="reference_assembly")
        self.assembly = assembly
        self._lock = RLock()
        self._closed = False
        self._fasta = None
        self._references = frozenset()
        self._cache = OrderedDict()
        try:
            if not fasta.is_file() or not Path(f"{fasta}.fai").is_file():
                raise OSError
            self._fasta = pysam.FastaFile(str(fasta))
            self._references = frozenset(self._fasta.references)
        except (OSError, ValueError):
            if self._fasta is not None:
                self._fasta.close()
            self._fasta = None
        self._right = hgvs.normalizer.Normalizer(self, shuffle_direction=3, validate=False)
        self._left = hgvs.normalizer.Normalizer(self, shuffle_direction=5, validate=False)

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    @property
    def closed(self):
        with self._lock:
            return self._closed

    def get_seq(self, accession, start=None, end=None):
        """Provide synchronized sequence access for ``hgvs.normalizer``."""
        with self._lock:
            if self._closed:
                raise MendelPrimeError("reference_genome_closed")
            if self._fasta is None:
                raise MendelPrimeError("reference_fasta_unavailable")
            contig = self._contig_for(accession)
            if contig is None:
                raise MendelPrimeError("reference_contig_unavailable")
            try:
                length = self._fasta.get_reference_length(contig)
                lower = max(0, start or 0)
                upper = min(length, length if end is None else end)
                return self._fasta.fetch(contig, lower, max(lower, upper)).upper()
            except (KeyError, OSError, ValueError) as exc:
                raise MendelPrimeError("reference_read_failed") from exc

    def resolve(self, key):
        key = _key(key)
        with self._lock:
            if self._closed:
                return _unresolved(key, "reference_genome_closed")
            cached = self._cache.get(key)
            if cached is not None:
                self._cache.move_to_end(key)
                return cached
            result = self._resolve(key)
            if result["status"] != "unresolved":
                self._cache[key] = result
                if len(self._cache) > _CACHE_SIZE:
                    self._cache.popitem(last=False)
            return result

    def equivalent(self, first, second):
        first, second = _key(first), _key(second)
        if first == second:
            return True
        left, right = self.resolve(first), self.resolve(second)
        return (
            left["status"] == right["status"] == "sequence_verified"
            and left["normalized_hgvs"] == right["normalized_hgvs"]
        )

    def close(self):
        with self._lock:
            if self._closed:
                return
            self._closed = True
            self._cache.clear()
            if self._fasta is not None:
                self._fasta.close()
                self._fasta = None
            self._references = frozenset()

    def _resolve(self, key):
        assembly, chromosome, position, reference, alternate = key
        if self.assembly != "GRCh37" or assembly != "GRCh37":
            return _unresolved(key, "unsupported_assembly")
        if "N" in reference or "N" in alternate:
            return _unresolved(key, "ambiguous_sequence")
        if self._fasta is None:
            return _unresolved(key, "reference_fasta_unavailable")
        accession = self._mitochondrial[0] if chromosome == "MT" else _ACCESSIONS_37.get(chromosome)
        if accession is None:
            return _unresolved(key, "reference_accession_unconfirmed")
        contig = self._contig_for(accession)
        if contig is None:
            return _unresolved(key, "reference_contig_unavailable", accession=accession)
        start, stop = position - 1, position - 1 + len(reference)
        try:
            length = self._fasta.get_reference_length(contig)
        except (KeyError, OSError, ValueError):
            return _unresolved(key, "reference_read_failed", accession=accession)
        if stop > length:
            return _unresolved(key, "reference_interval_out_of_bounds", accession=accession)
        try:
            observed = self.get_seq(accession, start, stop)
        except MendelPrimeError as exc:
            return _unresolved(key, exc.code, accession=accession)
        if observed != reference:
            return _resolved_record(
                key,
                status="reference_mismatch",
                region=_raw_region(key),
                accession=accession,
                observed=observed,
                issues=("reference_mismatch",),
            )
        if accession == "NC_001807.4":
            return _unresolved(key, "mitochondrial_reference_unsupported", accession=accession, observed=observed)
        change = hgvs.sequencevariant.SequenceVariant(
            ac=accession,
            type="m" if chromosome == "MT" else "g",
            posedit=hgvs.posedit.PosEdit(
                hgvs.location.Interval(
                    hgvs.location.SimplePosition(position),
                    hgvs.location.SimplePosition(position + len(reference) - 1),
                ),
                hgvs.edit.NARefAlt(ref=reference, alt=alternate),
            ),
        )
        try:
            right = self._right.normalize(change)
            left = self._left.normalize(change)
        except MendelPrimeError as exc:
            return _unresolved(key, exc.code, accession=accession, observed=observed)
        except Exception:
            return _unresolved(key, "hgvs_normalization_failed", accession=accession, observed=observed)
        return _resolved_record(
            key,
            status="sequence_verified",
            region=(
                min(position, left.posedit.pos.start.base - 1),
                max(position + len(reference), right.posedit.pos.end.base + 1),
            ),
            accession=accession,
            observed=observed,
            normalized_hgvs=str(right),
            normalization_method=NORMALIZATION_METHOD,
        )

    @cached_property
    def _mitochondrial(self):
        matches = {
            contig: _MT_REFERENCES.get(hashlib.sha256(self._fasta.fetch(contig).upper().encode()).hexdigest())
            for contig in ("MT", "chrMT", "chrM", "M", "NC_012920.1", "NC_001807.4")
            if contig in self._references
        }
        accessions = set(matches.values())
        if not matches or None in accessions or len(accessions) != 1:
            return None, None
        return next(iter(accessions)), next(iter(matches))

    def _contig_for(self, accession):
        if accession in _MT_REFERENCES.values():
            return self._mitochondrial[1] if self._mitochondrial[0] == accession else None
        chromosome = _BY_ACCESSION_37.get(accession, accession)
        return next((value for value in (chromosome, f"chr{chromosome}") if value in self._references), None)


def _key(value):
    if type(value) is not tuple or len(value) != 5:
        raise MendelPrimeError("variant_key_required", field="variant")
    return variant_key(*value)


def _raw_region(key):
    return key[2], key[2] + len(key[3]) - 1


def _resolved_record(
    key,
    *,
    status,
    region,
    accession=None,
    observed=None,
    normalized_hgvs=None,
    normalization_method=None,
    issues=(),
):
    return {
        "raw_key": key,
        "status": status,
        "equivalent_region": region,
        "reference_accession": accession,
        "reference_observed": observed,
        "normalized_hgvs": normalized_hgvs,
        "normalization_method": normalization_method,
        "issues": tuple(issues),
        "hgvs_package_version": HGVS_PACKAGE_VERSION,
        "hgvs_nomenclature_version": HGVS_NOMENCLATURE_VERSION,
    }


def _unresolved(key, issue, *, accession=None, observed=None):
    return _resolved_record(
        key,
        status="unresolved",
        region=_raw_region(key),
        accession=accession,
        observed=observed,
        issues=(issue,),
    )


__all__ = ["ReferenceGenome"]
