from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator, TextIO


SCHEMA = """
CREATE TABLE metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE variants (
    variation_id TEXT NOT NULL,
    allele_id TEXT NOT NULL,
    assembly TEXT NOT NULL,
    chromosome TEXT NOT NULL,
    position INTEGER NOT NULL,
    reference TEXT NOT NULL,
    alternate TEXT NOT NULL,
    gene_symbol TEXT NOT NULL,
    rs_id TEXT NOT NULL,
    clinical_significance TEXT NOT NULL,
    review_status TEXT NOT NULL,
    name TEXT NOT NULL,
    PRIMARY KEY (
        variation_id, allele_id, assembly, chromosome, position, reference, alternate
    )
);
CREATE TABLE identifiers (
    identifier TEXT NOT NULL,
    variation_id TEXT NOT NULL,
    allele_id TEXT NOT NULL,
    assembly TEXT NOT NULL,
    PRIMARY KEY (identifier, variation_id, allele_id, assembly)
);
CREATE TABLE citations (
    variation_id TEXT NOT NULL,
    allele_id TEXT NOT NULL,
    pmid TEXT NOT NULL,
    PRIMARY KEY (variation_id, allele_id, pmid)
);
"""

INDEXES = """
CREATE INDEX variants_coordinate_idx
    ON variants (assembly, chromosome, position, reference, alternate);
CREATE INDEX variants_variation_id_idx ON variants (variation_id);
CREATE INDEX variants_rsid_idx ON variants (rs_id);
CREATE INDEX identifiers_value_idx ON identifiers (identifier);
CREATE INDEX citations_variant_idx ON citations (variation_id, allele_id);
"""

_HGVS = re.compile(r"(?<![A-Za-z0-9_])([cgmnpr]\.[^\s,;()]+)", re.IGNORECASE)


def build_index(
    variant_summary: Path,
    citations: Path,
    output: Path,
    *,
    release: str,
) -> dict[str, int]:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    working = output.with_name(f"{output.name}.building")
    for path in (working, Path(f"{working}-wal"), Path(f"{working}-shm")):
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    connection = sqlite3.connect(working)
    try:
        connection.executescript(SCHEMA)
        connection.execute("PRAGMA journal_mode=DELETE")
        connection.execute("PRAGMA synchronous=OFF")
        variant_count, identifier_count = _load_variants(connection, variant_summary)
        citation_count = _load_citations(connection, citations)
        connection.executescript(INDEXES)
        metadata = {
            "schema_version": "mgagent_clinvar_index_v1",
            "clinvar_release": release,
            "built_at": datetime.now(timezone.utc).isoformat(),
            "variant_summary_path": variant_summary.name,
            "variant_summary_size": str(variant_summary.stat().st_size),
            "var_citations_path": citations.name,
            "var_citations_size": str(citations.stat().st_size),
            "variant_row_count": str(variant_count),
            "identifier_row_count": str(identifier_count),
            "citation_row_count": str(citation_count),
        }
        connection.executemany(
            "INSERT INTO metadata (key, value) VALUES (?, ?)", metadata.items()
        )
        connection.commit()
        connection.execute("PRAGMA optimize")
    finally:
        connection.close()
    os.replace(working, output)
    return {
        "variants": variant_count,
        "identifiers": identifier_count,
        "citations": citation_count,
    }


def _load_variants(connection: sqlite3.Connection, path: Path) -> tuple[int, int]:
    variants: list[tuple[object, ...]] = []
    identifiers: list[tuple[str, str, str, str]] = []
    variant_count = identifier_count = 0
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            row = {_key(key): value for key, value in raw.items()}
            assembly = _assembly(row.get("Assembly", ""))
            position = _integer(row.get("PositionVCF", ""))
            reference = _allele(row.get("ReferenceAlleleVCF", ""))
            alternate = _allele(row.get("AlternateAlleleVCF", ""))
            chromosome = str(row.get("Chromosome") or "").removeprefix("chr")
            variation_id = str(row.get("VariationID") or "").strip()
            allele_id = str(row.get("AlleleID") or "").strip()
            if not (
                assembly
                and chromosome
                and position > 0
                and reference
                and alternate
                and variation_id.isdigit()
                and allele_id.isdigit()
            ):
                continue
            name = str(row.get("Name") or "").strip()
            variants.append(
                (
                    variation_id,
                    allele_id,
                    assembly,
                    chromosome,
                    position,
                    reference,
                    alternate,
                    str(row.get("GeneSymbol") or "").strip().upper(),
                    _digits(str(row.get("RS# (dbSNP)") or "")),
                    str(row.get("ClinicalSignificance") or "").strip(),
                    str(row.get("ReviewStatus") or "").strip(),
                    name,
                )
            )
            for identifier in _hgvs_identifiers(name):
                identifiers.append((identifier, variation_id, allele_id, assembly))
            if len(variants) >= 10_000:
                variant_count += _insert_variants(connection, variants)
                identifier_count += _insert_identifiers(connection, identifiers)
                variants.clear()
                identifiers.clear()
                connection.commit()
    variant_count += _insert_variants(connection, variants)
    identifier_count += _insert_identifiers(connection, identifiers)
    connection.commit()
    return variant_count, identifier_count


def _load_citations(connection: sqlite3.Connection, path: Path) -> int:
    batch: list[tuple[str, str, str]] = []
    count = 0
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            row = {_key(key): value for key, value in raw.items()}
            if str(row.get("citation_source") or "").strip() != "PubMed":
                continue
            variation_id = str(row.get("VariationID") or "").strip()
            allele_id = str(row.get("AlleleID") or "").strip()
            pmid = str(row.get("citation_id") or "").strip()
            if variation_id.isdigit() and allele_id.isdigit() and pmid.isdigit():
                batch.append((variation_id, allele_id, pmid))
            if len(batch) >= 10_000:
                count += _insert_citations(connection, batch)
                batch.clear()
                connection.commit()
    count += _insert_citations(connection, batch)
    connection.commit()
    return count


def _insert_variants(connection: sqlite3.Connection, rows: list[tuple[object, ...]]) -> int:
    before = connection.total_changes
    connection.executemany(
        """
        INSERT OR IGNORE INTO variants VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )
    return connection.total_changes - before


def _insert_identifiers(
    connection: sqlite3.Connection, rows: list[tuple[str, str, str, str]]
) -> int:
    before = connection.total_changes
    connection.executemany(
        "INSERT OR IGNORE INTO identifiers VALUES (?, ?, ?, ?)", rows
    )
    return connection.total_changes - before


def _insert_citations(
    connection: sqlite3.Connection, rows: list[tuple[str, str, str]]
) -> int:
    before = connection.total_changes
    connection.executemany(
        "INSERT OR IGNORE INTO citations VALUES (?, ?, ?)", rows
    )
    return connection.total_changes - before


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def _key(value: str | None) -> str:
    return str(value or "").removeprefix("#")


def _assembly(value: str) -> str:
    normalized = value.strip().lower()
    if normalized in {"grch37", "hg19"}:
        return "GRCh37"
    if normalized in {"grch38", "hg38"}:
        return "GRCh38"
    return ""


def _integer(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def _allele(value: str) -> str:
    allele = value.strip().upper()
    return allele if allele and allele not in {"NA", "-"} else ""


def _digits(value: str) -> str:
    match = re.search(r"\d+", value)
    return match.group(0) if match else ""


def _hgvs_identifiers(name: str) -> Iterator[str]:
    seen: set[str] = set()
    for match in _HGVS.finditer(name):
        identifier = match.group(1).rstrip(".")
        if identifier not in seen:
            seen.add(identifier)
            yield identifier


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the MendelPRIME ClinVar SQLite index")
    parser.add_argument("--variant-summary", type=Path, required=True)
    parser.add_argument("--var-citations", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--release", required=True)
    args = parser.parse_args()
    counts = build_index(
        args.variant_summary,
        args.var_citations,
        args.output,
        release=args.release,
    )
    print(
        f"built {args.output}: variants={counts['variants']} "
        f"identifiers={counts['identifiers']} citations={counts['citations']}"
    )


if __name__ == "__main__":
    main()
