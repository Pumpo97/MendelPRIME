from __future__ import annotations

import argparse
import csv
import gzip
import re
import sqlite3
from pathlib import Path

_TRANSCRIPT = re.compile(r"^(N[MR]_\d+(?:\.\d+)?)")
_CODING = re.compile(r":(c\.[^ ]+)")
_PROTEIN = re.compile(r"\((p\.[A-Za-z*?=]+\d+[A-Za-z*?=]+)\)")
_CODON = re.compile(r"p\.[A-Za-z*?=]+(\d+)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--variant-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)

    source = sqlite3.connect(args.source.resolve().as_uri() + "?mode=ro", uri=True)
    output = sqlite3.connect(args.output)
    output.executescript(
        """
        PRAGMA journal_mode=WAL;
        CREATE TABLE protein_variants (
            gene TEXT NOT NULL,
            transcript TEXT NOT NULL,
            coding_hgvs TEXT NOT NULL,
            protein_hgvs TEXT NOT NULL,
            codon INTEGER NOT NULL,
            variation_id TEXT NOT NULL,
            clinical_significance TEXT NOT NULL,
            review_status TEXT NOT NULL,
            phenotype_ids TEXT NOT NULL,
            phenotype_list TEXT NOT NULL,
            PRIMARY KEY (gene, transcript, coding_hgvs, protein_hgvs, variation_id)
        );
        CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        """
    )
    release = source.execute(
        "SELECT value FROM metadata WHERE key='clinvar_release'"
    ).fetchone()
    output.execute(
        "INSERT INTO metadata VALUES (?, ?)",
        ("clinvar_release", str(release[0]) if release else ""),
    )
    rows = _source_rows(source, args.variant_summary)
    batch = []
    for gene, name, variation_id, significance, review, phenotype_ids, phenotype_list in rows:
        transcript = _match(_TRANSCRIPT, name)
        coding = _match(_CODING, name)
        protein = _match(_PROTEIN, name)
        codon_text = _match(_CODON, protein)
        if not gene or not transcript or not coding or not protein or not codon_text:
            continue
        batch.append(
            (
                str(gene).upper(), transcript, coding, protein, int(codon_text),
                str(variation_id), str(significance), str(review or ""),
                str(phenotype_ids or ""), str(phenotype_list or ""),
            )
        )
        if len(batch) >= 10_000:
            output.executemany("INSERT OR IGNORE INTO protein_variants VALUES (?,?,?,?,?,?,?,?,?,?)", batch)
            output.commit()
            batch.clear()
    if batch:
        output.executemany("INSERT OR IGNORE INTO protein_variants VALUES (?,?,?,?,?,?,?,?,?,?)", batch)
    output.executescript(
        """
        CREATE INDEX protein_exact_idx ON protein_variants(gene, protein_hgvs);
        CREATE INDEX protein_codon_idx ON protein_variants(gene, codon);
        """
    )
    output.commit()
    print(output.execute("SELECT COUNT(*) FROM protein_variants").fetchone()[0])
    output.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    output.close()
    source.close()


def _match(pattern: re.Pattern[str], value: object) -> str:
    match = pattern.search(str(value or ""))
    return match.group(1) if match else ""


def _source_rows(source: sqlite3.Connection, summary: Path | None):
    if summary is None:
        return (
            (*row, "", "")
            for row in source.execute(
                "SELECT DISTINCT gene_symbol,name,variation_id,clinical_significance,review_status "
                "FROM variants WHERE clinical_significance IN "
                "('Pathogenic','Likely pathogenic','Pathogenic/Likely pathogenic')"
            )
        )

    def rows():
        with gzip.open(summary, "rt", encoding="utf-8", errors="ignore", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("Assembly") != "GRCh37" or row.get("ClinicalSignificance") not in {
                    "Pathogenic",
                    "Likely pathogenic",
                    "Pathogenic/Likely pathogenic",
                }:
                    continue
                yield (
                    row.get("GeneSymbol"), row.get("Name"), row.get("VariationID"),
                    row.get("ClinicalSignificance"), row.get("ReviewStatus"),
                    row.get("PhenotypeIDS"), row.get("PhenotypeList"),
                )
    return rows()


if __name__ == "__main__":
    main()
