from __future__ import annotations

import argparse
import hashlib
import json
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    tables = _tables(args.workbook)
    output = {
        "schema_version": "ai_cura_knowledgebase_v1",
        "source": {
            "doi": "10.1126/scitranslmed.adz4172",
            "file": args.workbook.name,
            "sha256": hashlib.sha256(args.workbook.read_bytes()).hexdigest(),
        },
        "gene_disease_consistency": _records(tables["ST8_Gene-disease consistency"]),
        "pp1_recommendations": _records(tables["ST9_modified PP1"]),
        "pp4_knowledgebase": _records(tables["ST10_PP4 knowledgebase"]),
        "ps4_thresholds": {
            "supporting": 1,
            "moderate": 2,
            "strong": 4,
            "source": "AI-CURA main text sensitivity-prioritized universal threshold",
        },
        "ps2_pm6_strength_thresholds": {
            "supporting": 0.5,
            "moderate": 1.0,
            "strong": 2.0,
            "very_strong": 4.0,
            "source": "ClinGen PS2/PM6 recommendation v1.1 used by AI-CURA",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    print({key: len(value) for key, value in output.items() if isinstance(value, list)})


def _records(rows: list[list[str]]) -> list[dict[str, str]]:
    header_index = next(index for index, row in enumerate(rows) if sum(bool(value) for value in row) > 1)
    header = [_key(value) for value in rows[header_index]]
    return [
        {header[index]: value for index, value in enumerate(row) if index < len(header) and value}
        for row in rows[header_index + 1 :]
        if any(row)
    ]


def _key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def _tables(path: Path) -> dict[str, list[list[str]]]:
    with zipfile.ZipFile(path) as archive:
        shared = _shared_strings(archive)
        rels = {
            item.attrib["Id"]: item.attrib["Target"]
            for item in ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        }
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        output = {}
        for sheet in workbook.findall(f".//{{{MAIN}}}sheet"):
            name = sheet.attrib["name"]
            if name not in {
                "ST8_Gene-disease consistency",
                "ST9_modified PP1",
                "ST10_PP4 knowledgebase",
            }:
                continue
            relationship = sheet.attrib[f"{{{REL}}}id"]
            target = "xl/" + rels[relationship].lstrip("/").removeprefix("xl/")
            output[name] = _sheet_rows(archive.read(target), shared)
        return output


def _shared_strings(archive: zipfile.ZipFile) -> list[str]:
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return [
        "".join(item.text or "" for item in value.iter(f"{{{MAIN}}}t"))
        for value in root
    ]


def _sheet_rows(raw: bytes, shared: list[str]) -> list[list[str]]:
    root = ET.fromstring(raw)
    output = []
    for row in root.findall(f".//{{{MAIN}}}sheetData/{{{MAIN}}}row"):
        values: dict[int, str] = {}
        for cell in row.findall(f"{{{MAIN}}}c"):
            reference = cell.attrib.get("r", "A1")
            column = _column_index(reference)
            node = cell.find(f"{{{MAIN}}}v")
            value = "" if node is None else node.text or ""
            if cell.attrib.get("t") == "s" and value:
                value = shared[int(value)]
            values[column] = value.strip()
        output.append([values.get(index, "") for index in range(max(values, default=-1) + 1)])
    return output


def _column_index(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference).group(0)  # type: ignore[union-attr]
    value = 0
    for letter in letters:
        value = value * 26 + ord(letter) - 64
    return value - 1


if __name__ == "__main__":
    main()
