#!/usr/bin/env python3
"""Prepare NeuroFate-only probe-to-gene mappings from GEO platform annotation tables."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
from pathlib import Path


PROBE_COLUMNS = ["ID", "ID_REF", "probe_id", "ProbeID", "Probe Set ID", "Reporter Identifier"]
GENE_COLUMNS = ["gene_symbol", "Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_assignment", "Gene Assignment", "GENE", "Gene symbol"]


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") if path.name.endswith(".gz") else path.open("r", encoding="utf-8", errors="replace", newline="")


def delimiter_for(path: Path) -> str:
    return "\t" if path.name.endswith((".txt", ".txt.gz", ".tsv", ".tsv.gz", ".annot", ".annot.gz")) else ","


def clean(value: str | None) -> str:
    return str(value or "").strip().strip('"')


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_axis_genes(axis_rows: list[dict[str, str]]) -> set[str]:
    genes: set[str] = set()
    for row in axis_rows:
        for gene in row.get("gene_members", "").replace(",", ";").split(";"):
            if gene.strip():
                genes.add(gene.strip().upper())
    return genes


def parse_aliases(path: Path) -> dict[str, str]:
    aliases: dict[str, str] = {}
    if not path.exists():
        return aliases
    for row in read_tsv(path):
        symbol = clean(row.get("gene_symbol", "")).upper()
        if not symbol:
            continue
        for key in ["ensembl_gene_id", "alias"]:
            value = clean(row.get(key, ""))
            if value:
                aliases[value.upper().split(".", 1)[0]] = symbol
    return aliases


def find_column(fieldnames: list[str], candidates: list[str]) -> str:
    lowered = {field.casefold(): field for field in fieldnames}
    for candidate in candidates:
        if candidate.casefold() in lowered:
            return lowered[candidate.casefold()]
    for field in fieldnames:
        if any(candidate.casefold() in field.casefold() for candidate in candidates):
            return field
    return ""


def split_gene_symbols(value: str, aliases: dict[str, str], axis_genes: set[str]) -> list[str]:
    tokens = re.split(r"\s*///\s*|[;,|]", clean(value))
    found: set[str] = set()
    for token in tokens:
        token = token.strip()
        if not token:
            continue
        for piece in re.split(r"\s+", token):
            upper = piece.strip().upper().split(".", 1)[0]
            if upper in axis_genes:
                found.add(upper)
            elif upper in aliases and aliases[upper] in axis_genes:
                found.add(aliases[upper])
    return sorted(found)


def read_platform(path: Path, axis_genes: set[str], aliases: dict[str, str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with open_text(path) as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#") and not line.startswith("!")), delimiter=delimiter_for(path))
        fieldnames = reader.fieldnames or []
        probe_col = find_column(fieldnames, PROBE_COLUMNS)
        gene_col = find_column(fieldnames, GENE_COLUMNS)
        if not probe_col or not gene_col:
            raise SystemExit(f"Could not identify probe/gene columns in platform file. Columns={fieldnames}")
        for row in reader:
            probe = clean(row.get(probe_col, ""))
            if not probe:
                continue
            for symbol in split_gene_symbols(row.get(gene_col, ""), aliases, axis_genes):
                rows.append(
                    {
                        "probe_id": probe,
                        "gene_symbol": symbol,
                        "source_gene_field": gene_col,
                        "mapping_status": "neurofate_axis_probe",
                    }
                )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare GEO platform probe mapping for NeuroFate axis genes only.")
    parser.add_argument("--platform-file", type=Path, required=True)
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--alias-table", type=Path, default=Path("metadata/neurofate_axis_gene_aliases.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--log-file", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = read_platform(args.platform_file, parse_axis_genes(read_tsv(args.axis_registry)), parse_aliases(args.alias_table))
    write_tsv(args.output, rows, ["probe_id", "gene_symbol", "source_gene_field", "mapping_status"])
    logging.info("Wrote %d NeuroFate-only probe mappings", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
