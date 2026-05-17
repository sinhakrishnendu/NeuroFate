#!/usr/bin/env python3
"""Compare Mathys 2019 var genes with the NeuroFate target panel."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file, mode="w", encoding="utf-8"),
        ],
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %s rows: %d", path, len(rows))


def build_overlap(panel_rows: list[dict[str, str]], var_rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    symbol_lookup: dict[str, dict[str, str]] = {}
    id_lookup: dict[str, dict[str, str]] = {}
    for row in var_rows:
        symbol = row.get("gene_symbol", "")
        gene_id = row.get("gene_id", "")
        if symbol:
            symbol_lookup[symbol.upper()] = row
        if gene_id:
            id_lookup[gene_id.upper()] = row

    overlap_rows: list[dict[str, str]] = []
    missing_rows: list[dict[str, str]] = []
    for panel in panel_rows:
        gene = panel["gene_symbol"]
        match = symbol_lookup.get(gene.upper()) or id_lookup.get(gene.upper())
        status = "present" if match else "missing"
        row = {
            "gene_symbol": gene,
            "priority_tier": panel.get("priority_tier", ""),
            "biological_role": panel.get("biological_role", ""),
            "mathys_status": status,
            "mathys_gene_symbol": match.get("gene_symbol", "") if match else "",
            "mathys_gene_id": match.get("gene_id", "") if match else "",
            "mathys_var_index": match.get("var_index", "") if match else "",
            "notes": "metadata-only panel overlap; no expression values loaded",
        }
        overlap_rows.append(row)
        if status == "missing":
            missing_rows.append(row)
    return overlap_rows, missing_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare Mathys 2019 target-gene overlap tables.")
    parser.add_argument(
        "--var",
        type=Path,
        default=Path("data/interim/external/mathys_2019/mathys_var_genes.tsv"),
    )
    parser.add_argument("--panel", type=Path, default=Path("metadata/target_gene_panel_v1.tsv"))
    parser.add_argument(
        "--overlap-output",
        type=Path,
        default=Path("results/tables/mathys_gene_overlap.tsv"),
    )
    parser.add_argument(
        "--missing-output",
        type=Path,
        default=Path("results/tables/mathys_missing_target_genes.tsv"),
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/37_prepare_mathys_gene_panel_overlap.log"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    panel_rows = read_tsv(args.panel)
    var_rows = read_tsv(args.var)
    overlap_rows, missing_rows = build_overlap(panel_rows, var_rows)
    fieldnames = [
        "gene_symbol",
        "priority_tier",
        "biological_role",
        "mathys_status",
        "mathys_gene_symbol",
        "mathys_gene_id",
        "mathys_var_index",
        "notes",
    ]
    write_tsv(args.overlap_output, overlap_rows, fieldnames)
    write_tsv(args.missing_output, missing_rows, fieldnames)
    logging.info("Panel genes: %d", len(panel_rows))
    logging.info("Mathys var genes: %d", len(var_rows))
    logging.info("Target genes present: %d", len(overlap_rows) - len(missing_rows))
    logging.info("Target genes missing: %d", len(missing_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
