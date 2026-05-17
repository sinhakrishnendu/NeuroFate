#!/usr/bin/env python3
"""Extract NeuroFate target genes from Mathys 2019 CSV counts as sparse-like TSV."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from pathlib import Path
from typing import TextIO


MAX_TARGET_GENES = 64
OUTPUT_COLUMNS = ["cell_id", "gene_symbol", "expression_value"]


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


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def open_output(path: Path) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return path.open("w", encoding="utf-8", newline="")


def read_panel_genes(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        genes = [row["gene_symbol"] for row in csv.DictReader(handle, delimiter="\t")]
    if len(genes) > MAX_TARGET_GENES:
        raise ValueError(f"Target panel has {len(genes)} genes; maximum allowed is {MAX_TARGET_GENES}.")
    return genes


def infer_orientation(counts_path: Path, panel_genes: set[str], max_rows: int = 500) -> str:
    with open_text(counts_path) as handle:
        reader = csv.reader(handle)
        header = next(reader)
        header_hits = sum(1 for item in header[1:] if item.upper() in panel_genes)
        first_column_hits = 0
        for index, row in enumerate(reader):
            if index >= max_rows:
                break
            if row and row[0].upper() in panel_genes:
                first_column_hits += 1
    if first_column_hits > header_hits:
        return "genes_as_rows"
    if header_hits > first_column_hits:
        return "cells_as_rows"
    raise RuntimeError("Could not infer Mathys count orientation; pass --orientation explicitly.")


def is_nonzero(value: str) -> bool:
    try:
        return float(value) != 0.0
    except ValueError:
        return False


def extract_genes_as_rows(counts_path: Path, panel_genes: set[str], output_path: Path) -> int:
    written = 0
    with open_text(counts_path) as source, open_output(output_path) as target:
        reader = csv.reader(source)
        header = next(reader)
        cell_ids = header[1:]
        writer = csv.DictWriter(target, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in reader:
            if not row:
                continue
            gene = row[0]
            if gene.upper() not in panel_genes:
                continue
            for cell_id, value in zip(cell_ids, row[1:], strict=False):
                if not is_nonzero(value):
                    continue
                writer.writerow({"cell_id": cell_id, "gene_symbol": gene, "expression_value": value})
                written += 1
    return written


def extract_cells_as_rows(counts_path: Path, panel_genes: set[str], output_path: Path) -> int:
    written = 0
    with open_text(counts_path) as source, open_output(output_path) as target:
        reader = csv.reader(source)
        header = next(reader)
        target_columns = [
            (index, gene)
            for index, gene in enumerate(header)
            if index > 0 and gene.upper() in panel_genes
        ]
        writer = csv.DictWriter(target, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in reader:
            if not row:
                continue
            cell_id = row[0]
            for index, gene in target_columns:
                if index >= len(row) or not is_nonzero(row[index]):
                    continue
                writer.writerow({"cell_id": cell_id, "gene_symbol": gene, "expression_value": row[index]})
                written += 1
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract Mathys target genes from CSV counts.")
    parser.add_argument("--counts", type=Path, default=Path("data/raw/external/mathys_2019/GSE138852_counts.csv.gz"))
    parser.add_argument("--panel", type=Path, default=Path("metadata/target_gene_panel_v1.tsv"))
    parser.add_argument("--output", type=Path, default=Path("data/interim/external/mathys_2019/mathys_sparse_gene_panel_expression.tsv.gz"))
    parser.add_argument("--orientation", choices=["auto", "genes_as_rows", "cells_as_rows"], default="auto")
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/41_extract_mathys_target_gene_panel.log"))
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    genes = read_panel_genes(args.panel)
    panel_gene_set = {gene.upper() for gene in genes}
    orientation = infer_orientation(args.counts, panel_gene_set) if args.orientation == "auto" else args.orientation
    logging.info("Mathys count orientation: %s", orientation)
    logging.info("Target genes requested: %d", len(genes))
    logging.info("Output: %s", args.output)
    if args.dry_run:
        logging.info("Dry run only. No sparse-like expression output was written.")
        return 0
    if orientation == "genes_as_rows":
        written = extract_genes_as_rows(args.counts, panel_gene_set, args.output)
    else:
        written = extract_cells_as_rows(args.counts, panel_gene_set, args.output)
    logging.info("Wrote Mathys sparse-like target-gene rows: %d", written)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
