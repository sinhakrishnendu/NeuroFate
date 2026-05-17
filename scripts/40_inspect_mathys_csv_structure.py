#!/usr/bin/env python3
"""Inspect Mathys 2019 CSV count/covariate structure without single-cell workflows."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from pathlib import Path
from typing import TextIO


MATHYS_COVARIATE_COLUMNS = [
    "oupSample.batchCond",
    "oupSample.cellType",
    "oupSample.cellType_batchCond",
    "oupSample.subclustID",
    "oupSample.subclustCond",
]


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


def read_panel_genes(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {row["gene_symbol"].upper() for row in csv.DictReader(handle, delimiter="\t")}


def count_rows(path: Path) -> int:
    with open_text(path) as handle:
        reader = csv.reader(handle)
        next(reader, None)
        return sum(1 for _ in reader)


def inspect_counts(counts_path: Path, panel_genes: set[str], max_rows: int) -> dict[str, str]:
    with open_text(counts_path) as handle:
        reader = csv.reader(handle)
        header = next(reader)
        data_rows = []
        for index, row in enumerate(reader):
            if index >= max_rows:
                break
            data_rows.append(row)

    header_genes = sum(1 for item in header[1:] if item.upper() in panel_genes)
    first_column_genes = sum(1 for row in data_rows if row and row[0].upper() in panel_genes)
    if first_column_genes > header_genes:
        orientation = "genes_as_rows"
        inferred_cell_count = str(max(0, len(header) - 1))
        inferred_gene_axis = "first_column"
    elif header_genes > first_column_genes:
        orientation = "cells_as_rows"
        inferred_cell_count = str(count_rows(counts_path))
        inferred_gene_axis = "header"
    else:
        orientation = "ambiguous"
        inferred_cell_count = "unknown"
        inferred_gene_axis = "unknown"

    return {
        "counts_file": str(counts_path),
        "header_column_count": str(len(header)),
        "first_header": header[0] if header else "",
        "header_target_gene_hits": str(header_genes),
        "first_column_target_gene_hits_in_preview": str(first_column_genes),
        "preview_rows_examined": str(len(data_rows)),
        "inferred_orientation": orientation,
        "inferred_gene_axis": inferred_gene_axis,
        "inferred_cell_count": inferred_cell_count,
        "notes": "CSV-only inspection; no dense matrix or single-cell workflow",
    }


def inspect_covariates(covariates_path: Path, preview_output: Path, preview_rows: int) -> dict[str, str]:
    preview_output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(covariates_path) as source:
        reader = csv.DictReader(source)
        if reader.fieldnames is None:
            raise RuntimeError(f"Covariates CSV has no header: {covariates_path}")
        fieldnames = reader.fieldnames
        present = [column for column in MATHYS_COVARIATE_COLUMNS if column in fieldnames]
        likely_cell_id = fieldnames[0]
        row_count = 0
        with preview_output.open("w", encoding="utf-8", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in reader:
                row_count += 1
                if row_count <= preview_rows:
                    writer.writerow(row)

    return {
        "covariates_file": str(covariates_path),
        "covariate_column_count": str(len(fieldnames)),
        "covariate_row_count": str(row_count),
        "available_mathys_metadata_fields": ";".join(present),
        "likely_cell_identifier_column": likely_cell_id,
        "preview_output": str(preview_output),
        "notes": "Covariates preview only; no expression processing",
    }


def write_summary(path: Path, counts_info: dict[str, str], covariate_info: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {"section": "counts", "key": key, "value": value}
        for key, value in counts_info.items()
    ] + [
        {"section": "covariates", "key": key, "value": value}
        for key, value in covariate_info.items()
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["section", "key", "value"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote Mathys CSV structure summary: %s", path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Inspect Mathys GSE138852 CSV structure.")
    parser.add_argument("--counts", type=Path, default=Path("data/raw/external/mathys_2019/GSE138852_counts.csv.gz"))
    parser.add_argument("--covariates", type=Path, default=Path("data/raw/external/mathys_2019/GSE138852_covariates.csv.gz"))
    parser.add_argument("--panel", type=Path, default=Path("metadata/target_gene_panel_v1.tsv"))
    parser.add_argument("--summary-output", type=Path, default=Path("results/tables/mathys_csv_structure_summary.tsv"))
    parser.add_argument("--preview-output", type=Path, default=Path("data/interim/external/mathys_2019/mathys_covariates_preview.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/40_inspect_mathys_csv_structure.log"))
    parser.add_argument("--max-count-preview-rows", type=int, default=500)
    parser.add_argument("--covariate-preview-rows", type=int, default=25)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    logging.info("Starting Mathys CSV structure inspection.")
    panel_genes = read_panel_genes(args.panel)
    counts_info = inspect_counts(args.counts, panel_genes, args.max_count_preview_rows)
    covariate_info = inspect_covariates(args.covariates, args.preview_output, args.covariate_preview_rows)
    write_summary(args.summary_output, counts_info, covariate_info)
    logging.info("Mathys CSV inspection complete. No Scanpy, H5AD, or clustering workflow was used.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
