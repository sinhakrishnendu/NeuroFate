#!/usr/bin/env python3
"""Audit Mathys target-gene extraction outputs without raw data access."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import TextIO


TARGET_PANEL = Path("metadata/target_gene_panel_v1.tsv")
SPARSE_EXPRESSION = Path("data/interim/external/mathys_2019/mathys_sparse_gene_panel_expression.tsv.gz")
FEATURE_TABLE = Path("results/tables/mathys_2019_phase5_donor_feature_table.tsv")
LABEL_SUMMARY = Path("results/tables/mathys_2019_label_summary.tsv")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def target_genes(path: Path) -> set[str]:
    return {
        row.get("gene_symbol", "").strip()
        for row in read_tsv(path)
        if row.get("gene_symbol", "").strip()
    }


def sparse_stats(path: Path) -> tuple[set[str], int]:
    if not path.exists():
        return set(), 0
    genes: set[str] = set()
    rows = 0
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows += 1
            gene = row.get("gene_symbol", "").strip()
            if gene:
                genes.add(gene)
    return genes, rows


def diagnosis_groups(label_summary: Path) -> str:
    rows = read_tsv(label_summary)
    groups = [
        f"{row.get('label')}:{row.get('sample_count')}"
        for row in rows
        if row.get("label_field") == "label__diagnosis"
    ]
    return ";".join(groups) if groups else "unavailable"


def build_audit_row() -> dict[str, str]:
    requested = target_genes(TARGET_PANEL)
    extracted, sparse_rows = sparse_stats(SPARSE_EXPRESSION)
    feature_rows = read_tsv(FEATURE_TABLE)
    missing = sorted(requested - extracted) if requested else []
    if not SPARSE_EXPRESSION.exists():
        status = "missing_sparse_expression"
        warning = "Mathys sparse target-gene expression file is unavailable."
    elif not requested:
        status = "missing_target_panel"
        warning = "Target gene panel is unavailable."
    elif missing:
        status = "partial_extraction"
        warning = "Some requested target genes were not detected in extracted Mathys expression."
    else:
        status = "ok"
        warning = "External feasibility only; six sample units are insufficient for definitive validation."
    return {
        "dataset_id": "mathys_2019_gse138852",
        "requested_target_genes": str(len(requested)),
        "extracted_target_genes": str(len(extracted)),
        "missing_target_genes": ",".join(missing) if missing else "none",
        "sparse_expression_rows": str(sparse_rows),
        "feature_table_rows": str(len(feature_rows)),
        "sample_units": str(len(feature_rows)),
        "diagnosis_groups": diagnosis_groups(LABEL_SUMMARY),
        "status": status,
        "warning": warning,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit Mathys target-gene extraction outputs.")
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase13_mathys_gene_extraction_audit.tsv"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    row = build_audit_row()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "dataset_id",
                "requested_target_genes",
                "extracted_target_genes",
                "missing_target_genes",
                "sparse_expression_rows",
                "feature_table_rows",
                "sample_units",
                "diagnosis_groups",
                "status",
                "warning",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(row)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
