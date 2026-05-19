#!/usr/bin/env python3
"""Audit Phase 19 safe annotation map schema before Phase 20 feature building."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from pathlib import Path
from typing import TextIO


ID_COLUMN_HINTS = ["cell_id_expression", "cell_id", "normalized_cell_id", "cell_id_annotation_original", "cell_id_annotation_normalized"]
REQUIRED_COLUMNS = [
    "cell_id_expression",
    "cell_id_annotation_original",
    "cell_id_annotation_normalized",
    "normalized_cell_id",
    "barcode_core",
    "sample_id",
    "cell_type",
    "cluster_id",
    "annotation_source_sheet",
    "annotation_column_used",
    "match_status",
    "normalization_rule",
]


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_expression_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cell_id = row.get("cell_id", "")
            if cell_id:
                ids.add(cell_id)
    return ids


def read_safe_map(path: Path, preview_limit: int = 20) -> tuple[list[str], list[dict[str, str]], dict[str, set[str]]]:
    value_sets: dict[str, set[str]] = {}
    preview: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        for field in fieldnames:
            value_sets[field] = set()
        for index, row in enumerate(reader):
            if index < preview_limit:
                preview.append(row)
            for field in fieldnames:
                value = row.get(field, "")
                if value:
                    value_sets[field].add(value)
    return fieldnames, preview, value_sets


def candidate_columns(fieldnames: list[str]) -> list[str]:
    candidates = [field for field in ID_COLUMN_HINTS if field in fieldnames]
    candidates.extend(
        field
        for field in fieldnames
        if field not in candidates and ("cell" in field.lower() or field.lower().endswith("_id"))
    )
    return candidates


def write_preview(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_audit(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "audit_item",
        "value",
        "candidate_join_column",
        "overlap_count",
        "overlap_rate",
        "recommended_join_column",
        "notes",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_audit(safe_map: Path, expression: Path) -> tuple[list[dict[str, str]], list[dict[str, str]], list[str]]:
    expression_ids = read_expression_ids(expression)
    fieldnames, preview, values = read_safe_map(safe_map)
    candidates = candidate_columns(fieldnames)
    overlap_rows: list[dict[str, str]] = []
    best_column = ""
    best_overlap = -1
    denominator = max(1, len(expression_ids))
    for column in candidates:
        overlap = len(expression_ids & values.get(column, set()))
        if overlap > best_overlap:
            best_overlap = overlap
            best_column = column
        overlap_rows.append(
            {
                "audit_item": "candidate_join_column_overlap",
                "value": column,
                "candidate_join_column": column,
                "overlap_count": str(overlap),
                "overlap_rate": f"{overlap / denominator:.8g}",
                "recommended_join_column": "",
                "notes": "Direct overlap with expression cell_id values.",
            }
        )
    rows = [
        {
            "audit_item": "safe_map_columns",
            "value": ",".join(fieldnames),
            "candidate_join_column": "",
            "overlap_count": "",
            "overlap_rate": "",
            "recommended_join_column": best_column,
            "notes": "Schema columns present in safe map.",
        },
        {
            "audit_item": "missing_required_columns",
            "value": ",".join(column for column in REQUIRED_COLUMNS if column not in fieldnames),
            "candidate_join_column": "",
            "overlap_count": "",
            "overlap_rate": "",
            "recommended_join_column": best_column,
            "notes": "Should be empty for a Phase 20-ready safe map.",
        },
        {
            "audit_item": "expression_cell_id_column_name",
            "value": "cell_id",
            "candidate_join_column": "",
            "overlap_count": str(len(expression_ids)),
            "overlap_rate": "1",
            "recommended_join_column": best_column,
            "notes": "Expression sparse table join key.",
        },
        *overlap_rows,
        {
            "audit_item": "has_cell_type_or_cluster",
            "value": str("cell_type" in fieldnames or "cluster_id" in fieldnames).lower(),
            "candidate_join_column": "",
            "overlap_count": "",
            "overlap_rate": "",
            "recommended_join_column": best_column,
            "notes": "Feature builder requires cell_type or cluster_id.",
        },
        {
            "audit_item": "has_sample_id",
            "value": str("sample_id" in fieldnames).lower(),
            "candidate_join_column": "",
            "overlap_count": "",
            "overlap_rate": "",
            "recommended_join_column": best_column,
            "notes": "Sample-level aggregation requires sample_id.",
        },
        {
            "audit_item": "has_match_status",
            "value": str("match_status" in fieldnames).lower(),
            "candidate_join_column": "",
            "overlap_count": "",
            "overlap_rate": "",
            "recommended_join_column": best_column,
            "notes": "Safe map should mark linked rows as matched.",
        },
    ]
    return rows, preview, fieldnames


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit Phase 19 safe annotation map schema.")
    parser.add_argument("--safe-map", type=Path, default=Path("data/interim/external/gse243639_pd_snpc/gse243639_safe_cell_annotation_map.tsv"))
    parser.add_argument("--expression", type=Path, default=Path("data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase20_safe_annotation_map_schema_audit.tsv"))
    parser.add_argument("--preview-output", type=Path, default=Path("results/reports/phase20_safe_annotation_map_preview.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/100_audit_phase19_safe_annotation_map_schema.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows, preview, fieldnames = build_audit(args.safe_map, args.expression)
    write_audit(args.output, rows)
    write_preview(args.preview_output, preview, fieldnames)
    logging.info("Wrote Phase 20 safe-map schema audit.")
    print(f"Wrote {args.output}")
    print(f"Wrote {args.preview_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
