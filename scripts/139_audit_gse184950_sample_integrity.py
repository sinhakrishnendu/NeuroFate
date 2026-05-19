#!/usr/bin/env python3
"""Audit GSE184950 sample identifiers in axis outputs against series-matrix metadata."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


OUTPUT_COLUMNS = [
    "expected_samples",
    "observed_axis_score_samples",
    "valid_axis_score_samples",
    "invalid_axis_score_samples",
    "invalid_sample_ids",
    "axis_gene_audit_samples",
    "invalid_axis_gene_audit_samples",
    "invalid_axis_gene_audit_ids",
]


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def audit(axis_scores: list[dict[str, str]], metadata: list[dict[str, str]], axis_gene_audit: list[dict[str, str]]) -> dict[str, str]:
    expected = {row.get("sample_name", "") for row in metadata if row.get("sample_name")}
    observed_axis = {row.get("sample_id", "") for row in axis_scores if row.get("sample_id")}
    observed_audit = {row.get("sample_id", "") for row in axis_gene_audit if row.get("sample_id")}
    invalid_axis = sorted(observed_axis - expected)
    invalid_audit = sorted(observed_audit - expected)
    valid_axis = sorted(observed_axis & expected)
    return {
        "expected_samples": str(len(expected)),
        "observed_axis_score_samples": str(len(observed_axis)),
        "valid_axis_score_samples": str(len(valid_axis)),
        "invalid_axis_score_samples": str(len(invalid_axis)),
        "invalid_sample_ids": ";".join(invalid_axis),
        "axis_gene_audit_samples": str(len(observed_audit)),
        "invalid_axis_gene_audit_samples": str(len(invalid_audit)),
        "invalid_axis_gene_audit_ids": ";".join(invalid_audit),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit GSE184950 sample integrity.")
    parser.add_argument("--axis-scores", type=Path, default=Path("results/tables/phase25_gse184950_axis_scores.tsv"))
    parser.add_argument("--sample-metadata", type=Path, default=Path("results/tables/phase25_gse184950_series_sample_metadata.tsv"))
    parser.add_argument("--axis-gene-audit", type=Path, default=Path("results/tables/phase26_gse184950_axis_gene_extraction_audit.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase27_gse184950_sample_integrity_audit.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/139_audit_gse184950_sample_integrity.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    row = audit(read_tsv(args.axis_scores), read_tsv(args.sample_metadata), read_tsv(args.axis_gene_audit))
    write_tsv(args.output, [row])
    logging.info("GSE184950 sample integrity audit: %s", row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
