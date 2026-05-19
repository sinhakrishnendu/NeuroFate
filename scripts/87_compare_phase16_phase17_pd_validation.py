#!/usr/bin/env python3
"""Compare Phase 16 global and Phase 17 cell-type-aware PD validation outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def feature_count(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), [])
    return sum(
        1
        for field in header
        if field.startswith(("gene_mean__", "gene_detection__", "cell_fraction__", "celltype_gene_mean__", "celltype_gene_detection__", "index__"))
    )


def select_metric(rows: list[dict[str, str]], preferred_mode: str) -> dict[str, str]:
    for row in rows:
        if row.get("validation_mode") == preferred_mode and row.get("model", "logistic_regression") == "logistic_regression":
            return row
    return rows[0] if rows else {}


def build_comparison(args: argparse.Namespace) -> list[dict[str, str]]:
    phase16 = select_metric(read_tsv(args.phase16_metrics), "gse243639_pd_internal")
    phase17 = select_metric(read_tsv(args.phase17_metrics), "repeated_stratified_split")
    return [
        {
            "phase": "phase16_global_sample_features",
            "feature_table": str(args.phase16_features),
            "feature_count": str(feature_count(args.phase16_features)),
            "auroc": phase16.get("auroc", "unavailable"),
            "auprc": phase16.get("auprc", "unavailable"),
            "balanced_accuracy": phase16.get("balanced_accuracy", "unavailable"),
            "brier_score": phase16.get("brier_score", "unavailable"),
            "reliability_flag": phase16.get("reliability_flag", "unavailable"),
        },
        {
            "phase": "phase17_celltype_aware_features",
            "feature_table": str(args.phase17_features),
            "feature_count": str(feature_count(args.phase17_features)),
            "auroc": phase17.get("auroc", "unavailable"),
            "auprc": phase17.get("auprc", "unavailable"),
            "balanced_accuracy": phase17.get("balanced_accuracy", "unavailable"),
            "brier_score": phase17.get("brier_score", "unavailable"),
            "reliability_flag": phase17.get("reliability_flag", "unavailable"),
        },
    ]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["phase", "feature_table", "feature_count", "auroc", "auprc", "balanced_accuracy", "brier_score", "reliability_flag"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    phase16 = rows[0]
    phase17 = rows[1]
    lines = [
        "# Phase 17 PD Improvement Summary",
        "",
        "This comparison is descriptive and should not be used as a clinical or cross-disease claim.",
        "",
        f"- Phase 16 AUROC: {phase16['auroc']} with reliability `{phase16['reliability_flag']}`.",
        f"- Phase 17 AUROC: {phase17['auroc']} with reliability `{phase17['reliability_flag']}`.",
        f"- Phase 16 feature count: {phase16['feature_count']}.",
        f"- Phase 17 feature count: {phase17['feature_count']}.",
        "",
        "Interpretation: cell-type-aware features are considered an improvement only if AUROC, stability, permutation support, and reliability category improve together.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare Phase 16 and Phase 17 GSE243639 validation.")
    parser.add_argument("--phase16-metrics", type=Path, default=Path("results/tables/phase16_gse243639_external_validation_metrics.tsv"))
    parser.add_argument("--phase17-metrics", type=Path, default=Path("results/tables/phase17_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--phase16-features", type=Path, default=Path("results/tables/phase16_gse243639_feature_table.tsv"))
    parser.add_argument("--phase17-features", type=Path, default=Path("results/tables/phase17_gse243639_celltype_feature_table.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase17_pd_validation_comparison.tsv"))
    parser.add_argument("--summary-output", type=Path, default=Path("results/reports/phase17_pd_improvement_summary.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = build_comparison(args)
    write_tsv(args.output, rows)
    write_summary(args.summary_output, rows)
    print(f"Wrote {args.output}")
    print(f"Wrote {args.summary_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
