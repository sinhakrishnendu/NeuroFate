#!/usr/bin/env python3
"""Compare Phase 16, Phase 17, and repaired Phase 18 GSE243639 validation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FEATURE_PREFIXES = (
    "gene_mean__",
    "gene_detection__",
    "cell_fraction__",
    "celltype_gene_mean__",
    "celltype_gene_detection__",
    "index__",
)


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
    return sum(1 for field in header if field.startswith(FEATURE_PREFIXES))


def metric_row(path: Path, preferred_mode: str) -> dict[str, str]:
    rows = read_tsv(path)
    for row in rows:
        if row.get("validation_mode") == preferred_mode and row.get("model", "logistic_regression") in {"", "logistic_regression"}:
            return row
    return rows[0] if rows else {}


def match_rate(path: Path) -> str:
    rows = read_tsv(path)
    return rows[0].get("match_rate", "unavailable") if rows else "unavailable"


def interpretation(flag: str, phase: str) -> str:
    if phase == "phase17_defective_celltype_attempt":
        return "technical diagnostic failure; do not use for biological interpretation"
    if flag == "technical_failure_annotation_join":
        return "repaired annotation join still failed technical thresholds"
    if flag == "moderate_pd_internal_validation":
        return "repaired cell-type-aware signal may support moderate internal PD evidence"
    if flag == "preliminary_pd_internal_signal":
        return "repaired cell-type-aware signal remains preliminary"
    return "weak or unavailable PD signal"


def build_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    phases = [
        (
            "phase16_global_features",
            args.phase16_features,
            args.phase16_metrics,
            "gse243639_pd_internal",
            "unavailable",
        ),
        (
            "phase17_defective_celltype_attempt",
            args.phase17_features,
            args.phase17_metrics,
            "repeated_stratified_split",
            "technical_failure_reported_by_phase18",
        ),
        (
            "phase18_repaired_celltype_features",
            args.phase18_features,
            args.phase18_metrics,
            "repeated_stratified_split",
            match_rate(args.annotation_summary),
        ),
    ]
    rows = []
    for phase, feature_path, metric_path, preferred_mode, annotation_rate in phases:
        metric = metric_row(metric_path, preferred_mode)
        flag = metric.get("reliability_flag", "unavailable")
        rows.append(
            {
                "phase": phase,
                "feature_count": str(feature_count(feature_path)),
                "annotation_match_rate": annotation_rate,
                "auroc": metric.get("auroc", "unavailable"),
                "auprc": metric.get("auprc", "unavailable"),
                "balanced_accuracy": metric.get("balanced_accuracy", "unavailable"),
                "empirical_pvalue": metric.get("empirical_permutation_pvalue", "unavailable"),
                "reliability_flag": flag,
                "interpretation": interpretation(flag, phase),
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "phase",
        "feature_count",
        "annotation_match_rate",
        "auroc",
        "auprc",
        "balanced_accuracy",
        "empirical_pvalue",
        "reliability_flag",
        "interpretation",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    lines = [
        "# Phase 18 PD Repair Summary",
        "",
        "Phase 17 is treated as a technical diagnostic failure when annotation matching collapses cell-type features.",
        "Final PD interpretation should use Phase 18 when the repaired annotation join passes technical checks.",
        "",
    ]
    for row in rows:
        lines.append(
            f"- {row['phase']}: features={row['feature_count']}, AUROC={row['auroc']}, match_rate={row['annotation_match_rate']}, reliability={row['reliability_flag']}."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare Phase 16/17/18 GSE243639 PD validation.")
    parser.add_argument("--phase16-features", type=Path, default=Path("results/tables/phase16_gse243639_feature_table.tsv"))
    parser.add_argument("--phase17-features", type=Path, default=Path("results/tables/phase17_gse243639_celltype_feature_table.tsv"))
    parser.add_argument("--phase18-features", type=Path, default=Path("results/tables/phase18_gse243639_celltype_feature_table.tsv"))
    parser.add_argument("--phase16-metrics", type=Path, default=Path("results/tables/phase16_gse243639_external_validation_metrics.tsv"))
    parser.add_argument("--phase17-metrics", type=Path, default=Path("results/tables/phase17_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--phase18-metrics", type=Path, default=Path("results/tables/phase18_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--annotation-summary", type=Path, default=Path("results/tables/phase18_gse243639_annotation_match_summary.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase18_pd_validation_comparison.tsv"))
    parser.add_argument("--summary-output", type=Path, default=Path("results/reports/phase18_pd_repair_summary.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = build_rows(args)
    write_tsv(args.output, rows)
    write_summary(args.summary_output, rows)
    print(f"Wrote {args.output}")
    print(f"Wrote {args.summary_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
