#!/usr/bin/env python3
"""Compare Phase 16, 17, 18, and 20 GSE243639 PD validation outputs."""

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


def feature_group_match_rate(path: Path) -> str:
    rows = read_tsv(path)
    for row in rows:
        if row.get("annotation_match_rate"):
            return row["annotation_match_rate"]
    return "unavailable"


def simple_match_rate(path: Path) -> str:
    rows = read_tsv(path)
    return rows[0].get("match_rate", "unavailable") if rows else "unavailable"


def interpretation(phase: str, flag: str) -> str:
    if phase == "phase16_global_features":
        return "valid global sample-level PD extension"
    if phase == "phase17_defective_celltype_attempt":
        return "technically defective annotation attempt; do not use biologically"
    if phase == "phase18_annotation_join_failure":
        return "annotation join still failed before safe-map consumption repair"
    if flag == "moderate_pd_internal_validation":
        return "corrected safe-map cell-type-aware analysis may support moderate internal PD evidence"
    if flag == "preliminary_pd_internal_signal":
        return "corrected safe-map analysis remains preliminary"
    if flag == "technical_failure_annotation_join":
        return "corrected safe-map path still failed technical thresholds"
    return "weak or unavailable corrected PD signal"


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
            "technical_failure_reported",
        ),
        (
            "phase18_annotation_join_failure",
            args.phase18_features,
            args.phase18_metrics,
            "repeated_stratified_split",
            simple_match_rate(args.phase18_annotation_summary),
        ),
        (
            "phase20_safe_map_celltype_features",
            args.phase20_features,
            args.phase20_metrics,
            "repeated_stratified_split",
            feature_group_match_rate(args.phase20_feature_groups),
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
                "interpretation": interpretation(phase, flag),
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
        "# Phase 20 PD Repair Summary",
        "",
        "- Phase 16 is valid as the global sample-level PD extension.",
        "- Phase 17 was technically defective and should not be interpreted biologically.",
        "- Phase 18 still failed annotation joining before safe-map consumption was repaired.",
        "- Phase 20 is the corrected safe-map-based cell-type-aware analysis if annotation match rate succeeds.",
        "",
    ]
    for row in rows:
        lines.append(
            f"- {row['phase']}: features={row['feature_count']}, match_rate={row['annotation_match_rate']}, AUROC={row['auroc']}, reliability={row['reliability_flag']}."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare Phase 16/17/18/20 GSE243639 PD validation.")
    parser.add_argument("--phase16-features", type=Path, default=Path("results/tables/phase16_gse243639_feature_table.tsv"))
    parser.add_argument("--phase17-features", type=Path, default=Path("results/tables/phase17_gse243639_celltype_feature_table.tsv"))
    parser.add_argument("--phase18-features", type=Path, default=Path("results/tables/phase18_gse243639_celltype_feature_table.tsv"))
    parser.add_argument("--phase20-features", type=Path, default=Path("results/tables/phase20_gse243639_celltype_feature_table.tsv"))
    parser.add_argument("--phase16-metrics", type=Path, default=Path("results/tables/phase16_gse243639_external_validation_metrics.tsv"))
    parser.add_argument("--phase17-metrics", type=Path, default=Path("results/tables/phase17_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--phase18-metrics", type=Path, default=Path("results/tables/phase18_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--phase20-metrics", type=Path, default=Path("results/tables/phase20_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--phase18-annotation-summary", type=Path, default=Path("results/tables/phase18_gse243639_annotation_match_summary.tsv"))
    parser.add_argument("--phase20-feature-groups", type=Path, default=Path("results/tables/phase20_gse243639_feature_group_counts.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase20_pd_validation_comparison.tsv"))
    parser.add_argument("--summary-output", type=Path, default=Path("results/reports/phase20_pd_repair_summary.md"))
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
