#!/usr/bin/env python3
"""Generate Phase 7 cross-cohort validation figures from summary TSVs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def save_bar(labels: list[str], values: list[float], title: str, ylabel: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10, 5.8))
    positions = range(len(labels))
    plt.bar(positions, values, color="#4D6A8A")
    plt.xticks(positions, labels, rotation=45, ha="right", fontsize=8)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def figure16_crosscohort_generalization(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "phase7_generalization_summary.tsv")
    save_bar(
        [row["validation_mode"] for row in rows],
        [to_float(row["mean_auroc"]) for row in rows],
        "Cross-Cohort NeuroFate Generalization",
        "Mean AUROC",
        figures_dir / "figure16_crosscohort_generalization.png",
    )


def figure17_cohort_transfer(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "phase7_crosscohort_metrics.tsv")
    transfer = [row for row in rows if row["validation_mode"] == "train_sea_ad_test_external"]
    save_bar(
        [row["test_cohort"] for row in transfer],
        [to_float(row["auroc"]) for row in transfer],
        "SEA-AD to External Cohort Transfer",
        "AUROC",
        figures_dir / "figure17_cohort_transfer_performance.png",
    )


def figure18_feature_stability(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "crosscohort_feature_overlap.tsv")
    shared = sum(1 for row in rows if row["status"] == "shared")
    missing = max(0, len(rows) - shared)
    save_bar(
        ["shared", "cohort_specific_or_missing"],
        [shared, missing],
        "Cross-Cohort Feature Stability",
        "Feature count",
        figures_dir / "figure18_feature_stability.png",
    )


def figure19_multicohort_scores(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "crosscohort_donor_feature_table.tsv")
    values = [to_float(row.get("index__NVI", "0")) for row in rows]
    figures_dir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8.5, 5.5))
    plt.hist(values, bins=min(20, max(5, len(values) // 4 if values else 5)), color="#8A5A44")
    plt.xlabel("NVI donor feature")
    plt.ylabel("Donor count")
    plt.title("Multi-Cohort NeuroFate Score Proxy")
    plt.tight_layout()
    plt.savefig(figures_dir / "figure19_multicohort_neurofate_scores.png", dpi=300)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 7 cross-cohort figures.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    figure16_crosscohort_generalization(args.tables_dir, args.figures_dir)
    figure17_cohort_transfer(args.tables_dir, args.figures_dir)
    figure18_feature_stability(args.tables_dir, args.figures_dir)
    figure19_multicohort_scores(args.tables_dir, args.figures_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
