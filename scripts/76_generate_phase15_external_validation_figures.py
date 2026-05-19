#!/usr/bin/env python3
"""Generate Phase 15 external validation planning figures from existing TSVs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def simple_bar(path: Path, title: str, labels: list[str], values: list[float], ylabel: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4))
    if labels:
        ax.bar(range(len(labels)), values)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
    else:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center")
        ax.axis("off")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Phase 15 external validation figures.")
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    args = parser.parse_args()
    triage = read_tsv(Path("results/reports/phase15_external_dataset_triage.tsv"))
    counts: dict[str, int] = {}
    for row in triage:
        counts[row.get("readiness_category", "unknown")] = counts.get(row.get("readiness_category", "unknown"), 0) + 1
    simple_bar(args.figures_dir / "figure31_external_dataset_triage.png", "External Dataset Triage", list(counts), list(counts.values()), "Datasets")

    registry = read_tsv(Path("metadata/phase15_external_validation_candidates.tsv"))
    priorities: dict[str, int] = {}
    for row in registry:
        priorities[row.get("priority", "unknown")] = priorities.get(row.get("priority", "unknown"), 0) + 1
    simple_bar(args.figures_dir / "figure32_external_feature_overlap.png", "External Cohort Priorities", list(priorities), list(priorities.values()), "Datasets")

    metrics = read_tsv(Path("results/tables/phase15_multi_external_validation_metrics.tsv"))
    overlap_labels = [row.get("dataset_id", "unknown") for row in metrics]
    overlap_values = [float(row.get("feature_overlap_count", 0) or 0) for row in metrics]
    simple_bar(args.figures_dir / "figure33_multi_external_validation.png", "Multi-External Feature Overlap", overlap_labels, overlap_values, "Shared features")

    reliability = read_tsv(Path("results/reports/phase15_external_validation_reliability.tsv"))
    reliability_counts: dict[str, int] = {}
    for row in reliability:
        key = row.get("reliability_category", "not_run")
        reliability_counts[key] = reliability_counts.get(key, 0) + 1
    simple_bar(args.figures_dir / "figure34_external_reliability_matrix.png", "External Reliability Matrix", list(reliability_counts), list(reliability_counts.values()), "Cohorts")
    print(f"Wrote Phase 15 figures under {args.figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
