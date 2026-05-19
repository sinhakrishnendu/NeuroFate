#!/usr/bin/env python3
"""Generate PNAS-oriented NeuroFate axis figures from donor/sample-level tables."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: str | None, default: float = math.nan) -> float:
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def save_placeholder(path: Path, title: str, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.axis("off")
    ax.text(0.5, 0.62, title, ha="center", va="center", fontsize=16, fontweight="bold")
    ax.text(0.5, 0.40, message, ha="center", va="center", fontsize=11, wrap=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_concept(path: Path) -> None:
    labels = ["SEA-AD\nAD anchor", "Axis scores\nDonor level", "GSE243639\nPD cohort"]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    ax.axis("off")
    for x, label, color in zip([1.5, 5.0, 8.5], labels, ["#4c78a8", "#72b7b2", "#f58518"], strict=False):
        ax.add_patch(plt.Rectangle((x - 1.1, 1.55), 2.2, 0.9, facecolor=color, alpha=0.85, edgecolor="black"))
        ax.text(x, 2.0, label, ha="center", va="center", fontsize=11, color="white", fontweight="bold")
    ax.annotate("", xy=(3.8, 2.0), xytext=(2.7, 2.0), arrowprops={"arrowstyle": "->", "lw": 2})
    ax.annotate("", xy=(7.3, 2.0), xytext=(6.2, 2.0), arrowprops={"arrowstyle": "->", "lw": 2})
    ax.text(5.0, 0.75, "Candidate shared and disease-specific neurodegeneration axes", ha="center", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def axis_score_means(rows: list[dict[str, str]]) -> dict[str, float]:
    if not rows:
        return {}
    axes = [column for column in rows[0] if column.startswith("axis__")]
    means: dict[str, float] = {}
    for axis in axes:
        values = [to_float(row.get(axis)) for row in rows]
        values = [value for value in values if not math.isnan(value)]
        if values:
            means[axis.replace("axis__", "")] = sum(values) / len(values)
    return means


def plot_axis_scores(path: Path, sea_rows: list[dict[str, str]], pd_rows: list[dict[str, str]]) -> None:
    sea = axis_score_means(sea_rows)
    pd = axis_score_means(pd_rows)
    axes = sorted(set(sea) | set(pd))
    if not axes:
        save_placeholder(path, "Axis Scores By Disease", "Run scripts/106_build_neurofate_axis_scores.py first.")
        return
    fig, ax = plt.subplots(figsize=(10, 5.5))
    positions = range(len(axes))
    width = 0.38
    ax.bar([x - width / 2 for x in positions], [sea.get(axis, 0.0) for axis in axes], width=width, label="SEA-AD")
    ax.bar([x + width / 2 for x in positions], [pd.get(axis, 0.0) for axis in axes], width=width, label="GSE243639")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(list(positions))
    ax.set_xticklabels(axes, rotation=45, ha="right")
    ax.set_ylabel("Mean standardized axis score")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_similarity(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        save_placeholder(path, "AD/PD Axis Similarity", "Run scripts/108_compare_ad_pd_axis_patterns.py first.")
        return
    axes = [row.get("axis_id", "") for row in rows]
    values = [1.0 if row.get("same_direction") == "true" else 0.0 for row in rows]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(axes, values, color="#72b7b2")
    ax.set_ylim(-0.1, 1.1)
    ax.set_ylabel("Direction agreement")
    ax.set_xticklabels(axes, rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_classification(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        save_placeholder(path, "Shared Vs Disease-Specific Axes", "Run scripts/108_compare_ad_pd_axis_patterns.py first.")
        return
    counts: dict[str, int] = {}
    for row in rows:
        category = row.get("classification", "unknown")
        counts[category] = counts.get(category, 0) + 1
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(list(counts), list(counts.values()), color="#4c78a8")
    ax.set_ylabel("Axis count")
    ax.set_xticklabels(list(counts), rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_random_controls(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        save_placeholder(path, "Random Axis Controls", "Run scripts/109_run_axis_randomization_controls.py first.")
        return
    axes = [f"{row.get('cohort', '')}:{row.get('axis_id', '')}" for row in rows]
    values = [to_float(row.get("empirical_pvalue"), math.nan) for row in rows]
    colors = ["#54a24b" if value < 0.10 else "#e45756" for value in values]
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.bar(range(len(axes)), [0 if math.isnan(value) else value for value in values], color=colors)
    ax.axhline(0.10, color="black", linewidth=1, linestyle="--")
    ax.set_ylabel("Empirical p-value")
    ax.set_xticks(range(len(axes)))
    ax.set_xticklabels(axes, rotation=60, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 21 PNAS-style axis figures.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.figures_dir.mkdir(parents=True, exist_ok=True)
    plot_concept(args.figures_dir / "figure51_neurofate_pnas_concept.png")
    plot_axis_scores(
        args.figures_dir / "figure52_axis_scores_by_disease.png",
        read_tsv(args.tables_dir / "phase21_sea_ad_axis_scores.tsv"),
        read_tsv(args.tables_dir / "phase21_gse243639_axis_scores.tsv"),
    )
    similarity = read_tsv(args.tables_dir / "phase21_ad_pd_axis_similarity.tsv")
    shared = read_tsv(args.tables_dir / "phase21_shared_vs_disease_specific_axes.tsv")
    plot_similarity(args.figures_dir / "figure53_ad_pd_axis_similarity.png", similarity)
    plot_classification(args.figures_dir / "figure54_shared_vs_specific_axes.png", shared)
    plot_random_controls(
        args.figures_dir / "figure55_random_axis_controls.png",
        read_tsv(args.tables_dir / "phase21_axis_empirical_pvalues.tsv"),
    )
    print(f"Wrote Phase 21 figures under {args.figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
