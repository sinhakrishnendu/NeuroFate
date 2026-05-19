#!/usr/bin/env python3
"""Generate endpoint-locked Phase 22 PNAS-style figures with matplotlib only."""

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


def placeholder(path: Path, title: str, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.axis("off")
    ax.text(0.5, 0.62, title, ha="center", fontsize=15, fontweight="bold")
    ax.text(0.5, 0.40, message, ha="center", fontsize=11, wrap=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_axis_effects(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        placeholder(path, "Endpoint-Locked Axis Effects", "Run scripts/112_test_axis_associations_endpoint_locked.py first.")
        return
    rows = rows[:30]
    labels = [f"{row['endpoint_id']}:{row['axis_id']}" for row in rows]
    values = [to_float(row.get("effect_size"), 0.0) for row in rows]
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.bar(range(len(rows)), values, color=["#4c78a8" if value >= 0 else "#e45756" for value in values])
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Endpoint-locked effect size")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, rotation=70, ha="right", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_primary_comparison(path: Path, rows: list[dict[str, str]]) -> None:
    primary = [row for row in rows if row.get("comparison_type", "").startswith("primary")]
    if not primary:
        placeholder(path, "Primary AD/PD Axis Comparison", "Run scripts/113_compare_ad_pd_axes_endpoint_locked.py first.")
        return
    fig, ax = plt.subplots(figsize=(6, 6))
    x = [to_float(row.get("ad_effect_size"), 0.0) for row in primary]
    y = [to_float(row.get("pd_effect_size"), 0.0) for row in primary]
    ax.scatter(x, y, s=55, color="#72b7b2", edgecolor="black")
    for row, xval, yval in zip(primary, x, y, strict=False):
        ax.text(xval, yval, row.get("axis_id", "")[:16], fontsize=7)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("SEA-AD cognitive endpoint effect")
    ax.set_ylabel("GSE243639 PD diagnosis effect")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_random_controls(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        placeholder(path, "Endpoint-Locked Random Controls", "Run scripts/114_run_endpoint_locked_random_axis_controls.py first.")
        return
    rows = rows[:40]
    labels = [f"{row['endpoint_id']}:{row['axis_id']}" for row in rows]
    values = [to_float(row.get("empirical_pvalue"), 1.0) for row in rows]
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.bar(range(len(rows)), values, color=["#54a24b" if value < 0.05 else "#f58518" for value in values])
    ax.axhline(0.05, color="black", linestyle="--", linewidth=1)
    ax.set_ylabel("Matched random-control empirical p-value")
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels(labels, rotation=70, ha="right", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_claim_strength(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        placeholder(path, "Endpoint-Locked Claim Strength", "Run scripts/115_build_endpoint_locked_axis_evidence_table.py first.")
        return
    counts: dict[str, int] = {}
    for row in rows:
        key = row.get("axis_claim_class", "unknown")
        counts[key] = counts.get(key, 0) + 1
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar(list(counts), list(counts.values()), color="#4c78a8")
    ax.set_ylabel("Evidence rows")
    ax.set_xticklabels(list(counts), rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 22 endpoint-locked figures.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.figures_dir.mkdir(parents=True, exist_ok=True)
    plot_axis_effects(args.figures_dir / "figure56_endpoint_locked_axis_effects.png", read_tsv(args.tables_dir / "phase22_endpoint_locked_axis_statistics.tsv"))
    plot_primary_comparison(args.figures_dir / "figure57_primary_ad_pd_axis_comparison.png", read_tsv(args.tables_dir / "phase22_endpoint_locked_ad_pd_axis_similarity.tsv"))
    plot_random_controls(args.figures_dir / "figure58_endpoint_locked_random_controls.png", read_tsv(args.tables_dir / "phase22_endpoint_locked_axis_empirical_pvalues.tsv"))
    plot_claim_strength(args.figures_dir / "figure59_endpoint_locked_claim_strength.png", read_tsv(args.tables_dir / "phase22_endpoint_locked_axis_evidence_table.tsv"))
    print(f"Wrote Phase 22 endpoint-locked figures under {args.figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
