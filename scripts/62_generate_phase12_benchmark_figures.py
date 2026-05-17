#!/usr/bin/env python3
"""Generate Phase 12 benchmark figures from existing tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def table_or_empty(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def placeholder(path: Path, title: str, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_repeated(summary: pd.DataFrame, path: Path) -> None:
    if summary.empty or "auroc_mean" not in summary:
        placeholder(path, "Repeated Benchmark Stability", "Repeated benchmark summary is unavailable.")
        return
    top = summary.sort_values("auroc_mean", ascending=False).head(20)
    fig, ax = plt.subplots(figsize=(9, 5))
    labels = top["task"].astype(str) + "\n" + top["model"].astype(str)
    ax.bar(range(len(top)), top["auroc_mean"], yerr=top.get("auroc_sd"))
    ax.set_xticks(range(len(top)))
    ax.set_xticklabels(labels, rotation=60, ha="right")
    ax.set_ylabel("Mean AUROC")
    ax.set_ylim(0, 1)
    ax.set_title("Repeated Benchmark Stability")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_permutation(pvalues: pd.DataFrame, path: Path) -> None:
    if pvalues.empty or "empirical_pvalue" not in pvalues:
        placeholder(path, "Permutation Controls", "Permutation-control p-values are unavailable.")
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(pvalues["task"].astype(str), pvalues["empirical_pvalue"].astype(float))
    ax.axhline(0.05, color="black", linestyle="--", linewidth=1)
    ax.set_ylabel("Empirical p-value")
    ax.set_title("Permutation Controls")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_ablation(ablation: pd.DataFrame, path: Path) -> None:
    if ablation.empty or "delta_auroc_when_removed" not in ablation:
        placeholder(path, "Feature Ablation", "Feature-ablation results are unavailable.")
        return
    top = ablation.sort_values("delta_auroc_when_removed", ascending=False).head(20)
    fig, ax = plt.subplots(figsize=(9, 5))
    labels = top["task"].astype(str) + "\n" + top["feature_group"].astype(str)
    ax.bar(range(len(top)), top["delta_auroc_when_removed"].astype(float))
    ax.set_xticks(range(len(top)))
    ax.set_xticklabels(labels, rotation=60, ha="right")
    ax.set_ylabel("Delta AUROC when removed")
    ax.set_title("Feature Ablation")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_claim_strength(evidence: pd.DataFrame, path: Path) -> None:
    if evidence.empty or "evidence_category" not in evidence:
        placeholder(path, "Claim Strength Matrix", "Evidence strength matrix is unavailable.")
        return
    categories = evidence["evidence_category"].value_counts()
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(categories.index.astype(str), categories.values)
    ax.set_ylabel("Tasks")
    ax.set_title("Claim Strength Matrix")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 12 benchmark figures.")
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.figures_dir.mkdir(parents=True, exist_ok=True)
    plot_repeated(
        table_or_empty(Path("results/tables/phase12_repeated_benchmark_summary.tsv")),
        args.figures_dir / "figure27_repeated_benchmark_stability.png",
    )
    plot_permutation(
        table_or_empty(Path("results/tables/phase12_empirical_pvalues.tsv")),
        args.figures_dir / "figure28_permutation_controls.png",
    )
    plot_ablation(
        table_or_empty(Path("results/tables/phase12_feature_group_importance.tsv")),
        args.figures_dir / "figure29_feature_ablation.png",
    )
    plot_claim_strength(
        table_or_empty(Path("results/reports/evidence_strength_matrix.tsv")),
        args.figures_dir / "figure30_claim_strength_matrix.png",
    )
    print(f"Wrote Phase 12 figures under {args.figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
