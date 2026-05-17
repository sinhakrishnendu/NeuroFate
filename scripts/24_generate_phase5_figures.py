#!/usr/bin/env python3
"""Generate Phase 5 predictive modeling figures from Phase 5 TSV outputs."""

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


def save_bar(
    labels: list[str],
    values: list[float],
    title: str,
    ylabel: str,
    out_path: Path,
    color: str = "#476A6F",
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10.5, 5.8))
    positions = range(len(labels))
    plt.bar(positions, values, color=color)
    plt.xticks(positions, labels, rotation=45, ha="right", fontsize=8)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def figure9_model_performance(tables_dir: Path, figures_dir: Path) -> None:
    rows = [
        row for row in read_tsv(tables_dir / "phase5_model_metrics.tsv")
        if row.get("model_name") != "not_run"
    ]
    rows.sort(key=lambda row: (row["task_id"], row["model_name"]))
    labels = [f"{row['task_id']}\n{row['model_name']}" for row in rows]
    values = [to_float(row["auroc"]) for row in rows]
    save_bar(
        labels,
        values,
        "Phase 5 Donor-Level Model Performance",
        "AUROC",
        figures_dir / "figure9_model_performance.png",
        "#4D6A8A",
    )


def figure10_feature_importance(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "phase5_feature_importance.tsv")
    rows.sort(key=lambda row: to_float(row["absolute_importance"]), reverse=True)
    top = rows[:20]
    labels = [row["feature"][:45] for row in top]
    values = [to_float(row["absolute_importance"]) for row in top]
    save_bar(
        labels,
        values,
        "Top NeuroFate Predictive Features",
        "Absolute importance",
        figures_dir / "figure10_feature_importance.png",
        "#6E5E4F",
    )


def figure11_score_distribution(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "phase5_neurofate_scores.tsv")
    values = [
        to_float(row["neurofate_neurodegeneration_risk_score"])
        for row in rows
        if row.get("neurofate_neurodegeneration_risk_score") != "nan"
    ]
    figures_dir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8.5, 5.5))
    plt.hist(values, bins=min(20, max(5, len(values) // 3 if values else 5)), color="#8A5A44")
    plt.xlabel("NeuroFate Neurodegeneration Risk Score")
    plt.ylabel("Donor count")
    plt.title("Donor-Level NeuroFate Score Distribution")
    plt.tight_layout()
    plt.savefig(figures_dir / "figure11_neurofate_score_distribution.png", dpi=300)
    plt.close()


def figure12_donor_risk_heatmap(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "phase5_neurofate_scores.tsv")
    rows = sorted(
        rows,
        key=lambda row: to_float(row.get("neurofate_neurodegeneration_risk_score", "0")),
        reverse=True,
    )[:60]
    score_columns = [
        "dementia_vs_reference_oof_probability",
        "high_vs_low_ad_neuropathology_oof_probability",
        "apoe_risk_prediction_oof_probability",
        "mixed_pathology_burden_oof_probability",
        "neurofate_neurodegeneration_risk_score",
    ]
    grid = [
        [to_float(row.get(column, "0")) for column in score_columns]
        for row in rows
    ]
    labels = [row["donor_id"] for row in rows]

    figures_dir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8, max(5, 0.18 * len(labels))))
    image = plt.imshow(grid, aspect="auto", cmap="magma", vmin=0, vmax=1)
    plt.colorbar(image, label="Risk probability")
    plt.xticks(range(len(score_columns)), score_columns, rotation=45, ha="right", fontsize=8)
    plt.yticks(range(len(labels)), labels, fontsize=6)
    plt.title("Donor Risk Heatmap")
    plt.tight_layout()
    plt.savefig(figures_dir / "figure12_donor_risk_heatmap.png", dpi=300)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 5 predictive modeling figures.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    figure9_model_performance(args.tables_dir, args.figures_dir)
    figure10_feature_importance(args.tables_dir, args.figures_dir)
    figure11_score_distribution(args.tables_dir, args.figures_dir)
    figure12_donor_risk_heatmap(args.tables_dir, args.figures_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
