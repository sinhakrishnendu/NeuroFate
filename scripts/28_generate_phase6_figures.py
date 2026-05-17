#!/usr/bin/env python3
"""Generate Phase 6 MPS neural-model figures from Phase 6 TSV outputs."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
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
    plt.figure(figsize=(9.5, 5.5))
    positions = range(len(labels))
    plt.bar(positions, values, color="#4D6A8A")
    plt.xticks(positions, labels, rotation=45, ha="right", fontsize=8)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def figure13_model_performance(tables_dir: Path, figures_dir: Path) -> None:
    rows = [
        row for row in read_tsv(tables_dir / "phase6_mps_model_metrics.tsv")
        if row.get("model_name") != "not_run"
    ]
    rows.sort(key=lambda row: row["task_id"])
    save_bar(
        [row["task_id"] for row in rows],
        [to_float(row["auroc"]) for row in rows],
        "NeuroFate MPS MLP Performance",
        "AUROC",
        figures_dir / "figure13_mps_model_performance.png",
    )


def figure14_training_curves(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "phase6_mps_training_log.tsv")
    by_task: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_task[row["task_id"]].append(row)

    figures_dir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(9.5, 5.5))
    for task_id, task_rows in sorted(by_task.items()):
        task_rows.sort(key=lambda row: int(row["epoch"]))
        epochs = [int(row["epoch"]) for row in task_rows]
        values = [to_float(row["validation_loss"]) for row in task_rows]
        plt.plot(epochs, values, label=task_id)
    plt.xlabel("Epoch")
    plt.ylabel("Validation loss")
    plt.title("NeuroFate MPS Training Curves")
    if by_task:
        plt.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(figures_dir / "figure14_mps_training_curves.png", dpi=300)
    plt.close()


def figure15_prediction_distribution(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "phase6_mps_predictions.tsv")
    values = [to_float(row["predicted_probability"]) for row in rows]
    figures_dir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8.5, 5.5))
    plt.hist(values, bins=min(20, max(5, len(values) // 3 if values else 5)), color="#8A5A44")
    plt.xlabel("Predicted disease probability")
    plt.ylabel("Donor count")
    plt.title("NeuroFate MPS Prediction Distribution")
    plt.tight_layout()
    plt.savefig(figures_dir / "figure15_mps_prediction_distribution.png", dpi=300)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 6 MPS model figures.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    figure13_model_performance(args.tables_dir, args.figures_dir)
    figure14_training_curves(args.tables_dir, args.figures_dir)
    figure15_prediction_distribution(args.tables_dir, args.figures_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
