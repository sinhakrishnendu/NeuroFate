#!/usr/bin/env python3
"""Draft Phase 5 predictive modeling results text from model TSV outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


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


def best_metrics(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    usable = [row for row in rows if row.get("model_name") != "not_run"]
    by_task: dict[str, dict[str, str]] = {}
    for row in usable:
        task = row["task_id"]
        if task not in by_task or to_float(row["auroc"]) > to_float(by_task[task]["auroc"]):
            by_task[task] = row
    return [by_task[task] for task in sorted(by_task)]


def top_features(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = sorted(rows, key=lambda row: to_float(row.get("absolute_importance", "0")), reverse=True)
    return rows[:15]


def score_summary(rows: list[dict[str, str]]) -> str:
    values = [
        to_float(row["neurofate_neurodegeneration_risk_score"])
        for row in rows
        if row.get("neurofate_neurodegeneration_risk_score") != "nan"
    ]
    if not values:
        return "No out-of-fold NeuroFate scores were available; labeled donor counts may be too small."
    values.sort()
    median = values[len(values) // 2]
    return (
        f"NeuroFate scores were available for {len(values)} donors; "
        f"range={values[0]:.4g}-{values[-1]:.4g}, median={median:.4g}."
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Draft Phase 5 predictive modeling results text.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/tables/phase5_results_summary.txt"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metrics = read_tsv(args.tables_dir / "phase5_model_metrics.tsv")
    importance = read_tsv(args.tables_dir / "phase5_feature_importance.tsv")
    scores = read_tsv(args.tables_dir / "phase5_neurofate_scores.tsv")

    lines: list[str] = [
        "Phase 5 Predictive Neurodegeneration Modeling Summary",
        "",
        "This text is generated from donor-level Phase 5 model outputs and should be reviewed before manuscript use.",
        "",
        "Best model per task by AUROC:",
    ]
    for row in best_metrics(metrics):
        lines.append(
            "- {task}: {model}, AUROC={auroc}, AUPRC={auprc}, balanced_accuracy={bal}, Brier={brier}".format(
                task=row.get("task_id", ""),
                model=row.get("model_name", ""),
                auroc=row.get("auroc", ""),
                auprc=row.get("auprc", ""),
                bal=row.get("balanced_accuracy", ""),
                brier=row.get("brier_score", ""),
            )
        )
    lines.extend(["", "Top predictive features:"])
    for row in top_features(importance):
        lines.append(
            "- {feature} ({task}, {model}): importance={importance}, direction={direction}".format(
                feature=row.get("feature", ""),
                task=row.get("task_id", ""),
                model=row.get("model_name", ""),
                importance=row.get("absolute_importance", ""),
                direction=row.get("direction", ""),
            )
        )
    lines.extend(["", score_summary(scores), ""])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
