#!/usr/bin/env python3
"""Draft Phase 6 MPS neural-model results text from Phase 6 TSV outputs."""

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
    return sorted(
        [row for row in rows if row.get("model_name") != "not_run"],
        key=lambda row: to_float(row.get("auroc", "0")),
        reverse=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Draft Phase 6 MPS neural-model results text.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/tables/phase6_results_summary.txt"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metrics = read_tsv(args.tables_dir / "phase6_mps_model_metrics.tsv")
    predictions = read_tsv(args.tables_dir / "phase6_mps_predictions.tsv")
    training_log = read_tsv(args.tables_dir / "phase6_mps_training_log.tsv")

    lines = [
        "Phase 6 Apple Silicon Metal NeuroFate Neural Model Summary",
        "",
        "This text is generated from Phase 6 donor-level MPS model outputs and should be reviewed before manuscript use.",
        "",
        "Model performance by task:",
    ]
    for row in best_metrics(metrics):
        lines.append(
            "- {task}: device={device}, best_epoch={epoch}, AUROC={auroc}, AUPRC={auprc}, balanced_accuracy={bal}, Brier={brier}, model={model_path}".format(
                task=row.get("task_id", ""),
                device=row.get("selected_device", ""),
                epoch=row.get("best_epoch", ""),
                auroc=row.get("auroc", ""),
                auprc=row.get("auprc", ""),
                bal=row.get("balanced_accuracy", ""),
                brier=row.get("brier_score", ""),
                model_path=row.get("model_path", ""),
            )
        )
    skipped = [row for row in metrics if row.get("model_name") == "not_run"]
    if skipped:
        lines.extend(["", "Tasks skipped because donor labels were insufficient:"])
        for row in skipped:
            lines.append(f"- {row.get('task_id', '')}: {row.get('notes', '')}")

    lines.extend(
        [
            "",
            f"Prediction rows: {len(predictions)}",
            f"Training-log rows: {len(training_log)}",
            "",
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
