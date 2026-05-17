#!/usr/bin/env python3
"""Draft Phase 7 cross-cohort validation results text."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Draft Phase 7 cross-cohort validation summary text.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/tables/phase7_results_summary.txt"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metrics = read_tsv(args.tables_dir / "phase7_crosscohort_metrics.tsv")
    summary = read_tsv(args.tables_dir / "phase7_generalization_summary.tsv")
    feature_overlap = read_tsv(args.tables_dir / "crosscohort_feature_overlap.tsv")

    lines = [
        "Phase 7 Cross-Cohort External Validation Summary",
        "",
        "This text is generated from donor-level cross-cohort validation TSVs and should be reviewed before manuscript use.",
        "",
        "Generalization summary:",
    ]
    for row in summary:
        lines.append(
            "- {mode}: completed={n}, mean_AUROC={auroc}, mean_AUPRC={auprc}".format(
                mode=row.get("validation_mode", ""),
                n=row.get("n_completed_comparisons", ""),
                auroc=row.get("mean_auroc", ""),
                auprc=row.get("mean_auprc", ""),
            )
        )
    lines.extend(["", "Cohort-specific validation rows:"])
    for row in metrics[:20]:
        lines.append(
            "- {mode}: train={train}, test={test}, AUROC={auroc}, AUPRC={auprc}, note={note}".format(
                mode=row.get("validation_mode", ""),
                train=row.get("train_cohort", ""),
                test=row.get("test_cohort", ""),
                auroc=row.get("auroc", ""),
                auprc=row.get("auprc", ""),
                note=row.get("notes", ""),
            )
        )
    shared = sum(1 for row in feature_overlap if row.get("status") == "shared")
    lines.extend(["", f"Shared feature count across available cohorts: {shared}", ""])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
