#!/usr/bin/env python3
"""Run lightweight multi-external validation on donor/sample feature tables only."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


RELIABILITY_CATEGORIES = [
    "reliable_external_validation",
    "preliminary_external_feasibility",
    "insufficient_sample_size",
    "insufficient_feature_overlap",
    "failed_label_mapping",
    "moderate_pd_internal_validation",
    "preliminary_cross_disease_feature_transfer",
    "insufficient_cross_disease_validation",
]
MIN_RELIABLE_N = 20
MIN_FEATURE_OVERLAP = 10


def setup_logging(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=path, level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def feature_columns(rows: list[dict[str, str]]) -> set[str]:
    if not rows:
        return set()
    return {col for col in rows[0] if col.startswith(("gene_mean__", "gene_detection__", "index__", "cell_fraction__", "celltype_index__"))}


def label_counts(rows: list[dict[str, str]]) -> tuple[int, int]:
    positives = 0
    negatives = 0
    for row in rows:
        label = (row.get("label__diagnosis") or row.get("label__Cognitive_Status") or "").lower()
        if any(token in label for token in ["pd", "ad", "dementia", "case", "disease", "high"]):
            positives += 1
        elif any(token in label for token in ["control", "reference", "none", "low", "no dementia"]):
            negatives += 1
    return positives, negatives


def reliability(n_test: int, positives: int, negatives: int, overlap: int, dataset_id: str = "") -> str:
    if positives == 0 or negatives == 0:
        return "failed_label_mapping"
    if "gse243639" in dataset_id or "_pd_" in dataset_id:
        if n_test >= MIN_RELIABLE_N and overlap >= MIN_FEATURE_OVERLAP:
            return "preliminary_cross_disease_feature_transfer"
        return "insufficient_cross_disease_validation"
    if n_test < MIN_RELIABLE_N:
        return "preliminary_external_feasibility" if n_test >= 6 else "insufficient_sample_size"
    if overlap < MIN_FEATURE_OVERLAP:
        return "insufficient_feature_overlap"
    return "reliable_external_validation"


def parse_external_arg(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Use dataset_id=path for --external-feature-table.")
    dataset_id, path = value.split("=", 1)
    return dataset_id, Path(path)


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run donor-level multi-external validation.")
    parser.add_argument("--sea-ad-features", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--external-feature-table", action="append", type=parse_external_arg, default=[])
    parser.add_argument("--pd-feature-table", type=Path, default=Path("results/tables/phase16_gse243639_feature_table.tsv"))
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--reports-dir", type=Path, default=Path("results/reports"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/75_run_multi_external_validation.log"))
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log_file)
    external_feature_tables = list(args.external_feature_table)
    if not external_feature_tables and args.pd_feature_table.exists():
        external_feature_tables.append(("gse243639_pd_snpc", args.pd_feature_table))
    if args.dry_run or not external_feature_tables:
        print("Manual multi-external validation plan only. Provide --external-feature-table dataset_id=path to run lightweight donor-level validation.")
        return 0
    sea_rows = read_tsv(args.sea_ad_features)
    sea_features = feature_columns(sea_rows)
    metric_rows: list[dict[str, str]] = []
    prediction_rows: list[dict[str, str]] = []
    reliability_rows: list[dict[str, str]] = []
    for dataset_id, path in external_feature_tables:
        external_rows = read_tsv(path)
        external_features = feature_columns(external_rows)
        overlap = len(sea_features & external_features)
        positives, negatives = label_counts(external_rows)
        flag = reliability(len(external_rows), positives, negatives, overlap, dataset_id)
        metric_rows.append(
            {
                "validation_mode": "train_sea_ad_test_external",
                "dataset_id": dataset_id,
                "task_id": "external_case_control_or_disease_status",
                "n_train": str(len(sea_rows)),
                "n_test": str(len(external_rows)),
                "auroc": "not_run_model_placeholder",
                "auprc": "not_run_model_placeholder",
                "balanced_accuracy": "not_run_model_placeholder",
                "brier_score": "not_run_model_placeholder",
                "calibration_intercept": "not_feasible",
                "calibration_slope": "not_feasible",
                "positive_class_count": str(positives),
                "negative_class_count": str(negatives),
                "feature_overlap_count": str(overlap),
                "reliability_flag": flag,
            }
        )
        reliability_rows.append(
            {
                "dataset_id": dataset_id,
                "n_test": str(len(external_rows)),
                "positive_class_count": str(positives),
                "negative_class_count": str(negatives),
                "feature_overlap_count": str(overlap),
                "reliability_category": flag,
                "note": "Model execution is intentionally lightweight and must be reviewed before manuscript claims.",
            }
        )
        for row in external_rows:
            prediction_rows.append({"dataset_id": dataset_id, "unit_id": row.get("dataset_unit_id", ""), "prediction": "not_run_model_placeholder"})
    write_tsv(args.tables_dir / "phase15_multi_external_validation_metrics.tsv", metric_rows, list(metric_rows[0]))
    write_tsv(args.tables_dir / "phase15_multi_external_predictions.tsv", prediction_rows, ["dataset_id", "unit_id", "prediction"])
    write_tsv(args.reports_dir / "phase15_external_validation_reliability.tsv", reliability_rows, list(reliability_rows[0]))
    print(f"Wrote {args.tables_dir / 'phase15_multi_external_validation_metrics.tsv'}")
    print(f"Wrote {args.reports_dir / 'phase15_external_validation_reliability.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
