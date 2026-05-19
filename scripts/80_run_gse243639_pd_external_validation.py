#!/usr/bin/env python3
"""Run conservative donor/sample-level GSE243639 PD validation."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from collections import Counter
from pathlib import Path

import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


RANDOM_SEED = 163
FEATURE_PREFIXES = ("gene_mean__", "gene_detection__", "index__", "cell_fraction__", "celltype_index__")


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file, mode="w", encoding="utf-8"),
        ],
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %s rows: %d", path, len(rows))


def to_float(value: str | None) -> float:
    try:
        return float(value) if value not in (None, "") else 0.0
    except ValueError:
        return 0.0


def to_text(value: float) -> str:
    return "nan" if math.isnan(value) else f"{value:.8g}"


def feature_columns(*tables: list[dict[str, str]]) -> list[str]:
    feature_sets = []
    for rows in tables:
        if rows:
            feature_sets.append({field for field in rows[0] if field.startswith(FEATURE_PREFIXES)})
    if not feature_sets:
        return []
    return sorted(set.intersection(*feature_sets))


def pd_label(row: dict[str, str]) -> int | None:
    value = (row.get("label__diagnosis_binary") or "").strip()
    if value in {"0", "1"}:
        return int(value)
    diagnosis = (row.get("diagnosis") or "").lower()
    if "parkinson" in diagnosis:
        return 1
    if "control" in diagnosis:
        return 0
    return None


def filtered_pd_rows(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], np.ndarray]:
    kept = []
    labels = []
    for row in rows:
        label = pd_label(row)
        if label is None:
            continue
        kept.append(row)
        labels.append(label)
    return kept, np.asarray(labels, dtype=int)


def matrix(rows: list[dict[str, str]], features: list[str]) -> np.ndarray:
    return np.asarray([[to_float(row.get(feature)) for feature in features] for row in rows], dtype=float)


def model_pipeline() -> Pipeline:
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(max_iter=1000, class_weight="balanced", random_state=RANDOM_SEED)),
        ]
    )


def reliability_flag(n_samples: int, positives: int, negatives: int, feature_overlap: int, auroc_sd: float) -> str:
    if positives < 2 or negatives < 2:
        return "insufficient_cross_disease_validation"
    if n_samples >= 20 and positives >= 10 and negatives >= 10 and feature_overlap >= 10 and auroc_sd <= 0.15:
        return "moderate_pd_internal_validation"
    if n_samples >= 20 and feature_overlap >= 10:
        return "preliminary_cross_disease_feature_transfer"
    return "insufficient_cross_disease_validation"


def run_pd_internal_validation(
    pd_rows: list[dict[str, str]],
    features: list[str],
    repeats: int = 10,
) -> tuple[dict[str, str], list[dict[str, str]]]:
    rows, y = filtered_pd_rows(pd_rows)
    counts = Counter(y.tolist())
    if len(rows) < 8 or len(counts) < 2 or min(counts.values()) < 2 or not features:
        return skipped_metric("gse243639_pd_internal", len(rows), counts, len(features), "insufficient sample units, labels, or features"), []
    X = matrix(rows, features)
    splitter = RepeatedStratifiedKFold(n_splits=min(5, min(counts.values())), n_repeats=repeats, random_state=RANDOM_SEED)
    metrics: list[dict[str, float]] = []
    prediction_rows: list[dict[str, str]] = []
    for split_id, (train_idx, test_idx) in enumerate(splitter.split(X, y), start=1):
        model = model_pipeline()
        model.fit(X[train_idx], y[train_idx])
        probabilities = model.predict_proba(X[test_idx])[:, 1]
        predictions = (probabilities >= 0.5).astype(int)
        y_test = y[test_idx]
        metrics.append(
            {
                "auroc": float(roc_auc_score(y_test, probabilities)),
                "auprc": float(average_precision_score(y_test, probabilities)),
                "balanced_accuracy": float(balanced_accuracy_score(y_test, predictions)),
                "brier_score": float(brier_score_loss(y_test, probabilities)),
            }
        )
        for row_index, label, probability, prediction in zip(test_idx, y_test, probabilities, predictions, strict=False):
            prediction_rows.append(
                {
                    "validation_mode": "gse243639_pd_internal",
                    "split_id": str(split_id),
                    "sample_id": rows[int(row_index)].get("sample_id", ""),
                    "true_label": str(int(label)),
                    "predicted_probability": to_text(float(probability)),
                    "predicted_label": str(int(prediction)),
                }
            )
    auroc_values = [item["auroc"] for item in metrics]
    auroc_sd = float(np.std(auroc_values, ddof=1)) if len(auroc_values) > 1 else 0.0
    metric_row = {
        "validation_mode": "gse243639_pd_internal",
        "task_id": "parkinsons_vs_control",
        "auroc": to_text(float(np.mean(auroc_values))),
        "auroc_sd": to_text(auroc_sd),
        "auprc": to_text(float(np.mean([item["auprc"] for item in metrics]))),
        "balanced_accuracy": to_text(float(np.mean([item["balanced_accuracy"] for item in metrics]))),
        "brier_score": to_text(float(np.mean([item["brier_score"] for item in metrics]))),
        "n_samples": str(len(rows)),
        "positive_count": str(counts.get(1, 0)),
        "negative_count": str(counts.get(0, 0)),
        "feature_overlap_count": str(len(features)),
        "reliability_flag": reliability_flag(len(rows), counts.get(1, 0), counts.get(0, 0), len(features), auroc_sd),
        "notes": "Independent PD cohort internal validation using sample-level NeuroFate features only.",
    }
    return metric_row, prediction_rows


def skipped_metric(mode: str, n_samples: int, counts: Counter[int], feature_count: int, reason: str) -> dict[str, str]:
    return {
        "validation_mode": mode,
        "task_id": "parkinsons_vs_control",
        "auroc": "nan",
        "auroc_sd": "nan",
        "auprc": "nan",
        "balanced_accuracy": "nan",
        "brier_score": "nan",
        "n_samples": str(n_samples),
        "positive_count": str(counts.get(1, 0)),
        "negative_count": str(counts.get(0, 0)),
        "feature_overlap_count": str(feature_count),
        "reliability_flag": "insufficient_cross_disease_validation",
        "notes": reason,
    }


def cross_disease_feature_row(sea_rows: list[dict[str, str]], pd_rows: list[dict[str, str]], features: list[str]) -> dict[str, str]:
    pd_labeled, y = filtered_pd_rows(pd_rows)
    counts = Counter(y.tolist())
    return {
        "validation_mode": "cross_disease_feature_transfer_feasibility",
        "task_id": "ad_feature_space_to_pd_context",
        "auroc": "not_run_not_biologically_equivalent",
        "auroc_sd": "not_run_not_biologically_equivalent",
        "auprc": "not_run_not_biologically_equivalent",
        "balanced_accuracy": "not_run_not_biologically_equivalent",
        "brier_score": "not_run_not_biologically_equivalent",
        "n_samples": str(len(pd_labeled)),
        "positive_count": str(counts.get(1, 0)),
        "negative_count": str(counts.get(0, 0)),
        "feature_overlap_count": str(len(features)),
        "reliability_flag": "preliminary_cross_disease_feature_transfer" if sea_rows and len(pd_labeled) >= 20 and len(features) >= 10 else "insufficient_cross_disease_validation",
        "notes": "Reports structural feature overlap only; does not claim AD-trained diagnostic transfer to PD.",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run GSE243639 PD sample-level validation.")
    parser.add_argument("--sea-ad-features", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--pd-features", type=Path, default=Path("results/tables/phase16_gse243639_feature_table.tsv"))
    parser.add_argument("--metrics-output", type=Path, default=Path("results/tables/phase16_gse243639_external_validation_metrics.tsv"))
    parser.add_argument("--predictions-output", type=Path, default=Path("results/tables/phase16_gse243639_external_predictions.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/80_run_gse243639_pd_external_validation.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    sea_rows = read_tsv(args.sea_ad_features) if args.sea_ad_features.exists() else []
    pd_rows = read_tsv(args.pd_features)
    features = feature_columns(sea_rows, pd_rows) if sea_rows else feature_columns(pd_rows)
    metric, predictions = run_pd_internal_validation(pd_rows, features)
    metrics = [metric, cross_disease_feature_row(sea_rows, pd_rows, features)]
    metric_fields = [
        "validation_mode",
        "task_id",
        "auroc",
        "auroc_sd",
        "auprc",
        "balanced_accuracy",
        "brier_score",
        "n_samples",
        "positive_count",
        "negative_count",
        "feature_overlap_count",
        "reliability_flag",
        "notes",
    ]
    write_tsv(args.metrics_output, metrics, metric_fields)
    write_tsv(args.predictions_output, predictions, ["validation_mode", "split_id", "sample_id", "true_label", "predicted_probability", "predicted_label"])
    logging.info("GSE243639 validation used sample-level tables only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
