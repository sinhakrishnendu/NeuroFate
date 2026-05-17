#!/usr/bin/env python3
"""Run Mathys 2019 external validation from donor/sample-level feature tables."""

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
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


RANDOM_SEED = 91
FEATURE_PREFIXES = (
    "gene_mean__",
    "gene_detection__",
    "index__",
    "cell_fraction__",
    "celltype_index__",
)


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
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %s rows: %d", path, len(rows))


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def to_text(value: float) -> str:
    if isinstance(value, float) and math.isnan(value):
        return "nan"
    return f"{value:.8g}"


def target_label(row: dict[str, str]) -> int | None:
    for field in ["label__Overall_AD_neuropathological_Change", "label__diagnosis"]:
        value = (row.get(field) or "").strip().lower()
        if "high" in value or value == "ad":
            return 1
        if "none" in value or "low" in value or "control" in value or value == "ct":
            return 0
    return None


def feature_columns(*tables: list[dict[str, str]]) -> list[str]:
    feature_sets = []
    for rows in tables:
        if rows:
            feature_sets.append({field for field in rows[0] if field.startswith(FEATURE_PREFIXES)})
    if not feature_sets:
        return []
    return sorted(set.intersection(*feature_sets))


def filtered(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], np.ndarray]:
    kept: list[dict[str, str]] = []
    labels: list[int] = []
    for row in rows:
        label = target_label(row)
        if label is None:
            continue
        kept.append(row)
        labels.append(label)
    return kept, np.asarray(labels, dtype=int)


def matrix(rows: list[dict[str, str]], features: list[str]) -> np.ndarray:
    return np.asarray([[to_float(row.get(feature, "0")) for feature in features] for row in rows], dtype=float)


def model_pipeline() -> Pipeline:
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(max_iter=1000, class_weight="balanced", random_state=RANDOM_SEED)),
        ]
    )


def can_evaluate(y_train: np.ndarray, y_test: np.ndarray) -> tuple[bool, str]:
    if len(y_train) < 4 or len(y_test) < 2:
        return False, "insufficient sample units"
    if len(Counter(y_train.tolist())) < 2:
        return False, "training table has fewer than two classes"
    if len(Counter(y_test.tolist())) < 2:
        return False, "test table has fewer than two classes"
    return True, "ok"


def evaluate(
    mode: str,
    model: Pipeline,
    X_test: np.ndarray,
    y_test: np.ndarray,
    test_rows: list[dict[str, str]],
    n_train: int,
) -> tuple[dict[str, str], list[dict[str, str]]]:
    probabilities = model.predict_proba(X_test)[:, 1]
    predictions = (probabilities >= 0.5).astype(int)
    metric_row = {
        "validation_mode": mode,
        "task_id": "high_vs_low_ad_or_diagnosis",
        "n_train": str(n_train),
        "n_test": str(len(y_test)),
        "auroc": to_text(float(roc_auc_score(y_test, probabilities))),
        "auprc": to_text(float(average_precision_score(y_test, probabilities))),
        "balanced_accuracy": to_text(float(balanced_accuracy_score(y_test, predictions))),
        "brier_score": to_text(float(brier_score_loss(y_test, probabilities))),
        "notes": "logistic regression using shared donor-level features",
    }
    prediction_rows = [
        {
            "validation_mode": mode,
            "donor_id": row.get("donor_id", ""),
            "true_label": str(int(label)),
            "predicted_probability": to_text(float(probability)),
            "predicted_label": str(int(prediction)),
        }
        for row, label, probability, prediction in zip(test_rows, y_test, probabilities, predictions, strict=False)
    ]
    return metric_row, prediction_rows


def skipped(mode: str, reason: str) -> dict[str, str]:
    return {
        "validation_mode": mode,
        "task_id": "high_vs_low_ad_or_diagnosis",
        "n_train": "0",
        "n_test": "0",
        "auroc": "nan",
        "auprc": "nan",
        "balanced_accuracy": "nan",
        "brier_score": "nan",
        "notes": reason,
    }


def run_external(sea_rows: list[dict[str, str]], mathys_rows: list[dict[str, str]], features: list[str]) -> tuple[dict[str, str], list[dict[str, str]]]:
    sea, y_sea = filtered(sea_rows)
    mathys, y_mathys = filtered(mathys_rows)
    ok, reason = can_evaluate(y_sea, y_mathys)
    if not ok:
        return skipped("train_sea_ad_test_mathys", reason), []
    model = model_pipeline()
    model.fit(matrix(sea, features), y_sea)
    return evaluate("train_sea_ad_test_mathys", model, matrix(mathys, features), y_mathys, mathys, len(y_sea))


def run_internal(mathys_rows: list[dict[str, str]], features: list[str]) -> tuple[dict[str, str], list[dict[str, str]]]:
    mathys, y = filtered(mathys_rows)
    counts = Counter(y.tolist())
    if len(mathys) < 8 or len(counts) < 2 or min(counts.values()) < 2:
        return skipped("mathys_internal_diagnostic", "insufficient Mathys sample units"), []
    indices = np.arange(len(mathys))
    train_idx, test_idx = train_test_split(indices, test_size=0.30, random_state=RANDOM_SEED, stratify=y)
    X = matrix(mathys, features)
    model = model_pipeline()
    model.fit(X[train_idx], y[train_idx])
    test_rows = [mathys[int(index)] for index in test_idx]
    return evaluate("mathys_internal_diagnostic", model, X[test_idx], y[test_idx], test_rows, len(train_idx))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Mathys CSV external validation.")
    parser.add_argument("--sea-ad-features", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--mathys-features", type=Path, default=Path("results/tables/mathys_2019_phase5_donor_feature_table.tsv"))
    parser.add_argument("--metrics-output", type=Path, default=Path("results/tables/phase9_mathys_external_validation_metrics.tsv"))
    parser.add_argument("--predictions-output", type=Path, default=Path("results/tables/phase9_mathys_external_predictions.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/43_run_mathys_external_validation.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    sea_rows = read_tsv(args.sea_ad_features)
    mathys_rows = read_tsv(args.mathys_features)
    features = feature_columns(sea_rows, mathys_rows)
    if not features:
        raise RuntimeError("No shared SEA-AD/Mathys feature columns found.")
    logging.info("Shared feature columns: %d", len(features))
    metric_rows: list[dict[str, str]] = []
    prediction_rows: list[dict[str, str]] = []
    for metric, predictions in [run_external(sea_rows, mathys_rows, features), run_internal(mathys_rows, features)]:
        metric_rows.append(metric)
        prediction_rows.extend(predictions)
    write_tsv(
        args.metrics_output,
        metric_rows,
        ["validation_mode", "task_id", "n_train", "n_test", "auroc", "auprc", "balanced_accuracy", "brier_score", "notes"],
    )
    write_tsv(
        args.predictions_output,
        prediction_rows,
        ["validation_mode", "donor_id", "true_label", "predicted_probability", "predicted_label"],
    )
    logging.info("Mathys external validation complete using donor/sample-level tables only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
