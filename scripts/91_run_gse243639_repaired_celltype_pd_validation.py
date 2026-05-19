#!/usr/bin/env python3
"""Run repaired sample-level GSE243639 cell-type-aware PD validation."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from collections import Counter
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import LeaveOneOut, RepeatedStratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


FEATURE_PREFIXES = (
    "gene_mean__",
    "gene_detection__",
    "cell_fraction__",
    "celltype_gene_mean__",
    "celltype_gene_detection__",
    "index__",
)
MIN_REPAIRED_FEATURES = 20
MIN_MATCH_RATE = 0.50


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
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


def feature_columns(rows: list[dict[str, str]]) -> list[str]:
    if not rows:
        return []
    return sorted(field for field in rows[0] if field.startswith(FEATURE_PREFIXES))


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


def filtered_table(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], np.ndarray]:
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


def logistic_model(seed: int) -> Pipeline:
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(max_iter=1000, class_weight="balanced", random_state=seed)),
        ]
    )


def random_forest_model(seed: int) -> Pipeline:
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("model", RandomForestClassifier(n_estimators=200, max_depth=3, class_weight="balanced", random_state=seed)),
        ]
    )


def evaluate_predictions(y_true: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    predictions = (probabilities >= 0.5).astype(int)
    return {
        "auroc": float(roc_auc_score(y_true, probabilities)),
        "auprc": float(average_precision_score(y_true, probabilities)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, predictions)),
        "brier_score": float(brier_score_loss(y_true, probabilities)),
    }


def repeated_split_predictions(
    X: np.ndarray,
    y: np.ndarray,
    model_name: str,
    seed: int,
    repeats: int,
) -> tuple[list[dict[str, float]], list[dict[str, str]]]:
    counts = Counter(y.tolist())
    splitter = RepeatedStratifiedKFold(n_splits=min(5, min(counts.values())), n_repeats=repeats, random_state=seed)
    metrics = []
    predictions = []
    for split_id, (train_idx, test_idx) in enumerate(splitter.split(X, y), start=1):
        model = logistic_model(seed + split_id) if model_name == "logistic_regression" else random_forest_model(seed + split_id)
        model.fit(X[train_idx], y[train_idx])
        probabilities = model.predict_proba(X[test_idx])[:, 1]
        metrics.append(evaluate_predictions(y[test_idx], probabilities))
        for row_index, label, probability in zip(test_idx, y[test_idx], probabilities, strict=False):
            predictions.append(
                {
                    "model": model_name,
                    "validation_mode": "repeated_stratified_split",
                    "split_id": str(split_id),
                    "row_index": str(int(row_index)),
                    "true_label": str(int(label)),
                    "predicted_probability": to_text(float(probability)),
                    "predicted_label": str(int(probability >= 0.5)),
                }
            )
    return metrics, predictions


def leave_one_out_predictions(X: np.ndarray, y: np.ndarray, seed: int) -> tuple[dict[str, float], list[dict[str, str]]]:
    probabilities = np.zeros(len(y), dtype=float)
    for split_id, (train_idx, test_idx) in enumerate(LeaveOneOut().split(X), start=1):
        model = logistic_model(seed + split_id)
        model.fit(X[train_idx], y[train_idx])
        probabilities[test_idx] = model.predict_proba(X[test_idx])[:, 1]
    metrics = evaluate_predictions(y, probabilities)
    predictions = [
        {
            "model": "logistic_regression",
            "validation_mode": "leave_one_out",
            "split_id": "loo",
            "row_index": str(index),
            "true_label": str(int(label)),
            "predicted_probability": to_text(float(probability)),
            "predicted_label": str(int(probability >= 0.5)),
        }
        for index, (label, probability) in enumerate(zip(y, probabilities, strict=False))
    ]
    return metrics, predictions


def permutation_control(X: np.ndarray, y: np.ndarray, seed: int, observed_auroc: float, n_permutations: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    null_aurocs = []
    for index in range(n_permutations):
        shuffled = rng.permutation(y)
        metrics, _predictions = repeated_split_predictions(X, shuffled, "logistic_regression", seed + 1000 + index, repeats=1)
        if metrics:
            null_aurocs.append(float(np.mean([row["auroc"] for row in metrics])))
    if not null_aurocs:
        return math.nan, math.nan
    empirical_p = (sum(value >= observed_auroc for value in null_aurocs) + 1) / (len(null_aurocs) + 1)
    return float(empirical_p), float(np.mean(null_aurocs))


def bootstrap_ci(y: np.ndarray, probabilities: np.ndarray, seed: int, n_bootstrap: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values = []
    indices = np.arange(len(y))
    for _ in range(n_bootstrap):
        sampled = rng.choice(indices, size=len(indices), replace=True)
        if len(set(y[sampled].tolist())) < 2:
            continue
        values.append(float(roc_auc_score(y[sampled], probabilities[sampled])))
    if not values:
        return math.nan, math.nan
    return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def read_match_rate(path: Path) -> float:
    rows = read_tsv(path)
    if not rows:
        return 0.0
    return to_float(rows[0].get("match_rate") or rows[0].get("annotation_match_rate"))


def failure_metric(reason: str, n_samples: int, positives: int, negatives: int, feature_count: int, match_rate: float) -> list[dict[str, str]]:
    return [
        {
            "model": "not_run",
            "validation_mode": "technical_audit",
            "task_id": "parkinsons_vs_control",
            "auroc": "nan",
            "auroc_sd": "nan",
            "auroc_ci_low": "nan",
            "auroc_ci_high": "nan",
            "auprc": "nan",
            "balanced_accuracy": "nan",
            "brier_score": "nan",
            "empirical_permutation_pvalue": "nan",
            "null_auroc_mean": "nan",
            "n_samples": str(n_samples),
            "positive_count": str(positives),
            "negative_count": str(negatives),
            "feature_count": str(feature_count),
            "annotation_match_rate": f"{match_rate:.6g}",
            "reliability_flag": "technical_failure_annotation_join",
            "notes": reason,
        }
    ]


def reliability(observed_auroc: float, auroc_sd: float, empirical_p: float, match_rate: float, feature_count: int) -> str:
    if match_rate < MIN_MATCH_RATE or feature_count < MIN_REPAIRED_FEATURES:
        return "technical_failure_annotation_join"
    if observed_auroc >= 0.70 and auroc_sd <= 0.15 and empirical_p <= 0.05:
        return "moderate_pd_internal_validation"
    if observed_auroc >= 0.60 and empirical_p <= 0.20:
        return "preliminary_pd_internal_signal"
    return "weak_pd_signal"


def feature_importance(features: list[str], X: np.ndarray, y: np.ndarray, seed: int) -> list[dict[str, str]]:
    model = random_forest_model(seed)
    model.fit(X, y)
    importances = getattr(model.named_steps["model"], "feature_importances_", np.zeros(len(features)))
    ranked = sorted(zip(features, importances, strict=False), key=lambda item: float(item[1]), reverse=True)
    return [
        {"rank": str(index), "feature": feature, "importance": to_text(float(value)), "model": "random_forest_baseline"}
        for index, (feature, value) in enumerate(ranked[:50], start=1)
    ]


def summarize_metrics(
    model_name: str,
    mode: str,
    metrics: list[dict[str, float]],
    n_samples: int,
    counts: Counter[int],
    feature_count: int,
    match_rate: float,
    empirical_p: float = math.nan,
    null_auroc: float = math.nan,
    ci_low: float = math.nan,
    ci_high: float = math.nan,
) -> dict[str, str]:
    aurocs = [row["auroc"] for row in metrics]
    observed = float(np.mean(aurocs)) if aurocs else math.nan
    auroc_sd = float(np.std(aurocs, ddof=1)) if len(aurocs) > 1 else 0.0
    return {
        "model": model_name,
        "validation_mode": mode,
        "task_id": "parkinsons_vs_control",
        "auroc": to_text(observed),
        "auroc_sd": to_text(auroc_sd),
        "auroc_ci_low": to_text(ci_low),
        "auroc_ci_high": to_text(ci_high),
        "auprc": to_text(float(np.mean([row["auprc"] for row in metrics])) if metrics else math.nan),
        "balanced_accuracy": to_text(float(np.mean([row["balanced_accuracy"] for row in metrics])) if metrics else math.nan),
        "brier_score": to_text(float(np.mean([row["brier_score"] for row in metrics])) if metrics else math.nan),
        "empirical_permutation_pvalue": to_text(empirical_p),
        "null_auroc_mean": to_text(null_auroc),
        "n_samples": str(n_samples),
        "positive_count": str(counts.get(1, 0)),
        "negative_count": str(counts.get(0, 0)),
        "feature_count": str(feature_count),
        "annotation_match_rate": f"{match_rate:.6g}",
        "reliability_flag": reliability(observed, auroc_sd, empirical_p, match_rate, feature_count) if model_name == "logistic_regression" else "supporting_baseline_only",
        "notes": "Sample-level repaired cell-type-aware PD/control validation; no medical interpretation.",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run repaired GSE243639 cell-type-aware PD validation.")
    parser.add_argument("--features", type=Path, default=Path("results/tables/phase18_gse243639_celltype_feature_table.tsv"))
    parser.add_argument("--annotation-summary", type=Path, default=Path("results/tables/phase18_gse243639_annotation_match_summary.tsv"))
    parser.add_argument("--metrics-output", type=Path, default=Path("results/tables/phase18_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--predictions-output", type=Path, default=Path("results/tables/phase18_gse243639_celltype_predictions.tsv"))
    parser.add_argument("--importance-output", type=Path, default=Path("results/tables/phase18_gse243639_celltype_feature_importance.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/91_run_gse243639_repaired_celltype_pd_validation.log"))
    parser.add_argument("--n-permutations", type=int, default=100)
    parser.add_argument("--n-bootstrap", type=int, default=200)
    parser.add_argument("--n-repeats", type=int, default=10)
    parser.add_argument("--seed", type=int, default=181)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows, y = filtered_table(read_tsv(args.features))
    counts = Counter(y.tolist())
    features = feature_columns(rows)
    match_rate = read_match_rate(args.annotation_summary)
    metric_fields = [
        "model",
        "validation_mode",
        "task_id",
        "auroc",
        "auroc_sd",
        "auroc_ci_low",
        "auroc_ci_high",
        "auprc",
        "balanced_accuracy",
        "brier_score",
        "empirical_permutation_pvalue",
        "null_auroc_mean",
        "n_samples",
        "positive_count",
        "negative_count",
        "feature_count",
        "annotation_match_rate",
        "reliability_flag",
        "notes",
    ]
    if len(rows) < 8 or len(counts) < 2 or min(counts.values()) < 2 or len(features) < MIN_REPAIRED_FEATURES or match_rate < MIN_MATCH_RATE:
        reason = "annotation join or repaired feature table did not meet minimum sample, class, feature, or match-rate requirements"
        write_tsv(args.metrics_output, failure_metric(reason, len(rows), counts.get(1, 0), counts.get(0, 0), len(features), match_rate), metric_fields)
        write_tsv(args.predictions_output, [], ["model", "validation_mode", "split_id", "row_index", "true_label", "predicted_probability", "predicted_label"])
        write_tsv(args.importance_output, [], ["rank", "feature", "importance", "model"])
        return 0
    X = matrix(rows, features)
    logistic_metrics, split_predictions = repeated_split_predictions(X, y, "logistic_regression", args.seed, args.n_repeats)
    rf_metrics, rf_predictions = repeated_split_predictions(X, y, "random_forest_baseline", args.seed + 500, args.n_repeats)
    loo_metrics, loo_predictions = leave_one_out_predictions(X, y, args.seed + 900)
    observed_auroc = float(np.mean([row["auroc"] for row in logistic_metrics]))
    empirical_p, null_auroc = permutation_control(X, y, args.seed + 1200, observed_auroc, args.n_permutations)
    loo_probabilities = np.asarray([float(row["predicted_probability"]) for row in loo_predictions], dtype=float)
    ci_low, ci_high = bootstrap_ci(y, loo_probabilities, args.seed + 1400, args.n_bootstrap)
    metric_rows = [
        summarize_metrics("logistic_regression", "repeated_stratified_split", logistic_metrics, len(rows), counts, len(features), match_rate, empirical_p, null_auroc, ci_low, ci_high),
        summarize_metrics("logistic_regression", "leave_one_out", [loo_metrics], len(rows), counts, len(features), match_rate, empirical_p, null_auroc, ci_low, ci_high),
        summarize_metrics("random_forest_baseline", "repeated_stratified_split", rf_metrics, len(rows), counts, len(features), match_rate),
    ]
    write_tsv(args.metrics_output, metric_rows, metric_fields)
    write_tsv(args.predictions_output, [*split_predictions, *rf_predictions, *loo_predictions], ["model", "validation_mode", "split_id", "row_index", "true_label", "predicted_probability", "predicted_label"])
    write_tsv(args.importance_output, feature_importance(features, X, y, args.seed + 1600), ["rank", "feature", "importance", "model"])
    logging.info("Phase 18 validation completed on repaired sample-level feature table only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
