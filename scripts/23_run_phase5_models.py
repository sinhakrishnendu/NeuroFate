#!/usr/bin/env python3
"""Run Phase 5 donor-level predictive modeling from the donor feature table.

Models are lightweight scikit-learn baselines trained at donor level only.
Labels are derived from donor metadata columns and are excluded from features.
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


RANDOM_SEED = 17
TEST_SIZE = 0.30

LABEL_PREFIX = "label__"
FEATURE_PREFIXES = (
    "gene_mean__",
    "gene_detection__",
    "index__",
    "cell_fraction__",
    "celltype_index__",
)

TASK_LABELS = {
    "dementia_vs_reference": "Dementia vs Reference",
    "high_vs_low_ad_neuropathology": "High vs Low AD neuropathology",
    "apoe_risk_prediction": "APOE-risk prediction",
    "mixed_pathology_burden": "Mixed-pathology burden prediction",
}


@dataclass(frozen=True)
class ModelSpec:
    name: str
    builder: Callable[[], object]


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


def normalize_label(value: str | None) -> str:
    return (value or "").strip().lower()


def is_missing(label: str) -> bool:
    lowered = normalize_label(label)
    return lowered in {"", "missing", "nan", "none", "unknown", "not available"}


def target_dementia_vs_reference(row: dict[str, str]) -> int | None:
    cognitive = normalize_label(row.get("label__Cognitive_Status"))
    reference = normalize_label(row.get("label__Neurotypical_reference"))
    if "dementia" in cognitive and "no dementia" not in cognitive and "without dementia" not in cognitive:
        return 1
    if "reference" in reference and "false" not in reference and "no" not in reference:
        return 0
    if any(term in cognitive for term in ["no dementia", "cognitively normal", "control", "reference"]):
        return 0
    return None


def target_high_vs_low_ad(row: dict[str, str]) -> int | None:
    pathology = normalize_label(row.get("label__Overall_AD_neuropathological_Change"))
    if "high" in pathology:
        return 1
    if any(term in pathology for term in ["not", "none", "low"]):
        return 0
    return None


def target_apoe_risk(row: dict[str, str]) -> int | None:
    genotype = normalize_label(row.get("label__APOE_Genotype"))
    if is_missing(genotype):
        return None
    if "4" in genotype or "e4" in genotype or "epsilon4" in genotype:
        return 1
    if "2" in genotype or "3" in genotype or "e2" in genotype or "e3" in genotype:
        return 0
    return None


def burden_positive(label: str) -> bool | None:
    lowered = normalize_label(label)
    if is_missing(lowered):
        return None
    negative_terms = ["none", "not", "absent", "no ", "low", "0"]
    positive_terms = ["high", "intermediate", "moderate", "severe", "limbic", "neocortical", "present", "stage"]
    if any(term in lowered for term in positive_terms):
        return True
    if any(term in lowered for term in negative_terms):
        return False
    return None


def target_mixed_pathology(row: dict[str, str]) -> int | None:
    labels = [
        row.get("label__Highest_Lewy_Body_Disease", ""),
        row.get("label__LATE", ""),
        row.get("label__Overall_CAA_Score", ""),
    ]
    calls = [burden_positive(label) for label in labels]
    if any(call is True for call in calls):
        return 1
    if calls and all(call is False for call in calls):
        return 0
    return None


TARGET_BUILDERS: dict[str, Callable[[dict[str, str]], int | None]] = {
    "dementia_vs_reference": target_dementia_vs_reference,
    "high_vs_low_ad_neuropathology": target_high_vs_low_ad,
    "apoe_risk_prediction": target_apoe_risk,
    "mixed_pathology_burden": target_mixed_pathology,
}


def model_specs() -> list[ModelSpec]:
    return [
        ModelSpec(
            "logistic_regression",
            lambda: LogisticRegression(
                max_iter=1000,
                class_weight="balanced",
                random_state=RANDOM_SEED,
            ),
        ),
        ModelSpec(
            "elastic_net",
            lambda: LogisticRegression(
                penalty="elasticnet",
                solver="saga",
                l1_ratio=0.50,
                max_iter=2500,
                class_weight="balanced",
                random_state=RANDOM_SEED,
            ),
        ),
        ModelSpec(
            "random_forest_baseline",
            lambda: RandomForestClassifier(
                n_estimators=200,
                max_depth=4,
                min_samples_leaf=2,
                class_weight="balanced",
                random_state=RANDOM_SEED,
                n_jobs=1,
            ),
        ),
        ModelSpec(
            "gradient_boosting_baseline",
            lambda: GradientBoostingClassifier(random_state=RANDOM_SEED),
        ),
    ]


def build_pipeline(model: object) -> Pipeline:
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", model),
        ]
    )


def feature_columns(fieldnames: list[str]) -> list[str]:
    return [
        field
        for field in fieldnames
        if field.startswith(FEATURE_PREFIXES)
    ]


def build_feature_matrix(
    rows: list[dict[str, str]],
    selected_features: list[str],
) -> np.ndarray:
    return np.asarray(
        [[to_float(row.get(feature, "0")) for feature in selected_features] for row in rows],
        dtype=float,
    )


def build_task_dataset(
    rows: list[dict[str, str]],
    selected_features: list[str],
    task_id: str,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    builder = TARGET_BUILDERS[task_id]
    donor_ids: list[str] = []
    task_rows: list[dict[str, str]] = []
    labels: list[int] = []
    for row in rows:
        label = builder(row)
        if label is None:
            continue
        donor_ids.append(row["donor_id"])
        task_rows.append(row)
        labels.append(label)
    return donor_ids, build_feature_matrix(task_rows, selected_features), np.asarray(labels, dtype=int)


def can_model(y: np.ndarray) -> tuple[bool, str]:
    counts = Counter(y.tolist())
    if len(counts) < 2:
        return False, "skipped: fewer than two target classes"
    if min(counts.values()) < 2:
        return False, "skipped: at least one class has fewer than two donors"
    if len(y) < 6:
        return False, "skipped: fewer than six labeled donors"
    return True, "ok"


def train_test_indices(y: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray]:
    indices = np.arange(len(y))
    stratify = y if min(Counter(y.tolist()).values()) >= 2 else None
    train_idx, test_idx = train_test_split(
        indices,
        test_size=TEST_SIZE,
        random_state=seed,
        stratify=stratify,
    )
    return np.asarray(train_idx, dtype=int), np.asarray(test_idx, dtype=int)


def safe_metric(metric: Callable[..., float], *args: object) -> float:
    try:
        return float(metric(*args))
    except ValueError:
        return float("nan")


def evaluate_model(
    task_id: str,
    model_name: str,
    pipeline: Pipeline,
    donor_ids: list[str],
    X: np.ndarray,
    y: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
) -> dict[str, str]:
    X_train = X[train_idx]
    X_test = X[test_idx]
    y_train = y[train_idx]
    y_test = y[test_idx]
    pipeline.fit(X_train, y_train)
    probabilities = pipeline.predict_proba(X_test)[:, 1]
    predictions = (probabilities >= 0.5).astype(int)
    positive_rate = float(np.mean(y))
    auroc = safe_metric(roc_auc_score, y_test, probabilities) if len(set(y_test.tolist())) == 2 else float("nan")
    auprc = safe_metric(average_precision_score, y_test, probabilities) if len(set(y_test.tolist())) == 2 else float("nan")
    balanced_accuracy = safe_metric(balanced_accuracy_score, y_test, predictions)
    brier = safe_metric(brier_score_loss, y_test, probabilities)
    return {
        "task_id": task_id,
        "target_label": TASK_LABELS[task_id],
        "model_name": model_name,
        "n_donors": str(len(donor_ids)),
        "n_train": str(len(train_idx)),
        "n_test": str(len(test_idx)),
        "positive_rate": to_text(positive_rate),
        "auroc": to_text(auroc),
        "auprc": to_text(auprc),
        "balanced_accuracy": to_text(balanced_accuracy),
        "brier_score": to_text(brier),
        "calibration_metric": "brier_score_lower_is_better",
        "random_seed": str(RANDOM_SEED),
        "notes": "donor-level train/test split; labels excluded from features",
    }


def feature_importance_rows(
    task_id: str,
    model_name: str,
    pipeline: Pipeline,
    selected_features: list[str],
) -> list[dict[str, str]]:
    model = pipeline.named_steps["model"]
    if hasattr(model, "coef_"):
        values = np.ravel(model.coef_)
    elif hasattr(model, "feature_importances_"):
        values = np.ravel(model.feature_importances_)
    else:
        return []
    rows: list[dict[str, str]] = []
    for feature, importance in zip(selected_features, values, strict=False):
        rows.append(
            {
                "task_id": task_id,
                "model_name": model_name,
                "feature": feature,
                "importance": to_text(float(importance)),
                "absolute_importance": to_text(abs(float(importance))),
                "direction": "positive" if float(importance) >= 0 else "negative",
            }
        )
    rows.sort(key=lambda row: float(row["absolute_importance"]), reverse=True)
    return rows[:50]


def out_of_fold_probabilities(
    X: np.ndarray,
    y: np.ndarray,
    seed: int,
) -> np.ndarray:
    probabilities = np.full(len(y), np.nan, dtype=float)
    counts = Counter(y.tolist())
    if len(counts) < 2 or min(counts.values()) < 2:
        return probabilities
    n_splits = min(5, min(counts.values()))
    splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for train_idx, holdout_idx in splitter.split(X, y):
        pipeline = build_pipeline(
            LogisticRegression(
                max_iter=1000,
                class_weight="balanced",
                random_state=seed,
            )
        )
        pipeline.fit(X[train_idx], y[train_idx])
        probabilities[holdout_idx] = pipeline.predict_proba(X[holdout_idx])[:, 1]
    return probabilities


def run_models(
    rows: list[dict[str, str]],
    selected_features: list[str],
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    metrics_rows: list[dict[str, str]] = []
    importance_rows: list[dict[str, str]] = []
    score_lookup: dict[str, dict[str, float]] = defaultdict(dict)

    for task_id in TARGET_BUILDERS:
        donor_ids, X, y = build_task_dataset(rows, selected_features, task_id)
        ok, reason = can_model(y)
        logging.info("%s labeled donors: %d (%s)", task_id, len(y), reason)
        if not ok:
            metrics_rows.append(
                {
                    "task_id": task_id,
                    "target_label": TASK_LABELS[task_id],
                    "model_name": "not_run",
                    "n_donors": str(len(y)),
                    "n_train": "0",
                    "n_test": "0",
                    "positive_rate": to_text(float(np.mean(y)) if len(y) else float("nan")),
                    "auroc": "nan",
                    "auprc": "nan",
                    "balanced_accuracy": "nan",
                    "brier_score": "nan",
                    "calibration_metric": "not_available",
                    "random_seed": str(RANDOM_SEED),
                    "notes": reason,
                }
            )
            continue

        train_idx, test_idx = train_test_indices(y, RANDOM_SEED)
        for spec in model_specs():
            pipeline = build_pipeline(spec.builder())
            metrics_rows.append(
                evaluate_model(task_id, spec.name, pipeline, donor_ids, X, y, train_idx, test_idx)
            )
            importance_rows.extend(feature_importance_rows(task_id, spec.name, pipeline, selected_features))

        oof = out_of_fold_probabilities(X, y, RANDOM_SEED)
        for donor_id, probability in zip(donor_ids, oof, strict=False):
            score_lookup[donor_id][task_id] = float(probability)

    score_rows = build_score_rows(rows, score_lookup)
    return metrics_rows, importance_rows, score_rows


def build_score_rows(
    rows: list[dict[str, str]],
    score_lookup: dict[str, dict[str, float]],
) -> list[dict[str, str]]:
    score_rows: list[dict[str, str]] = []
    task_columns = [f"{task_id}_oof_probability" for task_id in TARGET_BUILDERS]
    for row in rows:
        donor_id = row["donor_id"]
        scores = score_lookup.get(donor_id, {})
        numeric_scores = [
            value
            for value in scores.values()
            if isinstance(value, float) and not math.isnan(value)
        ]
        neurofate_score = sum(numeric_scores) / len(numeric_scores) if numeric_scores else float("nan")
        out_row = {
            "donor_id": donor_id,
            "neurofate_neurodegeneration_risk_score": to_text(neurofate_score),
            "available_task_scores": str(len(numeric_scores)),
            "score_definition": "mean out-of-fold probability across available donor-level Phase 5 tasks",
        }
        for task_id, column in zip(TARGET_BUILDERS, task_columns, strict=False):
            out_row[column] = to_text(scores.get(task_id, float("nan")))
        score_rows.append(out_row)
    return score_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 5 donor-level interpretable models.")
    parser.add_argument(
        "--features",
        type=Path,
        default=Path("results/tables/phase5_donor_feature_table.tsv"),
    )
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/23_run_phase5_models.log"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    logging.info("Starting Phase 5 donor-level predictive modeling.")
    rows = read_tsv(args.features)
    if not rows:
        raise RuntimeError(f"No donor feature rows found: {args.features}")
    selected_features = feature_columns(list(rows[0].keys()))
    if not selected_features:
        raise RuntimeError("No Phase 5 feature columns found.")
    logging.info("Donor rows: %d", len(rows))
    logging.info("Feature columns used: %d", len(selected_features))

    metrics_rows, importance_rows, score_rows = run_models(rows, selected_features)
    write_tsv(
        args.tables_dir / "phase5_model_metrics.tsv",
        metrics_rows,
        [
            "task_id",
            "target_label",
            "model_name",
            "n_donors",
            "n_train",
            "n_test",
            "positive_rate",
            "auroc",
            "auprc",
            "balanced_accuracy",
            "brier_score",
            "calibration_metric",
            "random_seed",
            "notes",
        ],
    )
    write_tsv(
        args.tables_dir / "phase5_feature_importance.tsv",
        importance_rows,
        ["task_id", "model_name", "feature", "importance", "absolute_importance", "direction"],
    )
    write_tsv(
        args.tables_dir / "phase5_neurofate_scores.tsv",
        score_rows,
        [
            "donor_id",
            "neurofate_neurodegeneration_risk_score",
            "available_task_scores",
            "dementia_vs_reference_oof_probability",
            "high_vs_low_ad_neuropathology_oof_probability",
            "apoe_risk_prediction_oof_probability",
            "mixed_pathology_burden_oof_probability",
            "score_definition",
        ],
    )
    logging.info("Phase 5 modeling complete.")
    logging.info("No H5AD file, single-cell workflow, or deep model was used.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
