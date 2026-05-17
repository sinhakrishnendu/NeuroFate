#!/usr/bin/env python3
"""Run Mathys robustness diagnostics without deep models or single-cell workflows."""

from __future__ import annotations

import argparse
import csv
import logging
import math
import random
from collections import Counter
from pathlib import Path

import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, brier_score_loss, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


RANDOM_SEED = 103
DEFAULT_PERMUTATIONS = 100
DEFAULT_BOOTSTRAPS = 500
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
    if not path.exists():
        return []
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
        return float("nan")


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


def feature_columns(rows: list[dict[str, str]]) -> list[str]:
    if not rows:
        return []
    return [field for field in rows[0] if field.startswith(FEATURE_PREFIXES)]


def labeled_matrix(rows: list[dict[str, str]]) -> tuple[list[str], np.ndarray, np.ndarray]:
    features = feature_columns(rows)
    donor_ids: list[str] = []
    x_rows: list[list[float]] = []
    labels: list[int] = []
    for row in rows:
        label = target_label(row)
        if label is None:
            continue
        donor_ids.append(row.get("donor_id", "missing"))
        x_rows.append([0.0 if math.isnan(to_float(row.get(feature, "0"))) else to_float(row.get(feature, "0")) for feature in features])
        labels.append(label)
    return donor_ids, np.asarray(x_rows, dtype=float), np.asarray(labels, dtype=int)


def model_pipeline() -> Pipeline:
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(max_iter=1000, class_weight="balanced", random_state=RANDOM_SEED)),
        ]
    )


def leave_one_sample_out(rows: list[dict[str, str]]) -> tuple[dict[str, str], list[float], np.ndarray]:
    donor_ids, X, y = labeled_matrix(rows)
    counts = Counter(y.tolist())
    if len(y) < 6 or len(counts) < 2:
        return {
            "check_id": "leave_one_sample_out",
            "n_units": str(len(y)),
            "auroc": "nan",
            "auprc": "nan",
            "balanced_accuracy": "nan",
            "brier_score": "nan",
            "reliability": "not_feasible",
            "notes": "Too few labeled sample units or only one class.",
        }, [], y
    probabilities = np.full(len(y), np.nan, dtype=float)
    for holdout in range(len(y)):
        train_idx = [index for index in range(len(y)) if index != holdout]
        if len(set(y[train_idx].tolist())) < 2:
            continue
        model = model_pipeline()
        model.fit(X[train_idx], y[train_idx])
        probabilities[holdout] = model.predict_proba(X[[holdout]])[:, 1][0]
    valid = ~np.isnan(probabilities)
    if valid.sum() < 2 or len(set(y[valid].tolist())) < 2:
        reliability = "unreliable_small_n"
        auroc = auprc = balanced = brier = float("nan")
    else:
        predictions = (probabilities[valid] >= 0.5).astype(int)
        auroc = float(roc_auc_score(y[valid], probabilities[valid]))
        auprc = float(average_precision_score(y[valid], probabilities[valid]))
        balanced = float(balanced_accuracy_score(y[valid], predictions))
        brier = float(brier_score_loss(y[valid], probabilities[valid]))
        reliability = "preliminary_only" if valid.sum() < 20 else "diagnostic"
    return {
        "check_id": "leave_one_sample_out",
        "n_units": str(len(y)),
        "auroc": to_text(auroc),
        "auprc": to_text(auprc),
        "balanced_accuracy": to_text(balanced),
        "brier_score": to_text(brier),
        "reliability": reliability,
        "notes": "LOSO diagnostic from Mathys sample-level table; not a definitive external claim.",
    }, probabilities.tolist(), y


def predictions_from_phase9(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = read_tsv(path)
    labels: list[int] = []
    probabilities: list[float] = []
    for row in rows:
        label = row.get("true_label", "")
        probability = to_float(row.get("predicted_probability", "nan"))
        if label in {"0", "1"} and not math.isnan(probability):
            labels.append(int(label))
            probabilities.append(probability)
    return np.asarray(labels, dtype=int), np.asarray(probabilities, dtype=float)


def permutation_test(y: np.ndarray, probabilities: np.ndarray, n_permutations: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if len(y) < 4 or len(set(y.tolist())) < 2 or len(set(probabilities.tolist())) < 2:
        return [
            {
                "check_id": "permutation_label_test",
                "n_permutations": str(n_permutations),
                "observed_auroc": "nan",
                "permutation_p_value": "nan",
                "reliability": "not_feasible",
                "notes": "Insufficient labels or probability variation.",
            }
        ]
    observed = float(roc_auc_score(y, probabilities))
    rng = random.Random(RANDOM_SEED)
    extreme = 0
    for _ in range(n_permutations):
        permuted = list(y)
        rng.shuffle(permuted)
        if len(set(permuted)) < 2:
            continue
        permuted_auroc = float(roc_auc_score(permuted, probabilities))
        if permuted_auroc >= observed:
            extreme += 1
    p_value = (extreme + 1) / (n_permutations + 1)
    rows.append(
        {
            "check_id": "permutation_label_test",
            "n_permutations": str(n_permutations),
            "observed_auroc": to_text(observed),
            "permutation_p_value": to_text(p_value),
            "reliability": "preliminary_only" if len(y) < 20 else "diagnostic",
            "notes": "Permutation diagnostic for label/probability separation.",
        }
    )
    return rows


def bootstrap_ci(y: np.ndarray, probabilities: np.ndarray, n_bootstraps: int) -> list[dict[str, str]]:
    if len(y) < 4 or len(set(y.tolist())) < 2:
        return [
            {
                "metric": "auroc",
                "n_bootstraps": str(n_bootstraps),
                "estimate": "nan",
                "ci_low": "nan",
                "ci_high": "nan",
                "reliability": "not_feasible",
                "notes": "Insufficient labels for AUROC bootstrap.",
            }
        ]
    rng = np.random.default_rng(RANDOM_SEED)
    values: list[float] = []
    for _ in range(n_bootstraps):
        indices = rng.integers(0, len(y), len(y))
        if len(set(y[indices].tolist())) < 2:
            continue
        values.append(float(roc_auc_score(y[indices], probabilities[indices])))
    if not values:
        low = high = estimate = float("nan")
        reliability = "not_feasible"
    else:
        estimate = float(roc_auc_score(y, probabilities))
        low = float(np.percentile(values, 2.5))
        high = float(np.percentile(values, 97.5))
        reliability = "wide_or_small_n" if len(y) < 20 else "diagnostic"
    return [
        {
            "metric": "auroc",
            "n_bootstraps": str(n_bootstraps),
            "estimate": to_text(estimate),
            "ci_low": to_text(low),
            "ci_high": to_text(high),
            "reliability": reliability,
            "notes": "Bootstrap confidence interval from available Mathys predictions.",
        }
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 10 Mathys robustness diagnostics.")
    parser.add_argument("--mathys-features", type=Path, default=Path("results/tables/mathys_2019_phase5_donor_feature_table.tsv"))
    parser.add_argument("--predictions", type=Path, default=Path("results/tables/phase9_mathys_external_predictions.tsv"))
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    parser.add_argument("--bootstraps", type=int, default=DEFAULT_BOOTSTRAPS)
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/47_run_mathys_robustness_checks.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = read_tsv(args.mathys_features)
    loso_row, loso_probabilities, loso_labels = leave_one_sample_out(rows)
    phase9_labels, phase9_probabilities = predictions_from_phase9(args.predictions)
    labels = phase9_labels if len(phase9_labels) else loso_labels
    probabilities = phase9_probabilities if len(phase9_probabilities) else np.asarray(loso_probabilities, dtype=float)
    valid = ~np.isnan(probabilities)
    labels = labels[valid] if len(labels) == len(probabilities) else labels
    probabilities = probabilities[valid]
    write_tsv(args.tables_dir / "phase10_mathys_robustness_metrics.tsv", [loso_row], ["check_id", "n_units", "auroc", "auprc", "balanced_accuracy", "brier_score", "reliability", "notes"])
    write_tsv(args.tables_dir / "phase10_mathys_permutation_test.tsv", permutation_test(labels, probabilities, args.permutations), ["check_id", "n_permutations", "observed_auroc", "permutation_p_value", "reliability", "notes"])
    write_tsv(args.tables_dir / "phase10_mathys_bootstrap_ci.tsv", bootstrap_ci(labels, probabilities, args.bootstraps), ["metric", "n_bootstraps", "estimate", "ci_low", "ci_high", "reliability", "notes"])
    logging.info("Phase 10 Mathys robustness diagnostics complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
