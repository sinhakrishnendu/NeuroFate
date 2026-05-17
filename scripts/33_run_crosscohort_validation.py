#!/usr/bin/env python3
"""Run donor-level cross-cohort NeuroFate validation."""

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
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


RANDOM_SEED = 71
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


def normalize_label(value: str | None) -> str:
    return (value or "").strip().lower()


def target_high_vs_low_ad(row: dict[str, str]) -> int | None:
    pathology = normalize_label(row.get("label__Overall_AD_neuropathological_Change"))
    if "high" in pathology:
        return 1
    if any(term in pathology for term in ["none", "low", "not"]):
        return 0
    return None


def feature_columns(rows: list[dict[str, str]]) -> list[str]:
    if not rows:
        return []
    return [column for column in rows[0] if column.startswith(FEATURE_PREFIXES)]


def matrix(rows: list[dict[str, str]], features: list[str]) -> np.ndarray:
    return np.asarray([[to_float(row.get(feature, "0")) for feature in features] for row in rows], dtype=float)


def pipeline() -> Pipeline:
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "model",
                LogisticRegression(max_iter=1000, class_weight="balanced", random_state=RANDOM_SEED),
            ),
        ]
    )


def valid_rows(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], np.ndarray]:
    filtered: list[dict[str, str]] = []
    labels: list[int] = []
    for row in rows:
        label = target_high_vs_low_ad(row)
        if label is None:
            continue
        filtered.append(row)
        labels.append(label)
    return filtered, np.asarray(labels, dtype=int)


def can_fit(y_train: np.ndarray, y_test: np.ndarray) -> tuple[bool, str]:
    if len(y_train) < 4 or len(y_test) < 2:
        return False, "insufficient donor count"
    if len(Counter(y_train.tolist())) < 2:
        return False, "training set has fewer than two classes"
    if len(Counter(y_test.tolist())) < 2:
        return False, "test set has fewer than two classes"
    return True, "ok"


def metrics_row(
    mode: str,
    train_cohort: str,
    test_cohort: str,
    model: Pipeline,
    X_test: np.ndarray,
    y_test: np.ndarray,
    n_train: int,
) -> dict[str, str]:
    probabilities = model.predict_proba(X_test)[:, 1]
    predictions = (probabilities >= 0.5).astype(int)
    return {
        "validation_mode": mode,
        "train_cohort": train_cohort,
        "test_cohort": test_cohort,
        "task_id": "high_vs_low_ad_neuropathology",
        "n_train": str(n_train),
        "n_test": str(len(y_test)),
        "auroc": to_text(float(roc_auc_score(y_test, probabilities))),
        "auprc": to_text(float(average_precision_score(y_test, probabilities))),
        "balanced_accuracy": to_text(float(balanced_accuracy_score(y_test, predictions))),
        "brier_score": to_text(float(brier_score_loss(y_test, probabilities))),
        "notes": "donor-level cross-cohort logistic regression",
    }


def skipped_row(mode: str, train_cohort: str, test_cohort: str, reason: str) -> dict[str, str]:
    return {
        "validation_mode": mode,
        "train_cohort": train_cohort,
        "test_cohort": test_cohort,
        "task_id": "high_vs_low_ad_neuropathology",
        "n_train": "0",
        "n_test": "0",
        "auroc": "nan",
        "auprc": "nan",
        "balanced_accuracy": "nan",
        "brier_score": "nan",
        "notes": reason,
    }


def run_transfer(rows: list[dict[str, str]], features: list[str]) -> list[dict[str, str]]:
    cohorts = sorted({row.get("cohort_id", "missing") for row in rows})
    sea_ad_labels = {"sea_ad", "sea-ad", "seaad"}
    train_rows = [row for row in rows if row.get("cohort_id", "").lower() in sea_ad_labels]
    results: list[dict[str, str]] = []
    for cohort in cohorts:
        if cohort.lower() in sea_ad_labels:
            continue
        test_rows = [row for row in rows if row.get("cohort_id") == cohort]
        train_filtered, y_train = valid_rows(train_rows)
        test_filtered, y_test = valid_rows(test_rows)
        ok, reason = can_fit(y_train, y_test)
        if not ok:
            results.append(skipped_row("train_sea_ad_test_external", "sea_ad", cohort, reason))
            continue
        model = pipeline()
        model.fit(matrix(train_filtered, features), y_train)
        results.append(
            metrics_row(
                "train_sea_ad_test_external",
                "sea_ad",
                cohort,
                model,
                matrix(test_filtered, features),
                y_test,
                len(y_train),
            )
        )
    if not results:
        results.append(skipped_row("train_sea_ad_test_external", "sea_ad", "external", "no external cohort rows available"))
    return results


def run_leave_one_cohort_out(rows: list[dict[str, str]], features: list[str]) -> list[dict[str, str]]:
    cohorts = sorted({row.get("cohort_id", "missing") for row in rows})
    results: list[dict[str, str]] = []
    for cohort in cohorts:
        train_rows = [row for row in rows if row.get("cohort_id") != cohort]
        test_rows = [row for row in rows if row.get("cohort_id") == cohort]
        train_filtered, y_train = valid_rows(train_rows)
        test_filtered, y_test = valid_rows(test_rows)
        ok, reason = can_fit(y_train, y_test)
        if not ok:
            results.append(skipped_row("leave_one_cohort_out", "all_except_" + cohort, cohort, reason))
            continue
        model = pipeline()
        model.fit(matrix(train_filtered, features), y_train)
        results.append(
            metrics_row(
                "leave_one_cohort_out",
                "all_except_" + cohort,
                cohort,
                model,
                matrix(test_filtered, features),
                y_test,
                len(y_train),
            )
        )
    return results


def run_pooled(rows: list[dict[str, str]], features: list[str]) -> list[dict[str, str]]:
    filtered, y = valid_rows(rows)
    counts = Counter(y.tolist())
    if len(filtered) < 8 or len(counts) < 2 or min(counts.values()) < 2:
        return [skipped_row("pooled_multicohort_training", "pooled", "pooled_holdout", "insufficient labeled donors")]
    indices = np.arange(len(filtered))
    from sklearn.model_selection import train_test_split

    train_idx, test_idx = train_test_split(indices, test_size=0.30, random_state=RANDOM_SEED, stratify=y)
    model = pipeline()
    X = matrix(filtered, features)
    model.fit(X[train_idx], y[train_idx])
    return [
        metrics_row(
            "pooled_multicohort_training",
            "pooled",
            "pooled_holdout",
            model,
            X[test_idx],
            y[test_idx],
            len(train_idx),
        )
    ]


def summary_rows(metrics: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for mode in sorted({row["validation_mode"] for row in metrics}):
        usable = [row for row in metrics if row["validation_mode"] == mode and row["auroc"] != "nan"]
        rows.append(
            {
                "validation_mode": mode,
                "n_completed_comparisons": str(len(usable)),
                "mean_auroc": to_text(sum(to_float(row["auroc"]) for row in usable) / len(usable)) if usable else "nan",
                "mean_auprc": to_text(sum(to_float(row["auprc"]) for row in usable) / len(usable)) if usable else "nan",
                "notes": "summary of donor-level cross-cohort validation",
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 7 cross-cohort validation.")
    parser.add_argument(
        "--features",
        type=Path,
        default=Path("results/tables/crosscohort_donor_feature_table.tsv"),
    )
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/33_run_crosscohort_validation.log"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = read_tsv(args.features)
    features = feature_columns(rows)
    logging.info("Cross-cohort donor rows: %d", len(rows))
    logging.info("Feature columns: %d", len(features))
    if not rows or not features:
        raise RuntimeError("Cross-cohort donor table must contain rows and feature columns.")

    metrics = [
        *run_transfer(rows, features),
        *run_leave_one_cohort_out(rows, features),
        *run_pooled(rows, features),
    ]
    write_tsv(
        args.tables_dir / "phase7_crosscohort_metrics.tsv",
        metrics,
        [
            "validation_mode",
            "train_cohort",
            "test_cohort",
            "task_id",
            "n_train",
            "n_test",
            "auroc",
            "auprc",
            "balanced_accuracy",
            "brier_score",
            "notes",
        ],
    )
    write_tsv(
        args.tables_dir / "phase7_generalization_summary.tsv",
        summary_rows(metrics),
        ["validation_mode", "n_completed_comparisons", "mean_auroc", "mean_auprc", "notes"],
    )
    logging.info("Cross-cohort validation complete using donor-level feature table only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
