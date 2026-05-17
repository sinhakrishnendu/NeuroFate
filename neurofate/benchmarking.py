"""Lightweight donor-level benchmarking helpers for NeuroFate."""

from __future__ import annotations

from math import sqrt
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


FEATURE_PREFIXES = (
    "gene_mean__",
    "gene_detection__",
    "index__",
    "cell_fraction__",
    "celltype_index__",
    "mean_",
)
EXCLUDED_COLUMNS = {"donor_id", "sample_id", "cell_id", "cohort_id", "n_cells"}
LABEL_PREFIX = "label__"


def load_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def load_donor_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing donor-level feature table: {path}")
    return pd.read_csv(path, sep="\t")


def select_feature_columns(frame: pd.DataFrame) -> list[str]:
    columns = []
    for column in frame.columns:
        if column in EXCLUDED_COLUMNS or column.startswith(LABEL_PREFIX):
            continue
        if column.startswith(FEATURE_PREFIXES):
            columns.append(column)
    return columns


def feature_matrix(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    return frame[columns].apply(pd.to_numeric, errors="coerce").fillna(0.0)


def task_labels(frame: pd.DataFrame, task: str) -> pd.Series:
    if task == "dementia_vs_reference":
        values = frame.get("label__Cognitive_Status", pd.Series(index=frame.index, dtype=object))
        labels = values.map({"Dementia": 1, "Reference": 0, "No dementia": 0})
    elif task == "high_vs_low_ad_neuropathology":
        values = frame.get(
            "label__Overall_AD_neuropathological_Change",
            pd.Series(index=frame.index, dtype=object),
        )
        labels = values.map({"High": 1, "Intermediate": 1, "Low": 0, "Reference": 0})
    elif task == "apoe_risk_prediction":
        values = frame.get("label__APOE_Genotype", pd.Series(index=frame.index, dtype=object))
        labels = values.astype(str).map(lambda text: np.nan if text == "missing" else int("4" in text))
    elif task == "mixed_pathology_burden":
        lewy = frame.get("label__Highest_Lewy_Body_Disease", pd.Series("", index=frame.index)).astype(str)
        late = frame.get("label__LATE", pd.Series("", index=frame.index)).astype(str)
        caa = frame.get("label__Overall_CAA_Score", pd.Series("", index=frame.index)).astype(str)
        mixed = (
            lewy.str.contains("limbic|neocortical|brainstem", case=False, regex=True)
            | late.str.contains("stage [123]", case=False, regex=True)
            | caa.str.contains("moderate|severe", case=False, regex=True)
        )
        labels = mixed.astype(int)
    else:
        labels = pd.Series(index=frame.index, dtype=float)
    return labels


def make_model(name: str, seed: int):
    if name == "logistic_regression":
        return Pipeline(
            [
                ("scale", StandardScaler()),
                ("model", LogisticRegression(max_iter=1000, class_weight="balanced", random_state=seed)),
            ]
        )
    if name == "elastic_net":
        return Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        max_iter=1500,
                        penalty="elasticnet",
                        solver="saga",
                        l1_ratio=0.5,
                        class_weight="balanced",
                        random_state=seed,
                    ),
                ),
            ]
        )
    if name == "random_forest":
        return RandomForestClassifier(
            n_estimators=200,
            max_depth=4,
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=seed,
        )
    if name == "gradient_boosting":
        return GradientBoostingClassifier(random_state=seed)
    raise ValueError(f"Unsupported model: {name}")


def metric_dict(y_true: np.ndarray, score: np.ndarray) -> dict[str, float]:
    prediction = (score >= 0.5).astype(int)
    try:
        auroc = float(roc_auc_score(y_true, score))
    except ValueError:
        auroc = float("nan")
    try:
        auprc = float(average_precision_score(y_true, score))
    except ValueError:
        auprc = float("nan")
    return {
        "auroc": auroc,
        "auprc": auprc,
        "balanced_accuracy": float(balanced_accuracy_score(y_true, prediction)),
        "brier": float(brier_score_loss(y_true, score)),
    }


def train_evaluate_split(
    frame: pd.DataFrame,
    task: str,
    model_name: str,
    seed: int,
    test_size: float,
    min_class_count: int,
    feature_columns: list[str] | None = None,
) -> dict[str, Any]:
    features = feature_columns or select_feature_columns(frame)
    labels = task_labels(frame, task)
    valid = labels.notna()
    labels = labels[valid].astype(int)
    data = feature_matrix(frame.loc[valid], features)
    counts = labels.value_counts()
    if len(counts) < 2 or int(counts.min()) < min_class_count:
        return {
            "task": task,
            "model": model_name,
            "seed": seed,
            "status": "skipped_insufficient_class_count",
            "n_samples": int(len(labels)),
            "n_features": int(len(features)),
        }
    x_train, x_test, y_train, y_test = train_test_split(
        data,
        labels,
        test_size=test_size,
        random_state=seed,
        stratify=labels,
    )
    model = make_model(model_name, seed)
    model.fit(x_train, y_train)
    if hasattr(model, "predict_proba"):
        score = model.predict_proba(x_test)[:, 1]
    else:
        raw = model.decision_function(x_test)
        score = 1.0 / (1.0 + np.exp(-raw))
    metrics = metric_dict(y_test.to_numpy(), np.asarray(score))
    return {
        "task": task,
        "model": model_name,
        "seed": seed,
        "status": "ok",
        "n_samples": int(len(labels)),
        "n_features": int(len(features)),
        **metrics,
    }


def summarize_metric_rows(rows: list[dict[str, Any]], metric_names: list[str]) -> list[dict[str, Any]]:
    frame = pd.DataFrame([row for row in rows if row.get("status") == "ok"])
    if frame.empty:
        return []
    summaries = []
    for (task, model), group in frame.groupby(["task", "model"], dropna=False):
        row: dict[str, Any] = {
            "task": task,
            "model": model,
            "n_repeats": int(len(group)),
            "n_samples_min": int(group["n_samples"].min()),
            "n_features": int(group["n_features"].max()),
        }
        for metric in metric_names:
            values = pd.to_numeric(group[metric], errors="coerce").dropna()
            if values.empty:
                row[f"{metric}_mean"] = np.nan
                row[f"{metric}_sd"] = np.nan
                row[f"{metric}_ci95_low"] = np.nan
                row[f"{metric}_ci95_high"] = np.nan
                continue
            mean = float(values.mean())
            sd = float(values.std(ddof=1)) if len(values) > 1 else 0.0
            half_width = 1.96 * sd / sqrt(len(values)) if len(values) > 1 else 0.0
            row[f"{metric}_mean"] = mean
            row[f"{metric}_sd"] = sd
            row[f"{metric}_ci95_low"] = mean - half_width
            row[f"{metric}_ci95_high"] = mean + half_width
        summaries.append(row)
    return summaries


def feature_groups(columns: list[str]) -> dict[str, list[str]]:
    return {
        "gene_level_features": [
            column for column in columns if column.startswith(("gene_mean__", "gene_detection__", "mean_"))
        ],
        "cell_fraction_features": [column for column in columns if column.startswith("cell_fraction__")],
        "celltype_index_features": [column for column in columns if column.startswith("celltype_index__")],
        "inflammatory_signatures": [
            column
            for column in columns
            if "inflammatory" in column.lower()
            or any(gene in column for gene in ["TREM2", "TYROBP", "GPNMB", "HLA_DRA", "AIF1", "IL1B", "TNF"])
        ],
        "astrocyte_signatures": [
            column for column in columns if "ASI" in column or "GFAP" in column or "Astrocyte" in column
        ],
        "neuronal_signatures": [
            column
            for column in columns
            if "neuronal" in column.lower()
            or any(gene in column for gene in ["SLC17A7", "SST", "PVALB", "LAMP5"])
        ],
        "mitochondrial_neurodegeneration_signatures": [
            column
            for column in columns
            if "mitochondrial" in column.lower()
            or "NVI" in column
            or any(gene in column for gene in ["PINK1", "PRKN", "SNCA", "MAPT", "APOE"])
        ],
    }
