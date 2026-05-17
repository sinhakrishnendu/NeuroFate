#!/usr/bin/env python3
"""Train a small donor-level NeuroFate MLP with PyTorch MPS support."""

from __future__ import annotations

import argparse
import csv
import logging
import math
import random
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split

import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


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


@dataclass
class SplitData:
    train_idx: np.ndarray
    val_idx: np.ndarray
    test_idx: np.ndarray


class NeuroFateMLP(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int, dropout: float):
        super().__init__()
        second_hidden = max(4, hidden_dim // 2)
        self.network = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, second_hidden),
            nn.ReLU(),
            nn.Linear(second_hidden, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.network(x).squeeze(-1)


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


def load_config(path: Path) -> dict[str, float | int | str]:
    with path.open("r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle) or {}
    return loaded


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.backends.mps.is_available():
        torch.mps.manual_seed(seed)


def selected_device(preference: str) -> torch.device:
    if preference == "mps" and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


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
    positive_terms = ["high", "intermediate", "moderate", "severe", "limbic", "neocortical", "present", "stage"]
    negative_terms = ["none", "not", "absent", "no ", "low", "0"]
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


TARGET_BUILDERS = {
    "dementia_vs_reference": target_dementia_vs_reference,
    "high_vs_low_ad_neuropathology": target_high_vs_low_ad,
    "apoe_risk_prediction": target_apoe_risk,
    "mixed_pathology_burden": target_mixed_pathology,
}


def feature_columns(fieldnames: list[str]) -> list[str]:
    return [field for field in fieldnames if field.startswith(FEATURE_PREFIXES)]


def build_task_dataset(
    rows: list[dict[str, str]],
    features: list[str],
    task_id: str,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    donors: list[str] = []
    task_rows: list[dict[str, str]] = []
    labels: list[int] = []
    builder = TARGET_BUILDERS[task_id]
    for row in rows:
        label = builder(row)
        if label is None:
            continue
        donors.append(row["donor_id"])
        task_rows.append(row)
        labels.append(label)
    matrix = np.asarray([[to_float(row.get(feature, "0")) for feature in features] for row in task_rows], dtype=np.float32)
    return donors, matrix, np.asarray(labels, dtype=np.float32)


def can_train(labels: np.ndarray) -> tuple[bool, str]:
    counts = Counter(labels.astype(int).tolist())
    if len(labels) < 8:
        return False, "fewer than eight labeled donors"
    if len(counts) < 2:
        return False, "fewer than two target classes"
    if min(counts.values()) < 3:
        return False, "at least one target class has fewer than three donors"
    return True, "ok"


def split_indices(labels: np.ndarray, config: dict[str, float | int | str]) -> SplitData:
    seed = int(config["seed"])
    test_size = float(config["test_size"])
    validation_size = float(config["validation_size"])
    indices = np.arange(len(labels))
    train_val_idx, test_idx = train_test_split(
        indices,
        test_size=test_size,
        random_state=seed,
        stratify=labels,
    )
    relative_validation = validation_size / max(1e-8, 1.0 - test_size)
    train_idx, val_idx = train_test_split(
        train_val_idx,
        test_size=relative_validation,
        random_state=seed,
        stratify=labels[train_val_idx],
    )
    return SplitData(train_idx=np.asarray(train_idx), val_idx=np.asarray(val_idx), test_idx=np.asarray(test_idx))


def standardize_from_train(
    X: np.ndarray,
    train_idx: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = X[train_idx].mean(axis=0)
    std = X[train_idx].std(axis=0)
    std[std == 0] = 1.0
    return ((X - mean) / std).astype(np.float32), mean.astype(np.float32), std.astype(np.float32)


def make_loader(
    X: np.ndarray,
    y: np.ndarray,
    indices: np.ndarray,
    batch_size: int,
    device: torch.device,
    shuffle: bool,
) -> DataLoader:
    x_tensor = torch.tensor(X[indices], dtype=torch.float32, device=device)
    y_tensor = torch.tensor(y[indices], dtype=torch.float32, device=device)
    dataset = TensorDataset(x_tensor, y_tensor)
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle)


def loss_on_loader(
    model: nn.Module,
    loader: DataLoader,
    criterion: nn.Module,
) -> float:
    model.eval()
    losses: list[float] = []
    with torch.no_grad():
        for x_batch, y_batch in loader:
            logits = model(x_batch)
            loss = criterion(logits, y_batch)
            losses.append(float(loss.detach().cpu().item()))
    return float(np.mean(losses)) if losses else float("nan")


def predict_probabilities(
    model: nn.Module,
    X: np.ndarray,
    indices: np.ndarray,
    device: torch.device,
) -> np.ndarray:
    model.eval()
    x_tensor = torch.tensor(X[indices], dtype=torch.float32, device=device)
    with torch.no_grad():
        logits = model(x_tensor)
        probabilities = torch.sigmoid(logits).detach().cpu().numpy()
    return probabilities.astype(float)


def safe_metric(metric, *args) -> float:
    try:
        return float(metric(*args))
    except ValueError:
        return float("nan")


def train_task(
    task_id: str,
    donors: list[str],
    X: np.ndarray,
    y: np.ndarray,
    config: dict[str, float | int | str],
    device: torch.device,
    models_dir: Path,
) -> tuple[dict[str, str], list[dict[str, str]], list[dict[str, str]]]:
    split = split_indices(y, config)
    X_scaled, mean, std = standardize_from_train(X, split.train_idx)
    hidden_dim = int(config["hidden_dim"])
    dropout = float(config["dropout"])
    batch_size = int(config["batch_size"])
    epochs = int(config["epochs"])
    patience = int(config["early_stopping_patience"])

    model = NeuroFateMLP(input_dim=X_scaled.shape[1], hidden_dim=hidden_dim, dropout=dropout).to(device)
    positives = float(y[split.train_idx].sum())
    negatives = float(len(split.train_idx) - positives)
    pos_weight = torch.tensor([negatives / positives], dtype=torch.float32, device=device) if positives else None
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(config["learning_rate"]),
        weight_decay=float(config["weight_decay"]),
    )
    train_loader = make_loader(X_scaled, y, split.train_idx, batch_size, device, shuffle=True)
    val_loader = make_loader(X_scaled, y, split.val_idx, batch_size, device, shuffle=False)
    test_loader = make_loader(X_scaled, y, split.test_idx, batch_size, device, shuffle=False)

    best_state: dict[str, torch.Tensor] | None = None
    best_val_loss = float("inf")
    best_epoch = 0
    wait = 0
    log_rows: list[dict[str, str]] = []
    for epoch in range(1, epochs + 1):
        model.train()
        train_losses: list[float] = []
        for x_batch, y_batch in train_loader:
            optimizer.zero_grad(set_to_none=True)
            logits = model(x_batch)
            loss = criterion(logits, y_batch)
            loss.backward()
            optimizer.step()
            train_losses.append(float(loss.detach().cpu().item()))

        train_loss = float(np.mean(train_losses)) if train_losses else float("nan")
        val_loss = loss_on_loader(model, val_loader, criterion)
        log_rows.append(
            {
                "task_id": task_id,
                "epoch": str(epoch),
                "train_loss": to_text(train_loss),
                "validation_loss": to_text(val_loss),
                "selected_device": str(device),
            }
        )
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_epoch = epoch
            best_state = {name: tensor.detach().cpu().clone() for name, tensor in model.state_dict().items()}
            wait = 0
        else:
            wait += 1
            if wait >= patience:
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    test_loss = loss_on_loader(model, test_loader, criterion)
    test_probabilities = predict_probabilities(model, X_scaled, split.test_idx, device)
    test_labels = y[split.test_idx].astype(int)
    test_predictions = (test_probabilities >= 0.5).astype(int)

    auroc = safe_metric(roc_auc_score, test_labels, test_probabilities) if len(set(test_labels.tolist())) == 2 else float("nan")
    auprc = safe_metric(average_precision_score, test_labels, test_probabilities) if len(set(test_labels.tolist())) == 2 else float("nan")
    balanced_accuracy = safe_metric(balanced_accuracy_score, test_labels, test_predictions)
    brier = safe_metric(brier_score_loss, test_labels, test_probabilities)

    models_dir.mkdir(parents=True, exist_ok=True)
    model_path = models_dir / f"neurofate_mps_{task_id}.pt"
    torch.save(
        {
            "task_id": task_id,
            "model_state_dict": model.state_dict(),
            "feature_mean": mean,
            "feature_std": std,
            "config": config,
            "input_dim": int(X_scaled.shape[1]),
        },
        model_path,
    )

    prediction_rows: list[dict[str, str]] = []
    for donor_index, probability, label, prediction in zip(split.test_idx, test_probabilities, test_labels, test_predictions, strict=False):
        prediction_rows.append(
            {
                "task_id": task_id,
                "donor_id": donors[int(donor_index)],
                "split": "test",
                "true_label": str(int(label)),
                "predicted_probability": to_text(float(probability)),
                "predicted_label": str(int(prediction)),
            }
        )

    metric_row = {
        "task_id": task_id,
        "target_label": TASK_LABELS[task_id],
        "model_name": "neurofate_mps_mlp",
        "n_donors": str(len(donors)),
        "n_train": str(len(split.train_idx)),
        "n_validation": str(len(split.val_idx)),
        "n_test": str(len(split.test_idx)),
        "selected_device": str(device),
        "best_epoch": str(best_epoch),
        "validation_loss": to_text(best_val_loss),
        "test_loss": to_text(test_loss),
        "auroc": to_text(auroc),
        "auprc": to_text(auprc),
        "balanced_accuracy": to_text(balanced_accuracy),
        "brier_score": to_text(brier),
        "model_path": str(model_path),
        "notes": "small donor-level MLP; train-only standardization; donor-level split",
    }
    return metric_row, log_rows, prediction_rows


def skipped_metric_row(task_id: str, n_donors: int, reason: str) -> dict[str, str]:
    return {
        "task_id": task_id,
        "target_label": TASK_LABELS[task_id],
        "model_name": "not_run",
        "n_donors": str(n_donors),
        "n_train": "0",
        "n_validation": "0",
        "n_test": "0",
        "selected_device": "not_selected",
        "best_epoch": "0",
        "validation_loss": "nan",
        "test_loss": "nan",
        "auroc": "nan",
        "auprc": "nan",
        "balanced_accuracy": "nan",
        "brier_score": "nan",
        "model_path": "",
        "notes": reason,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train the Phase 6 NeuroFate donor-level MPS MLP.")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/neurofate_mps_model_config.yaml"),
    )
    parser.add_argument(
        "--features",
        type=Path,
        default=Path("results/tables/phase5_donor_feature_table.tsv"),
    )
    parser.add_argument("--task", choices=sorted(TARGET_BUILDERS), default="dementia_vs_reference")
    parser.add_argument("--all-tasks", action="store_true")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--models-dir", type=Path, default=Path("results/models"))
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/27_train_neurofate_mps_model.log"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    config = load_config(args.config)
    set_seed(int(config["seed"]))
    device = selected_device(str(config["device_preference"]))
    logging.info("Selected device: %s", device)
    logging.info("Feature table: %s", args.features)

    rows = read_tsv(args.features)
    if not rows:
        raise RuntimeError(f"No donor rows found in {args.features}")
    features = feature_columns(list(rows[0].keys()))
    if not features:
        raise RuntimeError("No donor-level feature columns found.")
    logging.info("Donor rows: %d", len(rows))
    logging.info("Feature columns: %d", len(features))

    tasks = sorted(TARGET_BUILDERS) if args.all_tasks else [args.task]
    metric_rows: list[dict[str, str]] = []
    training_log_rows: list[dict[str, str]] = []
    prediction_rows: list[dict[str, str]] = []
    for task_id in tasks:
        donors, X, y = build_task_dataset(rows, features, task_id)
        ok, reason = can_train(y)
        logging.info("%s: labeled donors=%d, status=%s", task_id, len(donors), reason)
        if not ok:
            metric_rows.append(skipped_metric_row(task_id, len(donors), reason))
            continue
        metric_row, log_rows, task_prediction_rows = train_task(
            task_id,
            donors,
            X,
            y,
            config,
            device,
            args.models_dir,
        )
        metric_rows.append(metric_row)
        training_log_rows.extend(log_rows)
        prediction_rows.extend(task_prediction_rows)

    write_tsv(
        args.tables_dir / "phase6_mps_model_metrics.tsv",
        metric_rows,
        [
            "task_id",
            "target_label",
            "model_name",
            "n_donors",
            "n_train",
            "n_validation",
            "n_test",
            "selected_device",
            "best_epoch",
            "validation_loss",
            "test_loss",
            "auroc",
            "auprc",
            "balanced_accuracy",
            "brier_score",
            "model_path",
            "notes",
        ],
    )
    write_tsv(
        args.tables_dir / "phase6_mps_training_log.tsv",
        training_log_rows,
        ["task_id", "epoch", "train_loss", "validation_loss", "selected_device"],
    )
    write_tsv(
        args.tables_dir / "phase6_mps_predictions.tsv",
        prediction_rows,
        ["task_id", "donor_id", "split", "true_label", "predicted_probability", "predicted_label"],
    )
    logging.info("Phase 6 donor-level MPS modeling script complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
