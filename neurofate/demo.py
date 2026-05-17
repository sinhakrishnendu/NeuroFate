"""Tiny bundled NeuroFate demo used for installation smoke tests."""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from importlib import resources
from pathlib import Path
from typing import Protocol


class ReadableTextResource(Protocol):
    name: str

    def is_file(self) -> bool: ...

    def open(self, mode: str = "r", *args, **kwargs): ...


RESOURCE_PACKAGE = "neurofate.resources.tiny_demo"


def demo_resource(name: str) -> ReadableTextResource:
    return resources.files(RESOURCE_PACKAGE).joinpath(name)


def read_tsv(resource: ReadableTextResource | Path) -> list[dict[str, str]]:
    with resource.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_expression(resource: ReadableTextResource | Path) -> list[dict[str, str]]:
    if resource.name.endswith(".gz"):
        with gzip.open(resource, "rt", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    return read_tsv(resource)


def require_resource(resource: ReadableTextResource | Path, label: str) -> None:
    if not resource.is_file():
        raise FileNotFoundError(f"Missing {label}: {resource}")


def rank_auc(labels: list[int], scores: list[float]) -> float:
    positives = [score for label, score in zip(labels, scores) if label == 1]
    negatives = [score for label, score in zip(labels, scores) if label == 0]
    if not positives or not negatives:
        return 0.0
    wins = 0.0
    for positive in positives:
        for negative in negatives:
            if positive > negative:
                wins += 1.0
            elif positive == negative:
                wins += 0.5
    return wins / (len(positives) * len(negatives))


def balanced_accuracy(labels: list[int], predictions: list[int]) -> float:
    positives = sum(1 for label in labels if label == 1)
    negatives = sum(1 for label in labels if label == 0)
    true_positive = sum(1 for label, pred in zip(labels, predictions) if label == 1 and pred == 1)
    true_negative = sum(1 for label, pred in zip(labels, predictions) if label == 0 and pred == 0)
    sensitivity = true_positive / positives if positives else 0.0
    specificity = true_negative / negatives if negatives else 0.0
    return (sensitivity + specificity) / 2.0


def build_demo_outputs(
    metadata_resource: ReadableTextResource | Path | None = None,
    panel_resource: ReadableTextResource | Path | None = None,
    expression_resource: ReadableTextResource | Path | None = None,
    outdir: Path = Path("results/demo"),
) -> None:
    metadata_resource = metadata_resource or demo_resource("tiny_metadata.tsv")
    panel_resource = panel_resource or demo_resource("tiny_gene_panel.tsv")
    expression_resource = expression_resource or demo_resource("tiny_sparse_expression.tsv")
    require_resource(metadata_resource, "tiny metadata")
    require_resource(panel_resource, "tiny gene panel")
    require_resource(expression_resource, "tiny sparse expression")

    metadata = read_tsv(metadata_resource)
    panel = read_tsv(panel_resource)
    expression = read_expression(expression_resource)
    cell_to_meta = {row["cell_id"]: row for row in metadata}
    genes = [row["gene_symbol"] for row in panel]

    donor_gene_values: dict[tuple[str, str], list[float]] = defaultdict(list)
    donor_labels: dict[str, str] = {}
    donor_pathology: dict[str, str] = {}
    for row in metadata:
        donor_labels[row["donor_id"]] = row["diagnosis"]
        donor_pathology[row["donor_id"]] = row["pathology_group"]
    for row in expression:
        cell_meta = cell_to_meta[row["cell_id"]]
        donor_gene_values[(cell_meta["donor_id"], row["gene_symbol"])].append(
            float(row["expression_value"])
        )

    donor_rows: list[dict[str, str]] = []
    for donor_id in sorted(donor_labels):
        feature_values: dict[str, float] = {}
        for gene in genes:
            values = donor_gene_values.get((donor_id, gene), [])
            feature_values[f"mean_{gene}"] = sum(values) / len(values) if values else 0.0
        risk_score = (
            feature_values["mean_APOE"]
            + feature_values["mean_TREM2"]
            + feature_values["mean_GFAP"]
            + feature_values["mean_SNCA"]
            - feature_values["mean_SLC17A7"]
        )
        row = {
            "donor_id": donor_id,
            "diagnosis": donor_labels[donor_id],
            "pathology_group": donor_pathology[donor_id],
            "demo_neurofate_score": f"{risk_score:.6f}",
        }
        row.update({key: f"{value:.6f}" for key, value in feature_values.items()})
        donor_rows.append(row)

    outdir.mkdir(parents=True, exist_ok=True)
    feature_path = outdir / "demo_donor_feature_table.tsv"
    fieldnames = ["donor_id", "diagnosis", "pathology_group", "demo_neurofate_score"] + [
        f"mean_{gene}" for gene in genes
    ]
    with feature_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(donor_rows)

    labels = [1 if row["diagnosis"] == "disease" else 0 for row in donor_rows]
    scores = [float(row["demo_neurofate_score"]) for row in donor_rows]
    threshold = sorted(scores)[len(scores) // 2]
    predictions = [1 if score >= threshold else 0 for score in scores]
    brier = sum((score - label) ** 2 for score, label in zip(scores, labels)) / len(scores)
    metrics = {
        "task": "tiny_demo_disease_vs_reference",
        "n_donors": str(len(donor_rows)),
        "n_genes": str(len(genes)),
        "auroc": f"{rank_auc(labels, scores):.6f}",
        "balanced_accuracy": f"{balanced_accuracy(labels, predictions):.6f}",
        "brier_score_unscaled": f"{brier:.6f}",
        "note": "Synthetic smoke-test metrics only; not biological evidence.",
    }
    metrics_path = outdir / "demo_model_metrics.tsv"
    with metrics_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics), delimiter="\t")
        writer.writeheader()
        writer.writerow(metrics)

    report_path = outdir / "demo_report.md"
    report_path.write_text(
        "\n".join(
            [
                "# NeuroFate Tiny Demo Report",
                "",
                "This demo used bundled synthetic toy data only.",
                f"- Pseudo-donors: {len(donor_rows)}",
                f"- Cell rows: {len(metadata)}",
                f"- Genes: {len(genes)}",
                f"- AUROC smoke-test value: {metrics['auroc']}",
                "",
                "The demo confirms CLI and tabular workflow plumbing; it is not a biological result.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {feature_path}")
    print(f"Wrote {metrics_path}")
    print(f"Wrote {report_path}")
