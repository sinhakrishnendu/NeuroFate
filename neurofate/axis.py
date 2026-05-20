"""Public axis-scoring utilities for the NeuroFate CLI.

The functions in this module intentionally operate on donor/sample-level tables.
They do not read single-cell containers, H5AD files, or genome-wide binary
formats. The goal is a small, reproducible public API that can be packaged on
PyPI and exercised by the bundled demo.
"""

from __future__ import annotations

import csv
import gzip
import math
from collections import Counter
from pathlib import Path
from typing import Iterable


RESEARCH_USE_NOTICE = (
    "NeuroFate is intended for research use only. It is not validated for "
    "clinical diagnosis, patient-level decision-making, or treatment selection."
)


def open_text(path: Path, mode: str = "rt"):
    """Open plain or gzip-compressed text."""

    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode.replace("t", ""), encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def normalize_token(value: str) -> str:
    return str(value).strip().casefold()


def parse_binary_label(value: str, positive_class: str, negative_class: str) -> int | None:
    normalized = normalize_token(value)
    positive_values = {normalize_token(positive_class), "1", "true", "case", "disease", "ad", "pd"}
    negative_values = {normalize_token(negative_class), "0", "false", "control", "reference"}
    if normalized in positive_values:
        return 1
    if normalized in negative_values:
        return 0
    return None


def load_axis_registry(path: Path) -> dict[str, list[str]]:
    axes: dict[str, list[str]] = {}
    for row in read_tsv(path):
        axis_id = row.get("axis_id", "").strip()
        members = row.get("gene_members", "").replace(",", ";")
        genes = [gene.strip().upper() for gene in members.split(";") if gene.strip()]
        if axis_id and genes:
            axes[axis_id] = genes
    if not axes:
        raise ValueError(f"No axes with gene_members found in {path}")
    return axes


def detect_expression_orientation(header: list[str], first_row: dict[str, str]) -> str:
    first_col = header[0].strip().casefold()
    if first_col in {"gene", "genes", "gene_symbol", "symbol", "feature", "feature_id", "id_ref"}:
        return "genes_rows"
    if first_col in {"sample", "sample_id", "donor", "donor_id"}:
        return "samples_rows"
    # Heuristic: expression matrices with genes as rows usually have many sample
    # columns and a first-column value that is not a sample identifier.
    if any(key.upper() in first_row for key in ["APOE", "SNCA", "MAPT", "GFAP"]):
        return "samples_rows"
    return "genes_rows"


def load_expression_matrix(
    path: Path,
    orientation: str = "auto",
    sample_id_column: str = "sample_id",
    gene_column: str = "gene_symbol",
) -> tuple[dict[str, dict[str, float]], list[str]]:
    """Load a compact sample-level expression table.

    Returns a mapping of sample_id -> gene_symbol -> expression value and a list
    of warnings. The intended input is already sample/donor-level or a compact
    axis-gene matrix, not a large raw matrix.
    """

    warnings: list[str] = []
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Expression table has no header: {path}")
        header = list(reader.fieldnames)
        rows = list(reader)

    if not rows:
        raise ValueError(f"Expression table has no data rows: {path}")
    resolved_orientation = (
        detect_expression_orientation(header, rows[0]) if orientation == "auto" else orientation
    )
    sample_gene: dict[str, dict[str, float]] = {}

    if resolved_orientation == "samples_rows":
        if sample_id_column not in header:
            raise ValueError(
                f"Sample-row expression requires sample id column {sample_id_column!r}; "
                f"available columns are {header}"
            )
        for row in rows:
            sample_id = row[sample_id_column].strip()
            if not sample_id:
                continue
            sample_gene.setdefault(sample_id, {})
            for key, value in row.items():
                if key == sample_id_column:
                    continue
                try:
                    sample_gene[sample_id][key.strip().upper()] = float(value)
                except (TypeError, ValueError):
                    continue
    elif resolved_orientation == "genes_rows":
        gene_col = gene_column if gene_column in header else header[0]
        sample_columns = [col for col in header if col != gene_col]
        if not sample_columns:
            raise ValueError(f"Gene-row expression has no sample columns: {path}")
        for sample_id in sample_columns:
            sample_gene.setdefault(sample_id, {})
        for row in rows:
            gene = row.get(gene_col, "").strip().upper()
            if not gene:
                continue
            for sample_id in sample_columns:
                try:
                    sample_gene[sample_id][gene] = float(row.get(sample_id, ""))
                except (TypeError, ValueError):
                    continue
    else:
        raise ValueError("orientation must be auto, genes_rows, or samples_rows")

    warnings.append(f"expression_orientation={resolved_orientation}")
    return sample_gene, warnings


def mean(values: Iterable[float]) -> float:
    values = list(values)
    return sum(values) / len(values) if values else math.nan


def build_axis_score_tables(
    expression: Path,
    metadata: Path,
    axis_registry: Path,
    outdir: Path,
    sample_id_column: str,
    endpoint_column: str,
    positive_class: str,
    negative_class: str,
    orientation: str = "auto",
    gene_column: str = "gene_symbol",
) -> dict[str, Path]:
    axes = load_axis_registry(axis_registry)
    sample_gene, warnings = load_expression_matrix(
        expression,
        orientation=orientation,
        sample_id_column=sample_id_column,
        gene_column=gene_column,
    )
    metadata_rows = read_tsv(metadata)
    if not metadata_rows:
        raise ValueError(f"Metadata table has no rows: {metadata}")
    if sample_id_column not in metadata_rows[0]:
        raise ValueError(f"Metadata does not contain sample id column {sample_id_column!r}")
    if endpoint_column not in metadata_rows[0]:
        raise ValueError(f"Metadata does not contain endpoint column {endpoint_column!r}")

    label_counts: Counter[str] = Counter()
    score_rows: list[dict[str, object]] = []
    warnings_rows: list[dict[str, object]] = [{"warning": warning} for warning in warnings]
    for row in metadata_rows:
        sample_id = row.get(sample_id_column, "").strip()
        if sample_id not in sample_gene:
            warnings_rows.append({"warning": f"metadata_sample_without_expression={sample_id}"})
            continue
        label = parse_binary_label(row.get(endpoint_column, ""), positive_class, negative_class)
        if label is None:
            warnings_rows.append({"warning": f"ambiguous_endpoint_sample={sample_id}"})
            continue
        label_counts[str(label)] += 1
        output_row: dict[str, object] = {
            "sample_id": sample_id,
            endpoint_column: row.get(endpoint_column, ""),
            "label": label,
            "research_use_only": "true",
        }
        for axis_id, genes in axes.items():
            available_values = [sample_gene[sample_id][gene] for gene in genes if gene in sample_gene[sample_id]]
            output_row[axis_id] = f"{mean(available_values):.6f}" if available_values else ""
        score_rows.append(output_row)

    if not score_rows:
        raise ValueError("No samples with unambiguous endpoint labels and expression values were found.")

    coverage_rows: list[dict[str, object]] = []
    universe = {gene for genes in sample_gene.values() for gene in genes}
    for axis_id, genes in axes.items():
        mapped = [gene for gene in genes if gene in universe]
        missing = [gene for gene in genes if gene not in universe]
        coverage_rows.append(
            {
                "axis_id": axis_id,
                "axis_gene_count": len(genes),
                "mapped_gene_count": len(mapped),
                "mapped_genes": ";".join(mapped),
                "missing_genes": ";".join(missing),
                "coverage_fraction": f"{len(mapped) / len(genes):.6f}",
            }
        )

    label_rows = [
        {
            "label": label,
            "count": count,
            "meaning": positive_class if label == "1" else negative_class,
        }
        for label, count in sorted(label_counts.items())
    ]

    config_rows = [
        {"key": "expression", "value": str(expression)},
        {"key": "metadata", "value": str(metadata)},
        {"key": "axis_registry", "value": str(axis_registry)},
        {"key": "sample_id_column", "value": sample_id_column},
        {"key": "endpoint_column", "value": endpoint_column},
        {"key": "positive_class", "value": positive_class},
        {"key": "negative_class", "value": negative_class},
        {"key": "research_use_only", "value": RESEARCH_USE_NOTICE},
    ]

    outdir.mkdir(parents=True, exist_ok=True)
    score_path = outdir / "axis_scores.tsv"
    coverage_path = outdir / "axis_feature_coverage.tsv"
    labels_path = outdir / "label_summary.tsv"
    config_path = outdir / "run_config.yaml"
    warnings_path = outdir / "warnings.tsv"
    score_fields = ["sample_id", endpoint_column, "label", "research_use_only"] + list(axes)
    write_tsv(score_path, score_rows, score_fields)
    write_tsv(
        coverage_path,
        coverage_rows,
        ["axis_id", "axis_gene_count", "mapped_gene_count", "mapped_genes", "missing_genes", "coverage_fraction"],
    )
    write_tsv(labels_path, label_rows, ["label", "count", "meaning"])
    write_tsv(warnings_path, warnings_rows, ["warning"])
    config_path.write_text(
        "\n".join(f"{row['key']}: {row['value']}" for row in config_rows) + "\n",
        encoding="utf-8",
    )
    return {
        "axis_scores": score_path,
        "axis_feature_coverage": coverage_path,
        "label_summary": labels_path,
        "run_config": config_path,
        "warnings": warnings_path,
    }


def score_research_risk(axis_scores: Path, outdir: Path) -> dict[str, Path]:
    rows = read_tsv(axis_scores)
    if not rows:
        raise ValueError(f"Axis score table has no rows: {axis_scores}")
    metadata_columns = {"sample_id", "label", "research_use_only"}
    axis_columns = [
        col
        for col in rows[0]
        if col not in metadata_columns and not col.startswith("label__") and col != "diagnosis"
    ]
    risk_rows: list[dict[str, object]] = []
    for row in rows:
        numeric = []
        for col in axis_columns:
            try:
                numeric.append(float(row.get(col, "")))
            except (TypeError, ValueError):
                continue
        raw_score = mean(numeric)
        risk_rows.append(
            {
                "sample_id": row.get("sample_id", ""),
                "neurofate_research_risk_score": f"{raw_score:.6f}" if not math.isnan(raw_score) else "",
                "axis_count_used": len(numeric),
                "research_use_only": "true",
            }
        )

    outdir.mkdir(parents=True, exist_ok=True)
    risk_path = outdir / "neurofate_risk_scores.tsv"
    report_path = outdir / "risk_score_report.md"
    write_tsv(
        risk_path,
        risk_rows,
        ["sample_id", "neurofate_research_risk_score", "axis_count_used", "research_use_only"],
    )
    report_path.write_text(
        "\n".join(
            [
                "# NeuroFate Research Risk Score Report",
                "",
                RESEARCH_USE_NOTICE,
                "",
                f"- Input axis score table: `{axis_scores}`",
                f"- Samples scored: {len(risk_rows)}",
                f"- Axis columns considered: {len(axis_columns)}",
                "",
                "The score is an exploratory donor/sample-level research score and must not be used",
                "for clinical diagnosis, patient-level decision-making, or treatment selection.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return {"risk_scores": risk_path, "risk_score_report": report_path}

