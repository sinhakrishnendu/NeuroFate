#!/usr/bin/env python3
"""Audit Mathys 2019 external validation for robustness and no-overclaiming."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import math
from collections import Counter
from pathlib import Path
from typing import TextIO


FEATURE_PREFIXES = (
    "gene_mean__",
    "gene_detection__",
    "index__",
    "cell_fraction__",
    "celltype_index__",
)
SMALL_N_WARNING_THRESHOLD = 20
MIN_RELIABLE_TEST_N = 20


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


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


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


def count_csv_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with open_text(path) as handle:
        reader = csv.reader(handle)
        next(reader, None)
        return sum(1 for _ in reader)


def count_sparse_rows_and_genes(path: Path) -> tuple[int, set[str]]:
    if not path.exists():
        return 0, set()
    genes: set[str] = set()
    row_count = 0
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row_count += 1
            gene = row.get("gene_symbol", "")
            if gene:
                genes.add(gene)
    return row_count, genes


def panel_genes(path: Path) -> list[str]:
    rows = read_tsv(path)
    return [row["gene_symbol"] for row in rows if row.get("gene_symbol")]


def feature_columns(rows: list[dict[str, str]]) -> set[str]:
    if not rows:
        return set()
    return {field for field in rows[0] if field.startswith(FEATURE_PREFIXES)}


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def diagnosis_counts(feature_rows: list[dict[str, str]]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for row in feature_rows:
        label = row.get("label__diagnosis") or row.get("label__Overall_AD_neuropathological_Change") or "missing"
        counts[label] += 1
    return counts


def celltype_counts_from_covariates(path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    if not path.exists():
        return counts
    with open_text(path) as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            counts[row.get("oupSample.cellType", "missing") or "missing"] += 1
    return counts


def gene_audit_rows(panel: list[str], overlap_rows: list[dict[str, str]], extracted_genes: set[str]) -> list[dict[str, str]]:
    overlap_status = {
        row.get("gene_symbol", ""): row.get("mathys_status", "")
        for row in overlap_rows
        if row.get("gene_symbol")
    }
    rows: list[dict[str, str]] = []
    for gene in panel:
        if gene in overlap_status:
            status = overlap_status[gene] or "unknown"
            source = "mathys_gene_overlap.tsv"
        elif gene in extracted_genes:
            status = "present"
            source = "extracted_sparse_expression"
        else:
            status = "missing"
            source = "panel_vs_extracted_sparse_expression"
        rows.append(
            {
                "gene_symbol": gene,
                "mathys_status": status,
                "evidence_source": source,
                "notes": "Phase 10 audit; present can be inferred from extracted sparse-like expression if overlap table is absent",
            }
        )
    return rows


def feature_overlap_rows(sea_rows: list[dict[str, str]], mathys_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    sea_features = feature_columns(sea_rows)
    mathys_features = feature_columns(mathys_rows)
    rows: list[dict[str, str]] = []
    for feature in sorted(sea_features | mathys_features):
        in_sea = feature in sea_features
        in_mathys = feature in mathys_features
        rows.append(
            {
                "feature": feature,
                "in_sea_ad": str(in_sea).lower(),
                "in_mathys": str(in_mathys).lower(),
                "status": "shared" if in_sea and in_mathys else "missing_from_mathys" if in_sea else "mathys_only",
                "notes": "shared features are eligible for external feasibility testing",
            }
        )
    return rows


def metric_reliability(metrics: list[dict[str, str]], predictions: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    warning_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
    prediction_counts = Counter(row.get("validation_mode", "missing") for row in predictions)
    for row in metrics:
        mode = row.get("validation_mode", "missing")
        n_test = int(float(row.get("n_test", "0") or 0))
        if n_test == 0:
            n_test = prediction_counts.get(mode, 0)
        auroc = to_float(row.get("auroc", "nan"))
        bal = to_float(row.get("balanced_accuracy", "nan"))
        reliable = n_test >= MIN_RELIABLE_TEST_N
        audit_rows.append(
            {
                "metric": f"{mode}_n_test",
                "value": str(n_test),
                "status": "reliable_size" if reliable else "small_n_unreliable",
                "notes": f"Minimum reviewer-facing reliability threshold is n_test >= {MIN_RELIABLE_TEST_N}.",
            }
        )
        if not reliable:
            warning_rows.append(
                {
                    "warning_id": f"{mode}_small_n",
                    "severity": "high",
                    "trigger": f"n_test={n_test}",
                    "recommendation": "Report Mathys as preliminary external feasibility evidence only.",
                }
            )
        contradiction = (
            not math.isnan(auroc)
            and not math.isnan(bal)
            and ((auroc >= 0.80 and bal <= 0.55) or (auroc <= 0.55 and bal >= 0.75))
        )
        audit_rows.append(
            {
                "metric": f"{mode}_auroc_balanced_accuracy_consistency",
                "value": f"AUROC={row.get('auroc', 'nan')}; balanced_accuracy={row.get('balanced_accuracy', 'nan')}",
                "status": "potential_contradiction" if contradiction else "consistent_or_not_testable",
                "notes": "Large AUROC/balanced-accuracy disagreement can indicate threshold instability or tiny sample effects.",
            }
        )
        if contradiction:
            warning_rows.append(
                {
                    "warning_id": f"{mode}_metric_contradiction",
                    "severity": "medium",
                    "trigger": f"AUROC={row.get('auroc', 'nan')}; balanced_accuracy={row.get('balanced_accuracy', 'nan')}",
                    "recommendation": "Discuss threshold-dependent instability and avoid binary validation claims.",
                }
            )
    return audit_rows, warning_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit Mathys external validation for reviewer-proof reporting.")
    parser.add_argument("--counts", type=Path, default=Path("data/raw/external/mathys_2019/GSE138852_counts.csv.gz"))
    parser.add_argument("--covariates", type=Path, default=Path("data/raw/external/mathys_2019/GSE138852_covariates.csv.gz"))
    parser.add_argument("--sparse-expression", type=Path, default=Path("data/interim/external/mathys_2019/mathys_sparse_gene_panel_expression.tsv.gz"))
    parser.add_argument("--mathys-features", type=Path, default=Path("results/tables/mathys_2019_phase5_donor_feature_table.tsv"))
    parser.add_argument("--sea-ad-features", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--panel", type=Path, default=Path("metadata/target_gene_panel_v1.tsv"))
    parser.add_argument("--gene-overlap", type=Path, default=Path("results/tables/mathys_gene_overlap.tsv"))
    parser.add_argument("--metrics", type=Path, default=Path("results/tables/phase9_mathys_external_validation_metrics.tsv"))
    parser.add_argument("--predictions", type=Path, default=Path("results/tables/phase9_mathys_external_predictions.tsv"))
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/46_audit_mathys_external_validation.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    mathys_features = read_tsv(args.mathys_features)
    sea_features = read_tsv(args.sea_ad_features)
    metrics = read_tsv(args.metrics)
    predictions = read_tsv(args.predictions)
    overlap = read_tsv(args.gene_overlap)
    panel = panel_genes(args.panel)
    sparse_rows, extracted_genes = count_sparse_rows_and_genes(args.sparse_expression)
    mathys_sample_units = len(mathys_features)
    mathys_cells = count_csv_rows(args.covariates)
    celltype_counts = celltype_counts_from_covariates(args.covariates)
    gene_rows = gene_audit_rows(panel, overlap, extracted_genes)
    feature_rows = feature_overlap_rows(sea_features, mathys_features)
    metric_rows, warning_rows = metric_reliability(metrics, predictions)

    found = sum(1 for row in gene_rows if row["mathys_status"] == "present")
    missing = sum(1 for row in gene_rows if row["mathys_status"] == "missing")
    shared_features = sum(1 for row in feature_rows if row["status"] == "shared")
    missing_features = sum(1 for row in feature_rows if row["status"] == "missing_from_mathys")
    audit_rows = [
        {"metric": "mathys_cells", "value": str(mathys_cells), "status": "observed", "notes": "Rows in covariates CSV."},
        {"metric": "mathys_inferred_sample_units", "value": str(mathys_sample_units), "status": "small_n_warning" if mathys_sample_units < SMALL_N_WARNING_THRESHOLD else "acceptable", "notes": "Feature-table rows used as sample-level units."},
        {"metric": "mathys_diagnosis_counts", "value": ";".join(f"{k}={v}" for k, v in sorted(diagnosis_counts(mathys_features).items())), "status": "observed", "notes": "Sample-level diagnosis/pathology labels."},
        {"metric": "mathys_cell_type_counts", "value": ";".join(f"{k}={v}" for k, v in sorted(celltype_counts.items())), "status": "observed", "notes": "Cell-type counts from covariates."},
        {"metric": "target_genes_requested", "value": str(len(panel)), "status": "observed", "notes": "NeuroFate target panel size."},
        {"metric": "target_genes_found", "value": str(found), "status": "observed", "notes": "From overlap table or extracted sparse expression fallback."},
        {"metric": "target_genes_missing", "value": str(missing), "status": "observed", "notes": "Panel genes not observed in Mathys overlap/extraction audit."},
        {"metric": "sparse_extracted_expression_rows", "value": str(sparse_rows), "status": "observed", "notes": "Rows in Mathys sparse-like target-gene expression table."},
        {"metric": "shared_feature_columns_with_sea_ad", "value": str(shared_features), "status": "observed", "notes": "Feature columns eligible for SEA-AD to Mathys transfer."},
        {"metric": "sea_ad_features_missing_from_mathys", "value": str(missing_features), "status": "observed", "notes": "Missing columns require omission or harmonization."},
        *metric_rows,
    ]
    if mathys_sample_units < SMALL_N_WARNING_THRESHOLD:
        warning_rows.append(
            {
                "warning_id": "mathys_sample_level_small_n",
                "severity": "high",
                "trigger": f"sample_units={mathys_sample_units}",
                "recommendation": "Treat Mathys as preliminary external feasibility evidence, not definitive external validation.",
            }
        )
    if found == 0 and sparse_rows > 0:
        warning_rows.append(
            {
                "warning_id": "gene_overlap_table_empty_but_expression_present",
                "severity": "medium",
                "trigger": f"sparse_rows={sparse_rows}; overlap_present={found}",
                "recommendation": "Regenerate gene overlap from extracted sparse expression or counts orientation audit.",
            }
        )

    write_tsv(args.tables_dir / "phase10_mathys_validation_audit.tsv", audit_rows, ["metric", "value", "status", "notes"])
    write_tsv(args.tables_dir / "phase10_mathys_gene_audit.tsv", gene_rows, ["gene_symbol", "mathys_status", "evidence_source", "notes"])
    write_tsv(args.tables_dir / "phase10_mathys_feature_overlap_audit.tsv", feature_rows, ["feature", "in_sea_ad", "in_mathys", "status", "notes"])
    write_tsv(args.tables_dir / "phase10_validation_warning_flags.tsv", warning_rows, ["warning_id", "severity", "trigger", "recommendation"])
    logging.info("Phase 10 Mathys validation audit complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
