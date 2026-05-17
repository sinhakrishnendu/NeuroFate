#!/usr/bin/env python3
"""Audit donor-level NeuroFate feature tables for predictor leakage risks."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FEATURE_PREFIXES = (
    "gene_mean__",
    "gene_detection__",
    "index__",
    "cell_fraction__",
    "celltype_index__",
    "mean_",
)
IDENTIFIER_COLUMNS = {"donor_id", "sample_id", "cell_id"}
LABEL_PREFIX = "label__"
LABEL_TERMS = (
    "diagnosis",
    "pathology",
    "cognitive",
    "apoe_genotype",
    "braak",
    "cerad",
    "neuropathological",
    "dementia",
    "disease_status",
    "neurotypical_reference",
    "highest_lewy_body_disease",
    "late",
    "overall_caa",
)


def read_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        return next(reader, [])


def classify_column(column: str) -> dict[str, str]:
    lowered = column.lower()
    is_feature = column.startswith(FEATURE_PREFIXES)
    if column in IDENTIFIER_COLUMNS:
        return {
            "column_name": column,
            "column_role": "identifier",
            "leakage_risk": "high",
            "reason": "identifier columns must not be predictors",
            "recommended_action": "exclude from predictor matrix",
        }
    if column == "cohort_id" or lowered.endswith("cohort_id"):
        return {
            "column_name": column,
            "column_role": "cohort",
            "leakage_risk": "medium",
            "reason": "cohort identifiers can leak external validation source labels",
            "recommended_action": "use only for explicit stratification or reporting",
        }
    if column.startswith(LABEL_PREFIX):
        return {
            "column_name": column,
            "column_role": "label",
            "leakage_risk": "high",
            "reason": "label columns are task targets, not predictors",
            "recommended_action": "derive task labels before feature selection, then drop",
        }
    if is_feature and any(term in lowered for term in LABEL_TERMS):
        return {
            "column_name": column,
            "column_role": "feature",
            "leakage_risk": "medium",
            "reason": "feature-like column contains label-associated wording",
            "recommended_action": "review manually before modeling",
        }
    if is_feature:
        return {
            "column_name": column,
            "column_role": "predictor",
            "leakage_risk": "low",
            "reason": "column matches approved feature prefixes",
            "recommended_action": "allowed as predictor if numeric",
        }
    if column == "n_cells":
        return {
            "column_name": column,
            "column_role": "quality_covariate",
            "leakage_risk": "medium",
            "reason": "cell count may capture sampling or batch effects",
            "recommended_action": "use only after sensitivity analysis",
        }
    return {
        "column_name": column,
        "column_role": "unknown",
        "leakage_risk": "medium",
        "reason": "column does not match approved predictor prefixes",
        "recommended_action": "exclude unless explicitly justified",
    }


def audit_feature_table(input_path: Path, output_path: Path) -> list[dict[str, str]]:
    if not input_path.exists():
        raise FileNotFoundError(f"Missing donor feature table: {input_path}")
    rows = [classify_column(column) for column in read_header(input_path)]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "column_name",
                "column_role",
                "leakage_risk",
                "reason",
                "recommended_action",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit NeuroFate donor-level predictors for leakage.")
    parser.add_argument("--input", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/feature_leakage_audit.tsv"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        rows = audit_feature_table(args.input, args.output)
    except FileNotFoundError as exc:
        print(f"Leakage audit could not start: {exc}")
        print("Run Phase 5 donor feature table generation first or provide --input.")
        return 2
    high = sum(1 for row in rows if row["leakage_risk"] == "high")
    medium = sum(1 for row in rows if row["leakage_risk"] == "medium")
    print(f"Wrote {args.output}")
    print(f"High-risk columns: {high}")
    print(f"Medium-risk columns: {medium}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
