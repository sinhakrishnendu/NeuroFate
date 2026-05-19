#!/usr/bin/env python3
"""Generate a conservative Phase 16 GSE243639 PD validation report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def table_lines(rows: list[dict[str, str]], fields: list[str]) -> list[str]:
    if not rows:
        return ["- Not available."]
    output = []
    for row in rows:
        output.append("- " + ", ".join(f"{field}={row.get(field, '')}" for field in fields))
    return output


def build_report(output: Path, tables_dir: Path) -> None:
    audit = read_tsv(tables_dir / "phase16_gse243639_gene_extraction_audit.tsv")
    labels = read_tsv(tables_dir / "phase16_gse243639_label_summary.tsv")
    schema = read_tsv(tables_dir / "phase16_gse243639_feature_schema_alignment.tsv")
    metrics = read_tsv(tables_dir / "phase16_gse243639_external_validation_metrics.tsv")
    reliability = next((row.get("reliability_flag", "") for row in metrics if row.get("validation_mode") == "gse243639_pd_internal"), "not_available")
    shared = sum(1 for row in schema if row.get("status") == "shared")
    lines = [
        "# Phase 16 GSE243639 PD Validation Report",
        "",
        "This report summarizes an independent PD cohort extension using sample-level NeuroFate features. It is not medical validation and does not establish cause-and-effect biology.",
        "",
        "## 1. Dataset Summary",
        "- Dataset: GSE243639 human substantia nigra pars compacta single-nucleus RNA-seq.",
        "- Planned sample units: 29 postmortem brains with Parkinson's/control labels.",
        *table_lines(labels, ["label_field", "label", "sample_count"]),
        "",
        "## 2. Clinical Metadata Parsing",
        "- Clinical metadata are parsed with a semicolon delimiter and a 1-based header line of 6.",
        "- Canonical fields include sample ID, diagnosis, age, sex, PMI, RIN, Lewy body fields, CERAD, and Braak.",
        "",
        "## 3. Count Matrix Orientation",
        "- The count table is treated as genes-as-rows and nuclei/cells-as-columns.",
        "- Sample IDs are derived from the cell ID prefix before the underscore.",
        "",
        "## 4. Target Gene Extraction",
        *table_lines(audit, ["requested_target_genes", "extracted_target_genes", "missing_target_genes", "cell_columns", "sparse_expression_rows"]),
        "",
        "## 5. Sample-Level Feature Construction",
        f"- Shared SEA-AD/GSE243639 schema fields: {shared}.",
        "- Features are aggregated to sample level before any validation step.",
        "",
        "## 6. PD/Control Validation",
        *table_lines(metrics, ["validation_mode", "task_id", "auroc", "auroc_sd", "balanced_accuracy", "n_samples", "reliability_flag"]),
        "",
        "## 7. Cross-Disease Interpretation",
        "- SEA-AD is an Alzheimer disease cohort and GSE243639 is a Parkinson disease cohort.",
        "- Cross-disease output is treated as feature-space transfer feasibility, not direct disease-label transfer.",
        "",
        "## 8. Reliability Category",
        f"- Current reliability category: {reliability}.",
        "",
        "## 9. Claim-Strength Implications",
        "- GSE243639 can strengthen the platform as an independent PD cohort extension when sample-level metrics are stable.",
        "- It should not be used to imply clinical deployment or direct AD-to-PD disease prediction.",
        "",
        "## 10. Limitations",
        "- Cohort-specific preprocessing, cell-type annotation compatibility, and target-gene overlap remain important constraints.",
        "- Larger independent cohorts are still needed before broad generalization language is appropriate.",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 16 GSE243639 report.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase16_gse243639_pd_validation_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    build_report(args.output, args.tables_dir)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
