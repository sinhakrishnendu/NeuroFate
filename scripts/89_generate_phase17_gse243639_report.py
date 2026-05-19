#!/usr/bin/env python3
"""Generate a conservative Phase 17 GSE243639 cell-type-aware PD report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def lines_for(rows: list[dict[str, str]], fields: list[str]) -> list[str]:
    if not rows:
        return ["- Not available."]
    output = []
    for row in rows[:10]:
        output.append("- " + ", ".join(f"{field}={row.get(field, '')}" for field in fields))
    return output


def build_report(output: Path, tables_dir: Path, reports_dir: Path) -> None:
    annotation = read_tsv(tables_dir / "phase17_gse243639_cell_annotation_summary.tsv")
    labels = read_tsv(tables_dir / "phase17_gse243639_celltype_label_summary.tsv")
    metrics = read_tsv(tables_dir / "phase17_gse243639_celltype_validation_metrics.tsv")
    comparison = read_tsv(tables_dir / "phase17_pd_validation_comparison.tsv")
    audit = read_tsv(reports_dir / "phase17_gse243639_umap_annotation_audit.tsv")
    reliability = next((row.get("reliability_flag", "") for row in metrics if row.get("model") == "logistic_regression"), "not_available")
    lines = [
        "# Phase 17 GSE243639 Cell-Type-Aware PD Report",
        "",
        "## 1. Why Phase 17 Was Needed",
        "Phase 16 established safe GSE243639 onboarding but the sample-level global target-gene signal was modest and unstable. Phase 17 tests whether existing cell annotations improve sample-level PD/control signal.",
        "",
        "## 2. UMAP/Annotation Inspection",
        "- The workbook is used only as an annotation and coordinate table.",
        "- No coordinate recomputation, clustering, or embedding analysis is performed.",
        *lines_for(audit, ["sheet_name", "column_name", "likely_role"]),
        "",
        "## 3. Cell-Type Annotation Mapping",
        "- Cell IDs are matched to the existing cell-to-sample map.",
        "- Candidate annotation columns are reported, and a conservative default cell-type column is selected.",
        *lines_for(annotation, ["sample_id", "cell_type", "cell_count"]),
        "",
        "## 4. Cell-Type-Aware Feature Construction",
        "- Features remain sample-level aggregates.",
        "- Feature groups include global gene means, global detection rates, cell fractions, cell-type gene means, cell-type detection rates, and target-gene indices.",
        *lines_for(labels, ["label_field", "label", "sample_count"]),
        "",
        "## 5. PD/Control Validation",
        *lines_for(metrics, ["model", "validation_mode", "auroc", "auroc_sd", "balanced_accuracy", "empirical_permutation_pvalue", "reliability_flag"]),
        "",
        "## 6. Permutation And Robustness Results",
        "- Permutation-label controls are used to compare observed AUROC against shuffled-label baselines.",
        "- Bootstrap intervals are reported when feasible from leave-one-out probabilities.",
        "",
        "## 7. Comparison With Phase 16",
        *lines_for(comparison, ["phase", "feature_count", "auroc", "balanced_accuracy", "reliability_flag"]),
        "",
        "## 8. Biological Interpretation",
        "- Cell-type-aware PD signal should be interpreted as exploratory unless it is stable, above permutation null, and biologically coherent.",
        "- Stronger signal would support the value of cell-state context in NeuroFate, not medical deployment.",
        "",
        "## 9. Reliability Category",
        f"- Current reliability category: {reliability}.",
        "",
        "## 10. Limitations And Next Steps",
        "- GSE243639 contains 29 sample-level units, so uncertainty remains substantial.",
        "- Additional independent PD cohorts and harmonized cell-type labels are needed before broad disease-transfer language is appropriate.",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 17 GSE243639 cell-type-aware report.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--reports-dir", type=Path, default=Path("results/reports"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase17_gse243639_celltype_pd_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    build_report(args.output, args.tables_dir, args.reports_dir)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
