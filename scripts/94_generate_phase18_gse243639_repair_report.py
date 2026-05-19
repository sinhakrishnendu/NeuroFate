#!/usr/bin/env python3
"""Generate the Phase 18 GSE243639 cell-annotation repair report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def lines_for(rows: list[dict[str, str]], fields: list[str], limit: int = 8) -> list[str]:
    if not rows:
        return ["- Not available."]
    return ["- " + ", ".join(f"{field}={row.get(field, '')}" for field in fields) for row in rows[:limit]]


def build_report(output: Path, tables_dir: Path, reports_dir: Path) -> None:
    id_audit = read_tsv(tables_dir / "phase18_gse243639_cell_id_matching_audit.tsv")
    match_summary = read_tsv(tables_dir / "phase18_gse243639_annotation_match_summary.tsv")
    candidates = read_tsv(reports_dir / "phase18_gse243639_annotation_column_candidates.tsv")
    feature_groups = read_tsv(tables_dir / "phase18_gse243639_feature_group_counts.tsv")
    metrics = read_tsv(tables_dir / "phase18_gse243639_celltype_validation_metrics.tsv")
    comparison = read_tsv(tables_dir / "phase18_pd_validation_comparison.tsv")
    phase19_decision = read_tsv(tables_dir / "phase19_gse243639_annotation_linkage_decision.tsv")
    phase19_category = phase19_decision[0].get("decision_category", "not_available") if phase19_decision else "not_available"
    phase19_safe = phase19_decision[0].get("safe_to_build_annotation_map", "false") if phase19_decision else "false"
    reliability = next((row.get("reliability_flag", "") for row in metrics if row.get("model") == "logistic_regression"), "not_available")
    phase19_lines = [
        f"- Phase 19 annotation-linkage status: {phase19_category}.",
        f"- Safe to build workbook-derived annotation map: {phase19_safe}.",
    ]
    if phase19_safe != "true" and phase19_category != "not_available":
        phase19_lines.append(
            "- Cell-type-aware GSE243639 validation is not currently supported from the available workbook because cell IDs cannot be safely linked to expression cells."
        )
    lines = [
        "# Phase 18 GSE243639 Cell-Type Repair Report",
        "",
        "## 1. Why Phase 18 Was Required",
        "Phase 17 produced too few cell-type-aware features and counted sparse expression rows as unmatched observations. Phase 18 treats that result as a technical audit failure rather than a biological finding.",
        "",
        "## 2. Cell-ID Mismatch Audit",
        *lines_for(id_audit, ["metric", "value", "recommended_normalization_rule"]),
        "",
        "## 3. Annotation Mapping Repair",
        "- The repaired mapper starts from expression/cell-map IDs and joins annotation rows back to those expression IDs.",
        "- It preserves original annotation IDs, normalized IDs, barcode cores, sample IDs, and match status.",
        *lines_for(match_summary, ["total_expression_cells", "matched_expression_cells", "unmatched_expression_cells", "match_rate", "warning"]),
        "",
        "## 4. Cell-Type/Cluster Label Confidence",
        "- Biological labels are preferred when available; numeric clusters are retained as cluster-labeled features with low confidence.",
        *lines_for(candidates, ["column_name", "annotation_kind", "biological_label_hits", "numeric_fraction", "confidence"]),
        "",
        "## 5. Repaired Feature Construction",
        "- Global Phase 16 gene mean and detection features are preserved.",
        "- Repaired cell fractions and cell-type/cluster gene features are added only after annotation matching.",
        *lines_for(feature_groups, ["feature_group", "feature_count", "annotation_match_rate", "warning"]),
        "",
        "## 6. Repaired PD/Control Validation",
        *lines_for(metrics, ["model", "validation_mode", "auroc", "auroc_sd", "balanced_accuracy", "empirical_permutation_pvalue", "reliability_flag"]),
        "",
        "## 7. Comparison With Phase 16 And Phase 17",
        *lines_for(comparison, ["phase", "feature_count", "annotation_match_rate", "auroc", "reliability_flag", "interpretation"]),
        "",
        "## 8. Biological Interpretation",
        "- Phase 18 should supersede Phase 17 only if annotation matching succeeds and feature counts recover.",
        "- Weak or failed repaired results should remain technical evidence, not a disease biology conclusion.",
        "",
        "## 9. Reliability Category",
        f"- Current repaired reliability category: {reliability}.",
        "",
        "## Phase 19 Annotation-Linkage Status",
        *phase19_lines,
        "",
        "## 10. Remaining Limitations",
        "- GSE243639 has 29 sample-level units, so uncertainty remains material.",
        "- Additional PD cohorts and clearer biological cell-type labels are needed before stronger generalization language is appropriate.",
        "- The workflow is research-only and does not establish cause-and-effect biology.",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 18 GSE243639 repair report.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--reports-dir", type=Path, default=Path("results/reports"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase18_gse243639_celltype_repair_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    build_report(args.output, args.tables_dir, args.reports_dir)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
