#!/usr/bin/env python3
"""Generate the Phase 20 GSE243639 safe cell-type PD report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bullets(rows: list[dict[str, str]], fields: list[str], limit: int = 8) -> list[str]:
    if not rows:
        return ["- Not available."]
    return ["- " + ", ".join(f"{field}={row.get(field, '')}" for field in fields) for row in rows[:limit]]


def build_report(args: argparse.Namespace) -> None:
    schema = read_tsv(args.schema_audit)
    groups = read_tsv(args.feature_groups)
    metrics = read_tsv(args.metrics)
    comparison = read_tsv(args.comparison)
    decision = read_tsv(args.phase19_decision)
    lines = [
        "# Phase 20 GSE243639 Safe Cell-Type PD Report",
        "",
        "## 1. Why Phase 20 Was Required",
        "Phase 19 identified a safe normalized ID linkage, but the feature builder still reported zero annotation matching. Phase 20 repairs the interface between the safe annotation map and the sample-level feature builder.",
        "",
        "## 2. Phase 19 Safe Linkage Status",
        *bullets(decision, ["decision_category", "direct_or_normalized_best_rule", "best_overlap_rate", "safe_to_build_annotation_map"]),
        "",
        "## 3. Safe Annotation Map Schema Audit",
        *bullets(schema, ["audit_item", "value", "recommended_join_column", "overlap_rate"], limit=12),
        "",
        "## 4. Corrected Feature Construction",
        *bullets(groups, ["feature_group", "feature_count", "annotation_match_rate", "unmatched_unique_expression_cells", "warning"], limit=12),
        "",
        "## 5. PD/Control Validation",
        *bullets(metrics, ["model", "validation_mode", "auroc", "auroc_sd", "balanced_accuracy", "empirical_permutation_pvalue", "reliability_flag"]),
        "",
        "## 6. Phase 16/17/18/20 Comparison",
        *bullets(comparison, ["phase", "feature_count", "annotation_match_rate", "auroc", "reliability_flag", "interpretation"]),
        "",
        "## 7. Interpretation",
        "- Phase 16 remains the valid global sample-level PD extension.",
        "- Phase 17 was technically defective and should not be used for biological interpretation.",
        "- Phase 18 showed the annotation join was still not being consumed correctly.",
        "- Phase 20 should be used for final cell-type-aware PD interpretation only if annotation match rate is high and validation reliability is not technical failure.",
        "",
        "## 8. Limitations",
        "- GSE243639 contains 29 sample-level units, so uncertainty remains material.",
        "- This is research software output, not medical validation.",
        "- The analysis does not establish cause-and-effect biology.",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 20 GSE243639 report.")
    parser.add_argument("--schema-audit", type=Path, default=Path("results/tables/phase20_safe_annotation_map_schema_audit.tsv"))
    parser.add_argument("--feature-groups", type=Path, default=Path("results/tables/phase20_gse243639_feature_group_counts.tsv"))
    parser.add_argument("--metrics", type=Path, default=Path("results/tables/phase20_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--comparison", type=Path, default=Path("results/tables/phase20_pd_validation_comparison.tsv"))
    parser.add_argument("--phase19-decision", type=Path, default=Path("results/tables/phase19_gse243639_annotation_linkage_decision.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase20_gse243639_safe_celltype_pd_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    build_report(args)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
