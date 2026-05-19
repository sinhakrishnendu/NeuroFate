#!/usr/bin/env python3
"""Generate Phase 26 GSE184950 axis-replication report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_by(rows: list[dict[str, str]], column: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = row.get(column, "") or "missing"
        counts[value] = counts.get(value, 0) + 1
    return counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate GSE184950 Phase 26 replication report.")
    parser.add_argument("--metadata", type=Path, default=Path("results/tables/phase25_gse184950_series_sample_metadata.tsv"))
    parser.add_argument("--nested-inventory", type=Path, default=Path("results/tables/phase26_gse184950_nested_archive_inventory.tsv"))
    parser.add_argument("--extraction-audit", type=Path, default=Path("results/tables/phase26_gse184950_selected_extraction_audit.tsv"))
    parser.add_argument("--axis-audit", type=Path, default=Path("results/tables/phase26_gse184950_axis_gene_extraction_audit.tsv"))
    parser.add_argument("--axis-scores", type=Path, default=Path("results/tables/phase25_gse184950_axis_scores.tsv"))
    parser.add_argument("--replication-stats", type=Path, default=Path("results/tables/phase25_gse184950_axis_replication_statistics.tsv"))
    parser.add_argument("--readiness", type=Path, default=Path("results/reports/phase23_pnas_readiness_matrix.tsv"))
    parser.add_argument("--phase27-clean-stats", type=Path, default=Path("results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase26_gse184950_axis_replication_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metadata = read_tsv(args.metadata)
    nested = read_tsv(args.nested_inventory)
    extraction = read_tsv(args.extraction_audit)
    axis_audit = read_tsv(args.axis_audit)
    scores = read_tsv(args.axis_scores)
    stats = read_tsv(args.replication_stats)
    phase27_stats = read_tsv(args.phase27_clean_stats)
    disease_counts = count_by(metadata, "disease_state")
    complete_samples = sorted({row.get("sample_id", "") for row in nested if row.get("complete_processed_matrix_set") == "true" and row.get("sample_id")})
    extracted_samples = sorted({row.get("sample_id", "") for row in extraction if row.get("status") == "extracted_processed_matrix_file_only"})
    scored_samples = sorted({row.get("sample_id", "") for row in scores if row.get("sample_id")})

    lines = [
        "# Phase 26 GSE184950 Axis Replication Report",
        "",
        "## 1. Dataset Summary",
        f"- Series-matrix samples: {len(metadata)}",
        "- Cohort role: independent PD/PDD-vs-control replication cohort.",
        "",
        "Disease-state counts:",
    ]
    lines.extend(f"- {label}: {count}" for label, count in sorted(disease_counts.items()))
    lines.extend(
        [
            "",
            "## 2. Series Matrix Metadata",
            "Phase 26 uses the Phase 25 series-matrix metadata rather than the incomplete add2 workbook.",
            "",
            "## 3. Nested Archive Contents",
            f"- Nested inventory rows: {len(nested)}",
            f"- Samples with complete processed matrix sets: {len(complete_samples)}",
            "",
            "## 4. Selective Extraction Audit",
            f"- Selected matrix-file extraction audit rows: {len(extraction)}",
            f"- Samples with extracted processed matrix files: {len(extracted_samples)}",
            "",
            "## 5. Axis Gene Coverage",
            f"- Axis extraction audit samples: {len(axis_audit)}",
            "",
            "## 6. Endpoint Definition",
            "`label__pd_pdd_vs_control` compares Parkinson's Disease plus Parkinson's Disease Dementia against Unaffected Control at sample level.",
            "",
            "## 7. PD/PDD vs Control Replication Statistics",
            f"- Axis-score samples: {len(scored_samples)}",
            f"- Replication statistic rows: {len(stats)}",
            "",
            "## 8. Comparison With Phase 22 Candidate Axes",
            "Use the `phase22_pd_direction` and `directional_replication_status` fields in the replication statistics table to assess direction consistency.",
            "",
            "## 9. PNAS Readiness Implication",
            "GSE184950 can strengthen the PD replication layer only if processed matrices are available, axis coverage is adequate, and endpoint-locked effects are directionally consistent. Claims remain candidate or preliminary until evidence is replicated and statistically stable.",
            "",
            "## 10. Remaining Limitations",
            "- This is not a clinical validation or diagnostic workflow.",
            "- No causal mechanism is established.",
            "- Raw sequence preprocessing is intentionally avoided.",
        ]
    )
    if phase27_stats:
        supported = sum(row.get("evidence_label") == "replicated_statistically_supported" for row in phase27_stats)
        lines.extend(
            [
                "",
                "## Phase 27 Clean Replication Update",
                f"- Clean Phase 27 replication statistic rows: {len(phase27_stats)}",
                f"- Statistically supported replicated axes: {supported}",
                "- Phase 27 clean tables supersede Phase 26/25 tables for GSE184950 interpretation.",
            ]
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
