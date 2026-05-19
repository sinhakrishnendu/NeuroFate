#!/usr/bin/env python3
"""Generate Phase 27 clean GSE184950 replication report with conservative interpretation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def value(row: dict[str, str], key: str, default: str = "0") -> str:
    return row.get(key, "") or default


def count_labels(rows: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        label = row.get("label__pd_pdd_vs_control", "") or "missing"
        counts[label] = counts.get(label, 0) + 1
    return counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 27 clean GSE184950 replication report.")
    parser.add_argument("--sample-integrity", type=Path, default=Path("results/tables/phase27_gse184950_sample_integrity_audit.tsv"))
    parser.add_argument("--axis-scores", type=Path, default=Path("results/tables/phase27_gse184950_axis_scores_clean.tsv"))
    parser.add_argument("--axis-coverage", type=Path, default=Path("results/tables/phase27_gse184950_axis_feature_coverage_clean.tsv"))
    parser.add_argument("--label-summary", type=Path, default=Path("results/tables/phase27_gse184950_axis_label_summary_clean.tsv"))
    parser.add_argument("--replication-stats", type=Path, default=Path("results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv"))
    parser.add_argument("--readiness", type=Path, default=Path("results/reports/phase23_pnas_readiness_matrix.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase27_gse184950_clean_replication_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    integrity = read_tsv(args.sample_integrity)
    integrity_row = integrity[0] if integrity else {}
    scores = read_tsv(args.axis_scores)
    coverage = read_tsv(args.axis_coverage)
    labels = read_tsv(args.label_summary)
    stats = read_tsv(args.replication_stats)
    readiness = read_tsv(args.readiness)
    label_counts = count_labels(scores)
    supported = [row for row in stats if row.get("evidence_label") == "replicated_statistically_supported"]
    directional = [row for row in stats if row.get("evidence_label") == "directionally_consistent_but_not_significant"]
    weak = [row for row in stats if row.get("evidence_label") in {"weak_or_no_replication", "opposite_direction"}]
    pnas_pd = next((row for row in readiness if row.get("criterion") == "independent_pd_replication"), {})

    lines = [
        "# Phase 27 GSE184950 Clean Replication Report",
        "",
        "## 1. Sample-Integrity Audit",
        f"- Expected samples: {value(integrity_row, 'expected_samples')}",
        f"- Observed clean axis-score samples: {len(scores)}",
        f"- Invalid axis-score samples before cleaning: {value(integrity_row, 'invalid_axis_score_samples')}",
        f"- Invalid sample IDs: {value(integrity_row, 'invalid_sample_ids', '')}",
        "",
        "## 2. Clean Sample Count",
        f"- Clean biological samples: {len(scores)}",
        "",
        "## 3. Axis Coverage",
        f"- Axis rows: {len(coverage)}",
        f"- Fully covered axes: {sum(row.get('status') == 'ok' for row in coverage)}",
        "",
        "## 4. Endpoint Definition",
        "`label__pd_pdd_vs_control` compares Parkinson's Disease plus Parkinson's Disease Dementia against Unaffected Control.",
        "",
        "Clean label counts:",
    ]
    lines.extend(f"- {label}: {count}" for label, count in sorted(label_counts.items()))
    lines.extend(
        [
            "",
            "## 5. Clean Replication Statistics",
            f"- Replication statistic rows: {len(stats)}",
            f"- Statistically supported replicated axes: {len(supported)}",
            f"- Directionally consistent but not significant axes: {len(directional)}",
            f"- Weak, absent, or opposite evidence axes: {len(weak)}",
            "",
            "## 6. Directional Consistency With Phase 22",
            "Direction-only agreement is reported as a preliminary signal and is not treated as replication unless p < 0.05 or FDR < 0.1.",
            "",
            "## 7. Statistical Support",
            "GSE184950 establishes the independent PD/PDD replication infrastructure. If FDR remains weak, it does not establish strong replicated axis biology.",
            "",
            "## 8. PNAS Implication",
            f"- Independent PD replication readiness: {pnas_pd.get('status', 'not_available')}",
            "PNAS-level shared-axis claims require statistically supported independent replication, not direction-only agreement.",
            "",
            "## 9. Conservative Interpretation",
            "GSE184950 currently supports independent PD replication feasibility and sample-level axis testing. Do not claim a validated AD/PD shared mechanism, clinical biomarker, diagnostic axis, or causal mechanism from weak FDR support.",
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
