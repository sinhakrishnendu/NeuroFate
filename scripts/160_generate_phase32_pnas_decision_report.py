#!/usr/bin/env python3
"""Generate a conservative Phase 32 PNAS decision report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_by(rows: list[dict[str, str]], column: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        key = row.get(column, "missing") or "missing"
        counts[key] = counts.get(key, 0) + 1
    return counts


def strongest_axis(rows: list[dict[str, str]]) -> dict[str, str]:
    for row in rows:
        if row.get("crosscohort_evidence_class") == "strong_ad_axis_with_nominal_external_replication":
            return row
    return rows[0] if rows else {}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 32 PNAS decision report.")
    parser.add_argument("--ranked", type=Path, default=Path("results/tables/phase32_axis_evidence_ranked.tsv"))
    parser.add_argument("--summary", type=Path, default=Path("results/tables/phase32_crosscohort_axis_evidence_summary.tsv"))
    parser.add_argument("--readiness", type=Path, default=Path("results/reports/phase23_pnas_readiness_matrix.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase32_pnas_decision_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ranked = read_tsv(args.ranked)
    summary = read_tsv(args.summary)
    readiness = read_tsv(args.readiness)
    strongest = strongest_axis(ranked)
    counts = count_by(summary, "crosscohort_evidence_class")
    readiness_status = {row.get("criterion", ""): row.get("status", "") for row in readiness}
    lines = [
        "# Phase 32 PNAS Decision Report",
        "",
        "## 1. Executive Summary",
        "NeuroFate-Axis now has independent AD replication evidence from GSE174367 bulk RNA. The support is nominal and endpoint-locked, not FDR-robust.",
        "",
        "## 2. SEA-AD Discovery Evidence",
        "SEA-AD remains the primary AD discovery anchor with endpoint-locked dementia-axis evidence from Phase 22.",
        "",
        "## 3. GSE174367 Independent AD Replication",
        "GSE174367 provides independent AD-vs-Control bulk RNA replication evidence. The neuronal vulnerability axis is the strongest replicated AD candidate.",
        f"Strongest axis: `{strongest.get('axis_id', 'unavailable')}`; GSE174367 p={strongest.get('gse174367_p', '')}; FDR={strongest.get('gse174367_fdr', '')}.",
        "",
        "## 4. GSE243639 PD Extension",
        "GSE243639 Phase 20 remains a preliminary PD extension with sample-level cell/cluster-aware signal. Phase 22 also provides endpoint-locked axis statistics for the same cohort; for the neuronal vulnerability axis these are same-direction with SEA-AD but not statistically strong after FDR correction.",
        f"GSE243639 axis effect for `{strongest.get('axis_id', 'unavailable')}`: effect={strongest.get('gse243639_axis_effect', '')}; p={strongest.get('gse243639_axis_p', '')}; FDR={strongest.get('gse243639_axis_fdr', '')}; empirical p={strongest.get('gse243639_axis_empirical_p', '')}.",
        "",
        "## 5. GSE184950 PD Replication Attempt",
        "GSE184950 provides clean independent PD/PDD replication infrastructure, but current axis-level effects are weak or direction-only.",
        "",
        "## 6. Cross-Cohort Axis Ranking",
        f"Evidence-class counts: {counts}.",
        "",
        "## 7. What Can Now Be Claimed",
        "NeuroFate-Axis supports nominal independent AD replication for the neuronal vulnerability axis and preliminary same-direction PD convergence where supported. This is an AD-replicated candidate axis with exploratory PD convergence, not a confirmed shared-disease mechanism.",
        "",
        "## 8. What Cannot Be Claimed",
        "Do not claim clinical use, diagnostic performance, causality, definitive shared AD/PD biology, or publication-ready mechanism-level evidence.",
        "",
        "## 9. Why PNAS Is Closer But Not Yet Complete",
        f"Independent AD replication status: {readiness_status.get('independent_ad_replication', 'unknown')}. Shared AD/PD axis claim status: {readiness_status.get('shared_ad_pd_axis_claim', 'unknown')}.",
        "PNAS-level shared AD/PD claims still require stronger PD axis replication and ideally FDR-robust support.",
        "",
        "## 10. Next Validation Priority",
        "The next priority is stronger independent PD axis replication, either through improved GSE184950 signal, another PD cohort, or a safer cell-type/pathology metadata route.",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
