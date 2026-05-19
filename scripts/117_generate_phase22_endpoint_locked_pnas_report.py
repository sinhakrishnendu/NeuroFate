#!/usr/bin/env python3
"""Generate the Phase 22 endpoint-locked PNAS axis report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize_registry(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "Endpoint registry is unavailable."
    return "\n".join(f"- `{row['endpoint_id']}`: `{row['source_column']}` ({row['endpoint_role']}, {row['endpoint_type']})." for row in rows)


def summarize_claims(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "Endpoint-locked evidence table has not yet been generated."
    candidates = [row for row in rows if "candidate" in row.get("axis_claim_class", "")]
    if not candidates:
        return "No endpoint-locked candidate axes are currently strong enough for PNAS-facing biological claims."
    return "\n".join(
        f"- `{row['axis_id']}` at `{row['endpoint_id']}`: {row['axis_claim_class']} "
        f"(effect={row['effect_size']}, FDR={row['fdr'] or 'NA'}, empirical p={row['empirical_pvalue'] or 'NA'})."
        for row in candidates[:20]
    )


def build_report(args: argparse.Namespace) -> str:
    endpoints = read_tsv(args.endpoint_registry)
    evidence = read_tsv(args.tables_dir / "phase22_endpoint_locked_axis_evidence_table.tsv")
    return f"""# Phase 22 Endpoint-Locked PNAS Axis Report

## 1. Why Phase 22 Was Required

Phase 21 created the NeuroFate-Axis framework, but its exploratory comparison selected the largest absolute effect across heterogeneous labels. Phase 22 supersedes that approach for PNAS-facing biological claims by locking endpoints before testing.

## 2. Endpoint Registry

{summarize_registry(endpoints)}

## 3. Primary AD Endpoint

The primary AD endpoint is `sea_ad_cognitive_dementia`, comparing Dementia versus No dementia and excluding Reference/missing values.

## 4. Primary PD Endpoint

The primary PD endpoint is `gse243639_pd_diagnosis`, comparing Parkinson's versus Control in the independent GSE243639 cohort.

## 5. Endpoint-Locked Axis Associations

Only endpoints declared in `metadata/neurofate_axis_endpoint_registry.tsv` are tested. Arbitrary labels are not scanned.

## 6. AD/PD Comparison

Primary AD/PD comparison uses `sea_ad_cognitive_dementia` versus `gse243639_pd_diagnosis`. Secondary comparisons are labelled explicitly and cannot be silently mixed into primary claims.

## 7. Matched Random-Axis Controls

Random-axis controls use the same endpoint and same association statistic as each curated axis, with random feature sets matched by feature count.

## 8. Candidate Shared Axes

{summarize_claims([row for row in evidence if "shared" in row.get("axis_claim_class", "")])}

## 9. Candidate Disease-Enriched Axes

{summarize_claims([row for row in evidence if "disease_enriched" in row.get("axis_claim_class", "")])}

## 10. What Can Be Claimed Now

Only endpoint-locked candidate or preliminary axis associations can be claimed, and only at donor/sample level.

## 11. What Remains Insufficient For PNAS

PNAS-grade claims still require independent replication, empirical random-control support, FDR-aware interpretation, and no-overclaiming audit clearance.

## 12. Next Validation Cohort Priority

The next priority is one additional donor/sample-level PD cohort to replicate GSE243639 axis patterns, followed by one larger AD external cohort beyond Mathys feasibility.
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 22 endpoint-locked PNAS report.")
    parser.add_argument("--endpoint-registry", type=Path, default=Path("metadata/neurofate_axis_endpoint_registry.tsv"))
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase22_endpoint_locked_pnas_axis_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(build_report(args), encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
