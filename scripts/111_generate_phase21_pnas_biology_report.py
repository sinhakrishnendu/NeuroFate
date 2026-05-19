#!/usr/bin/env python3
"""Generate a conservative PNAS-oriented NeuroFate biological discovery report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize_axis_claims(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "Phase 21 axis claim-strength results are not yet available."
    lines = []
    for row in rows[:12]:
        lines.append(
            f"- `{row.get('axis_id', 'axis')}`: {row.get('axis_classification', 'inconclusive_axis')} "
            f"with {row.get('claim_strength', 'axis_level_insufficient_validation')}."
        )
    return "\n".join(lines)


def summarize_coverage(rows: list[dict[str, str]], cohort: str) -> str:
    selected = [row for row in rows if row.get("cohort") == cohort]
    if not selected:
        return f"No {cohort} axis coverage table is available yet."
    ok = sum(1 for row in selected if row.get("status") == "ok")
    return f"{ok}/{len(selected)} axes have at least one available donor/sample-level feature in `{cohort}`."


def build_report(args: argparse.Namespace) -> str:
    coverage = read_tsv(args.tables_dir / "phase21_axis_feature_coverage.tsv")
    claims = read_tsv(args.tables_dir / "phase21_axis_claim_strength.tsv")
    random_controls = read_tsv(args.tables_dir / "phase21_axis_empirical_pvalues.tsv")
    random_status = (
        f"Random-axis controls are available for {len(random_controls)} cohort-axis pairs."
        if random_controls
        else "Random-axis controls have not yet been run."
    )
    return f"""# Phase 21 PNAS Biological Discovery Report

## 1. Biological Question

NeuroFate asks whether Alzheimer disease and Parkinson disease share donor-level neurodegeneration fate axes that capture conserved glial-inflammatory, myelin, and neuronal vulnerability programs while diverging in amyloid/tau- and synuclein-associated structure.

## 2. Cohorts Analyzed

- SEA-AD is the internal AD anchor cohort.
- Mathys 2019 remains preliminary AD external feasibility because the harmonized sample count is small.
- GSE243639 is an independent PD cohort. Phase 20 provides the corrected safe-map cell/cluster-aware sample-level PD feature table and remains preliminary.

## 3. NeuroFate Axis Definitions

Axes are curated donor/sample-level biological summaries. They are not clinical biomarkers and they do not establish disease causality.

## 4. SEA-AD AD Axis Findings

{summarize_coverage(coverage, "sea_ad")}

## 5. GSE243639 PD Axis Findings

{summarize_coverage(coverage, "gse243639_pd_snpc")}

## 6. Shared AD/PD Axes

Shared axes should be described as candidates only when AD and PD show same-direction donor/sample-level association and sufficient feature coverage.

{summarize_axis_claims([row for row in claims if row.get("axis_classification") == "shared_ad_pd_candidate"])}

## 7. Disease-Specific Axes

Disease-specific axes should be described as preliminary when one disease shows stronger association and the other remains weak or inconclusive.

{summarize_axis_claims([row for row in claims if "enriched" in row.get("axis_classification", "")])}

## 8. Random-Axis Controls

{random_status} Curated axes should not be interpreted as robust until they outperform random feature sets with conservative empirical support.

## 9. What Is Robust

The robust contribution is the reproducible, memory-safe, donor/sample-level framework for defining interpretable neurodegeneration axes from existing NeuroFate outputs.

## 10. What Remains Preliminary

The Phase 20 GSE243639 PD signal remains `preliminary_pd_internal_signal` because empirical permutation support is not significant at 0.05 and independent PD replication is still needed.

## 11. PNAS Readiness Assessment

The PNAS-oriented biological framing is plausible, but the project still needs replicated AD and PD axis-level validation, random-axis controls, and manuscript claims tied directly to evidence-strength tables.

## 12. Next Cohort Priorities

1. Add one larger AD external cohort beyond Mathys feasibility.
2. Add at least one additional PD donor/sample-level cohort.
3. Re-run Phase 21 axis comparisons after each new cohort is harmonized.
4. Keep all claims at candidate or preliminary level until replicated.
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate the Phase 21 PNAS biological discovery report.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase21_pnas_biological_discovery_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(build_report(args), encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
