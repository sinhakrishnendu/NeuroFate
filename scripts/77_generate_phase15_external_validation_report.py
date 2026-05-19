#!/usr/bin/env python3
"""Generate Phase 15 external validation expansion report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def lines_for_table(rows: list[dict[str, str]], fields: list[str]) -> list[str]:
    if not rows:
        return ["- Not available."]
    lines = []
    for row in rows:
        values = ", ".join(f"{field}={row.get(field, '')}" for field in fields)
        lines.append(f"- {values}")
    return lines


def build_report(output: Path) -> None:
    candidates = read_tsv(Path("metadata/phase15_external_validation_candidates.tsv"))
    triage = read_tsv(Path("results/reports/phase15_external_dataset_triage.tsv"))
    metrics = read_tsv(Path("results/tables/phase15_multi_external_validation_metrics.tsv"))
    reliability = read_tsv(Path("results/reports/phase15_external_validation_reliability.tsv"))
    deltas = read_tsv(Path("results/reports/phase15_claim_strength_delta.tsv"))
    next_dataset = next((row for row in triage if row.get("readiness_category") == "ready_for_manual_acquisition"), None)
    lines = [
        "# Phase 15 External Validation Expansion Report",
        "",
        "No datasets were downloaded and no external extraction or modeling was run by this report.",
        "",
        "## 1. Candidate Datasets",
        *lines_for_table(candidates, ["dataset_id", "disease", "priority", "accession_or_portal", "validation_role"]),
        "",
        "## 2. Acquisition Status",
        *lines_for_table(triage, ["dataset_id", "readiness_category", "recommended_next_step"]),
        "",
        "## 3. Format Status",
        "- Use `scripts/70_inspect_external_dataset_files.py` after manual acquisition to inventory local files and recommend format-specific plans.",
        "",
        "## 4. Metadata Harmonization Status",
        "- Use `scripts/71_inspect_external_metadata_safe.py` to map donor/sample/cell/diagnosis/cell-type fields without reading expression matrices.",
        "",
        "## 5. Gene-Panel Overlap",
        "- Use `scripts/72_plan_external_target_gene_overlap.py` with feature/gene files only; count matrices are not loaded.",
        "",
        "## 6. External Validation Readiness",
        *lines_for_table(reliability, ["dataset_id", "reliability_category", "n_test", "feature_overlap_count"]),
        "",
        "## 7. External Validation Results If Available",
        *lines_for_table(metrics, ["dataset_id", "validation_mode", "n_test", "reliability_flag"]),
        "",
        "## 8. Remaining Blockers",
        "- Manual acquisition of GSE243639 is the highest-priority next step for PD validation.",
        "- AD external cohorts require file-format triage and donor/sample metadata confirmation.",
        "- Controlled-access resources require user-managed approvals and provenance.",
        "",
        "## 9. Claim-Strength Implications",
        *lines_for_table(deltas, ["task", "old_claim_strength", "new_claim_strength", "phase15_external_status"]),
        "",
        "## 10. Next Recommended Dataset To Acquire",
        f"- {next_dataset['dataset_id']}: {next_dataset['recommended_next_step']}" if next_dataset else "- No ready public dataset found in triage table.",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 15 external validation report.")
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase15_external_validation_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    build_report(args.output)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
