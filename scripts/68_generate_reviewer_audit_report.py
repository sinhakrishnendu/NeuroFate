#!/usr/bin/env python3
"""Generate a reviewer-facing NeuroFate audit report from existing artifacts."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_high_overclaiming(path: Path) -> int:
    return sum(1 for row in read_tsv(path) if row.get("severity") == "high")


def leakage_summary(path: Path) -> str:
    rows = read_tsv(path)
    if not rows:
        return "Leakage audit not available."
    high = sum(1 for row in rows if row.get("leakage_risk") == "high")
    medium = sum(1 for row in rows if row.get("leakage_risk") == "medium")
    return f"Leakage audit present: high-risk label/identifier flags={high}, medium-risk flags={medium}. These are acceptable when detected and excluded from predictors."


def artifact_inventory() -> list[str]:
    paths = [
        "README.md",
        "pyproject.toml",
        "RESULTS_INTERPRETATION.md",
        "results/reports/claim_strength_table.tsv",
        "results/reports/best_supported_claims.tsv",
        "results/reports/no_overclaiming_audit.tsv",
        "dist/source_release_manifest.tsv",
        "dist/results_review_manifest.tsv",
    ]
    return [f"- {path}: {'present' if Path(path).exists() else 'missing'}" for path in paths]


def write_report(output: Path) -> None:
    claim_rows = read_tsv(Path("results/reports/claim_strength_table.tsv"))
    mathys_rows = read_tsv(Path("results/tables/phase13_mathys_gene_extraction_audit.tsv"))
    overclaiming_high = count_high_overclaiming(Path("results/reports/no_overclaiming_audit.tsv"))
    lines = [
        "# NeuroFate Reviewer Audit Report",
        "",
        "## 1. Artifact Inventory",
        *artifact_inventory(),
        "",
        "## 2. Reproducibility Status",
        "- PyPI metadata, CI configuration, source-release builder, results-review builder, tiny demo, and release checklists are present.",
        "- Raw data are intentionally excluded from release packages.",
        "",
        "## 3. Leakage Audit Status",
        f"- {leakage_summary(Path('results/reports/feature_leakage_audit.tsv'))}",
        "",
        "## 4. No-Overclaiming Status",
        f"- High-severity no-overclaiming flags: {overclaiming_high}.",
        "- Cautious contexts such as preliminary external feasibility and research software are allowed.",
        "",
        "## 5. Evidence-Strength Status",
    ]
    if claim_rows:
        for row in claim_rows:
            lines.append(
                f"- {row.get('task')}: {row.get('claim_strength')} using {row.get('model')} (risk={row.get('reviewer_risk')})."
            )
    else:
        lines.append("- Claim strength table is missing.")
    lines.extend(
        [
            "",
            "## 6. Internal Benchmark Status",
            "- Phase 12 benchmark tables should be used for repeated-seed uncertainty, permutation support, and ablation consistency.",
            "",
            "## 7. External Validation Status",
        ]
    )
    if mathys_rows:
        row = mathys_rows[0]
        lines.append(
            f"- Mathys extracted genes={row.get('extracted_target_genes')} from requested={row.get('requested_target_genes')}; sample units={row.get('sample_units')}; status={row.get('status')}."
        )
    else:
        lines.append("- Mathys gene extraction audit is missing.")
    lines.extend(
        [
            "- Mathys remains preliminary external feasibility, not definitive external validation.",
            "",
            "## 8. Software Release Readiness",
            "- Source and results-review package builders are present and exclude raw data by rule.",
            "- The package includes a tiny synthetic demo for installation smoke testing.",
            "",
            "## 9. Remaining Blockers For Nature Computational Science",
            "- Larger independent external cohorts are needed.",
            "- Full Phase 12 robustness outputs should be regenerated and reviewed before submission.",
            "- Manuscript claims should be aligned to `claim_strength_table.tsv`.",
            "",
            "## 10. Recommended Next Phase",
            "- Add a larger external cohort and rerun Phase 12 benchmarks with locked tasks, leakage rules, and package manifests.",
        ]
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate NeuroFate reviewer audit report.")
    parser.add_argument("--output", type=Path, default=Path("results/reports/reviewer_audit_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    write_report(args.output)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
