#!/usr/bin/env python3
"""Audit NeuroFate text outputs for claims that exceed the current evidence."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


SCAN_PATTERNS = [
    "results/tables/*summary*.txt",
    "results/reports/*.md",
    "results/reports/*.html",
    "results/reports/*.txt",
    "RESULTS_INTERPRETATION.md",
    "README.md",
    "manuscript/*.tex",
]

UNSAFE_PHRASES = [
    {
        "phrase": "validated across cohorts",
        "requires": "external_sample_n_at_least_20",
        "recommendation": "Use preliminary external feasibility evidence until a larger external cohort is available.",
    },
    {
        "phrase": "externally validated",
        "requires": "external_sample_n_at_least_20",
        "recommendation": "Use preliminary external feasibility unless a larger external validation cohort exists.",
    },
    {
        "phrase": "foundation model",
        "requires": "foundation_model_training",
        "recommendation": "Use small donor-level neural model unless foundation-model training exists.",
    },
    {
        "phrase": "causal",
        "requires": "causal_design",
        "recommendation": "Use association or statistical association unless a causal design is implemented.",
    },
    {
        "phrase": "clinical-grade",
        "requires": "clinical_validation",
        "recommendation": "Use research prototype unless clinical validation exists.",
    },
    {
        "phrase": "diagnostic",
        "requires": "clinical_validation",
        "recommendation": "Use research model or hypothesis generator unless clinical validation exists.",
    },
    {
        "phrase": "state-of-the-art",
        "requires": "independent_benchmark_superiority",
        "recommendation": "Use benchmarked or tested only when comparative benchmarks support it.",
    },
    {
        "phrase": "generalizable across diseases",
        "requires": "multi_disease_external_validation",
        "recommendation": "Use designed for future cross-disease testing unless multi-disease validation exists.",
    },
    {
        "phrase": "biomarker",
        "requires": "validated_biomarker_evidence",
        "recommendation": "Use candidate biomarker-like signal when phrased cautiously.",
    },
    {
        "phrase": "patient-level diagnosis",
        "requires": "clinical_validation",
        "recommendation": "Use donor-level research stratification, not patient-level diagnosis.",
    },
]

ALLOWED_CONTEXTS = [
    "preliminary external feasibility",
    "research software",
    "candidate biomarker-like signal",
    "not clinical",
    "not a clinical",
    "not clinical-grade",
    "not diagnostic",
    "not a diagnostic",
    "not causal",
    "not a foundation model",
    "not evidence of state-of-the-art",
    "does not establish patient-level diagnosis",
    "do not establish patient-level diagnosis",
    "use association language",
    "internal diagnostic",
    "internal_diagnostic",
    "diagnostic if enough sample units",
    "does not infer causality",
    "insufficient for",
    "do not claim",
    "forbidden_language",
    "cannot be claimed",
    "should not be described",
    "not yet",
]


def count_mathys_sample_units(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        return sum(1 for _ in reader)


def allowed_conditions(mathys_n: int) -> dict[str, bool]:
    return {
        "external_sample_n_at_least_20": mathys_n >= 20,
        "foundation_model_training": False,
        "causal_design": False,
        "clinical_validation": False,
        "independent_benchmark_superiority": False,
        "multi_disease_external_validation": False,
        "validated_biomarker_evidence": False,
    }


def iter_scan_files(patterns: list[str]) -> list[Path]:
    files: list[Path] = []
    for pattern in patterns:
        files.extend(path for path in sorted(Path(".").glob(pattern)) if path.is_file())
    return sorted(set(files))


def audit_file(path: Path, conditions: dict[str, bool]) -> list[dict[str, str]]:
    findings: list[dict[str, str]] = []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        return [
            {
                "source_path": str(path),
                "line_number": "",
                "phrase": "read_error",
                "severity": "warning",
                "allowed": "false",
                "reason": str(exc),
                "recommendation": "Inspect the file manually.",
                "line_text": "",
            }
        ]
    for line_number, line in enumerate(lines, start=1):
        lower_line = line.lower()
        for rule in UNSAFE_PHRASES:
            phrase = rule["phrase"]
            if phrase in lower_line:
                requirement = rule["requires"]
                context_allowed = any(context in lower_line for context in ALLOWED_CONTEXTS)
                allowed = conditions.get(requirement, False) or context_allowed
                severity = "ok" if allowed else "high"
                reason = (
                    f"context allows cautious usage of phrase requiring {requirement}"
                    if context_allowed
                    else f"requires {requirement}"
                )
                findings.append(
                    {
                        "source_path": str(path),
                        "line_number": str(line_number),
                        "phrase": phrase,
                        "severity": severity,
                        "allowed": str(allowed).lower(),
                        "reason": reason,
                        "recommendation": rule["recommendation"],
                        "line_text": line.strip()[:300],
                    }
                )
    return findings


def run_audit(output: Path) -> list[dict[str, str]]:
    mathys_n = count_mathys_sample_units(Path("results/tables/mathys_2019_phase5_donor_feature_table.tsv"))
    conditions = allowed_conditions(mathys_n)
    findings: list[dict[str, str]] = []
    for path in iter_scan_files(SCAN_PATTERNS):
        findings.extend(audit_file(path, conditions))
    if not findings:
        findings.append(
            {
                "source_path": "all_scanned_files",
                "line_number": "",
                "phrase": "none",
                "severity": "none",
                "allowed": "true",
                "reason": f"no configured unsafe phrases detected; mathys_sample_units={mathys_n}",
                "recommendation": "No action required.",
                "line_text": "",
            }
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "source_path",
                "line_number",
                "phrase",
                "severity",
                "allowed",
                "reason",
                "recommendation",
                "line_text",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(findings)
    return findings


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit NeuroFate outputs for overclaiming language.")
    parser.add_argument("--output", type=Path, default=Path("results/reports/no_overclaiming_audit.tsv"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    findings = run_audit(args.output)
    high = sum(1 for row in findings if row["severity"] == "high")
    print(f"Wrote {args.output}")
    print(f"High-severity overclaiming flags: {high}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
