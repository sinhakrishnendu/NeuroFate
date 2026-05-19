#!/usr/bin/env python3
"""Triage Phase 15 external validation candidate datasets without remote access."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REGISTRY = Path("metadata/phase15_external_validation_candidates.tsv")
TRIAGE_OUTPUT = Path("results/reports/phase15_external_dataset_triage.tsv")
SUMMARY_OUTPUT = Path("results/reports/phase15_external_dataset_priority_summary.md")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def readiness(row: dict[str, str]) -> tuple[str, str]:
    status = row.get("current_status", "").lower()
    access = row.get("access_type", "").lower()
    formats = row.get("expected_file_formats", "").lower()
    metadata_required = row.get("metadata_required", "").lower()
    if "controlled" in access or "controlled" in status:
        return "controlled_access_required", "Controlled or registered access is required before local validation."
    if "unknown" in formats or "unknown" in status:
        return "format_unknown", "Manual file inspection is needed after acquisition."
    if not metadata_required or "donor" not in metadata_required:
        return "insufficient_metadata", "Required donor/sample-level metadata are unclear."
    if row.get("priority") in {"A", "B", "C"}:
        return "ready_for_manual_acquisition", "Candidate is ready for user-led acquisition planning."
    return "deprioritized", "Lower priority optional validation layer."


def triage_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    output = []
    for row in rows:
        category, reason = readiness(row)
        output.append(
            {
                "dataset_id": row["dataset_id"],
                "disease": row["disease"],
                "priority": row["priority"],
                "accession_or_portal": row["accession_or_portal"],
                "readiness_category": category,
                "recommended_next_step": reason,
                "local_raw_dir": row["local_raw_dir"],
                "validation_role": row["validation_role"],
            }
        )
    return output


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "dataset_id",
                "disease",
                "priority",
                "accession_or_portal",
                "readiness_category",
                "recommended_next_step",
                "local_raw_dir",
                "validation_role",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["readiness_category"]] = counts.get(row["readiness_category"], 0) + 1
    lines = [
        "# Phase 15 External Dataset Priority Summary",
        "",
        "No remote access, downloads, or biological file processing were performed.",
        "",
        "## Readiness Counts",
    ]
    for category, count in sorted(counts.items()):
        lines.append(f"- {category}: {count}")
    lines.extend(["", "## Priority Order"])
    for row in sorted(rows, key=lambda item: item["priority"]):
        lines.append(
            f"- {row['priority']} / {row['dataset_id']}: {row['readiness_category']} - {row['recommended_next_step']}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Triage Phase 15 external validation candidates.")
    parser.add_argument("--registry", type=Path, default=REGISTRY)
    parser.add_argument("--output", type=Path, default=TRIAGE_OUTPUT)
    parser.add_argument("--summary", type=Path, default=SUMMARY_OUTPUT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = triage_rows(read_tsv(args.registry))
    write_tsv(args.output, rows)
    write_summary(args.summary, rows)
    print(f"Wrote {args.output}")
    print(f"Wrote {args.summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
