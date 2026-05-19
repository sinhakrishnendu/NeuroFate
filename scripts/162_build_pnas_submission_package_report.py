#!/usr/bin/env python3
"""Build a conservative PNAS submission package and bottleneck report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


COLUMNS = ["area", "current_status", "evidence", "submission_action", "stronger_claim_blocker"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def readiness_status(rows: list[dict[str, str]], criterion: str) -> str:
    for row in rows:
        if row.get("criterion") == criterion:
            return row.get("status", "missing")
    return "missing"


def strongest_axis(rows: list[dict[str, str]]) -> dict[str, str]:
    for row in rows:
        if row.get("axis_id") == "neuronal_vulnerability_axis":
            return row
    return rows[0] if rows else {}


def build_rows(readiness: list[dict[str, str]], evidence: list[dict[str, str]], manuscript: Path) -> list[dict[str, str]]:
    axis = strongest_axis(evidence)
    manuscript_status = "available" if manuscript.exists() and manuscript.stat().st_size > 0 else "missing"
    return [
        {
            "area": "software",
            "current_status": "submission_usable",
            "evidence": "Package scripts, registries, tests, safe extraction lanes, endpoint-locked reports",
            "submission_action": "Use as reproducible NeuroFate-Axis software framework",
            "stronger_claim_blocker": "None for software-resource framing; raw-data redistribution remains external",
        },
        {
            "area": "ad_discovery",
            "current_status": "strong_internal",
            "evidence": "SEA-AD endpoint-locked dementia axes with FDR support",
            "submission_action": "Present as internal AD discovery anchor",
            "stronger_claim_blocker": "External replication required for axis-specific generality",
        },
        {
            "area": "ad_replication",
            "current_status": readiness_status(readiness, "independent_ad_replication"),
            "evidence": (
                f"{axis.get('axis_id', 'axis')} GSE174367 p={axis.get('gse174367_p', '')}; "
                f"FDR={axis.get('gse174367_fdr', '')}; direction={axis.get('gse174367_direction_consistency', '')}"
            ),
            "submission_action": "Claim nominal independent AD replication where supported",
            "stronger_claim_blocker": "FDR-robust independent AD replication would strengthen the claim",
        },
        {
            "area": "pd_convergence",
            "current_status": readiness_status(readiness, "independent_pd_replication"),
            "evidence": (
                f"GSE243639 axis p={axis.get('gse243639_axis_p', '')}; "
                f"GSE184950 p={axis.get('gse184950_p', '')}; FDR={axis.get('gse184950_fdr', '')}"
            ),
            "submission_action": "Frame as preliminary PD convergence and replication infrastructure",
            "stronger_claim_blocker": "Statistically supported independent PD axis replication",
        },
        {
            "area": "shared_ad_pd_claim",
            "current_status": readiness_status(readiness, "shared_ad_pd_axis_claim"),
            "evidence": "AD side is stronger than PD side",
            "submission_action": "Do not claim a confirmed shared AD/PD mechanism",
            "stronger_claim_blocker": "At least one statistically supported PD replication axis aligned with AD evidence",
        },
        {
            "area": "manuscript",
            "current_status": manuscript_status,
            "evidence": str(manuscript),
            "submission_action": "Use the PNAS template draft with conservative claim framing",
            "stronger_claim_blocker": "Finalize author metadata, figures, and journal-specific administrative fields",
        },
    ]


def write_report(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# PNAS Submission Package Bottleneck Report",
        "",
        "This report separates blockers for a conservative NeuroFate-Axis submission from blockers for stronger shared AD/PD mechanism language.",
        "",
        "## Bottom Line",
        "The software and manuscript can now support a conservative submission centered on endpoint-locked AD axis discovery, nominal independent AD replication of neuronal vulnerability, and preliminary PD convergence. The stronger shared AD/PD mechanism claim remains blocked by weak PD axis replication.",
        "",
        "## Bottleneck Matrix",
        "",
        "| Area | Current status | Submission action | Stronger-claim blocker |",
        "|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['area']} | {row['current_status']} | {row['submission_action']} | {row['stronger_claim_blocker']} |"
        )
    lines.extend(
        [
            "",
            "## Required Before Submission",
            "- Replace placeholder author, affiliation, correspondence, DOI, and data-availability fields.",
            "- Add final figure files or keep the manuscript tables as the primary evidence display.",
            "- Re-run the no-overclaiming audit after any manual manuscript edits.",
            "",
            "## Required Before Stronger Shared AD/PD Claims",
            "- Obtain statistically supported independent PD axis replication.",
            "- Prefer a larger PD cohort or an endpoint/metadata route with stronger axis-level support.",
            "- Avoid clinical, diagnostic, causal, or definitive shared-mechanism language until that evidence exists.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build conservative PNAS submission package report.")
    parser.add_argument("--readiness", type=Path, default=Path("results/reports/phase23_pnas_readiness_matrix.tsv"))
    parser.add_argument("--evidence", type=Path, default=Path("results/tables/phase32_axis_evidence_ranked.tsv"))
    parser.add_argument("--manuscript", type=Path, default=Path("manuscript/Research_report.tex"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase33_pnas_submission_package_report.md"))
    parser.add_argument("--matrix-output", type=Path, default=Path("results/tables/phase33_pnas_bottleneck_matrix.tsv"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = build_rows(read_tsv(args.readiness), read_tsv(args.evidence), args.manuscript)
    write_tsv(args.matrix_output, rows)
    write_report(args.output, rows)
    print(f"Wrote {args.output}")
    print(f"Wrote {args.matrix_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
