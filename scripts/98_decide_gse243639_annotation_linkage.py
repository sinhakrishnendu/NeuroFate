#!/usr/bin/env python3
"""Decide whether GSE243639 workbook annotations can be safely linked."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


DECISION_COLUMNS = [
    "dataset_id",
    "decision_category",
    "direct_or_normalized_best_rule",
    "best_overlap_rate",
    "row_order_decision",
    "phase18_match_rate",
    "safe_to_build_annotation_map",
    "recommended_action",
    "notes",
]
SAFE_OVERLAP_THRESHOLD = 0.95


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: str | None, default: float = 0.0) -> float:
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def best_id_rule(rows: list[dict[str, str]]) -> dict[str, str]:
    id_rows = [row for row in rows if "candidate" not in row.get("rule_name", "")]
    if not id_rows:
        return {}
    return max(id_rows, key=lambda row: to_float(row.get("overlap_rate")))


def row_order_decision(rows: list[dict[str, str]]) -> str:
    decisions = {row.get("decision", "") for row in rows if row.get("decision")}
    if "safe_row_order_linkage" in decisions:
        return "safe_row_order_linkage"
    if "unsafe_row_order_linkage" in decisions:
        return "unsafe_row_order_linkage"
    if "inconclusive_row_order_linkage" in decisions:
        return "inconclusive_row_order_linkage"
    return "not_available"


def phase18_match_rate(rows: list[dict[str, str]]) -> float:
    if not rows:
        return 0.0
    row = rows[0]
    return to_float(row.get("match_rate") or row.get("value"))


def decide(normalization_rows: list[dict[str, str]], row_order_rows: list[dict[str, str]], phase18_rows: list[dict[str, str]]) -> dict[str, str]:
    best = best_id_rule(normalization_rows)
    best_rate = to_float(best.get("overlap_rate"))
    best_rule = best.get("rule_name", "unavailable")
    row_decision = row_order_decision(row_order_rows)
    match_rate = phase18_match_rate(phase18_rows)
    if best_rule == "raw_id" and best_rate >= SAFE_OVERLAP_THRESHOLD:
        category = "direct_id_linkage_safe"
        action = "Build annotation map by direct cell ID linkage."
        safe = "true"
    elif best_rate >= SAFE_OVERLAP_THRESHOLD:
        category = "normalized_id_linkage_safe"
        action = f"Build annotation map with audited normalization rule: {best_rule}."
        safe = "true"
    elif row_decision == "safe_row_order_linkage":
        category = "row_order_linkage_safe_with_caution"
        action = "Build a cautious row-order annotation map and label it as row-order linked."
        safe = "true"
    elif row_decision == "inconclusive_row_order_linkage":
        category = "annotation_linkage_inconclusive"
        action = "Do not build cell-type-aware features; use Phase 16 global sample-level PD extension."
        safe = "false"
    else:
        category = "annotation_linkage_unsafe"
        action = "Retire workbook annotation use for this dataset and use Phase 16 only."
        safe = "false"
    notes = (
        "Direct or normalized ID linkage requires >95% overlap. "
        "Row-order linkage requires exact row counts plus independent sample-consistency checks. "
        f"Observed Phase 18 match rate={match_rate:.6g}."
    )
    return {
        "dataset_id": "gse243639_pd_snpc",
        "decision_category": category,
        "direct_or_normalized_best_rule": best_rule,
        "best_overlap_rate": f"{best_rate:.8g}",
        "row_order_decision": row_decision,
        "phase18_match_rate": f"{match_rate:.8g}",
        "safe_to_build_annotation_map": safe,
        "recommended_action": action,
        "notes": notes,
    }


def write_tsv(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=DECISION_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def write_markdown(path: Path, row: dict[str, str]) -> None:
    lines = [
        "# Phase 19 GSE243639 Annotation-Linkage Decision",
        "",
        f"- Decision: `{row['decision_category']}`",
        f"- Best ID rule: `{row['direct_or_normalized_best_rule']}`",
        f"- Best ID overlap rate: {row['best_overlap_rate']}",
        f"- Row-order decision: `{row['row_order_decision']}`",
        f"- Phase 18 annotation match rate: {row['phase18_match_rate']}",
        f"- Safe to build annotation map: `{row['safe_to_build_annotation_map']}`",
        "",
        "## Recommended Action",
        row["recommended_action"],
        "",
        "## Scientific Interpretation",
    ]
    if row["safe_to_build_annotation_map"] != "true":
        lines.extend(
            [
                "Cell-type-aware GSE243639 validation is not currently supported from the available workbook because cell IDs cannot be safely linked to expression cells.",
                "Phase 16 remains the valid global sample-level PD extension for this dataset.",
                "Phase 17 and Phase 18 cell-type-aware outputs should be treated as technical diagnostics, not biological conclusions.",
            ]
        )
    else:
        lines.append("A cautious annotation map may be generated, but all linkage assumptions must remain visible in downstream reports.")
    lines.extend(["", "## Notes", row["notes"]])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Decide GSE243639 annotation-linkage safety.")
    parser.add_argument("--normalization-audit", type=Path, default=Path("results/tables/phase19_gse243639_deep_cell_id_overlap.tsv"))
    parser.add_argument("--row-order-audit", type=Path, default=Path("results/tables/phase19_gse243639_row_order_link_audit.tsv"))
    parser.add_argument("--phase18-match-summary", type=Path, default=Path("results/tables/phase18_gse243639_annotation_match_summary.tsv"))
    parser.add_argument("--output-md", type=Path, default=Path("results/reports/phase19_gse243639_annotation_linkage_decision.md"))
    parser.add_argument("--output-tsv", type=Path, default=Path("results/tables/phase19_gse243639_annotation_linkage_decision.tsv"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    row = decide(
        read_tsv(args.normalization_audit),
        read_tsv(args.row_order_audit),
        read_tsv(args.phase18_match_summary),
    )
    write_tsv(args.output_tsv, row)
    write_markdown(args.output_md, row)
    print(f"Wrote {args.output_tsv}")
    print(f"Wrote {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
