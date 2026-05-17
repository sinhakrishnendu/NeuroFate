#!/usr/bin/env python3
"""Draft Phase 4 manuscript-results text from statistical summary TSVs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def top_gene_associations(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(
        rows,
        key=lambda row: (
            to_float(row["fdr_p_value"]) if row.get("fdr_p_value") != "nan" else 1.0,
            -abs(to_float(row.get("effect_size_max_minus_min", "0"))),
        ),
    )[:10]


def top_vulnerable_celltypes(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    neuro_rows = [row for row in rows if row.get("signature_id") == "neurodegeneration_signature"]
    return sorted(neuro_rows, key=lambda row: int(row.get("rank", "999")))[:10]


def strongest_apoe_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(
        rows,
        key=lambda row: (
            row.get("index_id", ""),
            -to_float(row.get("donor_mean_index", "0")),
        ),
    )[:9]


def strongest_mixed_pathology(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(
        rows,
        key=lambda row: (
            to_float(row["fdr_p_value"]) if row.get("fdr_p_value") != "nan" else 1.0,
            -abs(to_float(row.get("effect_size_max_minus_min", "0"))),
        ),
    )[:10]


def render_gene_section(rows: list[dict[str, str]]) -> list[str]:
    lines = ["Top donor-aware disease-associated target-gene statistics:"]
    for row in top_gene_associations(rows):
        lines.append(
            "- {gene} vs {variable}: {method}, statistic={stat}, FDR={fdr}, effect={effect}, direction={direction}".format(
                gene=row.get("gene_symbol", ""),
                variable=row.get("test_variable", ""),
                method=row.get("test_method", ""),
                stat=row.get("rank_statistic", ""),
                fdr=row.get("fdr_p_value", ""),
                effect=row.get("effect_size_max_minus_min", ""),
                direction=row.get("direction", ""),
            )
        )
    return lines


def render_celltype_section(rows: list[dict[str, str]]) -> list[str]:
    lines = ["Top cell subclasses by neurodegeneration vulnerability index:"]
    for row in top_vulnerable_celltypes(rows):
        lines.append(
            "- rank {rank}: {celltype}, donor_mean_index={value}, CI95={low}..{high}".format(
                rank=row.get("rank", ""),
                celltype=row.get("cell_subclass", ""),
                value=row.get("donor_mean_index", ""),
                low=row.get("ci95_low", ""),
                high=row.get("ci95_high", ""),
            )
        )
    return lines


def render_apoe_section(rows: list[dict[str, str]]) -> list[str]:
    lines = ["APOE genotype-linked composite index patterns:"]
    for row in strongest_apoe_rows(rows):
        lines.append(
            "- {index_id} / {label}: donor_mean_index={value}, FDR={fdr}, effect={effect}".format(
                index_id=row.get("index_id", ""),
                label=row.get("group_label", ""),
                value=row.get("donor_mean_index", ""),
                fdr=row.get("fdr_p_value", ""),
                effect=row.get("effect_size_max_minus_min", ""),
            )
        )
    return lines


def render_mixed_section(rows: list[dict[str, str]]) -> list[str]:
    lines = ["Mixed-pathology composite index associations:"]
    for row in strongest_mixed_pathology(rows):
        lines.append(
            "- {index_id} vs {variable} / {label}: FDR={fdr}, effect={effect}, direction={direction}".format(
                index_id=row.get("index_id", ""),
                variable=row.get("grouping_variable", ""),
                label=row.get("group_label", ""),
                fdr=row.get("fdr_p_value", ""),
                effect=row.get("effect_size_max_minus_min", ""),
                direction=row.get("direction", ""),
            )
        )
    return lines


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Draft Phase 4 statistical results text.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/tables/phase4_results_summary.txt"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    gene_rows = read_tsv(args.tables_dir / "phase4_gene_statistics.tsv")
    vulnerability_rows = read_tsv(args.tables_dir / "phase4_celltype_vulnerability.tsv")
    apoe_rows = read_tsv(args.tables_dir / "phase4_apoe_analysis.tsv")
    mixed_rows = read_tsv(args.tables_dir / "phase4_mixed_pathology.tsv")

    lines: list[str] = [
        "Phase 4 Statistical Neurodegeneration Biology Summary",
        "",
        "This text is generated from Phase 4 statistical TSV outputs and should be reviewed before manuscript use.",
        "",
        *render_gene_section(gene_rows),
        "",
        *render_celltype_section(vulnerability_rows),
        "",
        *render_apoe_section(apoe_rows),
        "",
        *render_mixed_section(mixed_rows),
        "",
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
