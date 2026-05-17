#!/usr/bin/env python3
"""Generate concise Phase 3 manuscript results text from summary TSVs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return 0.0


def strongest_celltype_genes(tables_dir: Path) -> list[str]:
    rows = read_tsv(tables_dir / "gene_by_celltype_summary.tsv")
    best: dict[str, dict[str, str]] = {}
    for row in rows:
        gene = row["gene_symbol"]
        if gene not in best or to_float(row["mean_expression_all_cells"]) > to_float(
            best[gene]["mean_expression_all_cells"]
        ):
            best[gene] = row
    ranked = sorted(best.values(), key=lambda row: to_float(row["mean_expression_all_cells"]), reverse=True)
    return [
        f"{row['gene_symbol']} was highest in {row['cell_type']} "
        f"(mean={row['mean_expression_all_cells']}, detection={row['detection_rate']})."
        for row in ranked[:8]
    ]


def strongest_ad_trends(tables_dir: Path) -> list[str]:
    rows = read_tsv(tables_dir / "gene_by_ad_pathology.tsv")
    by_gene_variable: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        key = (row["gene_symbol"], row["pathology_variable"])
        by_gene_variable.setdefault(key, []).append(row)
    trends: list[tuple[float, str]] = []
    for (gene, variable), trend_rows in by_gene_variable.items():
        values = [to_float(row["mean_expression_all_cells"]) for row in trend_rows]
        if not values:
            continue
        trends.append((max(values) - min(values), f"{gene} varied across {variable} groups."))
    return [text for _, text in sorted(trends, reverse=True)[:8]]


def strongest_signature_rows(tables_dir: Path, filename: str) -> list[str]:
    rows = [
        row for row in read_tsv(tables_dir / filename)
        if row["grouping_variable"] in {"cell_type", "Overall AD neuropathological Change", "Cognitive Status"}
    ]
    ranked = sorted(rows, key=lambda row: to_float(row["mean_expression_per_gene_per_cell"]), reverse=True)
    return [
        f"{row['signature_name']} was strongest for {row['grouping_variable']}={row['group_label']} "
        f"(mean={row['mean_expression_per_gene_per_cell']})."
        for row in ranked[:6]
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 3 manuscript results text.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/tables/phase3_results_summary.txt"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    sections = [
        ("Strongest Cell-Type Enriched Genes", strongest_celltype_genes(args.tables_dir)),
        ("Strongest AD-Associated Trends", strongest_ad_trends(args.tables_dir)),
        (
            "Strongest Microglial Activation Signals",
            strongest_signature_rows(args.tables_dir, "microglial_activation_signature.tsv"),
        ),
        (
            "Neuronal Vulnerability Patterns",
            strongest_signature_rows(args.tables_dir, "neuronal_signature_summary.tsv")
            + strongest_signature_rows(args.tables_dir, "neurodegeneration_signature_summary.tsv"),
        ),
    ]

    with args.output.open("w", encoding="utf-8") as handle:
        handle.write("Phase 3 Sparse Expression Results Summary\n")
        handle.write("Generated from sparse target-gene summaries only.\n\n")
        for title, lines in sections:
            handle.write(f"{title}\n")
            if lines:
                for line in lines:
                    handle.write(f"- {line}\n")
            else:
                handle.write("- No rows available for this section.\n")
            handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
