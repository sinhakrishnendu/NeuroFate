#!/usr/bin/env python3
"""Generate the Phase 31 GSE174367 gene-mapping report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 31 GSE174367 gene identifier mapping report.")
    parser.add_argument("--gene-id-audit", type=Path, default=Path("results/tables/phase31_gse174367_bulk_gene_identifier_audit.tsv"))
    parser.add_argument("--gene-mapping", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_gene_mapping.tsv"))
    parser.add_argument("--coverage", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_gene_coverage.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase31_gse174367_gene_mapping_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    audit = read_tsv(args.gene_id_audit)
    mapping = read_tsv(args.gene_mapping)
    coverage = read_tsv(args.coverage)
    mapped_symbols = sorted({row.get("gene_symbol", "") for row in mapping if row.get("gene_symbol")})
    best = max((int(float(row.get("matched_axis_genes", "0") or 0)) for row in audit), default=0)
    lines = [
        "# Phase 31 GSE174367 Gene Identifier Mapping Report",
        "",
        "## 1. Why Phase 31 Was Needed",
        "Phase 30 fixed the sample map using the RDA-internal `targets` table, but conversion failed because `normExpr.reg` row identifiers did not directly match NeuroFate axis symbols.",
        "",
        "## 2. Gene Identifier Format In normExpr.reg",
        f"Gene-ID audit rows: {len(audit)}. Best matched axis genes in audit: {best}.",
        "",
        "## 3. NeuroFate Symbol-to-Ensembl Mapping",
        f"Mapped unique NeuroFate axis genes in the conversion table: {len(mapped_symbols)}.",
        "",
        "## 4. Unmapped Genes",
        f"Axis coverage rows: {len(coverage)}. Consult `phase31_gse174367_bulk_axis_gene_coverage.tsv` for missing genes by axis.",
        "",
        "## 5. Whether AD Replication Can Proceed",
        "AD replication can proceed only if the converter maps at least the conservative minimum number of axis genes and produces sample-level AD/Control axis scores.",
        "",
        "## 6. Conservative Interpretation",
        "Phase 31 repairs identifier mapping only. It does not establish deployment, decision-support, mechanism-proof, or definitive cross-disease evidence.",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
