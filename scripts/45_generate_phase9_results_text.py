#!/usr/bin/env python3
"""Draft Phase 9 Mathys CSV external-validation results text with no-overclaiming."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import TextIO


SMALL_N_WARNING_THRESHOLD = 20
DEFAULT_TARGET_PANEL = Path("metadata/target_gene_panel_v1.tsv")
DEFAULT_MATHYS_SPARSE = Path("data/interim/external/mathys_2019/mathys_sparse_gene_panel_expression.tsv.gz")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def sparse_extracted_genes(path: Path) -> set[str]:
    genes, _rows = sparse_expression_gene_stats(path)
    return genes


def sparse_expression_gene_stats(path: Path) -> tuple[set[str], int]:
    if not path.exists():
        return set(), 0
    genes: set[str] = set()
    rows = 0
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows += 1
            gene = row.get("gene_symbol", "")
            if gene:
                genes.add(gene)
    return genes, rows


def requested_target_genes(path: Path = DEFAULT_TARGET_PANEL) -> set[str]:
    rows = read_tsv(path)
    genes = {
        row.get("gene_symbol", "").strip()
        for row in rows
        if row.get("gene_symbol", "").strip()
    }
    return genes


def target_gene_overlap_counts(tables_dir: Path) -> tuple[int | None, int | None, str]:
    overlap = read_tsv(tables_dir / "mathys_gene_overlap.tsv")
    if overlap:
        present = sum(1 for row in overlap if row.get("mathys_status") == "present")
        missing = sum(1 for row in overlap if row.get("mathys_status") == "missing")
        if present or missing:
            return present, missing, "mathys_gene_overlap.tsv"

    audit_rows = read_tsv(tables_dir / "phase10_mathys_gene_audit.tsv")
    if audit_rows:
        present = sum(1 for row in audit_rows if row.get("mathys_status") == "present")
        missing = sum(1 for row in audit_rows if row.get("mathys_status") == "missing")
        if present or missing:
            return present, missing, "phase10_mathys_gene_audit.tsv"

    return None, None, "gene overlap table unavailable"


def mathys_gene_extraction_summary(tables_dir: Path) -> dict[str, object]:
    requested = requested_target_genes()
    extracted, sparse_rows = sparse_expression_gene_stats(DEFAULT_MATHYS_SPARSE)
    present, missing, overlap_source = target_gene_overlap_counts(tables_dir)
    missing_genes = sorted(requested - extracted) if requested and extracted else []
    return {
        "requested_target_genes": len(requested),
        "extracted_target_genes": len(extracted),
        "missing_target_genes": len(missing_genes),
        "sparse_expression_rows": sparse_rows,
        "overlap_present": present,
        "overlap_missing": missing,
        "overlap_source": overlap_source,
    }


def mathys_sample_count(tables_dir: Path) -> int:
    feature_rows = read_tsv(tables_dir / "mathys_2019_phase5_donor_feature_table.tsv")
    if feature_rows:
        return len(feature_rows)
    for row in read_tsv(tables_dir / "phase9_mathys_external_validation_metrics.tsv"):
        try:
            return int(float(row.get("n_test", "0") or 0))
        except ValueError:
            continue
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Draft Phase 9 Mathys CSV external-validation text.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase9_results_summary.txt"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    structure = read_tsv(args.tables_dir / "mathys_csv_structure_summary.tsv")
    labels = read_tsv(args.tables_dir / "mathys_2019_label_summary.tsv")
    metrics = read_tsv(args.tables_dir / "phase9_mathys_external_validation_metrics.tsv")
    gene_summary = mathys_gene_extraction_summary(args.tables_dir)
    n_samples = mathys_sample_count(args.tables_dir)
    evidence_label = (
        "preliminary external feasibility evidence"
        if n_samples < SMALL_N_WARNING_THRESHOLD
        else "external validation evidence requiring cohort-specific review"
    )

    lines = [
        "Phase 9 Mathys 2019 CSV External Validation Summary",
        "",
        "This text is generated from Phase 9 CSV-processing outputs and should be reviewed before manuscript use.",
        f"Evidence label: {evidence_label}.",
        "",
        f"Structure summary rows: {len(structure)}",
        f"Mathys sample-level units: {n_samples}",
        f"Requested target genes: {gene_summary['requested_target_genes']}",
        f"Extracted target genes: {gene_summary['extracted_target_genes']}",
        f"Missing target genes after extraction: {gene_summary['missing_target_genes']}",
        f"Sparse expression rows extracted: {gene_summary['sparse_expression_rows']}",
        f"Gene-overlap table status: {gene_summary['overlap_source']}",
    ]
    if gene_summary["overlap_present"] is not None:
        lines.append(
            "Gene-overlap table counts: present={present}, missing={missing}".format(
                present=gene_summary["overlap_present"],
                missing=gene_summary["overlap_missing"],
            )
        )
    if n_samples < SMALL_N_WARNING_THRESHOLD:
        lines.extend(
            [
                "",
                "Reviewer-facing warning: Mathys currently has small sample-level n; results should be described as preliminary external feasibility evidence and not as definitive cross-cohort validation.",
            ]
        )
    if gene_summary["overlap_source"] == "gene overlap table unavailable":
        lines.append(
            "Gene overlap table unavailable; extracted gene counts were derived from the sparse-like expression file instead of reporting present=0/missing=0."
        )

    lines.extend(["", "Mathys label summary:"])
    for row in labels:
        lines.append(f"- {row.get('label_field', '')} / {row.get('label', '')}: {row.get('sample_count', '')}")
    lines.extend(["", "External validation metrics:"])
    for row in metrics:
        lines.append(
            "- {mode}: AUROC={auroc}, AUPRC={auprc}, balanced_accuracy={bal}, Brier={brier}, note={note}".format(
                mode=row.get("validation_mode", ""),
                auroc=row.get("auroc", ""),
                auprc=row.get("auprc", ""),
                bal=row.get("balanced_accuracy", ""),
                brier=row.get("brier_score", ""),
                note=row.get("notes", ""),
            )
        )
    lines.append("")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
