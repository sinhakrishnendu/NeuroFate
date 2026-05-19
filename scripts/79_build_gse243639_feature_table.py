#!/usr/bin/env python3
"""Build sample-level GSE243639 NeuroFate features from target-gene sparse output."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import TextIO


FEATURE_PREFIXES = ("gene_mean__", "gene_detection__", "index__")
CLINICAL_FIELD_MAP = {
    "Sample ID": "sample_id",
    "Clinical diagnosis": "diagnosis",
    "Age": "age",
    "Sex": "sex",
    "PMI hours": "pmi",
    "RIN measure": "rin",
    "Braak stage for neurofibrillary tangles": "braak",
    "CERAD score for neuritic plaques": "cerad",
    "Lewy bodies presence in midbrain": "lewy_body_midbrain",
    "Lewy bodies presence in limbic regions (amygdala)": "lewy_body_limbic",
    "Lewy bodies presence in neocortical regions (frontal cortex)": "lewy_body_neocortical",
}
SIGNATURES = {
    "MAI": ["TREM2", "TYROBP", "GPNMB", "HLA-DRA", "AIF1"],
    "ASI": ["GFAP"],
    "NVI": ["PINK1", "PRKN", "SNCA", "MAPT", "APOE"],
    "inflammatory_index": ["TREM2", "TYROBP", "GPNMB", "HLA-DRA", "AIF1", "IL1B", "TNF", "NFKB1", "B2M"],
    "mitochondrial_index": ["PINK1", "PRKN", "LRRK2"],
    "neuronal_index": ["SLC17A7", "SST", "PVALB", "LAMP5"],
}


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file, mode="w", encoding="utf-8"),
        ],
    )


def open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %s rows: %d", path, len(rows))


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_") or "missing"


def diagnosis_binary(value: str) -> str:
    lowered = value.lower()
    if "parkinson" in lowered:
        return "1"
    if "control" in lowered:
        return "0"
    return ""


def read_clinical(path: Path, header_line: int = 6, delimiter: str = ";") -> dict[str, dict[str, str]]:
    with open_text(path) as handle:
        for _ in range(header_line - 1):
            next(handle, "")
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows: dict[str, dict[str, str]] = {}
        for raw in reader:
            if not raw:
                continue
            mapped = {canonical: (raw.get(source) or "").strip() for source, canonical in CLINICAL_FIELD_MAP.items()}
            sample_id = mapped.get("sample_id", "")
            if not sample_id:
                continue
            mapped["label__diagnosis_binary"] = diagnosis_binary(mapped.get("diagnosis", ""))
            rows[sample_id] = mapped
    return rows


def read_cell_counts(cell_map: Path) -> dict[str, int]:
    if not cell_map.exists():
        return {}
    counts: Counter[str] = Counter()
    with cell_map.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_id = row.get("sample_id", "")
            if sample_id:
                counts[sample_id] += 1
    return dict(counts)


def stream_expression(path: Path) -> tuple[set[str], dict[tuple[str, str], float], dict[tuple[str, str], int], dict[str, set[str]]]:
    genes: set[str] = set()
    sums: dict[tuple[str, str], float] = defaultdict(float)
    detections: dict[tuple[str, str], int] = defaultdict(int)
    cells_seen: dict[str, set[str]] = defaultdict(set)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_id = row["sample_id"]
            cell_id = row["cell_id"]
            gene = row["gene_symbol"]
            value = float(row["expression_value"])
            genes.add(gene)
            sums[(sample_id, gene)] += value
            detections[(sample_id, gene)] += 1
            cells_seen[sample_id].add(cell_id)
    return genes, sums, detections, cells_seen


def signature_denominator(signature: str, genes: set[str], sample_cells: int) -> int:
    present = len([gene for gene in SIGNATURES[signature] if gene in genes])
    return max(1, present * sample_cells)


def build_rows(
    clinical: dict[str, dict[str, str]],
    genes: set[str],
    sums: dict[tuple[str, str], float],
    detections: dict[tuple[str, str], int],
    cells_seen: dict[str, set[str]],
    cell_counts: dict[str, int],
) -> tuple[list[dict[str, str]], list[str]]:
    sorted_genes = sorted(genes)
    fieldnames = [
        "dataset_id",
        "dataset_unit_id",
        "unit_type",
        "sample_id",
        "donor_id",
        "diagnosis",
        "label__diagnosis_binary",
        "age",
        "sex",
        "pmi",
        "rin",
        "braak",
        "cerad",
        "lewy_body_midbrain",
        "lewy_body_limbic",
        "lewy_body_neocortical",
        "n_cells",
        *[f"gene_mean__{safe_name(gene)}" for gene in sorted_genes],
        *[f"gene_detection__{safe_name(gene)}" for gene in sorted_genes],
        *[f"index__{safe_name(signature)}" for signature in sorted(SIGNATURES)],
    ]
    rows: list[dict[str, str]] = []
    for sample_id in sorted(clinical):
        metadata = clinical[sample_id]
        n_cells = cell_counts.get(sample_id) or len(cells_seen.get(sample_id, set()))
        n_cells = max(1, n_cells)
        row = {
            "dataset_id": "gse243639_pd_snpc",
            "dataset_unit_id": sample_id,
            "unit_type": "sample",
            "sample_id": sample_id,
            "donor_id": sample_id,
            "diagnosis": metadata.get("diagnosis", ""),
            "label__diagnosis_binary": metadata.get("label__diagnosis_binary", ""),
            "age": metadata.get("age", ""),
            "sex": metadata.get("sex", ""),
            "pmi": metadata.get("pmi", ""),
            "rin": metadata.get("rin", ""),
            "braak": metadata.get("braak", ""),
            "cerad": metadata.get("cerad", ""),
            "lewy_body_midbrain": metadata.get("lewy_body_midbrain", ""),
            "lewy_body_limbic": metadata.get("lewy_body_limbic", ""),
            "lewy_body_neocortical": metadata.get("lewy_body_neocortical", ""),
            "n_cells": str(n_cells),
        }
        for gene in sorted_genes:
            row[f"gene_mean__{safe_name(gene)}"] = f"{sums.get((sample_id, gene), 0.0) / n_cells:.8g}"
            row[f"gene_detection__{safe_name(gene)}"] = f"{detections.get((sample_id, gene), 0) / n_cells:.8g}"
        for signature in sorted(SIGNATURES):
            total = sum(sums.get((sample_id, gene), 0.0) for gene in SIGNATURES[signature])
            row[f"index__{safe_name(signature)}"] = f"{total / signature_denominator(signature, genes, n_cells):.8g}"
        rows.append(row)
    return rows, fieldnames


def read_schema_fields(path: Path) -> set[str]:
    if not path.exists():
        return set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        return set(next(reader, []))


def schema_alignment(fieldnames: list[str], phase5_schema_path: Path) -> list[dict[str, str]]:
    pd_fields = set(fieldnames)
    phase5_fields = read_schema_fields(phase5_schema_path)
    rows = []
    for field in sorted(pd_fields | phase5_fields):
        rows.append(
            {
                "feature": field,
                "in_gse243639": str(field in pd_fields).lower(),
                "in_sea_ad_phase5": str(field in phase5_fields).lower(),
                "status": "shared" if field in pd_fields and field in phase5_fields else "schema_specific",
            }
        )
    return rows


def label_summary(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    counts: Counter[tuple[str, str]] = Counter()
    for row in rows:
        counts[("diagnosis", row.get("diagnosis", ""))] += 1
        counts[("label__diagnosis_binary", row.get("label__diagnosis_binary", ""))] += 1
    return [
        {"label_field": field, "label": label, "sample_count": str(count)}
        for (field, label), count in sorted(counts.items())
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build GSE243639 sample-level NeuroFate feature table.")
    parser.add_argument("--expression", type=Path, default=Path("data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz"))
    parser.add_argument("--clinical", type=Path, default=Path("data/raw/external/gse243639_pd_snpc/GSE243639_Clinical_data.csv.gz"))
    parser.add_argument("--cell-map", type=Path, default=Path("data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv"))
    parser.add_argument("--phase5-schema", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase16_gse243639_feature_table.tsv"))
    parser.add_argument("--schema-output", type=Path, default=Path("results/tables/phase16_gse243639_feature_schema_alignment.tsv"))
    parser.add_argument("--label-summary-output", type=Path, default=Path("results/tables/phase16_gse243639_label_summary.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/79_build_gse243639_feature_table.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    clinical = read_clinical(args.clinical, header_line=6, delimiter=";")
    genes, sums, detections, cells_seen = stream_expression(args.expression)
    cell_counts = read_cell_counts(args.cell_map)
    rows, fieldnames = build_rows(clinical, genes, sums, detections, cells_seen, cell_counts)
    write_tsv(args.output, rows, fieldnames)
    write_tsv(args.schema_output, schema_alignment(fieldnames, args.phase5_schema), ["feature", "in_gse243639", "in_sea_ad_phase5", "status"])
    write_tsv(args.label_summary_output, label_summary(rows), ["label_field", "label", "sample_count"])
    logging.info("GSE243639 sample-level feature rows: %d", len(rows))
    logging.info("No model training was performed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
