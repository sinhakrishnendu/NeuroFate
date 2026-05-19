#!/usr/bin/env python3
"""Build repaired sample-level, cell-type-aware GSE243639 NeuroFate features."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import TextIO


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
    "microglial_activation_index": ["TREM2", "TYROBP", "GPNMB", "HLA-DRA", "AIF1"],
    "astrocyte_stress_index": ["GFAP"],
    "neuronal_vulnerability_index": ["SLC17A7", "SST", "PVALB", "LAMP5", "NEFL", "NEFM"],
    "myelin_oligodendrocyte_index": ["MBP", "MOBP", "PLP1"],
    "synuclein_axis_index": ["SNCA", "LRRK2", "PINK1", "PRKN", "MAPT", "APOE"],
}
FEATURE_PREFIXES = (
    "gene_mean__",
    "gene_detection__",
    "cell_fraction__",
    "celltype_gene_mean__",
    "celltype_gene_detection__",
    "index__",
)


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_") or "missing"


def split_id(value: str) -> tuple[str, str]:
    stripped = value.strip().strip("\"'")
    if "_" in stripped:
        sample_id, barcode = stripped.split("_", 1)
        return sample_id, barcode
    return "", stripped


def remove_trailing(value: str) -> str:
    return re.sub(r"([.-])\d+$", "", value.strip())


def collapse_punctuation(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "", value.strip()).upper()


def normalize_by_rule(value: str, rule: str) -> str:
    if rule in {"", "raw_id", "direct_id_linkage_safe"}:
        return value.strip()
    if rule == "lowercase":
        return value.lower()
    if rule == "remove_quotes":
        return value.strip("\"'")
    if rule == "replace_dash_with_dot":
        return value.replace("-", ".")
    if rule == "replace_dot_with_dash":
        return value.replace(".", "-")
    if rule in {"remove_trailing_dot_or_dash_one", "remove_trailing_dot_suffix"}:
        return remove_trailing(value)
    if rule == "keep_barcode_only":
        return remove_trailing(split_id(value)[1])
    if rule == "keep_sample_prefix_only":
        return split_id(value)[0]
    if rule == "remove_sample_prefix":
        return split_id(value)[1]
    if rule == "collapse_punctuation":
        return collapse_punctuation(value)
    return value.strip()


def choose_join_column(fieldnames: list[str]) -> str:
    for candidate in ["cell_id_expression", "cell_id", "normalized_cell_id"]:
        if candidate in fieldnames:
            return candidate
    return ""


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
            mapped = {canonical: (raw.get(source) or "").strip() for source, canonical in CLINICAL_FIELD_MAP.items()}
            sample_id = mapped.get("sample_id", "")
            if sample_id:
                mapped["label__diagnosis_binary"] = diagnosis_binary(mapped.get("diagnosis", ""))
                rows[sample_id] = mapped
    return rows


def read_annotations(path: Path) -> tuple[dict[str, dict[str, str]], Counter[tuple[str, str]], Counter[str], str, str]:
    annotation_by_expression_id: dict[str, dict[str, str]] = {}
    celltype_counts: Counter[tuple[str, str]] = Counter()
    sample_counts: Counter[str] = Counter()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        join_column = choose_join_column(fieldnames)
        safe_map_mode = "normalization_rule" in fieldnames or "phase19_linkage_decision" in fieldnames
        normalization_rule = "raw_id"
        for row in reader:
            if "normalization_rule" in row and row.get("normalization_rule"):
                normalization_rule = row["normalization_rule"]
            if safe_map_mode and row.get("match_status", "") != "matched":
                continue
            if not safe_map_mode and row.get("match_status", "") == "unmatched":
                continue
            expression_id = row.get(join_column, "") if join_column else ""
            sample_id = row.get("sample_id", "")
            cell_type = row.get("cell_type", "") or (f"cluster_{safe_name(row.get('cluster_id', ''))}" if row.get("cluster_id") else "unmatched")
            if not expression_id or not sample_id:
                continue
            sample_counts[sample_id] += 1
            if cell_type != "unmatched":
                celltype_counts[(sample_id, cell_type)] += 1
            annotation_by_expression_id[expression_id] = {
                "sample_id": sample_id,
                "cell_type": cell_type,
                "match_status": "matched",
                "biological_celltype_confidence": row.get("biological_celltype_confidence", ""),
                "join_column": join_column,
                "normalization_rule": normalization_rule,
            }
    return annotation_by_expression_id, celltype_counts, sample_counts, join_column, normalization_rule


def stream_expression(
    expression_path: Path,
    annotations: dict[str, dict[str, str]],
    join_column: str,
    normalization_rule: str,
) -> tuple[
    set[str],
    set[str],
    dict[tuple[str, str], float],
    dict[tuple[str, str], int],
    dict[tuple[str, str, str], float],
    dict[tuple[str, str, str], int],
    set[str],
    set[str],
]:
    genes: set[str] = set()
    celltypes: set[str] = set()
    sample_gene_sum: dict[tuple[str, str], float] = defaultdict(float)
    sample_gene_detection: dict[tuple[str, str], int] = defaultdict(int)
    sample_celltype_gene_sum: dict[tuple[str, str, str], float] = defaultdict(float)
    sample_celltype_gene_detection: dict[tuple[str, str, str], int] = defaultdict(int)
    unmatched_cells: set[str] = set()
    matched_cells: set[str] = set()
    with open_text(expression_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cell_id = row["cell_id"]
            sample_id = row["sample_id"]
            gene = row["gene_symbol"]
            value = float(row["expression_value"])
            genes.add(gene)
            sample_gene_sum[(sample_id, gene)] += value
            sample_gene_detection[(sample_id, gene)] += 1
            lookup_id = normalize_by_rule(cell_id, normalization_rule) if join_column == "normalized_cell_id" else cell_id
            annotation = annotations.get(lookup_id)
            if annotation is None or annotation.get("match_status") == "unmatched":
                unmatched_cells.add(cell_id)
                continue
            matched_cells.add(cell_id)
            cell_type = annotation["cell_type"]
            celltypes.add(cell_type)
            sample_celltype_gene_sum[(sample_id, cell_type, gene)] += value
            sample_celltype_gene_detection[(sample_id, cell_type, gene)] += 1
    return genes, celltypes, sample_gene_sum, sample_gene_detection, sample_celltype_gene_sum, sample_celltype_gene_detection, unmatched_cells, matched_cells


def present_gene_count(signature: str, genes: set[str]) -> int:
    return max(1, len([gene for gene in SIGNATURES[signature] if gene in genes]))


def build_rows(
    clinical: dict[str, dict[str, str]],
    genes: set[str],
    celltypes: set[str],
    sample_gene_sum: dict[tuple[str, str], float],
    sample_gene_detection: dict[tuple[str, str], int],
    sample_celltype_gene_sum: dict[tuple[str, str, str], float],
    sample_celltype_gene_detection: dict[tuple[str, str, str], int],
    celltype_counts: Counter[tuple[str, str]],
    sample_counts: Counter[str],
) -> tuple[list[dict[str, str]], list[str]]:
    sorted_genes = sorted(genes)
    sorted_celltypes = sorted(celltypes)
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
        *[f"cell_fraction__{safe_name(cell_type)}" for cell_type in sorted_celltypes],
        *[f"celltype_gene_mean__{safe_name(cell_type)}__{safe_name(gene)}" for cell_type in sorted_celltypes for gene in sorted_genes],
        *[f"celltype_gene_detection__{safe_name(cell_type)}__{safe_name(gene)}" for cell_type in sorted_celltypes for gene in sorted_genes],
        *[f"index__{safe_name(signature)}" for signature in sorted(SIGNATURES)],
    ]
    rows: list[dict[str, str]] = []
    for sample_id in sorted(clinical):
        metadata = clinical[sample_id]
        n_cells = max(1, sample_counts.get(sample_id, 0))
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
            row[f"gene_mean__{safe_name(gene)}"] = f"{sample_gene_sum.get((sample_id, gene), 0.0) / n_cells:.8g}"
            row[f"gene_detection__{safe_name(gene)}"] = f"{sample_gene_detection.get((sample_id, gene), 0) / n_cells:.8g}"
        for cell_type in sorted_celltypes:
            celltype_count = celltype_counts.get((sample_id, cell_type), 0)
            row[f"cell_fraction__{safe_name(cell_type)}"] = f"{celltype_count / n_cells:.8g}"
            denominator = max(1, celltype_count)
            for gene in sorted_genes:
                row[f"celltype_gene_mean__{safe_name(cell_type)}__{safe_name(gene)}"] = f"{sample_celltype_gene_sum.get((sample_id, cell_type, gene), 0.0) / denominator:.8g}"
                row[f"celltype_gene_detection__{safe_name(cell_type)}__{safe_name(gene)}"] = f"{sample_celltype_gene_detection.get((sample_id, cell_type, gene), 0) / denominator:.8g}"
        for signature in sorted(SIGNATURES):
            total = sum(sample_gene_sum.get((sample_id, gene), 0.0) for gene in SIGNATURES[signature])
            row[f"index__{safe_name(signature)}"] = f"{total / (n_cells * present_gene_count(signature, genes)):.8g}"
        rows.append(row)
    return rows, fieldnames


def read_schema_fields(path: Path) -> set[str]:
    if not path.exists():
        return set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        return set(next(csv.reader(handle, delimiter="\t"), []))


def schema_alignment(fieldnames: list[str], phase5_schema_path: Path) -> list[dict[str, str]]:
    local = set(fieldnames)
    phase5 = read_schema_fields(phase5_schema_path)
    return [
        {
            "feature": field,
            "in_gse243639_celltype": str(field in local).lower(),
            "in_sea_ad_phase5": str(field in phase5).lower(),
            "status": "shared" if field in local and field in phase5 else "schema_specific",
        }
        for field in sorted(local | phase5)
    ]


def label_summary(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    counts: Counter[tuple[str, str]] = Counter()
    for row in rows:
        counts[("diagnosis", row.get("diagnosis", ""))] += 1
        counts[("label__diagnosis_binary", row.get("label__diagnosis_binary", ""))] += 1
    return [{"label_field": field, "label": label, "sample_count": str(count)} for (field, label), count in sorted(counts.items())]


def feature_group_counts(fieldnames: list[str], rows: list[dict[str, str]], match_rate: float, unmatched_cells: set[str]) -> list[dict[str, str]]:
    counts: Counter[str] = Counter()
    for field in fieldnames:
        group = next((prefix.rstrip("_") for prefix in FEATURE_PREFIXES if field.startswith(prefix)), "metadata")
        counts[group] += 1
    return [
        {
            "feature_group": group,
            "feature_count": str(count),
            "sample_rows": str(len(rows)),
            "annotation_match_rate": f"{match_rate:.6g}",
            "unmatched_unique_expression_cells": str(len(unmatched_cells)),
            "warning": "ok" if match_rate >= 0.90 and len(fieldnames) > 20 else "technical_failure_annotation_join_possible",
        }
        for group, count in sorted(counts.items())
    ]


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %s rows: %d", path, len(rows))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build repaired GSE243639 sample-level cell-type-aware feature table.")
    parser.add_argument("--expression", type=Path, default=Path("data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz"))
    parser.add_argument("--annotations", type=Path, default=Path("data/interim/external/gse243639_pd_snpc/gse243639_cell_annotation_map.tsv"))
    parser.add_argument("--clinical", type=Path, default=Path("data/raw/external/gse243639_pd_snpc/GSE243639_Clinical_data.csv.gz"))
    parser.add_argument("--phase5-schema", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase18_gse243639_celltype_feature_table.tsv"))
    parser.add_argument("--schema-output", type=Path, default=Path("results/tables/phase18_gse243639_celltype_schema_alignment.tsv"))
    parser.add_argument("--label-summary-output", type=Path, default=Path("results/tables/phase18_gse243639_celltype_label_summary.tsv"))
    parser.add_argument("--feature-group-output", type=Path, default=Path("results/tables/phase18_gse243639_feature_group_counts.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/85_build_gse243639_celltype_feature_table.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    clinical = read_clinical(args.clinical)
    annotations, celltype_counts, sample_counts, join_column, normalization_rule = read_annotations(args.annotations)
    (
        genes,
        celltypes,
        sample_gene_sum,
        sample_gene_detection,
        sample_celltype_gene_sum,
        sample_celltype_gene_detection,
        unmatched_cells,
        matched_cells,
    ) = stream_expression(args.expression, annotations, join_column, normalization_rule)
    unique_expression_cells = len(matched_cells | unmatched_cells)
    match_rate = len(matched_cells) / unique_expression_cells if unique_expression_cells else 0.0
    rows, fieldnames = build_rows(clinical, genes, celltypes, sample_gene_sum, sample_gene_detection, sample_celltype_gene_sum, sample_celltype_gene_detection, celltype_counts, sample_counts)
    write_tsv(args.output, rows, fieldnames)
    write_tsv(args.schema_output, schema_alignment(fieldnames, args.phase5_schema), ["feature", "in_gse243639_celltype", "in_sea_ad_phase5", "status"])
    write_tsv(args.label_summary_output, label_summary(rows), ["label_field", "label", "sample_count"])
    write_tsv(
        args.feature_group_output,
        feature_group_counts(fieldnames, rows, match_rate, unmatched_cells),
        ["feature_group", "feature_count", "sample_rows", "annotation_match_rate", "unmatched_unique_expression_cells", "warning"],
    )
    logging.info("Unmatched unique expression cells without annotation: %d", len(unmatched_cells))
    logging.info("Matched unique expression cells with annotation: %d", len(matched_cells))
    logging.info("Annotation join column: %s; normalization rule: %s", join_column, normalization_rule)
    logging.info("All repaired Phase 18 features are sample-level aggregates.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
