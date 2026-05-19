#!/usr/bin/env python3
"""Build donor/sample-level external feature tables from extracted target-gene TSVs."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from collections import defaultdict
from pathlib import Path
from statistics import mean


MIN_RELIABLE_UNITS = 20


def setup_logging(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=path, level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def canonical_mapping(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    mapping: dict[str, str] = {}
    for row in read_tsv(path):
        if row.get("canonical_field") != "unmapped":
            mapping[row["canonical_field"]] = row["source_field"]
    return mapping


def schema_features(path: Path) -> list[str]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
    return [col for col in header if col.startswith(("gene_mean__", "gene_detection__", "index__", "cell_fraction__", "celltype_index__"))]


def build_features(expression_rows: list[dict[str, str]], metadata_rows: list[dict[str, str]], mapping: dict[str, str]) -> tuple[list[dict[str, str]], str]:
    cell_field = mapping.get("cell_id", "cell_id")
    donor_field = mapping.get("donor_id") or mapping.get("sample_id") or "sample_id"
    diagnosis_field = mapping.get("diagnosis") or mapping.get("disease_status")
    cell_to_unit: dict[str, str] = {}
    unit_labels: dict[str, str] = {}
    for row in metadata_rows:
        cell_id = row.get(cell_field) or row.get("cell_id") or row.get("barcode")
        unit = row.get(donor_field) or row.get("sample_id") or row.get("donor_id") or "pseudo_donor"
        if cell_id:
            cell_to_unit[cell_id] = unit
        if diagnosis_field and row.get(diagnosis_field):
            unit_labels[unit] = row[diagnosis_field]
    values: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in expression_rows:
        cell_id = row.get("cell_id") or row.get("sample_id")
        unit = cell_to_unit.get(cell_id, cell_id or "bulk_sample")
        gene = row.get("gene_symbol", "")
        try:
            value = float(row.get("expression_value", "0") or 0)
        except ValueError:
            value = 0.0
        if gene:
            values[(unit, gene)].append(value)
    units = sorted({unit for unit, _gene in values})
    genes = sorted({gene for _unit, gene in values})
    unit_type = "donor" if "donor_id" in mapping else "sample" if "sample_id" in mapping else "pseudo-donor"
    output_rows: list[dict[str, str]] = []
    for unit in units:
        row = {
            "dataset_unit_id": unit,
            "unit_type": unit_type,
            "label__diagnosis": unit_labels.get(unit, "unmapped"),
        }
        for gene in genes:
            gene_values = values.get((unit, gene), [])
            row[f"gene_mean__{gene}"] = f"{mean(gene_values):.6g}" if gene_values else "0"
            row[f"gene_detection__{gene}"] = f"{sum(1 for value in gene_values if value > 0) / len(gene_values):.6g}" if gene_values else "0"
        output_rows.append(row)
    return output_rows, unit_type


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build generic external donor/sample feature table.")
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--expression", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--schema", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--schema-output", type=Path, required=True)
    parser.add_argument("--label-summary-output", type=Path, required=True)
    parser.add_argument("--log-file", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log_file)
    mapping = canonical_mapping(args.mapping)
    features, unit_type = build_features(read_tsv(args.expression), read_tsv(args.metadata), mapping)
    all_fields = sorted({field for row in features for field in row})
    fieldnames = ["dataset_unit_id", "unit_type", "label__diagnosis"] + [f for f in all_fields if f not in {"dataset_unit_id", "unit_type", "label__diagnosis"}]
    write_tsv(args.output, features, fieldnames)
    sea_ad_features = set(schema_features(args.schema))
    external_features = {field for field in fieldnames if field.startswith(("gene_mean__", "gene_detection__", "index__", "cell_fraction__", "celltype_index__"))}
    alignment_rows = [
        {"feature": feature, "status": "shared" if feature in sea_ad_features else "external_specific"}
        for feature in sorted(external_features | sea_ad_features)
    ]
    write_tsv(args.schema_output, alignment_rows, ["feature", "status"])
    label_counts: dict[str, int] = {}
    for row in features:
        label_counts[row["label__diagnosis"]] = label_counts.get(row["label__diagnosis"], 0) + 1
    label_rows = [{"label_field": "label__diagnosis", "label": label, "sample_count": str(count)} for label, count in sorted(label_counts.items())]
    write_tsv(args.label_summary_output, label_rows, ["label_field", "label", "sample_count"])
    if len(features) < MIN_RELIABLE_UNITS:
        logging.warning("External feature table has small n=%s and should be treated as feasibility.", len(features))
    logging.info("Built %s feature table rows using unit_type=%s", len(features), unit_type)
    print(f"Wrote {args.output}")
    print(f"Wrote {args.schema_output}")
    print(f"Wrote {args.label_summary_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
