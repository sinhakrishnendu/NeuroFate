#!/usr/bin/env python3
"""Parse GEO series-matrix sample metadata for AD replication cohorts."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from pathlib import Path
from urllib.parse import urlparse


METADATA_COLUMNS = [
    "cohort_id",
    "sample_title",
    "geo_accession",
    "source_name",
    "diagnosis",
    "disease_state",
    "condition",
    "brain_region",
    "age",
    "sex",
    "pmi",
    "braak",
    "ad_pathology",
    "supplementary_files",
    "inferred_ad_endpoint",
]

MANIFEST_COLUMNS = ["cohort_id", "sample_title", "geo_accession", "supplementary_file", "supplementary_file_name", "status"]


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.name.endswith(".gz") else path.open("r", encoding="utf-8", errors="replace")


def clean(value: str) -> str:
    value = str(value or "").strip()
    if len(value) >= 2 and value[0] == value[-1] == '"':
        value = value[1:-1]
    return value.strip()


def norm(value: str) -> str:
    return clean(value).lower().replace("_", " ").replace("-", " ")


def parse_series_rows(path: Path) -> dict[str, list[list[str]]]:
    rows: dict[str, list[list[str]]] = {}
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for raw in reader:
            if not raw:
                continue
            key = clean(raw[0])
            if key.startswith("!Sample_"):
                rows.setdefault(key, []).append([clean(value) for value in raw[1:]])
    return rows


def first_row(rows: dict[str, list[list[str]]], key: str) -> list[str]:
    return rows.get(key, [[]])[0]


def basename_from_url(value: str) -> str:
    if not value:
        return ""
    path = urlparse(value).path if "://" in value else value
    return Path(path).name


def canonical_characteristic_key(raw_key: str) -> str | None:
    key = norm(raw_key)
    if any(token in key for token in ["diagnosis", "diagnostic"]):
        return "diagnosis"
    if "disease state" in key or key == "disease":
        return "disease_state"
    if "condition" in key or "phenotype" in key:
        return "condition"
    if "brain region" in key or "region" in key or "tissue" in key:
        return "brain_region"
    if key == "age" or "age at death" in key:
        return "age"
    if key in {"sex", "gender"}:
        return "sex"
    if "postmortem" in key or key == "pmi" or "pmi" in key:
        return "pmi"
    if "braak" in key:
        return "braak"
    if "ad pathology" in key or "neuropath" in key or "cerad" in key:
        return "ad_pathology"
    return None


def parse_characteristics(rows: dict[str, list[list[str]]], sample_count: int) -> list[dict[str, str]]:
    parsed = [dict.fromkeys(["diagnosis", "disease_state", "condition", "brain_region", "age", "sex", "pmi", "braak", "ad_pathology"], "") for _ in range(sample_count)]
    for values in rows.get("!Sample_characteristics_ch1", []):
        for index, value in enumerate(values[:sample_count]):
            if ":" not in value:
                continue
            raw_key, raw_value = value.split(":", 1)
            key = canonical_characteristic_key(raw_key)
            if key and not parsed[index].get(key):
                parsed[index][key] = clean(raw_value)
    return parsed


def infer_ad_endpoint(row: dict[str, str]) -> str:
    combined = " ".join(row.get(key, "") for key in ["diagnosis", "disease_state", "condition", "ad_pathology"]).lower()
    if not combined.strip():
        return ""
    if any(token in combined for token in ["control", "normal", "non-demented", "nondemented", "unaffected"]):
        return "control"
    if "alzheimer" in combined or "ad" == combined.strip() or "dementia" in combined or "high pathology" in combined:
        return "ad"
    if "low pathology" in combined:
        return "control_or_low_pathology"
    return ""


def build_metadata(path: Path, cohort_id: str) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    rows = parse_series_rows(path)
    titles = first_row(rows, "!Sample_title")
    accessions = first_row(rows, "!Sample_geo_accession")
    sources = first_row(rows, "!Sample_source_name_ch1")
    supplementary_rows = [values for key, value_rows in rows.items() if key.startswith("!Sample_supplementary_file") for values in value_rows]
    sample_count = max([len(titles), len(accessions), len(sources), *[len(row) for row in supplementary_rows]] or [0])
    characteristics = parse_characteristics(rows, sample_count)
    metadata: list[dict[str, str]] = []
    manifest: list[dict[str, str]] = []
    for index in range(sample_count):
        title = titles[index] if index < len(titles) else ""
        geo = accessions[index] if index < len(accessions) else ""
        source = sources[index] if index < len(sources) else ""
        supplements = [row[index] for row in supplementary_rows if index < len(row) and row[index]]
        row = {
            "cohort_id": cohort_id,
            "sample_title": title,
            "geo_accession": geo,
            "source_name": source,
            **characteristics[index],
            "supplementary_files": ";".join(supplements),
        }
        row["inferred_ad_endpoint"] = infer_ad_endpoint(row)
        metadata.append(row)
        for supplement in supplements:
            manifest.append(
                {
                    "cohort_id": cohort_id,
                    "sample_title": title,
                    "geo_accession": geo,
                    "supplementary_file": supplement,
                    "supplementary_file_name": basename_from_url(supplement),
                    "status": "declared_in_series_matrix",
                }
            )
    counts: dict[str, int] = {}
    for row in metadata:
        label = row.get("inferred_ad_endpoint", "") or "unmapped"
        counts[label] = counts.get(label, 0) + 1
    label_summary = [{"cohort_id": cohort_id, "label": label, "count": str(count)} for label, count in sorted(counts.items())]
    return metadata, manifest, label_summary


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Parse GEO series-matrix metadata for an AD replication cohort.")
    parser.add_argument("--series-matrix", type=Path, required=True)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--file-manifest-output", type=Path)
    parser.add_argument("--label-summary-output", type=Path)
    parser.add_argument("--log-file", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_file = args.log_file or Path(f"results/logs/143_{args.cohort_id}_series_matrix.log")
    configure_logging(log_file)
    metadata, manifest, labels = build_metadata(args.series_matrix, args.cohort_id)
    output = args.output or Path(f"results/tables/phase28_{args.cohort_id}_sample_metadata.tsv")
    manifest_output = args.file_manifest_output or Path(f"results/tables/phase28_{args.cohort_id}_supplementary_file_manifest.tsv")
    label_output = args.label_summary_output or Path(f"results/tables/phase28_{args.cohort_id}_label_summary.tsv")
    write_tsv(output, metadata, METADATA_COLUMNS)
    write_tsv(manifest_output, manifest, MANIFEST_COLUMNS)
    write_tsv(label_output, labels, ["cohort_id", "label", "count"])
    logging.info("Parsed GEO series matrix cohort=%s samples=%d supplementary_files=%d", args.cohort_id, len(metadata), len(manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
