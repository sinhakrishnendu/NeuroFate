#!/usr/bin/env python3
"""Parse GSE184950 GEO series-matrix sample metadata without opening expression files."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from pathlib import Path
from urllib.parse import urlparse


METADATA_COLUMNS = [
    "sample_name",
    "title",
    "geo_accession",
    "source_name",
    "tissue",
    "disease_state",
    "donor_id",
    "age",
    "gender",
    "race",
    "ethnicity",
    "pmi_hours",
    "braak_stage",
    "supplementary_file_1",
    "processed_tar_name",
    "expected_archive_member",
    "label__pd_pdd_vs_control",
]

MANIFEST_COLUMNS = [
    "sample_name",
    "geo_accession",
    "disease_state",
    "supplementary_file_1",
    "processed_tar_name",
    "expected_archive_member",
    "status",
]


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
    value = (value or "").strip()
    if len(value) >= 2 and value[0] == value[-1] == '"':
        value = value[1:-1]
    return value.strip()


def norm_key(value: str) -> str:
    return clean(value).lower().replace("_", " ").replace("-", " ")


def parse_series_rows(path: Path) -> dict[str, list[list[str]]]:
    rows: dict[str, list[list[str]]] = {}
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for raw in reader:
            if not raw:
                continue
            key = clean(raw[0])
            if not key.startswith("!Sample_"):
                continue
            rows.setdefault(key, []).append([clean(value) for value in raw[1:]])
    return rows


def first_row(rows: dict[str, list[list[str]]], key: str) -> list[str]:
    return rows.get(key, [[]])[0]


def basename_from_url(value: str) -> str:
    if not value:
        return ""
    path = urlparse(value).path if "://" in value else value
    return Path(path).name


def characteristic_field_name(raw_key: str) -> str | None:
    mapping = {
        "tissue": "tissue",
        "disease state": "disease_state",
        "brain bank donor id": "donor_id",
        "age": "age",
        "gender": "gender",
        "race": "race",
        "ethnicity": "ethnicity",
        "postmortem interval hours": "pmi_hours",
        "braak stage": "braak_stage",
    }
    return mapping.get(norm_key(raw_key))


def disease_label(value: str) -> str:
    normalized = norm_key(value).replace("'", "")
    if normalized == "unaffected control" or "control" in normalized:
        return "0"
    if normalized in {"parkinsons disease", "parkinsons disease dementia"}:
        return "1"
    if "parkinson" in normalized:
        return "1"
    return ""


def parse_characteristics(rows: dict[str, list[list[str]]], sample_count: int) -> list[dict[str, str]]:
    parsed = [dict.fromkeys(["tissue", "disease_state", "donor_id", "age", "gender", "race", "ethnicity", "pmi_hours", "braak_stage"], "") for _ in range(sample_count)]
    for values in rows.get("!Sample_characteristics_ch1", []):
        for index, value in enumerate(values[:sample_count]):
            if ":" not in value:
                continue
            raw_key, raw_value = value.split(":", 1)
            field = characteristic_field_name(raw_key)
            if field:
                parsed[index][field] = clean(raw_value)
    return parsed


def build_metadata(path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    rows = parse_series_rows(path)
    titles = first_row(rows, "!Sample_title")
    accessions = first_row(rows, "!Sample_geo_accession")
    sources = first_row(rows, "!Sample_source_name_ch1")
    supplementary = first_row(rows, "!Sample_supplementary_file_1")
    sample_count = max(len(titles), len(accessions), len(sources), len(supplementary))
    characteristics = parse_characteristics(rows, sample_count)

    metadata_rows: list[dict[str, str]] = []
    manifest_rows: list[dict[str, str]] = []
    for index in range(sample_count):
        title = titles[index] if index < len(titles) else ""
        geo = accessions[index] if index < len(accessions) else ""
        source = sources[index] if index < len(sources) else ""
        link = supplementary[index] if index < len(supplementary) else ""
        processed_tar = basename_from_url(link)
        sample_name = title or processed_tar.replace(".tar.gz", "")
        disease = characteristics[index].get("disease_state", "")
        label = disease_label(disease)
        row = {
            "sample_name": sample_name,
            "title": title,
            "geo_accession": geo,
            "source_name": source,
            **characteristics[index],
            "supplementary_file_1": link,
            "processed_tar_name": processed_tar,
            "expected_archive_member": processed_tar,
            "label__pd_pdd_vs_control": label,
        }
        metadata_rows.append(row)
        manifest_rows.append(
            {
                "sample_name": sample_name,
                "geo_accession": geo,
                "disease_state": disease,
                "supplementary_file_1": link,
                "processed_tar_name": processed_tar,
                "expected_archive_member": processed_tar,
                "status": "expected_from_series_matrix" if processed_tar else "missing_supplementary_tar",
            }
        )

    label_counts: dict[str, int] = {}
    positive = 0
    negative = 0
    for row in metadata_rows:
        disease = row.get("disease_state", "") or "missing"
        label_counts[disease] = label_counts.get(disease, 0) + 1
        if row.get("label__pd_pdd_vs_control") == "1":
            positive += 1
        if row.get("label__pd_pdd_vs_control") == "0":
            negative += 1
    label_rows = [{"category": key, "count": str(value)} for key, value in sorted(label_counts.items())]
    label_rows.extend(
        [
            {"category": "combined_positive_pd_pdd", "count": str(positive)},
            {"category": "negative_unaffected_control", "count": str(negative)},
            {"category": "total_samples", "count": str(len(metadata_rows))},
        ]
    )
    return metadata_rows, manifest_rows, label_rows


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Parse GSE184950 GEO series-matrix metadata.")
    parser.add_argument("--series-matrix", type=Path, default=Path("data/raw/external/gse184950_pd_sn/GSE184950_series_matrix.txt.gz"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase25_gse184950_series_sample_metadata.tsv"))
    parser.add_argument("--processed-files-output", type=Path, default=Path("results/tables/phase25_gse184950_series_processed_file_manifest.tsv"))
    parser.add_argument("--label-summary-output", type=Path, default=Path("results/tables/phase25_gse184950_label_summary.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/131_parse_gse184950_series_matrix.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    metadata, manifest, labels = build_metadata(args.series_matrix)
    write_tsv(args.output, metadata, METADATA_COLUMNS)
    write_tsv(args.processed_files_output, manifest, MANIFEST_COLUMNS)
    write_tsv(args.label_summary_output, labels, ["category", "count"])
    logging.info("Parsed GSE184950 series matrix samples=%d", len(metadata))
    logging.info("No expression files, RAW archives, or matrices were opened.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
