#!/usr/bin/env python3
"""Parse sample metadata from a PD GEO series matrix without extracting expression."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from pathlib import Path


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") if path.name.endswith(".gz") else path.open("r", encoding="utf-8", errors="replace", newline="")


def clean(value: str | None) -> str:
    text = str(value or "").strip()
    if len(text) >= 2 and text[0] == text[-1] == '"':
        text = text[1:-1]
    return text.strip()


def norm(value: str | None) -> str:
    return clean(value).casefold().replace("_", " ").replace("-", " ")


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def collect_sample_rows(path: Path) -> dict[str, list[list[str]]]:
    rows: dict[str, list[list[str]]] = {}
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            key = clean(row[0])
            if key == "!series_matrix_table_begin":
                break
            if key.startswith("!Sample_"):
                rows.setdefault(key, []).append([clean(value) for value in row[1:]])
    return rows


def first(rows: dict[str, list[list[str]]], key: str) -> list[str]:
    return rows.get(key, [[]])[0]


def infer_fields(characteristics: list[str], title: str, source: str) -> tuple[str, str]:
    disease = ""
    region = ""
    for item in characteristics:
        if ":" in item:
            key, value = item.split(":", 1)
            nkey = norm(key)
            if any(token in nkey for token in ["disease", "diagnosis", "condition", "status"]):
                disease = clean(value)
            if any(token in nkey for token in ["region", "tissue", "brain"]):
                region = clean(value)
    blob = norm(" ".join([title, source, disease, *characteristics]))
    if not disease:
        if "control" in blob or "normal" in blob or "unaffected" in blob:
            disease = "Control"
        elif "parkinson" in blob or " pd" in f" {blob} ":
            disease = "Parkinson's disease"
    if not region:
        if "substantia nigra" in blob or "snpc" in blob:
            region = "substantia nigra"
        elif "frontal cortex" in blob:
            region = "frontal cortex"
    return disease, region


def label_for(disease: str, title: str, source: str) -> str:
    blob = norm(" ".join([disease, title, source]))
    if "control" in blob or "normal" in blob or "unaffected" in blob:
        return "0"
    if "parkinson" in blob or " pd" in f" {blob} ":
        return "1"
    return ""


def parse_metadata(path: Path, cohort_id: str) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    rows = collect_sample_rows(path)
    titles = first(rows, "!Sample_title")
    accessions = first(rows, "!Sample_geo_accession")
    sources = first(rows, "!Sample_source_name_ch1")
    platforms = first(rows, "!Sample_platform_id")
    supplement_rows = [values for key, values_list in rows.items() if key.startswith("!Sample_supplementary_file") for values in values_list]
    sample_count = max([len(titles), len(accessions), len(sources), len(platforms), *[len(row) for row in rows.get("!Sample_characteristics_ch1", [])], *[len(row) for row in supplement_rows]] or [0])
    metadata: list[dict[str, str]] = []
    platform_counts: dict[str, int] = {}
    label_counts: dict[str, int] = {}
    for index in range(sample_count):
        title = titles[index] if index < len(titles) else ""
        geo = accessions[index] if index < len(accessions) else ""
        source = sources[index] if index < len(sources) else ""
        platform = platforms[index] if index < len(platforms) else ""
        characteristics = [values[index] for values in rows.get("!Sample_characteristics_ch1", []) if index < len(values)]
        disease, region = infer_fields(characteristics, title, source)
        label = label_for(disease, title, source)
        supplements = [values[index] for values in supplement_rows if index < len(values) and values[index]]
        row = {
            "cohort_id": cohort_id,
            "sample_id": geo or title,
            "sample_title": title,
            "geo_accession": geo,
            "source_name": source,
            "disease_state": disease,
            "tissue_or_region": region,
            "platform_id": platform,
            "supplementary_files": ";".join(supplements),
            "label__pd_vs_control": label,
            "endpoint_status": "unambiguous" if label in {"0", "1"} else "ambiguous_or_missing",
        }
        metadata.append(row)
        platform_counts[platform or "missing"] = platform_counts.get(platform or "missing", 0) + 1
        label_counts[label or "missing"] = label_counts.get(label or "missing", 0) + 1
    label_summary = [{"cohort_id": cohort_id, "label__pd_vs_control": key, "count": str(value)} for key, value in sorted(label_counts.items())]
    platform_summary = [{"cohort_id": cohort_id, "platform_id": key, "sample_count": str(value)} for key, value in sorted(platform_counts.items())]
    return metadata, label_summary, platform_summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Parse Phase 34 PD GEO series matrix sample metadata.")
    parser.add_argument("--series-matrix", type=Path, required=True)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--label-summary-output", type=Path, required=True)
    parser.add_argument("--platform-output", type=Path, required=True)
    parser.add_argument("--log-file", type=Path, required=True)
    args = parser.parse_args()
    configure_logging(args.log_file)
    metadata, labels, platforms = parse_metadata(args.series_matrix, args.cohort_id)
    write_tsv(args.output, metadata, ["cohort_id", "sample_id", "geo_accession", "sample_title", "source_name", "disease_state", "tissue_or_region", "platform_id", "supplementary_files", "label__pd_vs_control", "endpoint_status"])
    write_tsv(args.label_summary_output, labels, ["cohort_id", "label__pd_vs_control", "count"])
    write_tsv(args.platform_output, platforms, ["cohort_id", "platform_id", "sample_count"])
    logging.info("Parsed Phase 34 metadata cohort=%s samples=%d", args.cohort_id, len(metadata))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
