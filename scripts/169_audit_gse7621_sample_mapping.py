#!/usr/bin/env python3
"""Audit GSE7621 series-matrix sample IDs against parsed metadata.

This audit reads only the GEO series-matrix table header. It does not read
genome-wide expression rows or write expression output.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from pathlib import Path


JOIN_KEYS = ["geo_accession", "sample_id", "sample_title", "source_name"]


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def clean(value: str | None) -> str:
    text = str(value or "").strip()
    if len(text) >= 2 and text[0] == text[-1] == '"':
        text = text[1:-1]
    return text.strip()


def normalize(value: str | None) -> str:
    return "".join(ch for ch in clean(value).casefold() if ch.isalnum())


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def expression_sample_ids(series_matrix: Path) -> list[str]:
    in_table = False
    with open_text(series_matrix) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            key = clean(row[0])
            if key == "!series_matrix_table_begin":
                in_table = True
                continue
            if not in_table:
                continue
            if key.lower() in {"id_ref", "id", "probe", "probe_id", "gene"}:
                return [clean(value) for value in row[1:]]
            raise SystemExit("Found series matrix table but did not find an ID_REF/header row.")
    raise SystemExit("No !series_matrix_table_begin block found in series matrix.")


def mapping_for(rows: list[dict[str, str]], key: str, normalized: bool = False) -> dict[str, dict[str, str]]:
    mapped: dict[str, dict[str, str]] = {}
    for row in rows:
        value = normalize(row.get(key, "")) if normalized else clean(row.get(key, ""))
        if value:
            mapped[value] = row
    return mapped


def count_labels(rows: list[dict[str, str]]) -> tuple[int, int]:
    pd = sum(1 for row in rows if clean(row.get("label__pd_vs_control")) == "1")
    control = sum(1 for row in rows if clean(row.get("label__pd_vs_control")) == "0")
    return pd, control


def audit(series_matrix: Path, metadata_path: Path) -> tuple[list[dict[str, str]], str, list[str]]:
    expr_samples = expression_sample_ids(series_matrix)
    metadata = read_tsv(metadata_path)
    rows: list[dict[str, str]] = []
    best_key = ""
    best_matched = -1
    for key in JOIN_KEYS:
        for normalized in [False, True]:
            label = f"normalized_{key}" if normalized else key
            metadata_map = mapping_for(metadata, key, normalized)
            expr_keys = [normalize(sample) if normalized else sample for sample in expr_samples]
            matched = sum(1 for sample in expr_keys if sample in metadata_map)
            if matched > best_matched:
                best_matched = matched
                best_key = label
            rows.append(
                {
                    "candidate_join_key": label,
                    "expression_sample_count": str(len(expr_samples)),
                    "metadata_sample_count": str(len(metadata)),
                    "matched_sample_count": str(matched),
                    "unmatched_expression_samples": str(len(expr_samples) - matched),
                    "unmatched_metadata_samples": str(max(0, len(metadata) - matched)),
                    "best_join_key": "",
                    "is_best_join_key": "pending",
                    "pd_count": "",
                    "control_count": "",
                }
            )
    pd_count, control_count = count_labels(metadata)
    for row in rows:
        row["best_join_key"] = best_key
        row["is_best_join_key"] = "true" if row["candidate_join_key"] == best_key else "false"
        row["pd_count"] = str(pd_count)
        row["control_count"] = str(control_count)
    preview = [
        "Expression sample IDs: " + ", ".join(expr_samples[:10]),
        "Metadata geo_accession: " + ", ".join(row.get("geo_accession", "") for row in metadata[:10]),
        f"Best join key: {best_key}",
        f"Matched samples: {best_matched}/{len(expr_samples)}",
        f"PD count: {pd_count}; Control count: {control_count}",
    ]
    return rows, best_key, preview


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit GSE7621 series-matrix-to-metadata sample mapping.")
    parser.add_argument("--series-matrix", type=Path, default=Path("data/raw/external/gse7621_pd_sn_bulk/GSE7621_series_matrix.txt.gz"))
    parser.add_argument("--metadata", type=Path, default=Path("results/tables/phase37_gse7621_pd_sn_bulk_sample_metadata.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase37_gse7621_sample_mapping_audit.tsv"))
    parser.add_argument("--preview-output", type=Path, default=Path("results/reports/phase37_gse7621_sample_mapping_preview.txt"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/169_gse7621_sample_mapping_audit.log"))
    args = parser.parse_args()
    configure_logging(args.log_file)
    rows, best_key, preview = audit(args.series_matrix, args.metadata)
    columns = ["candidate_join_key", "expression_sample_count", "metadata_sample_count", "matched_sample_count", "unmatched_expression_samples", "unmatched_metadata_samples", "best_join_key", "is_best_join_key", "pd_count", "control_count"]
    write_tsv(args.output, rows, columns)
    args.preview_output.parent.mkdir(parents=True, exist_ok=True)
    args.preview_output.write_text("\n".join(preview) + "\n", encoding="utf-8")
    logging.info("Best GSE7621 join key=%s", best_key)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
