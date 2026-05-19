#!/usr/bin/env python3
"""Inspect external metadata fields safely without reading expression matrices."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from pathlib import Path
from typing import TextIO


EXACT_CANONICAL_FIELDS = {
    "sample id": "sample_id",
    "clinical diagnosis": "diagnosis",
    "age": "age",
    "sex": "sex",
    "pmi hours": "pmi",
    "rin measure": "rin",
    "lewy bodies presence in midbrain": "lewy_body_midbrain",
    "lewy bodies presence in limbic regions (amygdala)": "lewy_body_limbic",
    "lewy bodies presence in neocortical regions (frontal cortex)": "lewy_body_neocortical",
    "cerad score for neuritic plaques": "cerad",
    "braak stage for neurofibrillary tangles": "braak",
}

CANONICAL_TERMS = {
    "donor_id": ["donor", "individual", "subject"],
    "sample_id": ["sample"],
    "cell_id": ["cell", "barcode"],
    "diagnosis": ["diagnosis", "condition", "disease"],
    "disease_status": ["status", "case_control"],
    "pathology": ["braak", "cerad", "pathology", "neuropath"],
    "brain_region": ["region", "brain"],
    "age": ["age"],
    "sex": ["sex", "gender"],
    "pmi": ["pmi", "postmortem"],
    "rin": ["rin"],
    "lewy_body_midbrain": ["lewy", "midbrain"],
    "lewy_body_limbic": ["lewy", "limbic", "amygdala"],
    "lewy_body_neocortical": ["lewy", "neocortical", "frontal"],
    "cerad": ["cerad"],
    "braak": ["braak"],
    "cell_type": ["celltype", "cell_type", "subclass", "cluster"],
    "batch": ["batch", "lane", "library"],
    "apoe_genotype": ["apoe"],
    "sequencing_platform": ["platform", "chemistry", "sequenc"],
}
DELIMITER_CHOICES = [",", "\t", ";"]


def setup_logging(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=path, level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")


def open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def normalize_delimiter(value: str | None) -> str | None:
    if value is None:
        return None
    if value in {"\\t", "tab", "TAB"}:
        return "\t"
    return value


def sniff_delimiter_from_lines(lines: list[str]) -> str:
    best_delimiter = ","
    best_score = -1
    for delimiter in DELIMITER_CHOICES:
        score = 0
        for line in lines:
            if line.strip():
                score = max(score, line.count(delimiter))
        if score > best_score:
            best_delimiter = delimiter
            best_score = score
    return best_delimiter


def parse_line(line: str, delimiter: str) -> list[str]:
    return next(csv.reader([line], delimiter=delimiter))


def detect_header_index(lines: list[str], delimiter: str) -> int:
    best_index = 0
    best_score = -1
    terms = set(EXACT_CANONICAL_FIELDS)
    for index, line in enumerate(lines):
        if not line.strip():
            continue
        fields = [field.strip() for field in parse_line(line, delimiter)]
        if len(fields) < 2:
            continue
        lowered = {field.lower() for field in fields}
        exact_hits = len(lowered & terms)
        keyword_hits = sum(
            1
            for field in lowered
            if any(term in field for values in CANONICAL_TERMS.values() for term in values)
        )
        score = exact_hits * 10 + keyword_hits + len(fields)
        if score > best_score:
            best_index = index
            best_score = score
    return best_index


def inspect_text_header(
    path: Path,
    delimiter: str | None,
    header_line: int | None = None,
    max_preview_rows: int = 5,
    max_scan_lines: int = 200,
) -> tuple[list[str], list[list[str]], str, int]:
    with open_text(path) as handle:
        lines = []
        for index, line in enumerate(handle):
            if index >= max_scan_lines:
                break
            lines.append(line.rstrip("\n\r"))
    if not lines:
        return [], [], delimiter or ",", 0
    chosen_delimiter = delimiter or sniff_delimiter_from_lines(lines)
    header_index = header_line - 1 if header_line is not None else detect_header_index(lines, chosen_delimiter)
    if header_index < 0 or header_index >= len(lines):
        raise RuntimeError(f"Header line {header_line} is outside the inspected range for {path}.")
    header = [field.strip() for field in parse_line(lines[header_index], chosen_delimiter)]
    preview = []
    for line in lines[header_index + 1 :]:
        if not line.strip():
            continue
        row = parse_line(line, chosen_delimiter)
        if row:
            preview.append(row)
        if len(preview) >= max_preview_rows:
            break
    return header, preview, chosen_delimiter, header_index


def inspect_hdf5_metadata(path: Path) -> list[str]:
    import h5py

    fields: list[str] = []
    with h5py.File(path, "r") as handle:
        for key in handle.keys():
            if key == "X":
                fields.append("FORBIDDEN_X_PRESENT_NOT_ACCESSED")
                continue
            if key in {"obs", "var"}:
                node = handle[key]
                if hasattr(node, "keys"):
                    fields.extend(f"{key}/{child}" for child in node.keys() if child != "X")
            else:
                fields.append(key)
    return fields


def suggest_mapping(field: str) -> tuple[str, str]:
    lowered = field.lower().replace(" ", "_")
    exact = EXACT_CANONICAL_FIELDS.get(field.strip().lower())
    if exact is not None:
        return exact, "exact_gse243639_or_canonical_match"
    for canonical, terms in CANONICAL_TERMS.items():
        if any(term in lowered for term in terms):
            return canonical, "name_similarity"
    return "unmapped", "no_simple_match"


def write_field_outputs(dataset_id: str, fields: list[str], output: Path, mapping_output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["dataset_id", "field_name", "field_role"], delimiter="\t")
        writer.writeheader()
        for field in fields:
            role, _reason = suggest_mapping(field)
            writer.writerow({"dataset_id": dataset_id, "field_name": field, "field_role": role})
    with mapping_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["dataset_id", "source_field", "canonical_field", "confidence", "reason"],
            delimiter="\t",
        )
        writer.writeheader()
        for field in fields:
            canonical, reason = suggest_mapping(field)
            confidence = "medium" if canonical != "unmapped" else "low"
            writer.writerow(
                {
                    "dataset_id": dataset_id,
                    "source_field": field,
                    "canonical_field": canonical,
                    "confidence": confidence,
                    "reason": reason,
                }
            )


def write_preview_output(path: Path, header: list[str], preview: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["preview_row_index", *header]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for index, row in enumerate(preview, start=1):
            values = {field: value for field, value in zip(header, row, strict=False)}
            values["preview_row_index"] = str(index)
            writer.writerow(values)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Safely inspect external metadata fields.")
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--metadata-file", type=Path, required=True)
    parser.add_argument("--format", choices=["auto", "csv", "tsv", "h5ad", "h5"], default="auto")
    parser.add_argument("--header-line", type=int, default=None, help="Optional 1-based header line for files with prose preambles.")
    parser.add_argument("--delimiter", default=None, help="Optional delimiter override: ',', ';', '\\t', or tab.")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mapping-output", type=Path, default=None)
    parser.add_argument("--preview-output", type=Path, default=None)
    parser.add_argument("--log-file", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log_file)
    file_format = args.format
    if file_format == "auto":
        suffix = args.metadata_file.name.lower()
        file_format = "h5ad" if suffix.endswith(".h5ad") else "h5" if suffix.endswith(".h5") else "tsv" if ".tsv" in suffix else "csv"
    if file_format in {"csv", "tsv"}:
        delimiter = normalize_delimiter(args.delimiter) or ("\t" if file_format == "tsv" else None)
        fields, preview, detected_delimiter, header_index = inspect_text_header(args.metadata_file, delimiter, args.header_line)
        preview_output = args.preview_output or Path(f"results/reports/phase15_{args.dataset_id}_metadata_preview.tsv")
        write_preview_output(preview_output, fields, preview)
        logging.info(
            "Text metadata delimiter=%r header_line=%s preview_rows=%s preview_output=%s",
            detected_delimiter,
            header_index + 1,
            len(preview),
            preview_output,
        )
    else:
        fields = inspect_hdf5_metadata(args.metadata_file)
    mapping_output = args.mapping_output or Path(f"results/reports/phase15_{args.dataset_id}_canonical_mapping_suggestions.tsv")
    write_field_outputs(args.dataset_id, fields, args.output, mapping_output)
    logging.info("Inspected %s metadata fields: %s", args.dataset_id, len(fields))
    print(f"Wrote {args.output}")
    print(f"Wrote {mapping_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
