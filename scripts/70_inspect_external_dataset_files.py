#!/usr/bin/env python3
"""Safely inventory external dataset files without opening biological matrices."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


FORMAT_SUFFIXES = {
    ".h5ad": "h5ad",
    ".h5": "h5",
    ".mtx": "mtx",
    ".gz": "compressed",
    ".csv": "csv",
    ".tsv": "tsv",
    ".txt": "txt",
    ".rds": "rds",
    ".loom": "loom",
}
LARGE_FILE_MB = 500


def setup_logging(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=path,
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def detect_format(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".mtx.gz"):
        return "mtx.gz"
    if name.endswith(".csv.gz"):
        return "csv.gz"
    if name.endswith(".tsv.gz"):
        return "tsv.gz"
    if name.endswith(".txt.gz"):
        return "txt.gz"
    return FORMAT_SUFFIXES.get(path.suffix.lower(), "unknown")


def likely_role(path: Path) -> str:
    name = path.name.lower()
    if any(token in name for token in ["meta", "covar", "clinical", "phenotype", "sample", "donor"]):
        return "metadata_or_covariates"
    if any(token in name for token in ["feature", "gene", "var"]):
        return "gene_or_feature_file"
    if any(token in name for token in ["count", "matrix", "expr", "expression", "mtx"]):
        return "count_matrix_candidate"
    if path.suffix.lower() in {".h5ad", ".h5", ".loom"}:
        return "container_candidate"
    return "unknown"


def inventory(dataset_id: str, input_dir: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if not input_dir.exists():
        return rows
    for path in sorted(input_dir.rglob("*")):
        if not path.is_file():
            continue
        size = path.stat().st_size
        rows.append(
            {
                "dataset_id": dataset_id,
                "path": str(path),
                "file_name": path.name,
                "size_bytes": str(size),
                "size_mb": f"{size / (1024 * 1024):.3f}",
                "detected_format": detect_format(path),
                "likely_role": likely_role(path),
                "large_file_flag": str(size >= LARGE_FILE_MB * 1024 * 1024).lower(),
            }
        )
    return rows


def recommendation_rows(dataset_id: str, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    formats = sorted({row["detected_format"] for row in rows})
    roles = sorted({row["likely_role"] for row in rows})
    if not rows:
        recommendation = "input_directory_missing_or_empty"
    elif "h5ad" in formats:
        recommendation = "use_h5ad_metadata_inspection_then_h5ad_sparse_plan"
    elif "mtx.gz" in formats or "mtx" in formats:
        recommendation = "use_mtx_features_barcodes_sparse_plan"
    elif "csv.gz" in formats or "csv" in formats or "tsv.gz" in formats or "tsv" in formats:
        recommendation = "inspect_headers_to_determine_matrix_orientation"
    elif "rds" in formats:
        recommendation = "manual_rds_conversion_required_outside_codex"
    else:
        recommendation = "format_unknown_manual_review_required"
    return [
        {
            "dataset_id": dataset_id,
            "detected_formats": ",".join(formats) if formats else "none",
            "detected_roles": ",".join(roles) if roles else "none",
            "format_recommendation": recommendation,
            "note": "No matrix contents were opened; this is filename and size inventory only.",
        }
    ]


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Inventory external dataset files safely.")
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-summary", type=Path, required=True)
    parser.add_argument("--format-output", type=Path, default=None)
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/70_inspect_external_dataset_files.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log_file)
    rows = inventory(args.dataset_id, args.input_dir)
    logging.info("Dataset %s inventory rows: %s", args.dataset_id, len(rows))
    write_tsv(
        args.output_summary,
        rows,
        ["dataset_id", "path", "file_name", "size_bytes", "size_mb", "detected_format", "likely_role", "large_file_flag"],
    )
    format_output = args.format_output or Path(f"results/reports/phase15_{args.dataset_id}_format_recommendation.tsv")
    write_tsv(
        format_output,
        recommendation_rows(args.dataset_id, rows),
        ["dataset_id", "detected_formats", "detected_roles", "format_recommendation", "note"],
    )
    print(f"Wrote {args.output_summary}")
    print(f"Wrote {format_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
