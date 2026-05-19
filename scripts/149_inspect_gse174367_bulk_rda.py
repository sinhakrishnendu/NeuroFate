#!/usr/bin/env python3
"""Safely inspect the GSE174367 processed bulk RNA RDA file."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import shutil
import subprocess
import tempfile
from pathlib import Path


SUMMARY_COLUMNS = [
    "object_name",
    "object_class",
    "nrow",
    "ncol",
    "row_names_preview",
    "col_names_preview",
    "expression_matrix_like",
    "metadata_table_like",
]


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_missing_runtime_outputs(output: Path, preview_output: Path, reason: str) -> None:
    write_rows(
        output,
        [
            {
                "object_name": "unavailable",
                "object_class": "missing_runtime",
                "nrow": "",
                "ncol": "",
                "row_names_preview": "",
                "col_names_preview": "",
                "expression_matrix_like": "unknown",
                "metadata_table_like": "unknown",
            }
        ],
    )
    preview_output.parent.mkdir(parents=True, exist_ok=True)
    preview_output.write_text(
        "GSE174367 bulk RDA inspection was not executed because no supported R runtime was available.\n"
        f"Reason: {reason}\n"
        "Install Rscript or pyreadr, then rerun the inspector. No conversion or analysis was attempted.\n",
        encoding="utf-8",
    )


def run_rscript(rda_gz: Path, output: Path, preview_output: Path) -> bool:
    rscript = shutil.which("Rscript")
    helper = Path(__file__).with_name("r_inspect_gse174367_bulk_rda.R")
    if not rscript or not helper.exists():
        return False
    output.parent.mkdir(parents=True, exist_ok=True)
    preview_output.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        rscript,
        str(helper),
        "--input",
        str(rda_gz),
        "--output",
        str(output),
        "--preview-output",
        str(preview_output),
    ]
    logging.info("Inspecting RDA with Rscript helper")
    subprocess.run(cmd, check=True)
    return True


def _preview_names(values: object) -> str:
    try:
        return ";".join(str(value) for value in list(values)[:8])
    except TypeError:
        return ""


def inspect_with_pyreadr(rda_gz: Path, output: Path, preview_output: Path) -> bool:
    try:
        import pyreadr  # type: ignore
    except ImportError:
        return False

    rows: list[dict[str, str]] = []
    preview_lines = ["# GSE174367 bulk RDA pyreadr preview", ""]
    with tempfile.TemporaryDirectory() as tmpdir:
        rda_path = Path(tmpdir) / rda_gz.name.removesuffix(".gz")
        with gzip.open(rda_gz, "rb") as src, rda_path.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        result = pyreadr.read_r(str(rda_path))
        for name, obj in result.items():
            shape = getattr(obj, "shape", ("", ""))
            nrow = str(shape[0]) if len(shape) >= 1 else ""
            ncol = str(shape[1]) if len(shape) >= 2 else ""
            row_names = _preview_names(getattr(obj, "index", []))
            col_names = _preview_names(getattr(obj, "columns", []))
            rows.append(
                {
                    "object_name": name,
                    "object_class": type(obj).__name__,
                    "nrow": nrow,
                    "ncol": ncol,
                    "row_names_preview": row_names,
                    "col_names_preview": col_names,
                    "expression_matrix_like": "true" if int(nrow or 0) > 100 and int(ncol or 0) > 5 else "false",
                    "metadata_table_like": "true" if 5 < int(nrow or 0) and 1 < int(ncol or 0) < 200 else "false",
                }
            )
            preview_lines.extend([f"## {name}", f"class: {type(obj).__name__}", f"shape: {shape}", f"columns: {col_names}", ""])
    write_rows(output, rows)
    preview_output.parent.mkdir(parents=True, exist_ok=True)
    preview_output.write_text("\n".join(preview_lines) + "\n", encoding="utf-8")
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Inspect GSE174367 bulk RNA RDA structure without conversion.")
    parser.add_argument("--rda-gz", type=Path, default=Path("data/raw/external/gse174367_ad_multiomics/GSE174367_bulkRNA_processed.rda.gz"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase29_gse174367_bulk_rda_structure.tsv"))
    parser.add_argument("--preview-output", type=Path, default=Path("results/reports/phase29_gse174367_bulk_rda_preview.txt"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/149_inspect_gse174367_bulk_rda.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    if not args.rda_gz.exists():
        raise SystemExit(f"Missing input RDA file: {args.rda_gz}")
    try:
        if run_rscript(args.rda_gz, args.output, args.preview_output):
            return 0
    except (subprocess.CalledProcessError, OSError) as exc:
        logging.warning("Rscript inspection failed: %s", exc)
    if inspect_with_pyreadr(args.rda_gz, args.output, args.preview_output):
        return 0
    write_missing_runtime_outputs(args.output, args.preview_output, "Rscript unavailable or failed, and pyreadr is not installed.")
    logging.warning("Wrote missing-runtime RDA inspection report")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
