#!/usr/bin/env python3
"""Summarize SEA-AD metadata-only TSV output.

This script reads only `sea_ad_obs_metadata_minimal.tsv`. It does not open H5AD,
has no h5py or Scanpy dependency, and does not access expression matrices.
"""

from __future__ import annotations

import argparse
import csv
import logging
from collections import Counter, defaultdict
from pathlib import Path


DONOR_SUMMARY_OUTPUT = "sea_ad_donor_summary.tsv"
CELLTYPE_BY_AD_OUTPUT = "sea_ad_celltype_by_ad_pathology.tsv"
CELLTYPE_BY_COG_OUTPUT = "sea_ad_celltype_by_cognitive_status.tsv"


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


def add_unique(mapping: dict[str, set[str]], key: str, value: str) -> None:
    if value:
        mapping[key].add(value)


def choose_cell_type(row: dict[str, str]) -> str:
    return row.get("Subclass") or row.get("Supertype") or row.get("Class") or "unannotated"


def summarize_metadata(metadata_path: Path) -> dict[str, object]:
    donor_cell_counts: Counter[str] = Counter()
    donor_values: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    celltype_by_ad: Counter[tuple[str, str]] = Counter()
    celltype_by_cog: Counter[tuple[str, str]] = Counter()
    donors_by_ad_celltype: dict[tuple[str, str], set[str]] = defaultdict(set)
    donors_by_cog_celltype: dict[tuple[str, str], set[str]] = defaultdict(set)
    total_rows = 0

    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total_rows += 1
            donor_id = row.get("Donor ID", "") or "unknown_donor"
            cell_type = choose_cell_type(row)
            ad_pathology = row.get("Overall AD neuropathological Change", "") or "unknown_ad_pathology"
            cognitive_status = row.get("Cognitive Status", "") or "unknown_cognitive_status"

            donor_cell_counts[donor_id] += 1
            for field in [
                "Sex",
                "Gender",
                "Age at Death",
                "Braak",
                "CERAD score",
                "Overall AD neuropathological Change",
                "Cognitive Status",
                "Brain Region",
                "APOE Genotype",
                "Highest Lewy Body Disease",
                "LATE",
                "Thal",
                "Overall CAA Score",
                "Neurotypical reference",
            ]:
                add_unique(donor_values[donor_id], field, row.get(field, ""))

            celltype_by_ad[(cell_type, ad_pathology)] += 1
            celltype_by_cog[(cell_type, cognitive_status)] += 1
            donors_by_ad_celltype[(cell_type, ad_pathology)].add(donor_id)
            donors_by_cog_celltype[(cell_type, cognitive_status)].add(donor_id)

    return {
        "total_rows": total_rows,
        "donor_cell_counts": donor_cell_counts,
        "donor_values": donor_values,
        "celltype_by_ad": celltype_by_ad,
        "celltype_by_cog": celltype_by_cog,
        "donors_by_ad_celltype": donors_by_ad_celltype,
        "donors_by_cog_celltype": donors_by_cog_celltype,
    }


def joined(values: set[str]) -> str:
    return "|".join(sorted(value for value in values if value))


def write_donor_summary(summary: dict[str, object], out_path: Path) -> None:
    donor_cell_counts: Counter[str] = summary["donor_cell_counts"]  # type: ignore[assignment]
    donor_values: dict[str, dict[str, set[str]]] = summary["donor_values"]  # type: ignore[assignment]
    fields = [
        "donor_id",
        "cell_count",
        "Sex",
        "Gender",
        "Age at Death",
        "Braak",
        "CERAD score",
        "Overall AD neuropathological Change",
        "Cognitive Status",
        "Brain Region",
        "APOE Genotype",
        "Highest Lewy Body Disease",
        "LATE",
        "Thal",
        "Overall CAA Score",
        "Neurotypical reference",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for donor_id, cell_count in sorted(donor_cell_counts.items()):
            row = {"donor_id": donor_id, "cell_count": cell_count}
            for field in fields[2:]:
                row[field] = joined(donor_values[donor_id].get(field, set()))
            writer.writerow(row)


def write_celltype_summary(
    counter: Counter[tuple[str, str]],
    donor_sets: dict[tuple[str, str], set[str]],
    group_column: str,
    out_path: Path,
) -> None:
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["cell_type", group_column, "cell_count", "donor_count"],
            delimiter="\t",
        )
        writer.writeheader()
        for (cell_type, group_value), cell_count in sorted(counter.items()):
            writer.writerow(
                {
                    "cell_type": cell_type,
                    group_column: group_value,
                    "cell_count": cell_count,
                    "donor_count": len(donor_sets[(cell_type, group_value)]),
                }
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize SEA-AD obs metadata TSV only.")
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path("data/interim/sea_ad/sea_ad_obs_metadata_minimal.tsv"),
    )
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/12_summarize_sea_ad_metadata.log"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)

    logging.info("Starting SEA-AD metadata-only summary.")
    logging.info("Metadata TSV: %s", args.metadata)
    logging.info("Tables directory: %s", args.tables_dir)

    if not args.metadata.exists():
        logging.error("Metadata TSV not found: %s", args.metadata)
        return 2

    args.tables_dir.mkdir(parents=True, exist_ok=True)
    summary = summarize_metadata(args.metadata)

    donor_summary_path = args.tables_dir / DONOR_SUMMARY_OUTPUT
    ad_path = args.tables_dir / CELLTYPE_BY_AD_OUTPUT
    cog_path = args.tables_dir / CELLTYPE_BY_COG_OUTPUT

    write_donor_summary(summary, donor_summary_path)
    write_celltype_summary(
        summary["celltype_by_ad"],  # type: ignore[arg-type]
        summary["donors_by_ad_celltype"],  # type: ignore[arg-type]
        "ad_pathology",
        ad_path,
    )
    write_celltype_summary(
        summary["celltype_by_cog"],  # type: ignore[arg-type]
        summary["donors_by_cog_celltype"],  # type: ignore[arg-type]
        "cognitive_status",
        cog_path,
    )

    logging.info("Rows summarized: %d", summary["total_rows"])
    logging.info("Donors summarized: %d", len(summary["donor_cell_counts"]))  # type: ignore[arg-type]
    logging.info("Wrote %s", donor_summary_path)
    logging.info("Wrote %s", ad_path)
    logging.info("Wrote %s", cog_path)
    logging.info("No H5AD/HDF5 files or expression matrices were opened.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
