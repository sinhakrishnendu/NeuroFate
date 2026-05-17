#!/usr/bin/env python3
"""Build harmonized cross-cohort donor-level NeuroFate feature tables."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


FEATURE_PREFIXES = (
    "gene_mean__",
    "gene_detection__",
    "index__",
    "cell_fraction__",
    "celltype_index__",
)
LABEL_PREFIX = "label__"

CELLTYPE_SYNONYMS = {
    "microglia": "Microglia",
    "mic": "Microglia",
    "astrocyte": "Astrocyte",
    "ast": "Astrocyte",
    "excitatory": "Excitatory_neuron",
    "inhibitory": "Inhibitory_neuron",
    "oligodendrocyte": "Oligodendrocyte",
}

PATHOLOGY_LABEL_MAP = {
    "not ad": "none_or_low",
    "none": "none_or_low",
    "low": "none_or_low",
    "intermediate": "intermediate",
    "moderate": "intermediate",
    "high": "high",
}


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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %s rows: %d", path, len(rows))


def parse_feature_table_arg(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Feature table arguments must use cohort_id=path.tsv")
    cohort_id, path = value.split("=", 1)
    return cohort_id, Path(path)


def feature_group(feature: str) -> str:
    for prefix in FEATURE_PREFIXES:
        if feature.startswith(prefix):
            return prefix.rstrip("_")
    return "non_feature"


def harmonize_label(value: str) -> str:
    lowered = value.strip().lower()
    for key, label in PATHOLOGY_LABEL_MAP.items():
        if key in lowered:
            return label
    return value


def harmonize_row(row: dict[str, str], cohort_id: str) -> dict[str, str]:
    harmonized = dict(row)
    harmonized["cohort_id"] = cohort_id
    for column, value in list(harmonized.items()):
        if column.startswith("label__"):
            harmonized[column] = harmonize_label(value)
        if column.startswith("cell_fraction__"):
            suffix = column.removeprefix("cell_fraction__").lower()
            for key, replacement in CELLTYPE_SYNONYMS.items():
                if key in suffix:
                    harmonized[f"cell_fraction__{replacement}"] = value
    return harmonized


def build_crosscohort_tables(
    inputs: list[tuple[str, Path]],
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[str]]:
    rows: list[dict[str, str]] = []
    features_by_cohort: dict[str, set[str]] = {}
    label_columns: set[str] = set()
    all_features: set[str] = set()

    for cohort_id, path in inputs:
        cohort_rows = read_tsv(path)
        if not cohort_rows:
            logging.warning("Skipping empty donor feature table: %s", path)
            continue
        features = {
            column
            for column in cohort_rows[0]
            if column.startswith(FEATURE_PREFIXES)
        }
        features_by_cohort[cohort_id] = features
        all_features.update(features)
        label_columns.update(column for column in cohort_rows[0] if column.startswith(LABEL_PREFIX))
        for row in cohort_rows:
            rows.append(harmonize_row(row, cohort_id))

    cohorts = sorted(features_by_cohort)
    fieldnames = ["cohort_id", "donor_id", "n_cells", *sorted(label_columns), *sorted(all_features)]
    for row in rows:
        for field in fieldnames:
            row.setdefault(field, "0" if field in all_features else "missing")

    overlap_rows: list[dict[str, str]] = []
    for feature in sorted(all_features):
        present = sorted(cohort for cohort, features in features_by_cohort.items() if feature in features)
        missing = sorted(set(cohorts) - set(present))
        overlap_rows.append(
            {
                "feature": feature,
                "feature_group": feature_group(feature),
                "present_in_cohorts": ",".join(present),
                "missing_in_cohorts": ",".join(missing),
                "n_present": str(len(present)),
                "n_missing": str(len(missing)),
                "status": "shared" if len(present) == len(cohorts) else "cohort_specific_or_missing",
            }
        )
    return rows, overlap_rows, fieldnames


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build harmonized cross-cohort donor feature tables.")
    parser.add_argument(
        "--feature-table",
        dest="feature_tables",
        action="append",
        type=parse_feature_table_arg,
        default=[],
        help="Cohort donor feature table as cohort_id=path.tsv. Repeat for each cohort.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/tables/crosscohort_donor_feature_table.tsv"),
    )
    parser.add_argument(
        "--overlap-output",
        type=Path,
        default=Path("results/tables/crosscohort_feature_overlap.tsv"),
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/32_build_crosscohort_feature_tables.log"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    inputs = args.feature_tables or [
        ("sea_ad", Path("results/tables/phase5_donor_feature_table.tsv")),
    ]
    logging.info("Feature tables requested: %d", len(inputs))
    rows, overlap_rows, fieldnames = build_crosscohort_tables(inputs)
    write_tsv(args.output, rows, fieldnames)
    write_tsv(
        args.overlap_output,
        overlap_rows,
        ["feature", "feature_group", "present_in_cohorts", "missing_in_cohorts", "n_present", "n_missing", "status"],
    )
    logging.info("Cross-cohort harmonization complete without reading expression matrices.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
