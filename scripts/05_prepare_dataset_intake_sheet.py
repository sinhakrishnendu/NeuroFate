#!/usr/bin/env python3
"""LIGHTWEIGHT dataset intake checklist preparation for NeuroFate Phase 1C.

This script reads small metadata TSV files only. It never downloads datasets, opens
remote URLs, reads h5ad/HDF5 files, computes checksums, or processes biological data.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


INTAKE_COLUMNS = [
    "dataset_id",
    "manuscript_module",
    "required_for_landmark_claim",
    "minimum_required_metadata",
    "expected_raw_format",
    "expected_processed_format",
    "manual_download_needed",
    "checksum_needed",
    "ethics_or_access_note",
    "ready_for_phase2",
    "blocking_issue",
    "notes",
]

PROVENANCE_TRACKING_FIELDS = [
    "source_url_or_accession",
    "date_accessed",
    "license_or_terms",
    "checksum_value",
    "file_size_expected",
    "file_size_observed",
    "verified",
]

DATASET_INTAKE_DEFAULTS = {
    "sea_ad_single_nucleus": {
        "manuscript_module": "single_cell_transcriptomics;alzheimers_disease",
        "minimum_required_metadata": (
            "donor_id,disease_status,brain_region,cell_or_nucleus_barcode,cell_type_or_cluster"
        ),
        "expected_processed_format": "validated_metadata_tsv_and_future_processed_h5ad",
        "ethics_or_access_note": "review SEA-AD terms and citation requirements",
        "blocking_issue": "manual_download_not_completed",
        "notes": "SEA-AD AD snRNA-seq placeholder; required for AD single-nucleus manuscript layer.",
    },
    "mathys_2019_ad_single_nucleus": {
        "manuscript_module": "single_cell_transcriptomics;alzheimers_disease",
        "minimum_required_metadata": (
            "donor_id,disease_status,cell_or_nucleus_barcode,cell_type_or_cluster,study_batch"
        ),
        "expected_processed_format": "validated_metadata_tsv_and_future_processed_h5ad",
        "ethics_or_access_note": "review publication/repository license and citation requirements",
        "blocking_issue": "manual_download_not_completed",
        "notes": "Mathys 2019 AD single-cell/single-nucleus placeholder.",
    },
    "rosmap_ad_transcriptomics": {
        "manuscript_module": "alzheimers_disease",
        "minimum_required_metadata": "subject_id,disease_status,brain_region,assay_type,batch_or_cohort",
        "expected_processed_format": "validated_transcriptomics_metadata_tsv",
        "ethics_or_access_note": "controlled-access or data-use terms may apply",
        "blocking_issue": "access_terms_not_reviewed",
        "notes": "ROSMAP-associated AD transcriptomics placeholder.",
    },
    "pd_single_cell_single_nucleus_placeholder": {
        "manuscript_module": "single_cell_transcriptomics;parkinsons_disease",
        "minimum_required_metadata": (
            "donor_id,disease_status,brain_region,cell_or_nucleus_barcode,cell_type_or_cluster"
        ),
        "expected_processed_format": "validated_metadata_tsv_and_future_processed_h5ad",
        "ethics_or_access_note": "review dataset-specific consent, license, and citation requirements",
        "blocking_issue": "specific_dataset_not_selected",
        "notes": "Parkinson single-cell/single-nucleus dataset placeholder.",
    },
    "gut_brain_microbiome_metabolite_placeholder": {
        "manuscript_module": "gut_brain_axis;microbiome_metabolite_layer",
        "minimum_required_metadata": (
            "compound_or_taxon_id,source_database,evidence_type,neuroactive_annotation"
        ),
        "expected_processed_format": "validated_signature_table_tsv",
        "ethics_or_access_note": "review source database licenses and attribution requirements",
        "blocking_issue": "manual_curation_not_completed",
        "notes": "Microbiome/metabolite layer placeholder for gut-brain axis features.",
    },
    "string_ppi_placeholder": {
        "manuscript_module": "protein_interaction_network_biology",
        "minimum_required_metadata": (
            "gene_or_protein_a,gene_or_protein_b,interaction_score,source_database"
        ),
        "expected_processed_format": "validated_edge_table_tsv",
        "ethics_or_access_note": "review STRING/PPI database license and citation requirements",
        "blocking_issue": "manual_download_not_completed",
        "notes": "STRING/PPI network layer placeholder; no graph processing yet.",
    },
    "evolutionary_ortholog_placeholder": {
        "manuscript_module": "evolutionary_conservation;positive_selection",
        "minimum_required_metadata": (
            "gene_id,ortholog_id,species,conservation_score_or_selection_metric"
        ),
        "expected_processed_format": "validated_evolutionary_feature_table_tsv",
        "ethics_or_access_note": "review ortholog/conservation database license and citation requirements",
        "blocking_issue": "source_database_not_selected",
        "notes": "Evolutionary ortholog/conservation layer placeholder.",
    },
}

MISSING_PROVENANCE_VALUES = {"", "pending", "not_accessed", "manual_accession_required"}


def configure_logging(log_file: Path | None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file, mode="w", encoding="utf-8"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
    )


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        return [], []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def generate_intake_rows(dataset_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    generated: list[dict[str, str]] = []
    for dataset in dataset_rows:
        dataset_id = dataset["dataset_id"]
        defaults = DATASET_INTAKE_DEFAULTS.get(dataset_id, {})
        generated.append(
            {
                "dataset_id": dataset_id,
                "manuscript_module": defaults.get("manuscript_module", "manual_mapping_required"),
                "required_for_landmark_claim": "true",
                "minimum_required_metadata": defaults.get(
                    "minimum_required_metadata",
                    "manual_metadata_requirements_required",
                ),
                "expected_raw_format": dataset.get("expected_file_type", "manual_format_required"),
                "expected_processed_format": defaults.get(
                    "expected_processed_format",
                    "validated_metadata_tsv",
                ),
                "manual_download_needed": dataset.get("heavy_to_download", "true"),
                "checksum_needed": "true",
                "ethics_or_access_note": defaults.get(
                    "ethics_or_access_note",
                    "manual license and access review required",
                ),
                "ready_for_phase2": "false",
                "blocking_issue": defaults.get("blocking_issue", "manual_intake_not_completed"),
                "notes": defaults.get("notes", dataset.get("notes", "manual intake required")),
            }
        )
    return generated


def compare_existing_rows(
    expected_rows: list[dict[str, str]],
    existing_rows: list[dict[str, str]],
) -> list[str]:
    issues: list[str] = []
    expected_ids = {row["dataset_id"] for row in expected_rows}
    existing_ids = {row.get("dataset_id", "") for row in existing_rows}
    missing = sorted(expected_ids - existing_ids)
    extra = sorted(existing_ids - expected_ids)
    if missing:
        issues.append(f"dataset_intake_checklist missing dataset_id rows: {', '.join(missing)}")
    if extra:
        issues.append(f"dataset_intake_checklist has extra dataset_id rows: {', '.join(extra)}")
    return issues


def summarize_missing_provenance(provenance_rows: list[dict[str, str]]) -> dict[str, list[str]]:
    missing_by_dataset: dict[str, list[str]] = {}
    for row in provenance_rows:
        dataset_id = row.get("dataset_id", "")
        missing = [
            field
            for field in PROVENANCE_TRACKING_FIELDS
            if row.get(field, "").strip() in MISSING_PROVENANCE_VALUES
        ]
        if missing:
            missing_by_dataset[dataset_id] = missing
    return missing_by_dataset


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LIGHTWEIGHT NeuroFate dataset intake prep.")
    parser.add_argument("--dataset-registry", type=Path, default=Path("metadata/dataset_registry.tsv"))
    parser.add_argument("--module-map", type=Path, default=Path("metadata/manuscript_module_map.tsv"))
    parser.add_argument(
        "--provenance-template",
        type=Path,
        default=Path("metadata/provenance_template.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("metadata/dataset_intake_checklist.tsv"),
    )
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/05_prepare_dataset_intake_sheet.log"))
    parser.add_argument("--write", action="store_true", help="Write the generated intake checklist.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview/validate only. This script never processes datasets in any mode.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    effective_dry_run = args.dry_run or not args.write

    logging.info("Starting LIGHTWEIGHT Phase 1C dataset intake preparation.")
    logging.info("Dry run: %s", effective_dry_run)
    logging.info("Dataset registry: %s", args.dataset_registry)
    logging.info("Manuscript module map: %s", args.module_map)
    logging.info("Provenance template: %s", args.provenance_template)
    logging.info("Checklist output: %s", args.output)

    _, dataset_rows = read_tsv(args.dataset_registry)
    module_columns, module_rows = read_tsv(args.module_map)
    intake_columns, intake_rows = read_tsv(args.output)
    _, provenance_rows = read_tsv(args.provenance_template)

    if not dataset_rows:
        logging.error("No dataset registry rows found.")
        return 2
    if not module_rows:
        logging.error("No manuscript module map rows found.")
        return 2

    expected_rows = generate_intake_rows(dataset_rows)
    issues: list[str] = []
    if intake_rows:
        missing_columns = [column for column in INTAKE_COLUMNS if column not in intake_columns]
        if missing_columns:
            issues.append(
                "dataset_intake_checklist missing required columns: "
                + ", ".join(missing_columns)
            )
        issues.extend(compare_existing_rows(expected_rows, intake_rows))
    else:
        logging.info("No existing checklist found; generated rows are ready for writing.")

    module_ids = {row["module_id"] for row in module_rows}
    for row in expected_rows:
        modules = row["manuscript_module"].split(";")
        unknown = [module for module in modules if module not in module_ids]
        if unknown:
            issues.append(
                f"{row['dataset_id']} references unknown manuscript modules: {', '.join(unknown)}"
            )

    missing_provenance = summarize_missing_provenance(provenance_rows)

    logging.info("Expected dataset intake rows: %d", len(expected_rows))
    for row in expected_rows:
        logging.info(
            "Expected dataset: %s | modules=%s | phase2_ready=%s",
            row["dataset_id"],
            row["manuscript_module"],
            row["ready_for_phase2"],
        )

    logging.info("Missing or pending provenance fields by dataset:")
    if missing_provenance:
        for dataset_id, fields in missing_provenance.items():
            logging.info("  %s: %s", dataset_id, ", ".join(fields))
    else:
        logging.info("  none")

    if issues:
        for issue in issues:
            logging.error(issue)
        logging.error("Dataset intake checklist validation failed.")
        logging.info("No datasets were opened, downloaded, checksummed, or processed.")
        return 2

    if effective_dry_run:
        logging.info("Dry run complete. No files were written.")
    else:
        write_tsv(args.output, expected_rows, INTAKE_COLUMNS)
        logging.info("Wrote dataset intake checklist: %s", args.output)

    logging.info("No datasets were opened, downloaded, checksummed, or processed.")
    logging.info("No HDF5/h5ad files were accessed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
