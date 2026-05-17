import csv
from pathlib import Path


REGISTRY = Path("metadata/external_validation_registry.tsv")
PLAN_SCRIPT = Path("scripts/30_prepare_external_validation_plan.py")
HARMONIZATION_SCRIPT = Path("scripts/32_build_crosscohort_feature_tables.py")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_external_validation_registry_columns_and_cohorts():
    rows = read_rows(REGISTRY)
    expected_columns = [
        "dataset_id",
        "disease",
        "modality",
        "species",
        "brain_region",
        "cells_or_donors",
        "accession",
        "download_status",
        "processing_status",
        "notes",
    ]
    assert rows
    assert list(rows[0].keys()) == expected_columns
    dataset_ids = {row["dataset_id"] for row in rows}
    for expected in [
        "mathys_2019_ad",
        "rosmap_ad_transcriptomics",
        "pd_sn_cortical_placeholder",
        "tabula_sapiens_brain_optional",
        "adkp_harmonized_optional",
    ]:
        assert expected in dataset_ids


def test_external_registry_marks_all_cohorts_not_downloaded():
    rows = read_rows(REGISTRY)
    assert all(row["download_status"] == "not_downloaded" for row in rows)
    assert all(row["processing_status"] in {"planned", "optional"} for row in rows)


def test_external_plan_mentions_metadata_and_gene_overlap():
    text = PLAN_SCRIPT.read_text(encoding="utf-8")
    assert "REQUIRED_METADATA_FIELDS" in text
    assert "build_gene_overlap_rows" in text
    assert "build_metadata_overlap_rows" in text
    assert "external_validation_gene_overlap.tsv" in text
    assert "external_validation_metadata_overlap.tsv" in text


def test_harmonization_strategy_declared():
    text = HARMONIZATION_SCRIPT.read_text(encoding="utf-8")
    for expected in [
        "CELLTYPE_SYNONYMS",
        "PATHOLOGY_LABEL_MAP",
        "cohort_id",
        "crosscohort_donor_feature_table.tsv",
        "crosscohort_feature_overlap.tsv",
        "cohort_specific_or_missing",
    ]:
        assert expected in text
