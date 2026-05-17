import csv
from pathlib import Path


REQUIRED_COLUMNS = {
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
}

EXPECTED_DATASETS = {
    "sea_ad_single_nucleus",
    "mathys_2019_ad_single_nucleus",
    "rosmap_ad_transcriptomics",
    "pd_single_cell_single_nucleus_placeholder",
    "gut_brain_microbiome_metabolite_placeholder",
    "string_ppi_placeholder",
    "evolutionary_ortholog_placeholder",
}


def read_rows():
    path = Path("metadata/dataset_intake_checklist.tsv")
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_dataset_intake_checklist_schema():
    rows = read_rows()
    assert rows
    assert REQUIRED_COLUMNS.issubset(rows[0].keys())


def test_dataset_intake_checklist_expected_rows():
    rows = read_rows()
    dataset_ids = {row["dataset_id"] for row in rows}
    assert EXPECTED_DATASETS == dataset_ids


def test_dataset_intake_checklist_phase2_blocked_by_default():
    rows = read_rows()
    assert all(row["manual_download_needed"] == "true" for row in rows)
    assert all(row["checksum_needed"] == "true" for row in rows)
    assert all(row["ready_for_phase2"] == "false" for row in rows)
    assert all(row["blocking_issue"] for row in rows)
