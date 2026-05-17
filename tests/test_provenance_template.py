import csv
from pathlib import Path


REQUIRED_COLUMNS = {
    "provenance_id",
    "dataset_id",
    "source_name",
    "source_url_or_accession",
    "source_database",
    "download_status",
    "download_command_manual_only",
    "date_accessed",
    "license_or_terms",
    "original_file_name",
    "local_expected_path",
    "checksum_algorithm",
    "checksum_value",
    "file_size_expected",
    "file_size_observed",
    "verified",
    "notes",
}


def read_rows():
    path = Path("metadata/provenance_template.tsv")
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_provenance_template_schema():
    rows = read_rows()
    assert rows
    assert REQUIRED_COLUMNS.issubset(rows[0].keys())


def test_provenance_template_covers_dataset_registry():
    rows = read_rows()
    provenance_dataset_ids = {row["dataset_id"] for row in rows}
    with Path("metadata/dataset_registry.tsv").open("r", encoding="utf-8", newline="") as handle:
        dataset_rows = list(csv.DictReader(handle, delimiter="\t"))
    registry_dataset_ids = {row["dataset_id"] for row in dataset_rows}
    assert registry_dataset_ids == provenance_dataset_ids


def test_provenance_template_is_manual_and_unverified_by_default():
    rows = read_rows()
    assert all(row["download_status"] == "not_started" for row in rows)
    assert all(row["download_command_manual_only"].startswith("MANUAL_ONLY:") for row in rows)
    assert all(row["checksum_algorithm"] == "sha256" for row in rows)
    assert all(row["checksum_value"] == "pending" for row in rows)
    assert all(row["verified"] == "false" for row in rows)
