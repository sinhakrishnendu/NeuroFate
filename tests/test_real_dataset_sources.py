import csv
from pathlib import Path


REQUIRED_COLUMNS = {
    "dataset_id",
    "real_source_name",
    "official_resource",
    "accession_or_resource_id",
    "access_type",
    "manual_download_required",
    "controlled_access",
    "recommended_priority",
    "expected_format",
    "expected_size_category",
    "notes",
}

EXPECTED_REAL_SOURCES = {
    "sea_ad_single_nucleus",
    "mathys_2019_ad_single_nucleus",
    "rosmap_ad_transcriptomics",
    "gut_brain_microbiome_metabolite_placeholder",
    "string_ppi_placeholder",
}

CONTROLLED_DATASETS = {
    "mathys_2019_ad_single_nucleus",
    "rosmap_ad_transcriptomics",
}


def read_rows(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_real_dataset_sources_schema_and_priority_sources():
    rows = read_rows("metadata/real_dataset_sources.tsv")
    assert rows
    assert REQUIRED_COLUMNS.issubset(rows[0].keys())
    dataset_ids = {row["dataset_id"] for row in rows}
    assert EXPECTED_REAL_SOURCES.issubset(dataset_ids)


def test_real_dataset_sources_reference_registry_ids():
    rows = read_rows("metadata/real_dataset_sources.tsv")
    registry_ids = {row["dataset_id"] for row in read_rows("metadata/dataset_registry.tsv")}
    assert {row["dataset_id"] for row in rows}.issubset(registry_ids)


def test_controlled_access_flags_are_correct():
    rows = read_rows("metadata/real_dataset_sources.tsv")
    by_id = {row["dataset_id"]: row for row in rows}
    for dataset_id in CONTROLLED_DATASETS:
        assert by_id[dataset_id]["controlled_access"] == "true"
    assert by_id["sea_ad_single_nucleus"]["controlled_access"] == "false"
    assert by_id["string_ppi_placeholder"]["controlled_access"] == "false"


def test_source_priorities_include_sea_ad_first():
    rows = read_rows("metadata/real_dataset_sources.tsv")
    by_id = {row["dataset_id"]: row for row in rows}
    assert by_id["sea_ad_single_nucleus"]["recommended_priority"] == "P0"
    assert by_id["sea_ad_single_nucleus"]["manual_download_required"] == "true"
