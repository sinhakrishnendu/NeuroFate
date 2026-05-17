import csv
from pathlib import Path


REQUIRED_COLUMNS = {
    "dataset_id",
    "disease_area",
    "modality",
    "source",
    "expected_file_type",
    "local_path",
    "status",
    "heavy_to_download",
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
    path = Path("metadata/dataset_registry.tsv")
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_dataset_registry_schema():
    rows = read_rows()
    assert rows
    assert REQUIRED_COLUMNS.issubset(rows[0].keys())


def test_dataset_registry_expected_placeholders():
    rows = read_rows()
    dataset_ids = {row["dataset_id"] for row in rows}
    assert EXPECTED_DATASETS.issubset(dataset_ids)
    assert all(row["status"] == "placeholder" for row in rows)


def test_dataset_registry_heavy_download_flags_are_explicit():
    rows = read_rows()
    assert all(row["heavy_to_download"] in {"true", "false"} for row in rows)
    assert any(row["heavy_to_download"] == "true" for row in rows)
