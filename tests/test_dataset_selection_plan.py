import csv
import importlib.util
from pathlib import Path


SELECTION_REQUIRED_COLUMNS = {
    "claim_id",
    "manuscript_claim",
    "disease_area",
    "required_dataset_id",
    "data_modality",
    "minimum_required_fields",
    "preferred_source",
    "access_type",
    "estimated_size_category",
    "manual_download_priority",
    "phase",
    "blocking_risk",
    "notes",
}

MVDP_REQUIRED_COLUMNS = {
    "phase",
    "dataset_id",
    "why_needed",
    "minimum_file_needed",
    "can_start_without_it",
    "first_analysis_enabled",
    "manual_action_needed",
    "notes",
}

EXPECTED_CLAIMS = {
    "claim_ad_single_cell_state_reconstruction",
    "claim_pd_single_cell_state_reconstruction",
    "claim_gut_brain_metabolite_network_layer",
    "claim_ppi_network_biology_layer",
    "claim_evolutionary_conservation_layer",
    "claim_positive_selection_layer",
    "claim_multimodal_neurofate_score",
    "claim_interpretability_reporting_layer",
}


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def registry_dataset_ids():
    return {row["dataset_id"] for row in read_tsv("metadata/dataset_registry.tsv")}


def split_dataset_ids(value):
    return [item.strip() for item in value.split(";") if item.strip()]


def test_dataset_selection_plan_schema_and_claims():
    rows = read_tsv("metadata/dataset_selection_plan.tsv")
    assert rows
    assert SELECTION_REQUIRED_COLUMNS.issubset(rows[0].keys())
    assert EXPECTED_CLAIMS == {row["claim_id"] for row in rows}


def test_dataset_selection_plan_references_registered_datasets():
    rows = read_tsv("metadata/dataset_selection_plan.tsv")
    known = registry_dataset_ids()
    for row in rows:
        for dataset_id in split_dataset_ids(row["required_dataset_id"]):
            assert dataset_id in known


def test_minimum_viable_dataset_plan_schema_and_dataset_ids():
    rows = read_tsv("metadata/minimum_viable_dataset_plan.tsv")
    assert rows
    assert MVDP_REQUIRED_COLUMNS.issubset(rows[0].keys())
    known = registry_dataset_ids()
    assert {row["dataset_id"] for row in rows}.issubset(known)
    assert any(row["can_start_without_it"] == "false" for row in rows)


def test_dataset_selection_validator_constants():
    path = Path("scripts/07_validate_dataset_selection_plan.py")
    spec = importlib.util.spec_from_file_location("validate_dataset_selection_plan", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    assert "claim_id" in module.SELECTION_COLUMNS
    assert "dataset_id" in module.MVDP_COLUMNS
    assert "claim_multimodal_neurofate_score" in module.EXPECTED_CLAIMS
