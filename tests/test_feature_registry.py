import csv
from pathlib import Path


REQUIRED_COLUMNS = {
    "feature_id",
    "feature_group",
    "biological_layer",
    "expected_input",
    "expected_output",
    "status",
    "notes",
}

EXPECTED_FEATURES = {
    "single_cell_expression_features",
    "cell_type_labels",
    "disease_state_labels",
    "pseudotime_trajectory_features",
    "microbiome_metabolite_features",
    "protein_interaction_features",
    "evolutionary_conservation_features",
    "positive_selection_features",
    "multimodal_neurofate_score",
}


def read_rows():
    path = Path("metadata/feature_registry.tsv")
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_feature_registry_schema():
    rows = read_rows()
    assert rows
    assert REQUIRED_COLUMNS.issubset(rows[0].keys())


def test_feature_registry_expected_features():
    rows = read_rows()
    feature_ids = {row["feature_id"] for row in rows}
    assert EXPECTED_FEATURES.issubset(feature_ids)
    assert all(row["status"] == "placeholder" for row in rows)


def test_feature_registry_has_unique_ids():
    rows = read_rows()
    feature_ids = [row["feature_id"] for row in rows]
    assert len(feature_ids) == len(set(feature_ids))
