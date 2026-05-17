import csv
from pathlib import Path


REQUIRED_COLUMNS = {
    "module_id",
    "manuscript_section",
    "scientific_layer",
    "expected_data_source",
    "expected_pipeline_script",
    "expected_result_table",
    "expected_result_figure",
    "status",
    "notes",
}

EXPECTED_MODULES = {
    "single_cell_transcriptomics",
    "alzheimers_disease",
    "parkinsons_disease",
    "gut_brain_axis",
    "microbiome_metabolite_layer",
    "protein_interaction_network_biology",
    "evolutionary_conservation",
    "positive_selection",
    "multimodal_ai_fate_prediction",
    "interpretability_reporting",
}


def read_rows():
    path = Path("metadata/manuscript_module_map.tsv")
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_manuscript_module_map_schema():
    rows = read_rows()
    assert rows
    assert REQUIRED_COLUMNS.issubset(rows[0].keys())


def test_manuscript_module_map_expected_modules():
    rows = read_rows()
    module_ids = {row["module_id"] for row in rows}
    assert EXPECTED_MODULES.issubset(module_ids)
    assert all(row["status"] == "placeholder" for row in rows)


def test_manuscript_module_map_pipeline_placeholders_exist():
    rows = read_rows()
    for row in rows:
        assert Path(row["expected_pipeline_script"]).exists()
