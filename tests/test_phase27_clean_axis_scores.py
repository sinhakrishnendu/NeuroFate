import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCORER = ROOT / "scripts/128_build_gse184950_axis_scores.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase27_scores", SCORER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_clean_axis_scores_exclude_non_metadata_and_missing_label_rows():
    module = load_module()
    metadata = [
        {"sample_name": "S01", "disease_state": "Unaffected Control", "label__pd_pdd_vs_control": "0"},
        {"sample_name": "S02", "disease_state": "Parkinson's Disease", "label__pd_pdd_vs_control": "1"},
        {"sample_name": "S03", "disease_state": "Parkinson's Disease Dementia", "label__pd_pdd_vs_control": "1"},
    ]
    expression = [
        {"sample_id": "S01", "gene_symbol": "SNCA", "mean_expression": "1"},
        {"sample_id": "S02", "gene_symbol": "SNCA", "mean_expression": "2"},
        {"sample_id": "S03", "gene_symbol": "SNCA", "mean_expression": "3"},
        {"sample_id": "processed_matrices", "gene_symbol": "SNCA", "mean_expression": "100"},
    ]
    rows, _coverage, labels = module.build_scores(expression, metadata, {"synuclein_axis": ["SNCA"]})
    assert [row["sample_id"] for row in rows] == ["S01", "S02", "S03"]
    assert "processed_matrices" not in {row["sample_id"] for row in rows}
    assert {row["label__pd_pdd_vs_control"] for row in rows} == {"0", "1"}
    assert labels == [{"label__pd_pdd_vs_control": "0", "count": "1"}, {"label__pd_pdd_vs_control": "1", "count": "2"}]


def test_clean_axis_score_defaults_are_phase27_outputs():
    text = SCORER.read_text(encoding="utf-8")
    assert "phase27_gse184950_axis_scores_clean.tsv" in text
    assert "phase27_gse184950_axis_label_summary_clean.tsv" in text
