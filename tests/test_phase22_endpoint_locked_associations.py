from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/112_test_axis_associations_endpoint_locked.py"


def test_endpoint_locked_script_uses_registry_not_arbitrary_labels():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "--endpoint-registry" in text
    assert "metadata/neurofate_axis_endpoint_registry.tsv" in text
    assert "source_column" in text
    assert "label_columns" not in text
    assert "LABEL_HINTS" not in text
    assert "best_axis_effect" not in text


def test_binary_direction_is_explicit():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "positive_class" in text
    assert "negative_class" in text
    assert "positive_minus_negative" in text
    assert "positive_n" in text
    assert "negative_n" in text
    assert "rank_biserial" in text
    assert "standardized_mean_difference" in text


def test_endpoint_locked_script_avoids_forbidden_single_cell_tools():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["scanpy", "read_h5ad", ".h5ad", "umap", "leiden", "cluster", "torch"]:
        assert forbidden not in text
