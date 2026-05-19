from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/106_build_neurofate_axis_scores.py"


def test_axis_score_script_excludes_label_columns():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "LABEL_PREFIXES" in text
    assert "METADATA_COLUMNS" in text
    assert "is_label_or_metadata_column" in text
    assert "feature_columns" in text
    assert "label__" in text
    assert "diagnosis" in text
    assert "apoe_genotype" in text


def test_axis_score_script_uses_donor_sample_tables_only():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    assert "phase5_donor_feature_table.tsv" in text
    assert "phase20_gse243639_celltype_feature_table.tsv" in text
    assert "scanpy" not in text
    assert "read_h5ad" not in text
    assert ".h5ad" not in text
    assert "umap" not in text
    assert "leiden" not in text
