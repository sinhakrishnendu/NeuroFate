import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/145_build_ad_replication_axis_scores_from_matrix.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase28_axis_builder", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_axis_builder_excludes_label_columns_from_gene_features():
    module = load_module()
    sample_values = {
        "S1": {"SNCA": 1.0, "MAPT": 2.0},
        "S2": {"SNCA": 2.0, "MAPT": 4.0},
        "bad": {"SNCA": 10.0},
    }
    metadata = {
        "S1": {"sample_id": "S1", "diagnosis": "Control", "label__ad_replication_binary": "0"},
        "S2": {"sample_id": "S2", "diagnosis": "AD", "label__ad_replication_binary": "1"},
    }
    axes = [{"axis_id": "amyloid_tau_axis", "gene_members": "SNCA;MAPT"}]
    rows, coverage, labels = module.build_scores(sample_values, metadata, axes, "diagnosis")
    assert [row["sample_id"] for row in rows] == ["S1", "S2"]
    assert "bad" not in {row["sample_id"] for row in rows}
    assert all("label__ad_replication_binary" in row for row in rows)
    assert coverage[0]["genes_found"] == "2"
    assert labels == [{"label__ad_replication_binary": "0", "count": "1"}, {"label__ad_replication_binary": "1", "count": "1"}]


def test_axis_builder_source_avoids_dense_and_single_cell_tools():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["scanpy", "anndata", "read_h5ad", "toarray(", "todense(", "umap", "leiden"]:
        assert forbidden not in text
    assert "is_metadata_column" in text
