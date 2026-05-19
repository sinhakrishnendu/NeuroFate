import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/127_extract_gse184950_axis_genes_from_10x.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase26_axis_extract", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_axis_extractor_aggregates_to_sample_level_from_sparse_matrix(tmp_path):
    module = load_module()
    sample_dir = tmp_path / "processed_matrices" / "A10"
    sample_dir.mkdir(parents=True)
    (sample_dir / "features.tsv").write_text("ENSG1\tSNCA\nENSG2\tMAPT\nENSG3\tGFAP\n", encoding="utf-8")
    (sample_dir / "barcodes.tsv").write_text("cell1\ncell2\ncell3\n", encoding="utf-8")
    (sample_dir / "matrix.mtx").write_text(
        "%%MatrixMarket matrix coordinate integer general\n"
        "3 3 4\n"
        "1 1 2\n"
        "1 3 1\n"
        "2 2 5\n"
        "3 1 1\n",
        encoding="utf-8",
    )
    rows, audit = module.extract_sample(sample_dir, {"SNCA", "GFAP"}, {"A10"})
    assert {row["gene_symbol"] for row in rows} == {"SNCA", "GFAP"}
    assert all(row["sample_id"] == "A10" for row in rows)
    assert all(set(row) >= {"sample_id", "gene_symbol", "mean_expression", "detection_rate", "n_cells"} for row in rows)
    assert audit["genes_found"] == "2"
    assert audit["n_cells"] == "3"


def test_axis_extractor_does_not_densify_or_write_cell_level_tables():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["toarray(", "todense(", "np.array(", "scanpy", "read_h5ad", "anndata"]:
        assert forbidden not in text
    assert "sample-level" in text
    assert "cell_id" not in text
    assert "mean_expression" in text
    assert "detection_rate" in text
