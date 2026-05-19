import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/150_convert_gse174367_bulk_rda_to_axis_matrix.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase29_bulk_converter", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_ad_endpoint_label_mapping():
    module = load_module()
    assert module.label_for_endpoint("AD") == "1"
    assert module.label_for_endpoint("Alzheimer's Disease") == "1"
    assert module.label_for_endpoint("Control") == "0"
    assert module.label_for_endpoint("unknown") == ""


def test_converter_is_axis_gene_only_and_has_sample_mapping():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    assert "axis_genes" in text
    assert "gene_symbol" in text
    assert "map_samples" in text
    assert "no expression sample ids matched" in text
    for forbidden in ["scanpy", "anndata", "read_h5ad", "toarray(", "todense(", "umap", "leiden"]:
        assert forbidden not in text


def test_coverage_rows_report_missing_axis_genes():
    module = load_module()
    rows = module.coverage_rows(
        [{"axis_id": "axis_a", "gene_members": "SNCA;MAPT;APP"}],
        {"SNCA", "APP"},
    )
    assert rows[0]["genes_found"] == "2"
    assert rows[0]["genes_missing"] == "1"
    assert rows[0]["missing_gene_members"] == "MAPT"
