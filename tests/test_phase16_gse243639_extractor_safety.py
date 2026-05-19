from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/78_extract_gse243639_target_gene_panel.py"


def load_extractor_module():
    spec = importlib.util.spec_from_file_location("gse243639_extractor", SCRIPT)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_extractor_does_not_import_scanpy_or_h5ad_tools() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for forbidden in ["scanpy", "anndata", "read_h5ad", "h5py", "AnnData"]:
        assert forbidden not in text


def test_extractor_does_not_densify_or_load_full_matrix() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for forbidden in ["toarray", "todense", "np.array", "pandas", "polars"]:
        assert forbidden not in text
    assert "for row in reader" in text
    assert "MAX_TARGET_GENES" in text


def test_sample_id_is_prefix_before_underscore() -> None:
    module = load_extractor_module()
    assert module.sample_id_from_cell_id("s.0096_AAACCCAAGTACGAGC.1") == "s.0096"


def test_extractor_outputs_required_columns_and_audit_fields() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for column in ["cell_id", "sample_id", "gene_symbol", "expression_value"]:
        assert column in text
    for field in ["requested_target_genes", "extracted_target_genes", "missing_target_genes", "cell_columns", "sparse_expression_rows"]:
        assert field in text
