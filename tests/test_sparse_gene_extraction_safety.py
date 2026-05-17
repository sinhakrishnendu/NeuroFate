from pathlib import Path


PLAN_SCRIPT = Path("scripts/14_plan_sparse_gene_extraction.py")
EXTRACT_SCRIPT = Path("scripts/15_sparse_gene_extraction_safe.py")


def test_planning_script_does_not_open_h5ad_or_expression_matrix():
    text = PLAN_SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "import h5py" not in lowered
    assert "read_h5ad" not in lowered
    assert "import scanpy" not in lowered
    assert "handle[\"X\"]" not in text
    assert "handle['X']" not in text


def test_sparse_extractor_has_hard_limits_and_manual_execute_gate():
    text = EXTRACT_SCRIPT.read_text(encoding="utf-8")
    assert "MAX_GENES_DEFAULT = 64" in text
    assert "CHUNK_SIZE_DEFAULT = 5000" in text
    assert "MAX_CHUNK_SIZE = 50000" in text
    assert "M5_MAX_HIGH_MEMORY_LIMIT_MB = 32768" in text
    assert "--m5-max-profile" in text
    assert "RUN_MANUAL_EXTRACTION" in text
    assert "--execute" in text
    assert "--dry-run" in text


def test_sparse_extractor_uses_csr_components_without_dense_conversion():
    text = EXTRACT_SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered
    assert "toarray(" not in lowered
    assert "todense(" not in lowered
    assert "np.array" not in lowered
    assert "matrix[\"indptr\"]" in text
    assert "matrix[\"indices\"]" in text
    assert "matrix[\"data\"]" in text


def test_sparse_extractor_declares_expected_output():
    text = EXTRACT_SCRIPT.read_text(encoding="utf-8")
    assert "sparse_gene_panel_expression.tsv.gz" in text
    assert "nonzero" in text.lower()
