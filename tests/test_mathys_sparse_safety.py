from pathlib import Path


PLAN_SCRIPT = Path("scripts/38_prepare_mathys_sparse_extraction.py")
FEATURE_SCRIPT = Path("scripts/39_build_mathys_donor_feature_table.py")
EXTERNAL_EXTRACTION_SCRIPT = Path("scripts/31_sparse_external_gene_extraction.py")


def test_mathys_sparse_planner_does_not_open_expression_files():
    text = PLAN_SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "import h5py" not in lowered
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered
    assert "toarray(" not in lowered
    assert "todense(" not in lowered
    assert "mathys_sparse_extraction_plan.tsv" in text
    assert "manual_mathys_sparse_extraction_template.sh" in text
    assert "RUN_MANUAL_EXTERNAL_EXTRACTION" in text


def test_external_sparse_extraction_has_mathys_compatible_safety_limits():
    text = EXTERNAL_EXTRACTION_SCRIPT.read_text(encoding="utf-8")
    assert "MAX_GENES = 64" in text
    assert "MAX_CHUNK_SIZE = 50000" in text
    assert "memory_limit_mb" in text
    assert "--dry-run" in text
    assert "--execute" in text
    assert "RUN_MANUAL_EXTERNAL_EXTRACTION" in text
    assert "No external expression file was opened" in text
    lowered = text.lower()
    assert "toarray(" not in lowered
    assert "todense(" not in lowered
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered


def test_mathys_donor_feature_builder_uses_sparse_tsvs_only():
    text = FEATURE_SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "import h5py" not in lowered
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered
    assert "toarray(" not in lowered
    assert "todense(" not in lowered
    assert "sparse_gene_panel_expression.tsv.gz" in text
    assert "mathys_2019_phase5_donor_feature_table.tsv" in text
    assert "mathys_2019_feature_schema_alignment.tsv" in text
    assert "phase5_donor_feature_table.tsv" in text
