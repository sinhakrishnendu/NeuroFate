from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/73_prepare_external_sparse_extraction_plan.py"


def test_extraction_planner_supports_required_formats() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for fmt in [
        "h5ad_csr",
        "csv_genes_as_rows",
        "csv_cells_as_rows",
        "mtx_features_barcodes",
        "bulk_matrix",
    ]:
        assert fmt in text


def test_extraction_planner_generates_manual_templates_only() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "RUN_MANUAL_EXTRACTION" in text
    assert "MANUAL_HEAVY" in text
    assert "Manual user execution only" in text
    assert "subprocess" not in text
    assert "os.system" not in text
    assert "check_call" not in text


def test_extraction_planner_does_not_import_scanpy_or_read_h5ad() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "scanpy" not in text
    assert "read_h5ad" not in text
    assert "anndata" not in text
