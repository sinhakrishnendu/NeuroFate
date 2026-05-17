from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "scripts/64_audit_mathys_gene_extraction.py"
PHASE9 = ROOT / "scripts/45_generate_phase9_results_text.py"


def test_mathys_gene_audit_script_exists_and_is_safe() -> None:
    text = AUDIT.read_text(encoding="utf-8").lower()
    assert AUDIT.exists()
    for forbidden in ["scanpy", "read_h5ad", "anndata", "import h5py", "pca", "umap", "clustering"]:
        assert forbidden not in text
    assert "phase13_mathys_gene_extraction_audit.tsv" in text


def test_phase9_text_reports_extraction_counts_not_zero_overlap_fallback() -> None:
    text = PHASE9.read_text(encoding="utf-8")
    for token in [
        "Requested target genes",
        "Extracted target genes",
        "Sparse expression rows extracted",
        "Gene-overlap table status",
        "gene overlap table unavailable",
    ]:
        assert token in text
    assert "present=0, missing=0" not in text
