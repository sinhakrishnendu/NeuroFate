from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PHASE26_SCRIPTS = [
    ROOT / "scripts/135_inspect_gse184950_nested_archives.py",
    ROOT / "scripts/136_extract_gse184950_processed_matrices_selective.py",
    ROOT / "scripts/137_generate_phase26_gse184950_replication_report.py",
    ROOT / "scripts/138_generate_phase26_gse184950_figures.py",
]


def test_phase26_scripts_do_not_use_disallowed_workflows():
    forbidden = [
        "scanpy",
        "anndata",
        "read_h5ad",
        "cellranger",
        "fasterq",
        "prefetch ",
        "umap(",
        "leiden",
        "toarray(",
        "todense(",
        "fit(",
    ]
    for path in PHASE26_SCRIPTS:
        text = path.read_text(encoding="utf-8").lower()
        for term in forbidden:
            assert term not in text, f"{term} found in {path}"


def test_phase26_extraction_is_matrix_component_limited():
    text = (ROOT / "scripts/136_extract_gse184950_processed_matrices_selective.py").read_text(encoding="utf-8")
    assert 'RUN_MANUAL_GSE184950_EXTRACTION") == "YES"' in text
    assert "ALLOWED_NAMES" in text
    assert "matrix.mtx.gz" in text
    assert "features.tsv.gz" in text
    assert "genes.tsv.gz" in text
    assert "barcodes.tsv.gz" in text
