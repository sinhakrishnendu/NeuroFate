from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PHASE29_SCRIPTS = [
    ROOT / "scripts/149_inspect_gse174367_bulk_rda.py",
    ROOT / "scripts/150_convert_gse174367_bulk_rda_to_axis_matrix.py",
    ROOT / "scripts/151_build_gse174367_bulk_axis_scores.py",
    ROOT / "scripts/152_test_gse174367_bulk_ad_axis_replication.py",
    ROOT / "scripts/153_generate_phase29_gse174367_ad_replication_report.py",
    ROOT / "scripts/154_generate_phase29_gse174367_figures.py",
]


def test_phase29_scripts_avoid_disallowed_single_cell_and_dense_workflows():
    forbidden = ["scanpy", "anndata", "read_h5ad", "toarray(", "todense(", "umap", "leiden", "louvain", "cellranger", "fasterq", "prefetch"]
    for script in PHASE29_SCRIPTS:
        text = script.read_text(encoding="utf-8").lower()
        for phrase in forbidden:
            assert phrase not in text


def test_phase29_documentation_keeps_claims_conservative():
    readme = (ROOT / "README.md").read_text(encoding="utf-8").lower()
    phase29 = readme.split("## phase 29:", 1)[1]
    assert "if the rda structure or sample mapping is ambiguous, it stops" in phase29
    assert "clinical, diagnostic, causal, or pnas-ready claims" in phase29
    assert "validates shared" not in phase29
    assert "clinical-grade" not in phase29
