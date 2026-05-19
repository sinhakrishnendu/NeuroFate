from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = [
    ROOT / "scripts/162_triage_phase34_pd_geo_files.py",
    ROOT / "scripts/163_parse_pd_geo_series_matrix_metadata.py",
    ROOT / "scripts/164_prepare_phase34_platform_axis_probe_mapping.py",
    ROOT / "scripts/165_build_pd_axis_scores_from_geo_expression.py",
    ROOT / "scripts/166_test_phase34_pd_axis_replication.py",
    ROOT / "scripts/167_generate_phase34_pd_replication_report.py",
]


def test_phase34_sources_avoid_single_cell_and_dense_tools():
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in SCRIPTS)
    for forbidden in ["scanpy", "anndata", "read_h5ad", "umap", "leiden", "louvain", "toarray(", "todense("]:
        assert forbidden not in combined


def test_phase34_sources_do_not_make_strong_shared_or_clinical_claims():
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in SCRIPTS)
    for phrase in ["clinical-grade", "diagnostic tool", "causal mechanism", "validated shared mechanism", "validated shared ad/pd"]:
        assert phrase not in combined
    assert "direction-only" in combined
    assert "does not establish confirmed cross-disease biology" in combined
