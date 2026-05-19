from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FILES = [
    ROOT / "scripts/168_audit_gse20141_sample_mapping.py",
    ROOT / "scripts/165_build_pd_axis_scores_from_geo_expression.py",
    ROOT / "scripts/166_test_phase34_pd_axis_replication.py",
    ROOT / "scripts/122_integrate_endpoint_locked_replication.py",
    ROOT / "scripts/123_build_pnas_readiness_matrix.py",
    ROOT / "scripts/167_generate_phase34_pd_replication_report.py",
]


def test_phase35_sources_avoid_disallowed_processing():
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in FILES)
    for forbidden in ["scanpy", "anndata", "read_h5ad", "umap", "leiden", "louvain", "cel processing", "chp processing", "toarray(", "todense("]:
        assert forbidden not in combined


def test_phase35_sources_do_not_make_strong_claims():
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in FILES)
    for phrase in ["clinical-grade", "diagnostic tool", "causal mechanism", "validated shared mechanism", "validated shared ad/pd"]:
        assert phrase not in combined
    assert "statistically_supported_pd_replication" in combined
