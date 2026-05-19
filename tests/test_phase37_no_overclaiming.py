from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FILES = [
    ROOT / "scripts/163_parse_pd_geo_series_matrix_metadata.py",
    ROOT / "scripts/169_audit_gse7621_sample_mapping.py",
    ROOT / "scripts/165_build_pd_axis_scores_from_geo_expression.py",
    ROOT / "scripts/166_test_phase34_pd_axis_replication.py",
    ROOT / "scripts/170_generate_phase37_gse7621_pd_replication_report.py",
    ROOT / "scripts/122_integrate_endpoint_locked_replication.py",
    ROOT / "scripts/123_build_pnas_readiness_matrix.py",
]


def test_phase37_sources_avoid_forbidden_processing_and_single_cell_routes():
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in FILES)
    for forbidden in [
        "scanpy",
        "anndata",
        "read_h5ad",
        "h5ad",
        "umap",
        "leiden",
        "louvain",
        "process fastq",
        "sra-tools",
        "toarray(",
        "todense(",
    ]:
        assert forbidden not in combined


def test_phase37_sources_do_not_make_affirmative_strong_claims():
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in FILES)
    for phrase in [
        "clinical-grade",
        "diagnostic tool",
        "causal mechanism",
        "validated shared mechanism",
        "validated shared ad/pd",
        "definitive shared mechanism",
    ]:
        assert phrase not in combined
    assert "direction-only support remains preliminary" in combined
