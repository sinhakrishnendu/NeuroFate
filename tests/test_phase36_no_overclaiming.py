from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FILES = [
    ROOT / "scripts/165_build_pd_axis_scores_from_geo_expression.py",
    ROOT / "scripts/166_test_phase34_pd_axis_replication.py",
    ROOT / "scripts/122_integrate_endpoint_locked_replication.py",
    ROOT / "scripts/123_build_pnas_readiness_matrix.py",
    ROOT / "scripts/167_generate_phase34_pd_replication_report.py",
]


def test_phase36_sources_avoid_forbidden_processing_paths():
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
        "sra",
        "cel/chp",
        "toarray(",
        "todense(",
    ]:
        assert forbidden not in combined


def test_phase36_sources_do_not_make_affirmative_strong_claims():
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
