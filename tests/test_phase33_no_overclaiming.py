from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FILES = [
    ROOT / "scripts/162_triage_pd_replication_geo_files.py",
    ROOT / "scripts/163_build_axis_scores_from_geo_series_matrix.py",
    ROOT / "scripts/164_prepare_geo_platform_gene_mapping.py",
    ROOT / "scripts/165_test_pd_axis_replication_microarray.py",
    ROOT / "scripts/166_generate_phase33_pd_replication_report.py",
]

DOC_FILES = [
    ROOT / "README.md",
    ROOT / "PNAS_DISCOVERY_STRATEGY.md",
    ROOT / "docs/pnas_validation_strategy.md",
    ROOT / "docs/external_validation_expansion.md",
]


def test_phase33_sources_do_not_use_disallowed_processing_tools():
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in FILES)
    for forbidden in ["scanpy", "anndata", "read_h5ad", "umap", "leiden", "louvain", "toarray(", "todense("]:
        assert forbidden not in combined


def test_phase33_sources_do_not_make_affirmative_overclaims():
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in FILES)
    docs = "\n".join(path.read_text(encoding="utf-8").lower() for path in DOC_FILES)
    affirmative_forbidden = [
        "diagnostic tool",
        "validated shared mechanism",
        "validated shared ad/pd",
        "definitive shared ad/pd mechanism claim",
    ]
    for phrase in affirmative_forbidden:
        assert phrase not in combined
    assert "direction-only" in combined + docs
    assert "does not download data automatically" in docs
    assert "clinical" in docs and "out of scope" in docs
