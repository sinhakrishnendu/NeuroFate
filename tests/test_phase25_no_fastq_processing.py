from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PHASE25_SCRIPTS = [
    ROOT / "scripts/131_parse_gse184950_series_matrix.py",
    ROOT / "scripts/132_reconcile_gse184950_archive_with_series_metadata.py",
    ROOT / "scripts/133_plan_gse184950_selective_tar_extraction.py",
    ROOT / "scripts/134_generate_phase25_gse184950_report.py",
]


def test_phase25_scripts_avoid_heavy_biology_workflows():
    forbidden = [
        "scanpy",
        "anndata",
        "read_h5ad",
        "cellranger",
        "fasterq",
        "prefetch ",
        "umap(",
        "leiden",
        "cluster_cells",
        "toarray(",
        "todense(",
        "fit(",
    ]
    for path in PHASE25_SCRIPTS:
        text = path.read_text(encoding="utf-8").lower()
        for term in forbidden:
            assert term not in text, f"{term} found in {path}"


def test_phase25_docs_are_conservative():
    combined = "\n".join(
        [
            (ROOT / "README.md").read_text(encoding="utf-8"),
            (ROOT / "PNAS_DISCOVERY_STRATEGY.md").read_text(encoding="utf-8"),
            (ROOT / "docs/pnas_validation_strategy.md").read_text(encoding="utf-8"),
            (ROOT / "docs/external_validation_expansion.md").read_text(encoding="utf-8"),
        ]
    ).lower()
    assert "phase 25" in combined
    assert "series matrix" in combined
    assert "pd/pdd" in combined
    assert "fastq/sra processing remains out of scope" in combined
    for overclaim in [
        "phase 25 proves",
        "phase 25 is clinical-grade",
        "phase 25 is a diagnostic tool",
        "phase 25 provides causal proof",
        "phase 25 is validated across diseases",
    ]:
        assert overclaim not in combined
