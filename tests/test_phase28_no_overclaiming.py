from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PHASE28_SCRIPTS = [
    ROOT / "scripts/143_parse_geo_series_matrix_generic.py",
    ROOT / "scripts/144_triage_ad_replication_files.py",
    ROOT / "scripts/145_build_ad_replication_axis_scores_from_matrix.py",
    ROOT / "scripts/146_plan_ad_snrna_replication_extraction.py",
    ROOT / "scripts/147_test_ad_replication_axis_associations.py",
    ROOT / "scripts/148_generate_phase28_ad_replication_report.py",
]


def test_phase28_scripts_avoid_disallowed_workflows():
    forbidden = ["scanpy", "anndata", "read_h5ad", "cellranger", "fasterq", "prefetch ", "umap(", "leiden", "toarray(", "todense("]
    for path in PHASE28_SCRIPTS:
        text = path.read_text(encoding="utf-8").lower()
        for term in forbidden:
            assert term not in text, f"{term} found in {path}"


def test_phase28_docs_do_not_overclaim_ad_replication():
    combined = "\n".join(
        [
            (ROOT / "README.md").read_text(encoding="utf-8"),
            (ROOT / "PNAS_DISCOVERY_STRATEGY.md").read_text(encoding="utf-8"),
            (ROOT / "docs/pnas_validation_strategy.md").read_text(encoding="utf-8"),
            (ROOT / "docs/external_validation_expansion.md").read_text(encoding="utf-8"),
        ]
    ).lower()
    assert "phase 28" in combined
    assert "independent ad replication" in combined
    for phrase in [
        "phase 28 proves",
        "phase 28 validates",
        "phase 28 is a clinical diagnostic",
        "phase 28 proves a causal ad mechanism",
        "phase 28 validates a shared mechanism",
    ]:
        assert phrase not in combined
