from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FILES = [
    ROOT / "scripts/159_build_crosscohort_axis_evidence_summary.py",
    ROOT / "scripts/160_generate_phase32_pnas_decision_report.py",
    ROOT / "scripts/161_generate_phase32_crosscohort_figures.py",
]


def test_phase32_sources_do_not_use_disallowed_processing_tools():
    forbidden = ["scanpy", "anndata", "read_h5ad", "umap", "leiden", "louvain", "todense(", "toarray("]
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in FILES)
    for phrase in forbidden:
        assert phrase not in combined


def test_phase32_sources_do_not_make_affirmative_overclaims():
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in FILES)
    affirmative_forbidden = [
        "clinical-grade",
        "diagnostic tool",
        "causal mechanism",
        "pnas-ready mechanism-level evidence",
        "definitive shared ad/pd claim",
    ]
    for phrase in affirmative_forbidden:
        assert phrase not in combined
    assert "do not claim" in combined
