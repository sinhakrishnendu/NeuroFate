from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FILES = [
    ROOT / "scripts/171_audit_gse7621_axis_direction_and_probe_mapping.py",
    ROOT / "scripts/172_summarize_gse7621_axis_score_distributions.py",
    ROOT / "scripts/173_generate_phase38_gse7621_interpretation_report.py",
    ROOT / "scripts/174_generate_phase38_gse7621_figures.py",
    ROOT / "scripts/159_build_crosscohort_axis_evidence_summary.py",
    ROOT / "scripts/123_build_pnas_readiness_matrix.py",
]


def test_phase38_sources_do_not_use_forbidden_processing_routes():
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in FILES)
    for forbidden in [
        "scanpy",
        "anndata",
        "read_h5ad",
        "h5ad",
        "leiden",
        "louvain",
        "process fastq",
        "sra-tools",
        "toarray(",
        "todense(",
    ]:
        assert forbidden not in combined


def test_phase38_sources_do_not_make_affirmative_clinical_or_shared_claims():
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in FILES)
    for phrase in [
        "clinical-grade",
        "diagnostic tool",
        "causal mechanism",
        "is a validated shared",
        "validated shared mechanism",
        "definitive shared mechanism",
    ]:
        assert phrase not in combined
    assert "not shared ad/pd replication" in combined
    assert "pd-divergent" in combined
