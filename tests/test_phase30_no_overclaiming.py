from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PHASE30_SCRIPTS = [
    ROOT / "scripts/155_audit_gse174367_bulk_sample_mapping.py",
    ROOT / "scripts/150_convert_gse174367_bulk_rda_to_axis_matrix.py",
    ROOT / "scripts/151_build_gse174367_bulk_axis_scores.py",
    ROOT / "scripts/156_test_gse174367_bulk_ad_axis_replication_clean.py",
]


def test_phase30_scripts_avoid_disallowed_single_cell_and_dense_workflows():
    forbidden = ["scanpy", "anndata", "read_h5ad", "umap", "leiden", "louvain", "cellranger", "fasterq", "prefetch", "todense(", "toarray("]
    for script in PHASE30_SCRIPTS:
        text = script.read_text(encoding="utf-8").lower()
        for phrase in forbidden:
            assert phrase not in text


def test_phase30_scripts_do_not_make_affirmative_clinical_or_causal_claims():
    forbidden_claims = [
        "clinical-grade",
        "diagnostic tool",
        "causal mechanism",
        "validated mechanism",
        "pnas-ready claim",
        "validated across diseases",
    ]
    combined = "\n".join(script.read_text(encoding="utf-8").lower() for script in PHASE30_SCRIPTS)
    for phrase in forbidden_claims:
        assert phrase not in combined
