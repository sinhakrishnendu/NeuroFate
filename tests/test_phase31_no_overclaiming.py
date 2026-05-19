from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PHASE31_FILES = [
    ROOT / "scripts/157_audit_gse174367_bulk_gene_identifiers.py",
    ROOT / "scripts/150_convert_gse174367_bulk_rda_to_axis_matrix.py",
    ROOT / "scripts/151_build_gse174367_bulk_axis_scores.py",
    ROOT / "scripts/156_test_gse174367_bulk_ad_axis_replication_clean.py",
    ROOT / "scripts/158_generate_phase31_gse174367_gene_mapping_report.py",
]


def test_phase31_files_avoid_disallowed_workflows():
    forbidden = ["scanpy", "anndata", "read_h5ad", "umap", "leiden", "louvain", "cellranger", "fasterq", "prefetch", "todense(", "toarray("]
    for path in PHASE31_FILES:
        text = path.read_text(encoding="utf-8").lower()
        for phrase in forbidden:
            assert phrase not in text


def test_phase31_files_avoid_affirmative_overclaims():
    forbidden_claims = [
        "clinical-grade",
        "diagnostic tool",
        "causal mechanism",
        "validated mechanism",
        "validated across diseases",
    ]
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in PHASE31_FILES)
    for phrase in forbidden_claims:
        assert phrase not in combined
