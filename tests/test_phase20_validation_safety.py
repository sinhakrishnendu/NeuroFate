from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "scripts/102_run_gse243639_phase20_celltype_pd_validation.py"
PHASE91 = ROOT / "scripts/91_run_gse243639_repaired_celltype_pd_validation.py"
COMPARISON = ROOT / "scripts/103_compare_phase16_17_18_20_pd_validation.py"


def test_phase20_validation_uses_sample_level_outputs_and_conservative_labels() -> None:
    text = VALIDATION.read_text(encoding="utf-8")
    assert "phase20_gse243639_celltype_feature_table.tsv" in text
    assert "phase20_gse243639_celltype_validation_metrics.tsv" in text
    for label in [
        "moderate_pd_internal_validation",
        "preliminary_pd_internal_signal",
        "weak_pd_signal",
        "technical_failure_annotation_join",
    ]:
        assert label in PHASE91.read_text(encoding="utf-8")


def test_phase20_validation_safety_forbidden_workloads_absent() -> None:
    combined = "\n".join([VALIDATION.read_text(encoding="utf-8"), COMPARISON.read_text(encoding="utf-8")]).lower()
    for forbidden in ["scanpy", "anndata", "read_h5ad", "umap.fit", "fit_transform", "leiden", "neighbors"]:
        assert forbidden not in combined


def test_phase20_comparison_marks_prior_phases_correctly() -> None:
    text = COMPARISON.read_text(encoding="utf-8")
    assert "Phase 16 is valid as the global sample-level PD extension" in text
    assert "Phase 17 was technically defective" in text
    assert "Phase 18 still failed annotation joining" in text
    assert "Phase 20 is the corrected safe-map-based cell-type-aware analysis" in text
