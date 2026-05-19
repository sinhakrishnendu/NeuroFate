import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "scripts/159_build_crosscohort_axis_evidence_summary.py"


def load_summary():
    spec = importlib.util.spec_from_file_location("summary", SUMMARY)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_opposite_direction_significant_pd_axis_is_divergent_not_shared():
    summary = load_summary()
    klass, safe, unsafe, next_step = summary.classify_axis(
        "synuclein_mitochondrial_axis",
        {"effect_size": "0.4", "pvalue": "0.01", "fdr": "0.02"},
        {},
        {},
        {},
        {"effect_size": "-0.76", "pvalue": "0.0018", "fdr": "0.018", "evidence_label": "opposite_direction"},
        "preliminary_pd_internal_signal",
    )
    assert klass == "pd_divergent_axis_candidate"
    assert "candidate PD-divergent axis" in safe
    assert "not as shared AD/PD replication" in safe
    assert "validated shared" in unsafe
    assert "another PD cohort" in next_step


def test_direction_only_neuronal_axis_does_not_upgrade_shared_claim():
    summary = load_summary()
    klass, safe, _unsafe, _next_step = summary.classify_axis(
        "neuronal_vulnerability_axis",
        {"effect_size": "-0.35", "pvalue": "0.006", "fdr": "0.015"},
        {"effect_size": "-0.26", "pvalue": "0.03", "fdr": "0.24"},
        {},
        {},
        {"effect_size": "-0.17", "pvalue": "0.497", "fdr": "0.71", "evidence_label": "directionally_consistent_but_not_significant"},
        "preliminary_pd_internal_signal",
    )
    assert klass == "strong_ad_axis_with_nominal_external_replication"
    assert "not as a shared mechanism" in safe or "not a definitive mechanism" in safe
