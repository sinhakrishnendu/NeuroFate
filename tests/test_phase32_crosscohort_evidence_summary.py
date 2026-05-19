import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/159_build_crosscohort_axis_evidence_summary.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase32_summary", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_neuronal_axis_nominal_ad_replication_is_ranked_top():
    module = load_module()
    phase22 = [
        {"axis_id": "neuronal_vulnerability_axis", "cohort": "sea_ad", "endpoint_id": "sea_ad_cognitive_dementia", "effect_size": "-0.34", "pvalue": "0.006", "fdr": "0.015"},
        {"axis_id": "astrocyte_stress_axis", "cohort": "sea_ad", "endpoint_id": "sea_ad_cognitive_dementia", "effect_size": "0.33", "pvalue": "0.008", "fdr": "0.017"},
    ]
    gse174367 = [
        {"axis_id": "neuronal_vulnerability_axis", "effect_size": "-0.2658", "pvalue": "0.0299", "fdr": "0.2438", "directional_consistency": "consistent"},
        {"axis_id": "astrocyte_stress_axis", "effect_size": "0.18", "pvalue": "0.135", "fdr": "0.338", "directional_consistency": "consistent"},
    ]
    gse184950 = [
        {"axis_id": "neuronal_vulnerability_axis", "effect_size": "-0.11", "pvalue": "0.59", "fdr": "0.91", "directional_consistency": "consistent"}
    ]
    gse243639 = [{"model": "logistic_regression", "validation_mode": "repeated_stratified_split", "reliability_flag": "preliminary_pd_internal_signal"}]
    rows = module.build_summary(phase22, gse174367, gse184950, gse243639)
    ranked = module.rank_rows(rows)
    assert ranked[0]["axis_id"] == "neuronal_vulnerability_axis"
    assert ranked[0]["crosscohort_evidence_class"] == "strong_ad_axis_with_nominal_external_replication"
    assert "nominal independent AD replication" in ranked[0]["safe_claim"]


def test_direction_only_pd_does_not_create_shared_mechanism_claim():
    module = load_module()
    klass, _safe, unsafe, _next = module.classify_axis(
        "axis",
        {"effect_size": "-0.3", "pvalue": "0.01", "fdr": "0.02"},
        {"effect_size": "-0.2", "pvalue": "0.03", "fdr": "0.24"},
        {"effect_size": "-0.1", "pvalue": "0.8", "fdr": "0.9", "directional_consistency": "consistent"},
        {"effect_size": "-0.1", "pvalue": "0.2", "fdr": "0.9", "empirical_pvalue": "0.01"},
        "preliminary_pd_internal_signal",
    )
    assert klass == "strong_ad_axis_with_nominal_external_replication"
    assert "Do not claim" in unsafe
