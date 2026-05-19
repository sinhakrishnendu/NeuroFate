from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_axis_comparison_uses_conservative_categories():
    text = (ROOT / "scripts/108_compare_ad_pd_axis_patterns.py").read_text(encoding="utf-8")
    for category in [
        "shared_ad_pd_candidate",
        "ad_enriched_axis",
        "pd_enriched_axis",
        "inconclusive_axis",
        "insufficient_coverage",
        "axis_level_preliminary_evidence",
        "axis_level_insufficient_validation",
    ]:
        assert category in text
    assert "validated across diseases" in text
    assert "Do not claim" in text


def test_claim_strength_table_accepts_phase21_axis_rows():
    text = (ROOT / "scripts/65_build_claim_strength_table.py").read_text(encoding="utf-8")
    assert "--phase21-axis-claims" in text
    assert "phase21_axis_claim_strength.tsv" in text
    assert "phase21_axis_claim_strength_delta.tsv" in text
    assert "phase21_axis_biology" in text
    assert "axis_level_preliminary_evidence" in text
