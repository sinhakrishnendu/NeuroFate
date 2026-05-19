from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_endpoint_locked_evidence_schema_and_rules_are_defined():
    text = (ROOT / "scripts/115_build_endpoint_locked_axis_evidence_table.py").read_text(encoding="utf-8")
    for column in [
        "axis_id",
        "endpoint_id",
        "cohort",
        "endpoint_role",
        "effect_size",
        "standardized_mean_difference",
        "pvalue",
        "fdr",
        "empirical_pvalue",
        "axis_claim_class",
        "allowed_claim",
        "disallowed_claim",
    ]:
        assert column in text
    assert "Strong claims are not allowed yet" in text
    assert "candidate_shared_axis_endpoint_locked" in text
    assert "preliminary_candidate" in text


def test_claim_strength_table_supports_phase22_and_supersedes_phase21():
    text = (ROOT / "scripts/65_build_claim_strength_table.py").read_text(encoding="utf-8")
    assert "--phase22-axis-evidence" in text
    assert "phase22_endpoint_locked_axis_evidence_table.tsv" in text
    assert "phase22_endpoint_locked_axis_claim_strength_delta.tsv" in text
    assert "phase22_endpoint_locked_axis_biology" in text
    assert "supersedes Phase 21" in text
