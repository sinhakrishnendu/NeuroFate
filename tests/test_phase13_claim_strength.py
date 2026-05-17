from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/65_build_claim_strength_table.py"


def test_claim_strength_schema_is_defined() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for column in [
        "claim_id",
        "task",
        "model",
        "evidence_layer",
        "internal_auroc",
        "internal_auroc_sd",
        "internal_auprc",
        "balanced_accuracy",
        "brier_score",
        "permutation_empirical_p",
        "feature_ablation_support",
        "external_validation_status",
        "external_sample_units",
        "leakage_status",
        "overclaiming_status",
        "claim_strength",
        "allowed_claim_text",
        "disallowed_claim_text",
        "reviewer_risk",
    ]:
        assert column in text


def test_claim_strength_rules_are_conservative() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "strong_internal" in text
    assert "preliminary_external_feasibility" in text
    assert "failed_or_unstable" in text
    assert "apoe_risk_prediction" in text
    assert "Do not claim clinical utility" in text
    assert "definitive cross-cohort validation" in text
