import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLAIM = ROOT / "scripts/65_build_claim_strength_table.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase32_claims", CLAIM)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_claim_strength_table_mentions_phase32_and_nominal_language():
    text = CLAIM.read_text(encoding="utf-8")
    assert "phase32_crosscohort_axis_consolidation" in text
    assert "nominal independent AD replication" in text
    assert "not as a definitive mechanism" in text


def test_phase32_claim_rows_do_not_turn_nominal_ad_into_shared_claim(tmp_path):
    module = load_module()
    evidence = tmp_path / "phase32.tsv"
    evidence.write_text(
        "axis_id\tcrosscohort_evidence_class\tgse174367_p\tsafe_claim\n"
        "neuronal_vulnerability_axis\tstrong_ad_axis_with_nominal_external_replication\t0.0299\tNominal AD support.\n",
        encoding="utf-8",
    )
    args = type(
        "Args",
        (),
        {
            "phase32_crosscohort_evidence": evidence,
            "leakage_audit": tmp_path / "missing.tsv",
            "overclaiming_audit": tmp_path / "missing2.tsv",
        },
    )()
    rows = []
    module.append_phase32_crosscohort_claims(rows, args)
    assert rows[0]["external_validation_status"] == "strong_ad_axis_with_nominal_external_replication"
    assert "AD replication candidate" in rows[0]["allowed_claim_text"]
    assert "definitive shared AD/PD mechanism" in rows[0]["disallowed_claim_text"]
