from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_replication_integration_does_not_upgrade_without_replication():
    text = (ROOT / "scripts/122_integrate_endpoint_locked_replication.py").read_text(encoding="utf-8")
    assert "replication_pending" in text
    assert "preliminary_replicated_candidate" in text
    assert '"claim_upgrade_allowed": "false"' in text
    assert "No independent replication statistics supplied" in text
    assert "Do not claim validated shared AD/PD mechanism" in text


def test_pnas_readiness_matrix_has_required_criteria():
    text = (ROOT / "scripts/123_build_pnas_readiness_matrix.py").read_text(encoding="utf-8")
    for criterion in [
        "software_reproducibility",
        "sea_ad_internal_evidence",
        "pd_gse243639_evidence",
        "endpoint_locked_axis_evidence",
        "matched_random_controls",
        "independent_ad_replication",
        "independent_pd_replication",
        "network_pathway_interpretation",
        "no_overclaiming_audit",
        "public_reproducibility_package",
    ]:
        assert criterion in text
