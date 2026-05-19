import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/165_test_pd_axis_replication_microarray.py"
INTEGRATION = ROOT / "scripts/122_integrate_endpoint_locked_replication.py"
READINESS = ROOT / "scripts/123_build_pnas_readiness_matrix.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_direction_only_pd_signal_does_not_become_supported_replication():
    module = load_module(SCRIPT, "phase33_pd_replication")
    label = module.evidence_label(effect=-0.25, reference_direction=-1, pvalue=0.4, fdr=0.9, n=30, positive_n=15, negative_n=15)
    assert label == "directionally_consistent_but_not_significant"
    supported = module.evidence_label(effect=-0.25, reference_direction=-1, pvalue=0.03, fdr=0.2, n=30, positive_n=15, negative_n=15)
    assert supported == "statistically_supported_pd_replication"


def test_phase33_integration_requires_statistical_support_for_upgrade():
    integration = load_module(INTEGRATION, "phase33_integration")
    discovery = [{"axis_id": "neuronal_vulnerability_axis", "endpoint_id": "gse243639_pd_diagnosis", "cohort": "gse243639_pd_snpc", "effect_size": "-0.2", "axis_claim_class": "axis_level_preliminary_evidence"}]
    replication = {
        "mock_pd": [
            {
                "axis_id": "neuronal_vulnerability_axis",
                "effect_size": "-0.3",
                "pvalue": "0.5",
                "fdr": "0.9",
                "n": "30",
                "positive_n": "15",
                "negative_n": "15",
                "evidence_label": "directionally_consistent_but_not_significant",
            }
        ]
    }
    rows = integration.integrate(discovery, replication)
    assert rows[0]["replication_status"] == "directionally_consistent_preliminary_signal"
    assert rows[0]["claim_upgrade_allowed"] == "false"


def test_readiness_source_knows_phase33_pd_replication_outputs():
    text = READINESS.read_text(encoding="utf-8")
    assert "phase33_*_pd_axis_replication_statistics.tsv" in text
    assert "statistically_supported_pd_replication" in text
