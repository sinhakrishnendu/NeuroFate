import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INTEGRATE = ROOT / "scripts/122_integrate_endpoint_locked_replication.py"
READINESS = ROOT / "scripts/123_build_pnas_readiness_matrix.py"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_phase29_direction_only_does_not_upgrade_claims():
    module = load(INTEGRATE, "phase29_integrate")
    discovery = [{"axis_id": "axis", "effect_size": "0.25", "endpoint_id": "sea_ad_cognitive_dementia", "cohort": "sea_ad"}]
    replication = {
        "gse174367_ad_multiomics_bulk": [
            {
                "axis_id": "axis",
                "effect_size": "0.2",
                "pvalue": "0.5",
                "fdr": "0.8",
                "n": "230",
                "positive_n": "118",
                "negative_n": "112",
                "evidence_label": "directionally_consistent_but_not_significant",
            }
        ]
    }
    rows = module.integrate(discovery, replication)
    assert rows[0]["replication_status"] == "directionally_consistent_preliminary_signal"
    assert rows[0]["claim_upgrade_allowed"] == "false"


def test_phase29_supported_ad_replication_can_upgrade_when_thresholds_pass():
    module = load(INTEGRATE, "phase29_integrate_supported")
    discovery = [{"axis_id": "axis", "effect_size": "0.25", "endpoint_id": "sea_ad_cognitive_dementia", "cohort": "sea_ad"}]
    replication = {
        "gse174367_ad_multiomics_bulk": [
            {
                "axis_id": "axis",
                "effect_size": "0.2",
                "pvalue": "0.01",
                "fdr": "0.05",
                "n": "230",
                "positive_n": "118",
                "negative_n": "112",
                "evidence_label": "statistically_supported_ad_replication",
            }
        ]
    }
    rows = module.integrate(discovery, replication)
    assert rows[0]["replication_status"] == "statistically_supported_replication"
    assert rows[0]["claim_upgrade_allowed"] == "true"


def test_readiness_recognizes_phase29_ad_replication(tmp_path, monkeypatch):
    module = load(READINESS, "phase29_readiness")
    monkeypatch.chdir(tmp_path)
    path = tmp_path / "results/tables/phase29_gse174367_bulk_axis_replication_statistics.tsv"
    path.parent.mkdir(parents=True)
    path.write_text("axis_id\tevidence_label\naxis\tdirectionally_consistent_but_not_significant\n", encoding="utf-8")
    assert module.phase28_ad_replication_available() is True
    assert module.phase28_ad_statistical_support() is False
    path.write_text("axis_id\tevidence_label\naxis\tstatistically_supported_ad_replication\n", encoding="utf-8")
    assert module.phase28_ad_statistical_support() is True
