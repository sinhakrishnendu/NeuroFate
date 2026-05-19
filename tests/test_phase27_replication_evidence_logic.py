import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLEAN = ROOT / "scripts/140_test_gse184950_axis_replication_clean.py"
INTEGRATE = ROOT / "scripts/122_integrate_endpoint_locked_replication.py"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_direction_only_is_not_statistically_supported_replication():
    module = load(CLEAN, "phase27_clean_rep")
    label, consistency = module.evidence_label(
        effect=0.2,
        pvalue=0.5,
        fdr=0.9,
        phase22="positive",
        n=34,
        positive_n=24,
        negative_n=10,
    )
    assert consistency == "consistent"
    assert label == "directionally_consistent_but_not_significant"


def test_statistical_support_required_for_upgrade():
    module = load(INTEGRATE, "phase27_integrate")
    discovery = [{"axis_id": "axis", "effect_size": "0.3", "endpoint_id": "gse243639_pd_diagnosis", "cohort": "gse243639"}]
    replication = {
        "gse184950_pd_sn": [
            {
                "axis_id": "axis",
                "effect_size": "0.2",
                "pvalue": "0.5",
                "fdr": "0.91",
                "n": "34",
                "positive_n": "24",
                "negative_n": "10",
                "evidence_label": "directionally_consistent_but_not_significant",
            }
        ]
    }
    rows = module.integrate(discovery, replication)
    assert rows[0]["replication_status"] == "directionally_consistent_preliminary_signal"
    assert rows[0]["claim_upgrade_allowed"] == "false"
