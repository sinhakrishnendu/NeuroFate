import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPLICATION = ROOT / "scripts/166_test_phase34_pd_axis_replication.py"
INTEGRATION = ROOT / "scripts/122_integrate_endpoint_locked_replication.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_direction_only_phase37_pd_evidence_does_not_upgrade_claim():
    integration = load_module(INTEGRATION, "integration")
    discovery = [
        {
            "axis_id": "neuronal_vulnerability_axis",
            "endpoint_id": "gse243639_pd_diagnosis",
            "cohort": "gse243639_pd_snpc",
            "effect_size": "-0.4",
            "axis_claim_class": "axis_level_preliminary_evidence",
        }
    ]
    replication = {
        "gse7621_pd_sn_bulk": [
            {
                "axis_id": "neuronal_vulnerability_axis",
                "effect_size": "-0.2",
                "pvalue": "0.25",
                "fdr": "0.8",
                "n": "25",
                "positive_n": "16",
                "negative_n": "9",
                "evidence_label": "directionally_consistent_but_not_significant",
            }
        ]
    }
    rows = integration.integrate(discovery, replication)
    assert rows[0]["replication_status"] == "directionally_consistent_preliminary_signal"
    assert rows[0]["claim_upgrade_allowed"] == "false"


def test_phase37_pd_replication_label_requires_p_or_fdr_support():
    repl = load_module(REPLICATION, "repl")
    assert (
        repl.evidence_label(
            effect=-0.2,
            ref_direction=-1,
            pvalue=0.25,
            fdr=0.8,
            n=25,
            positive_n=16,
            negative_n=9,
        )
        == "directionally_consistent_but_not_significant"
    )
    assert (
        repl.evidence_label(
            effect=-0.2,
            ref_direction=-1,
            pvalue=0.03,
            fdr=0.2,
            n=25,
            positive_n=16,
            negative_n=9,
        )
        == "statistically_supported_pd_replication"
    )
