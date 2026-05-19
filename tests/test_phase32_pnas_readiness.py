import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
READINESS = ROOT / "scripts/123_build_pnas_readiness_matrix.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase32_readiness", READINESS)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_phase32_readiness_marks_nominal_ad_but_not_shared_ready(tmp_path, monkeypatch):
    module = load_module()
    monkeypatch.chdir(tmp_path)
    table = tmp_path / "results/tables/phase32_crosscohort_axis_evidence_summary.tsv"
    table.parent.mkdir(parents=True)
    table.write_text(
        "axis_id\tcrosscohort_evidence_class\tgse174367_fdr\n"
        "neuronal_vulnerability_axis\tstrong_ad_axis_with_nominal_external_replication\t0.2438\n",
        encoding="utf-8",
    )
    pd = tmp_path / "results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv"
    pd.write_text("axis_id\tevidence_label\nneuronal_vulnerability_axis\tdirectionally_consistent_but_not_significant\n", encoding="utf-8")
    assert module.status_for("independent_ad_replication")[0] == "nominally_supported"
    assert module.status_for("independent_pd_replication")[0] == "available_but_preliminary"
    assert module.status_for("shared_ad_pd_axis_claim")[0] == "not_ready"
    assert module.status_for("pnas_biological_claim")[0] == "promising_but_not_ready"
