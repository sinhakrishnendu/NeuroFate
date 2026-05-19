import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TESTER = ROOT / "scripts/147_test_ad_replication_axis_associations.py"
READINESS = ROOT / "scripts/123_build_pnas_readiness_matrix.py"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_ad_replication_direction_only_is_not_supported():
    module = load(TESTER, "phase28_tester")
    label, consistency = module.label_for(
        effect=0.2,
        pvalue=0.4,
        fdr=0.8,
        phase22="positive",
        n=40,
        positive_n=20,
        negative_n=20,
    )
    assert consistency == "consistent"
    assert label == "directionally_consistent_but_not_significant"


def test_readiness_requires_statistically_supported_ad_replication(tmp_path, monkeypatch):
    module = load(READINESS, "phase28_readiness")
    monkeypatch.chdir(tmp_path)
    table = tmp_path / "results/tables/phase28_gse174367_ad_multiomics_axis_replication_statistics.tsv"
    table.parent.mkdir(parents=True)
    table.write_text("axis_id\tevidence_label\naxis\tdirectionally_consistent_but_not_significant\n", encoding="utf-8")
    assert module.phase28_ad_statistical_support() is False
    assert module.phase28_ad_replication_available() is True
    table.write_text("axis_id\tevidence_label\naxis\tstatistically_supported_ad_replication\n", encoding="utf-8")
    assert module.phase28_ad_statistical_support() is True
