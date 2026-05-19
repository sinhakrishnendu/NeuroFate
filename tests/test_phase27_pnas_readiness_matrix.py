import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/123_build_pnas_readiness_matrix.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase27_readiness", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_readiness_matrix_has_phase27_status_logic():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "phase27_gse184950_axis_replication_statistics_clean.tsv" in text
    assert "available_but_preliminary" in text
    assert "strong_pnas_ready" in text
    assert "satisfied_34_processed_matrices" in text


def test_gse184950_statistical_support_requires_supported_label(tmp_path, monkeypatch):
    module = load_module()
    monkeypatch.chdir(tmp_path)
    path = tmp_path / "results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv"
    path.parent.mkdir(parents=True)
    path.write_text("axis_id\tevidence_label\naxis\tdirectionally_consistent_but_not_significant\n", encoding="utf-8")
    assert module.gse184950_statistical_support() is False
    path.write_text("axis_id\tevidence_label\naxis\treplicated_statistically_supported\n", encoding="utf-8")
    assert module.gse184950_statistical_support() is True
