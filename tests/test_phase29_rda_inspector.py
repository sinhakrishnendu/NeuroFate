import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/149_inspect_gse174367_bulk_rda.py"
R_HELPER = ROOT / "scripts/r_inspect_gse174367_bulk_rda.R"


def load_module():
    spec = importlib.util.spec_from_file_location("phase29_rda_inspector", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_rda_inspector_writes_missing_runtime_report(tmp_path):
    module = load_module()
    output = tmp_path / "structure.tsv"
    preview = tmp_path / "preview.txt"
    module.write_missing_runtime_outputs(output, preview, "test runtime missing")
    assert "missing_runtime" in output.read_text(encoding="utf-8")
    assert "No conversion or analysis was attempted" in preview.read_text(encoding="utf-8")


def test_r_helper_exists_and_is_safe():
    text = R_HELPER.read_text(encoding="utf-8")
    assert "new.env" in text
    assert "load(" in text
    assert "write.table" in text
    assert "source(" not in text


def test_rda_inspector_source_uses_supported_runtimes_only():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "Rscript" in text
    assert "pyreadr" in text
    assert "write_missing_runtime_outputs" in text
