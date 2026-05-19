import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/133_plan_gse184950_selective_tar_extraction.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase25_plan", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_selective_extraction_plan_excludes_fastq_by_default(tmp_path):
    module = load_module()
    axis = tmp_path / "axis.tsv"
    axis.write_text("axis_id\tgene_members\naxis\tSNCA;MAPT\n", encoding="utf-8")
    inventory = [{"member_path": "GSE184950_RAW/A22.tar.gz"}]
    metadata = [{"sample_name": "A22", "processed_tar_name": "A22.tar.gz"}]
    rows = module.build_plan(inventory, metadata, axis)
    assert rows[0]["manual_action"] == "extract_processed_tar_only_then_run_127"
    assert rows[0]["fastq_handling"] == "skip_fastq_no_raw_preprocessing"
    assert "processed_matrices/A22" in rows[0]["output_directory"]


def test_selective_manual_script_is_guarded(tmp_path):
    module = load_module()
    script = tmp_path / "manual.sh"
    module.write_manual_script(script, [{"archive_member_path": "GSE184950_RAW/A22.tar.gz"}])
    text = script.read_text(encoding="utf-8")
    assert "RUN_MANUAL_GSE184950_EXTRACTION" in text
    assert '!= "YES"' in text
    assert "# tar -xf" in text
    assert "Do not run from Codex" in text
