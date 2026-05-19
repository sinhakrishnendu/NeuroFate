import importlib.util
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/136_extract_gse184950_processed_matrices_selective.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase26_selective", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_selective_extractor_requires_manual_environment_guard(tmp_path, monkeypatch):
    module = load_module()
    monkeypatch.delenv("RUN_MANUAL_GSE184950_EXTRACTION", raising=False)
    try:
        module.extract_selected(tmp_path / "missing.tar", [{"outer_member_path": "x", "nested_member_path": "matrix.mtx.gz"}], tmp_path / "out", execute=True)
    except SystemExit as exc:
        assert "RUN_MANUAL_GSE184950_EXTRACTION=YES" in str(exc)
    else:
        raise AssertionError("manual extraction guard did not stop execution")
    assert os.environ.get("RUN_MANUAL_GSE184950_EXTRACTION") is None


def test_selective_extractor_skips_fastq_and_unrelated_files():
    module = load_module()
    rows = [
        {"nested_member_path": "filtered_feature_bc_matrix/matrix.mtx.gz"},
        {"nested_member_path": "filtered_feature_bc_matrix/features.tsv.gz"},
        {"nested_member_path": "filtered_feature_bc_matrix/barcodes.tsv.gz"},
        {"nested_member_path": "sample_R1.fastq.gz"},
        {"nested_member_path": "metadata.txt"},
    ]
    selected = module.selected_inventory_rows(rows)
    assert [Path(row["nested_member_path"]).name for row in selected] == ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]


def test_selective_extractor_has_no_broad_archive_extraction():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "extractall" not in text
    assert ".extract(" not in text
    assert "ALLOWED_NAMES" in text
