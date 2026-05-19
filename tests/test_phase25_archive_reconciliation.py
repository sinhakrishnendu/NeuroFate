import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/132_reconcile_gse184950_archive_with_series_metadata.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase25_reconcile", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_archive_reconciliation_detects_expected_tar_members():
    module = load_module()
    inventory = [
        {"member_path": "GSE184950_RAW/A22.tar.gz", "likely_role": "per_sample_processed_archive"},
        {"member_path": "GSE184950_RAW/A22/filtered_feature_bc_matrix/matrix.mtx.gz", "likely_role": "tenx_matrix"},
        {"member_path": "GSE184950_RAW/A22/filtered_feature_bc_matrix/barcodes.tsv.gz", "likely_role": "tenx_barcodes"},
        {"member_path": "GSE184950_RAW/A22/filtered_feature_bc_matrix/features.tsv.gz", "likely_role": "tenx_features"},
        {"member_path": "GSE184950_RAW/A22_R1.fastq.gz", "likely_role": "raw_fastq_do_not_process_here"},
    ]
    manifest = [{"sample_name": "A22", "processed_tar_name": "A22.tar.gz", "expected_archive_member": "A22.tar.gz"}]
    rows = module.reconcile(inventory, manifest)
    assert rows[0]["found_in_archive"] == "true"
    assert rows[0]["contains_matrix_mtx"] == "true"
    assert rows[0]["contains_barcodes"] == "true"
    assert rows[0]["contains_features_or_genes"] == "true"
    assert rows[0]["contains_fastq"] == "true"
    assert rows[0]["processed_matrix_availability"] == "processed_10x_members_visible"


def test_archive_reconciliation_does_not_extract_files():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    assert "extractall" not in text
    assert ".extract(" not in text
    assert "tarfile.open" not in text
