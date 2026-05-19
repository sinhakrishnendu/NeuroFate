import importlib.util
import io
import tarfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/135_inspect_gse184950_nested_archives.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase26_nested", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def add_bytes(tar: tarfile.TarFile, name: str, data: bytes) -> None:
    info = tarfile.TarInfo(name)
    info.size = len(data)
    tar.addfile(info, io.BytesIO(data))


def test_nested_archive_inspector_lists_without_extracting(tmp_path):
    module = load_module()
    inner_bytes = io.BytesIO()
    with tarfile.open(fileobj=inner_bytes, mode="w:gz") as inner:
        add_bytes(inner, "filtered_feature_bc_matrix/matrix.mtx.gz", b"matrix")
        add_bytes(inner, "filtered_feature_bc_matrix/features.tsv.gz", b"features")
        add_bytes(inner, "filtered_feature_bc_matrix/barcodes.tsv.gz", b"barcodes")
    raw_tar = tmp_path / "GSE184950_RAW.tar"
    with tarfile.open(raw_tar, "w") as outer:
        add_bytes(outer, "GSM5602315_A10.tar.gz", inner_bytes.getvalue())
    rows = module.inspect_nested(raw_tar, [{"processed_tar_name": "GSM5602315_A10.tar.gz"}])
    assert len(rows) == 3
    assert {row["likely_role"] for row in rows} == {"tenx_matrix", "tenx_features", "tenx_barcodes"}
    assert all(row["sample_id"] == "A10" for row in rows)
    assert all(row["complete_processed_matrix_set"] == "true" for row in rows)


def test_nested_archive_inspector_source_does_not_extract_files():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "extractall" not in text
    assert ".extract(" not in text
    assert "extractfile" in text
