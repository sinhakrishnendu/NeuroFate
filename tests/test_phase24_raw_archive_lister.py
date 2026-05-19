from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/125_list_gse184950_raw_archive.py"
DOWNLOAD = ROOT / "scripts/manual_downloads/download_gse184950_raw_manual.sh"


def test_raw_archive_lister_lists_without_extracting():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "tarfile.open" in text
    assert "getmembers()" in text
    assert ".extract" not in text
    assert "extractall" not in text
    assert "Archive was listed only; no members were extracted." in text


def test_raw_download_template_is_guarded():
    text = DOWNLOAD.read_text(encoding="utf-8")
    assert "RUN_MANUAL_DOWNLOAD" in text
    assert 'RUN_MANUAL_DOWNLOAD}" != "YES"' in text
    assert "Do not run from Codex. Manual user execution only." in text
    assert "# curl" in text
    assert "# wget" in text
    assert "GSE184950_RAW.tar" in text
