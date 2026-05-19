from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = [
    ROOT / "scripts/manual_downloads/download_gse184950_manual.sh",
    ROOT / "scripts/manual_downloads/download_gse174367_manual.sh",
    ROOT / "scripts/manual_downloads/download_gse147528_manual.sh",
]


def test_manual_download_scripts_are_guarded():
    for script in SCRIPTS:
        text = script.read_text(encoding="utf-8")
        assert "RUN_MANUAL_DOWNLOAD" in text
        assert 'RUN_MANUAL_DOWNLOAD}" != "YES"' in text
        assert "Do not run from Codex. Manual user execution only." in text
        assert "MANUAL_HEAVY" in text


def test_manual_download_commands_remain_commented():
    for script in SCRIPTS:
        lines = script.read_text(encoding="utf-8").splitlines()
        risky = [line for line in lines if line.strip().startswith(("wget ", "curl ", "prefetch ", "fasterq-dump "))]
        assert risky == []
        assert any("# wget" in line or "# prefetch" in line for line in lines)
