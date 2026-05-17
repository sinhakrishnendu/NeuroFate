from pathlib import Path


SCRIPT_PATHS = [
    Path("scripts/manual_downloads/download_sea_ad_manual.sh"),
    Path("scripts/manual_downloads/download_string_manual.sh"),
    Path("scripts/manual_downloads/download_mathys_synapse_manual.sh"),
    Path("scripts/manual_downloads/download_rosmap_synapse_manual.sh"),
]


def test_manual_download_templates_exist():
    for path in SCRIPT_PATHS:
        assert path.exists()


def test_manual_download_templates_are_guarded():
    for path in SCRIPT_PATHS:
        text = path.read_text(encoding="utf-8")
        lines = text.splitlines()
        assert lines[0] == "set -euo pipefail"
        assert "DO NOT RUN FROM CODEX" in text
        assert "MANUAL_HEAVY" in text
        assert "RUN_MANUAL_DOWNLOAD" in text
        assert 'RUN_MANUAL_DOWNLOAD}" != "YES"' in text
        assert "Refusing to run" in text


def test_manual_download_templates_write_to_allowed_dirs_only():
    for path in SCRIPT_PATHS:
        text = path.read_text(encoding="utf-8")
        assert "data/raw/" in text or "data/external/" in text
        assert "/tmp" not in text
        assert "$HOME" not in text


def test_manual_download_templates_keep_heavy_commands_manual():
    for path in SCRIPT_PATHS:
        text = path.read_text(encoding="utf-8")
        heavy_lines = [line.strip() for line in text.splitlines() if "aws " in line or "curl " in line or "synapse " in line]
        assert heavy_lines
        assert all(line.startswith("#") for line in heavy_lines)
