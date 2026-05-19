from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "metadata/phase15_external_validation_candidates.tsv"


def read_registry() -> list[dict[str, str]]:
    with REGISTRY.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_phase15_candidate_registry_contains_priority_datasets() -> None:
    rows = read_registry()
    dataset_ids = {row["dataset_id"] for row in rows}
    assert "gse243639_pd_snpc" in dataset_ids
    assert "gse174367_ad_multiomics" in dataset_ids
    assert "gse147528_ad_snrna" in dataset_ids
    assert "rosmap_ad_bulk_rnaseq" in dataset_ids


def test_phase15_candidate_registry_required_columns() -> None:
    rows = read_registry()
    assert rows
    required = {
        "dataset_id",
        "disease",
        "priority",
        "modality",
        "accession_or_portal",
        "access_type",
        "expected_file_formats",
        "metadata_required",
        "local_raw_dir",
        "local_interim_dir",
        "validation_role",
    }
    assert required.issubset(rows[0])


def test_phase15_manual_download_templates_are_guarded() -> None:
    scripts = [
        "download_gse243639_manual.sh",
        "download_gse174367_manual.sh",
        "download_gse147528_manual.sh",
    ]
    for script_name in scripts:
        text = (ROOT / "scripts/manual_downloads" / script_name).read_text(encoding="utf-8")
        assert "RUN_MANUAL_DOWNLOAD" in text
        assert "Do not run from Codex. Manual user execution only." in text
        assert "MANUAL_HEAVY" in text
        assert "# wget" in text


def test_rosmap_access_template_is_manual_only() -> None:
    text = (ROOT / "scripts/manual_downloads/prepare_rosmap_access_manual.md").read_text(encoding="utf-8")
    assert "Manual" in text or "manual" in text
    assert "Do not run from Codex" in text
    assert "controlled" in text.lower()
