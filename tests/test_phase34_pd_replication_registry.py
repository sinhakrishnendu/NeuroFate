from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "metadata/phase34_pd_replication_registry.tsv"
MANUAL = [
    ROOT / "scripts/manual_downloads/download_gse20141_manual.sh",
    ROOT / "scripts/manual_downloads/download_gse7621_manual.sh",
    ROOT / "scripts/manual_downloads/download_gse8397_manual.sh",
    ROOT / "scripts/manual_downloads/download_gse20186_manual.sh",
]


def test_phase34_registry_contains_required_pd_cohorts():
    rows = list(csv.DictReader(REGISTRY.open("r", encoding="utf-8"), delimiter="\t"))
    accessions = {row["geo_accession"] for row in rows}
    assert {"GSE20141", "GSE7621", "GSE8397", "GSE20186"}.issubset(accessions)
    assert any(row["cohort_id"] == "gse20141_pd_snpc_lcm" and row["priority"] == "1" for row in rows)


def test_phase34_manual_download_templates_are_guarded():
    for path in MANUAL:
        text = path.read_text(encoding="utf-8")
        assert "RUN_MANUAL_DOWNLOAD" in text
        assert "Do not run from Codex. Manual user execution only." in text
        assert "series_matrix" in text
        assert "soft" in text.lower()
        assert "miniml" in text.lower()
