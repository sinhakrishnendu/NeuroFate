from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "metadata/phase33_pd_replication_registry.tsv"
MANUAL = [
    ROOT / "scripts/manual_downloads/download_gse20141_manual.sh",
    ROOT / "scripts/manual_downloads/download_gse20186_manual.sh",
    ROOT / "scripts/manual_downloads/download_pd_microarray_replication_manual.sh",
    ROOT / "scripts/manual_downloads/download_gse157783_manual.sh",
]


def test_phase33_registry_contains_priority_pd_targets():
    with REGISTRY.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    cohorts = {row["cohort_id"] for row in rows}
    accessions = {row["geo_accession"] for row in rows}
    assert "gse20141_pd_snpc_lcm" in cohorts
    assert "gse20186_pd_superseries" in cohorts
    assert {"GSE20141", "GSE20186"}.issubset(accessions)
    assert any(row["cohort_id"] == "gse20141_pd_snpc_lcm" and row["priority"] == "1" for row in rows)


def test_phase33_manual_download_scripts_are_guarded():
    for path in MANUAL:
        text = path.read_text(encoding="utf-8")
        assert "RUN_MANUAL_DOWNLOAD" in text
        assert "Do not run from Codex. Manual user execution only." in text
        assert "curl" in text or "wget" in text
