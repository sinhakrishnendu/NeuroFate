import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "metadata/phase28_ad_replication_registry.tsv"


def test_phase28_registry_contains_required_ad_cohorts():
    with REGISTRY.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    cohort_ids = {row["cohort_id"] for row in rows}
    accessions = {row["geo_accession"] for row in rows}
    assert {"gse174367_ad_multiomics", "gse147528_ad_progression", "gse157827_ad_snuc_optional"}.issubset(cohort_ids)
    assert {"GSE174367", "GSE147528", "GSE157827"}.issubset(accessions)
    first = next(row for row in rows if row["cohort_id"] == "gse174367_ad_multiomics")
    assert "bulk" in first["first_use_strategy"]


def test_phase28_manual_download_templates_are_guarded():
    for name in ["download_gse174367_manual.sh", "download_gse147528_manual.sh", "download_gse157827_manual.sh"]:
        text = (ROOT / "scripts/manual_downloads" / name).read_text(encoding="utf-8")
        assert "RUN_MANUAL_DOWNLOAD" in text
        assert '!= "YES"' in text
        assert "Do not run from Codex" in text
        assert "# curl" in text or "# wget" in text
