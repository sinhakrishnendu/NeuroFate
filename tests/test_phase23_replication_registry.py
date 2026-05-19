from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parents[1]


def read_registry():
    with (ROOT / "metadata/phase23_replication_cohort_registry.tsv").open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_phase23_registry_contains_target_replication_cohorts():
    rows = read_registry()
    accessions = {row["geo_accession"] for row in rows}
    cohort_ids = {row["cohort_id"] for row in rows}
    assert {"GSE184950", "GSE174367", "GSE147528"} <= accessions
    assert {"gse184950_pd_sn", "gse174367_ad_multiomics", "gse147528_ad_progression"} <= cohort_ids


def test_phase23_registry_has_local_dirs_and_endpoint_roles():
    rows = read_registry()
    for row in rows:
        assert row["local_raw_dir"].startswith("data/raw/external/")
        assert row["local_interim_dir"].startswith("data/interim/external/")
        assert row["endpoint_role"]
        assert row["status"] == "manual_acquisition_required"
