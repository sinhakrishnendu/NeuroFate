from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parents[1]


def read_registry():
    with (ROOT / "metadata/neurofate_axis_endpoint_registry.tsv").open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_endpoint_registry_contains_required_primary_endpoints():
    rows = read_registry()
    ids = {row["endpoint_id"] for row in rows}
    assert "sea_ad_cognitive_dementia" in ids
    assert "gse243639_pd_diagnosis" in ids
    sea = next(row for row in rows if row["endpoint_id"] == "sea_ad_cognitive_dementia")
    pd = next(row for row in rows if row["endpoint_id"] == "gse243639_pd_diagnosis")
    assert sea["source_column"] == "label__Cognitive_Status"
    assert sea["positive_class"] == "Dementia"
    assert sea["negative_class"] == "No dementia"
    assert pd["source_column"] == "diagnosis"
    assert pd["positive_class"] == "Parkinson's"
    assert pd["negative_class"] == "Control"
    assert sea["allowed_for_cross_disease_comparison"] == "true"
    assert pd["allowed_for_cross_disease_comparison"] == "true"


def test_endpoint_registry_has_secondary_pathology_endpoints():
    ids = {row["endpoint_id"] for row in read_registry()}
    assert "sea_ad_ad_pathology_ordinal" in ids
    assert "sea_ad_braak_ordinal" in ids
    assert "sea_ad_cerad_ordinal" in ids
    assert "gse243639_lewy_neocortical" in ids
