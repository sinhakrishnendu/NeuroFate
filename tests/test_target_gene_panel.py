import csv
from pathlib import Path


EXPECTED_GENES = {
    "SNCA",
    "MAPT",
    "APP",
    "PSEN1",
    "PSEN2",
    "APOE",
    "TREM2",
    "TYROBP",
    "GFAP",
    "AIF1",
    "P2RY12",
    "GPNMB",
    "LAMP5",
    "SST",
    "PVALB",
    "SLC17A7",
    "MBP",
    "MOBP",
    "PLP1",
    "NEFL",
    "NEFM",
    "B2M",
    "HLA-DRA",
    "CX3CR1",
    "IL1B",
    "TNF",
    "NFKB1",
    "LRRK2",
    "PINK1",
    "PRKN",
}

REQUIRED_COLUMNS = {
    "gene_symbol",
    "biological_role",
    "manuscript_relevance",
    "priority_tier",
    "expected_cell_types",
    "notes",
}


def read_panel():
    with Path("metadata/target_gene_panel_v1.tsv").open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_target_gene_panel_schema_and_size():
    rows = read_panel()
    assert rows
    assert REQUIRED_COLUMNS.issubset(rows[0].keys())
    assert len(rows) == 30


def test_target_gene_panel_expected_genes():
    rows = read_panel()
    assert {row["gene_symbol"] for row in rows} == EXPECTED_GENES


def test_target_gene_panel_has_priority_tiers_and_roles():
    rows = read_panel()
    assert all(row["priority_tier"].startswith("tier") for row in rows)
    assert all(row["biological_role"] for row in rows)
    assert all(row["manuscript_relevance"] for row in rows)
