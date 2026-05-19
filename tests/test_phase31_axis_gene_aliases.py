import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ALIAS = ROOT / "metadata/neurofate_axis_gene_aliases.tsv"


def read_aliases():
    with ALIAS.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_alias_table_has_required_schema_and_conservative_types():
    rows = read_aliases()
    assert rows
    assert set(rows[0]) == {"gene_symbol", "ensembl_gene_id", "alias_type", "source_note"}
    alias_types = {row["alias_type"] for row in rows}
    assert "ensembl_gene_id" in alias_types


def test_alias_table_contains_required_neurofate_ensembl_mappings():
    rows = read_aliases()
    by_gene = {row["gene_symbol"]: row for row in rows}
    required = {
        "SNCA": "ENSG00000145335",
        "MAPT": "ENSG00000186868",
        "APP": "ENSG00000142192",
        "APOE": "ENSG00000130203",
        "TREM2": "ENSG00000095970",
        "PRKN": "ENSG00000185345",
    }
    for gene, ensembl in required.items():
        assert by_gene[gene]["ensembl_gene_id"] == ensembl
        assert by_gene[gene]["source_note"] == "curated_human_ensembl_symbol_mapping_for_neurofate_axis_panel"


def test_alias_table_has_no_blank_ensembl_ids():
    rows = read_aliases()
    assert all(row["ensembl_gene_id"].startswith("ENSG") for row in rows)
    assert len(rows) >= 30
