import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONVERTER = ROOT / "scripts/150_convert_gse174367_bulk_rda_to_axis_matrix.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase31_converter", CONVERTER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_converter_builds_alias_lookup_conservatively():
    module = load_module()
    lookup = module.build_alias_lookup(
        [
            {"gene_symbol": "SNCA", "ensembl_gene_id": "ENSG00000145335", "alias_type": "ensembl_gene_id"},
            {"gene_symbol": "MAPT", "ensembl_gene_id": "", "alias_type": "ensembl_gene_id"},
        ]
    )
    assert lookup == {"ENSG00000145335": "SNCA"}


def test_converter_stops_when_too_few_genes_map():
    module = load_module()
    r_code = module.r_axis_extraction_code(["SNCA", "MAPT"], min_mapped_genes=10)
    assert "mapped_gene_count < min_mapped_genes" in r_code
    assert "stop rather than producing undercovered replication matrix" in r_code
    assert "ensembl_alias_match" in r_code
    assert "ensembl_alias_version_stripped" in r_code


def test_converter_writes_axis_only_matrix_not_full_expression():
    text = CONVERTER.read_text(encoding="utf-8")
    assert "axis_df" in text
    assert "write.table(axis_df" in text
    assert "write.table(expr" not in text
    assert "gene_mapping_output" in text
    assert "phase31_gse174367_bulk_axis_gene_mapping.tsv" in text
