import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "scripts/157_audit_gse174367_bulk_gene_identifiers.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase31_gene_audit", AUDIT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_gene_id_audit_declares_identifier_rules():
    text = AUDIT.read_text(encoding="utf-8")
    assert "gene_symbol" in text
    assert "ensembl_gene_id" in text
    assert "ensembl_gene_id_versioned" in text
    assert "entrez_numeric" in text
    assert "strip_version_suffix" in text
    assert "ensembl_alias_match" in text
    assert "ensembl_alias_version_stripped" in text


def test_gene_id_audit_axis_gene_parser():
    module = load_module()
    rows = [{"axis_id": "axis", "gene_members": "SNCA; MAPT, APP"}]
    assert module.axis_genes(rows) == ["APP", "MAPT", "SNCA"]
