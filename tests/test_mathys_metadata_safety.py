from pathlib import Path


MANUAL_DOWNLOAD = Path("scripts/manual_downloads/download_mathys2019_geo_manual.sh")
METADATA_SCRIPT = Path("scripts/36_inspect_external_h5ad_metadata.py")
OVERLAP_SCRIPT = Path("scripts/37_prepare_mathys_gene_panel_overlap.py")


def test_mathys_manual_download_template_is_guarded():
    text = MANUAL_DOWNLOAD.read_text(encoding="utf-8")
    assert "RUN_MANUAL_DOWNLOAD" in text
    assert "RUN_MANUAL_DOWNLOAD is not YES. Exiting without download." in text
    assert "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138852" in text
    assert "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/" in text
    assert "# wget" in text
    assert "# curl" in text
    assert "md5" in text
    assert "shasum -a 256" in text


def test_mathys_metadata_inspector_has_no_expression_access():
    text = METADATA_SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered
    assert "import anndata" not in lowered
    assert "scvi" not in lowered
    assert "scvelo" not in lowered
    assert '["X"]' not in text
    assert "Forbidden expression-matrix access" in text
    assert "obs/__categories" in text
    assert "mathys2019_metadata_summary.tsv" in text
    assert "mathys_var_genes.tsv" in text


def test_mathys_gene_overlap_is_metadata_only():
    text = OVERLAP_SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "import h5py" not in lowered
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered
    assert "expression_value" not in text
    assert "mathys_gene_overlap.tsv" in text
    assert "mathys_missing_target_genes.tsv" in text
    assert "target_gene_panel_v1.tsv" in text
