from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_phase24_scripts_avoid_scanpy_h5ad_anndata_umap_clustering():
    for script_name in [
        "124_parse_gse184950_geo_metadata_workbook.py",
        "125_list_gse184950_raw_archive.py",
        "126_plan_gse184950_processed_matrix_extraction.py",
        "127_extract_gse184950_axis_genes_from_10x.py",
        "128_build_gse184950_axis_scores.py",
        "129_test_gse184950_axis_replication.py",
    ]:
        text = (ROOT / "scripts" / script_name).read_text(encoding="utf-8").lower()
        for forbidden in ["scanpy", "read_h5ad", "anndata", ".h5ad", "umap", "leiden", "cluster"]:
            assert forbidden not in text


def test_extractor_is_manual_guarded_and_no_dense_conversion():
    text = (ROOT / "scripts/127_extract_gse184950_axis_genes_from_10x.py").read_text(encoding="utf-8").lower()
    assert "--run-manual-extraction" in text
    assert 'choices=["yes", "no"]' in text
    assert "toarray" not in text
    assert "todense" not in text
    assert "np.array" not in text
    assert "fastq" not in text
