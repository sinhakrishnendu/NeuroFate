from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/124_parse_gse184950_geo_metadata_workbook.py"


def test_metadata_parser_detects_samples_section_and_expected_columns():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "METADATA TEMPLATE" in text
    assert "find_sample_header" in text
    assert '"samples"' in text
    for column in [
        "sample_name",
        "disease_state",
        "donor_id",
        "pmi_hours",
        "braak_stage",
        "processed_data_file",
        "raw_files",
    ]:
        assert column in text


def test_metadata_parser_does_not_open_archives_or_expression():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["tarfile", "matrix.mtx", "fastq", "scanpy", "read_h5ad", "anndata"]:
        assert forbidden not in text
