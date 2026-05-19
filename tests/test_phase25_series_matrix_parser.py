import gzip
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/131_parse_gse184950_series_matrix.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase25_series", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def make_mock_series(path: Path) -> None:
    titles = [f"S{i:02d}" for i in range(1, 35)]
    accessions = [f"GSM{i:07d}" for i in range(5602315, 5602315 + 34)]
    diseases = ["Unaffected Control"] * 10 + ["Parkinson's Disease"] * 6 + ["Parkinson's Disease Dementia"] * 18
    links = [f"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM/mock/suppl/{title}.tar.gz" for title in titles]
    rows = [
        ["!Sample_title", *titles],
        ["!Sample_geo_accession", *accessions],
        ["!Sample_source_name_ch1", *(["substantia nigra"] * 34)],
        ["!Sample_characteristics_ch1", *(["tissue: substantia nigra"] * 34)],
        ["!Sample_characteristics_ch1", *[f"disease state: {disease}" for disease in diseases]],
        ["!Sample_characteristics_ch1", *[f"brain bank donor id: donor_{i:02d}" for i in range(1, 35)]],
        ["!Sample_characteristics_ch1", *(["age: 80"] * 34)],
        ["!Sample_characteristics_ch1", *(["gender: female"] * 34)],
        ["!Sample_characteristics_ch1", *(["race: reported"] * 34)],
        ["!Sample_characteristics_ch1", *(["ethnicity: reported"] * 34)],
        ["!Sample_characteristics_ch1", *(["postmortem interval hours: 10"] * 34)],
        ["!Sample_characteristics_ch1", *(["Braak stage: 3"] * 34)],
        ["!Sample_supplementary_file_1", *links],
    ]
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for row in rows:
            handle.write("\t".join(f'"{value}"' for value in row) + "\n")


def test_series_matrix_parser_extracts_34_records(tmp_path):
    module = load_module()
    series = tmp_path / "GSE184950_series_matrix.txt.gz"
    make_mock_series(series)
    metadata, manifest, labels = module.build_metadata(series)
    assert len(metadata) == 34
    assert len(manifest) == 34
    assert sum(row["label__pd_pdd_vs_control"] == "1" for row in metadata) == 24
    assert sum(row["label__pd_pdd_vs_control"] == "0" for row in metadata) == 10
    assert any(row["category"] == "Parkinson's Disease Dementia" and row["count"] == "18" for row in labels)
    assert metadata[0]["processed_tar_name"] == "S01.tar.gz"


def test_series_matrix_parser_is_metadata_only():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["scanpy", "anndata", "read_h5ad", "tarfile.open", "extractall", "fastq", "umap", "leiden"]:
        assert forbidden not in text
