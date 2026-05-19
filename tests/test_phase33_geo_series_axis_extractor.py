from pathlib import Path
import csv
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/163_build_axis_scores_from_geo_series_matrix.py"


def write_mock_series(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "!Sample_title\tControl_1\tPD_1",
                "!Sample_geo_accession\tGSM1\tGSM2",
                "!Sample_characteristics_ch1\t\"disease state: Control\"\t\"disease state: Parkinson's disease\"",
                "!series_matrix_table_begin",
                "ID_REF\tGSM1\tGSM2",
                "SNCA\t1\t3",
                "NEFL\t2\t4",
                "RANDOMGENE\t100\t200",
                "!series_matrix_table_end",
                "",
            ]
        ),
        encoding="utf-8",
    )


def test_series_matrix_extractor_writes_axis_scores_only(tmp_path):
    series = tmp_path / "mock_series_matrix.txt"
    axis_registry = tmp_path / "axis.tsv"
    output = tmp_path / "axis_scores.tsv"
    coverage = tmp_path / "coverage.tsv"
    labels = tmp_path / "labels.tsv"
    log_file = tmp_path / "run.log"
    write_mock_series(series)
    axis_registry.write_text(
        "axis_id\tgene_members\nneuronal_vulnerability_axis\tSNCA;NEFL\nsynuclein_mitochondrial_axis\tSNCA\n",
        encoding="utf-8",
    )
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--series-matrix",
            str(series),
            "--cohort-id",
            "mock_pd",
            "--axis-registry",
            str(axis_registry),
            "--output",
            str(output),
            "--coverage-output",
            str(coverage),
            "--label-summary-output",
            str(labels),
            "--log-file",
            str(log_file),
        ],
        cwd=ROOT,
        check=True,
    )
    rows = list(csv.DictReader(output.open("r", encoding="utf-8"), delimiter="\t"))
    assert [row["sample_id"] for row in rows] == ["GSM1", "GSM2"]
    assert {row["label__pd_vs_control"] for row in rows} == {"0", "1"}
    assert "RANDOMGENE" not in output.read_text(encoding="utf-8")
    assert "axis__neuronal_vulnerability_axis" in rows[0]


def test_series_matrix_extractor_source_avoids_single_cell_and_dense_tools():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["scanpy", "anndata", "read_h5ad", "umap", "leiden", "louvain", "toarray(", "todense("]:
        assert forbidden not in text
