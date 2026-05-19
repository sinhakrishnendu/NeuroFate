import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/163_parse_pd_geo_series_matrix_metadata.py"


def test_phase34_metadata_parser_extracts_pd_control_endpoint(tmp_path):
    series = tmp_path / "GSE_mock_series_matrix.txt"
    out = tmp_path / "metadata.tsv"
    labels = tmp_path / "labels.tsv"
    platforms = tmp_path / "platforms.tsv"
    log = tmp_path / "parser.log"
    series.write_text(
        "\n".join(
            [
                "!Sample_title\tcontrol SN\tPD SN",
                "!Sample_geo_accession\tGSM1\tGSM2",
                "!Sample_source_name_ch1\tsubstantia nigra\tsubstantia nigra",
                "!Sample_platform_id\tGPL1\tGPL1",
                "!Sample_characteristics_ch1\t\"disease state: Control\"\t\"disease state: Parkinson's disease\"",
                "!Sample_characteristics_ch1\t\"brain region: substantia nigra\"\t\"brain region: substantia nigra\"",
                "!series_matrix_table_begin",
                "ID_REF\tGSM1\tGSM2",
                "!series_matrix_table_end",
                "",
            ]
        ),
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
            "--output",
            str(out),
            "--label-summary-output",
            str(labels),
            "--platform-output",
            str(platforms),
            "--log-file",
            str(log),
        ],
        cwd=ROOT,
        check=True,
    )
    rows = list(csv.DictReader(out.open("r", encoding="utf-8"), delimiter="\t"))
    assert [row["label__pd_vs_control"] for row in rows] == ["0", "1"]
    assert rows[0]["endpoint_status"] == "unambiguous"
    summary = {row["label__pd_vs_control"]: row["count"] for row in csv.DictReader(labels.open("r", encoding="utf-8"), delimiter="\t")}
    assert summary == {"0": "1", "1": "1"}
