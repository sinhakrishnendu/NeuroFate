import csv
import gzip
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARSER = ROOT / "scripts/163_parse_pd_geo_series_matrix_metadata.py"


def test_phase37_gse7621_metadata_parser_emits_geo_ids_and_pd_labels(tmp_path):
    series = tmp_path / "GSE7621_series_matrix.txt.gz"
    output = tmp_path / "metadata.tsv"
    labels = tmp_path / "labels.tsv"
    platforms = tmp_path / "platforms.tsv"
    log = tmp_path / "parse.log"
    with gzip.open(series, "wt", encoding="utf-8", newline="") as handle:
        handle.write(
            "!Sample_title\t\"control SN 1\"\t\"PD SN 1\"\t\"PD SN 2\"\n"
            "!Sample_geo_accession\t\"GSMC1\"\t\"GSMP1\"\t\"GSMP2\"\n"
            "!Sample_source_name_ch1\t\"substantia nigra\"\t\"substantia nigra\"\t\"substantia nigra\"\n"
            "!Sample_platform_id\t\"GPL570\"\t\"GPL570\"\t\"GPL570\"\n"
            "!Sample_characteristics_ch1\t\"disease state: control\"\t\"disease state: Parkinson's disease\"\t\"disease state: PD\"\n"
            "!series_matrix_table_begin\n"
            "\"ID_REF\"\t\"GSMC1\"\t\"GSMP1\"\t\"GSMP2\"\n"
            "!series_matrix_table_end\n"
        )
    subprocess.run(
        [
            sys.executable,
            str(PARSER),
            "--series-matrix",
            str(series),
            "--cohort-id",
            "gse7621_pd_sn_bulk",
            "--output",
            str(output),
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
    rows = list(csv.DictReader(output.open("r", encoding="utf-8"), delimiter="\t"))
    assert [row["sample_id"] for row in rows] == ["GSMC1", "GSMP1", "GSMP2"]
    assert [row["geo_accession"] for row in rows] == ["GSMC1", "GSMP1", "GSMP2"]
    assert [row["label__pd_vs_control"] for row in rows] == ["0", "1", "1"]
    summary = {row["label__pd_vs_control"]: row["count"] for row in csv.DictReader(labels.open("r", encoding="utf-8"), delimiter="\t")}
    assert summary == {"0": "1", "1": "2"}
