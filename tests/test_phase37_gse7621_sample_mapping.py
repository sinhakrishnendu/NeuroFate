import csv
import gzip
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "scripts/169_audit_gse7621_sample_mapping.py"


def test_phase37_sample_mapping_reads_header_only_and_selects_geo_accession(tmp_path):
    series = tmp_path / "GSE7621_series_matrix.txt.gz"
    metadata = tmp_path / "metadata.tsv"
    output = tmp_path / "audit.tsv"
    preview = tmp_path / "preview.txt"
    log = tmp_path / "audit.log"
    with gzip.open(series, "wt", encoding="utf-8", newline="") as handle:
        handle.write(
            "!Series_title\t\"mock\"\n"
            "!series_matrix_table_begin\n"
            "\"ID_REF\"\t\"GSM1\"\t\"GSM2\"\t\"GSM3\"\n"
            "\"probe_should_not_be_read\"\t1\t2\t3\n"
            "!series_matrix_table_end\n"
        )
    metadata.write_text(
        "cohort_id\tsample_id\tgeo_accession\tsample_title\tsource_name\tlabel__pd_vs_control\n"
        "gse7621_pd_sn_bulk\tS1\tGSM1\tControl SN\tSN\t0\n"
        "gse7621_pd_sn_bulk\tS2\tGSM2\tPD SN\tSN\t1\n"
        "gse7621_pd_sn_bulk\tS3\tGSM3\tPD SN2\tSN\t1\n",
        encoding="utf-8",
    )
    subprocess.run(
        [
            sys.executable,
            str(AUDIT),
            "--series-matrix",
            str(series),
            "--metadata",
            str(metadata),
            "--output",
            str(output),
            "--preview-output",
            str(preview),
            "--log-file",
            str(log),
        ],
        cwd=ROOT,
        check=True,
    )
    rows = list(csv.DictReader(output.open("r", encoding="utf-8"), delimiter="\t"))
    best = next(row for row in rows if row["is_best_join_key"] == "true")
    assert best["candidate_join_key"] == "geo_accession"
    assert best["matched_sample_count"] == "3"
    assert best["pd_count"] == "2"
    assert best["control_count"] == "1"
    assert "probe_should_not_be_read" not in preview.read_text(encoding="utf-8")
