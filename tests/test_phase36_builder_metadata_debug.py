import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts/165_build_pd_axis_scores_from_geo_expression.py"


def test_metadata_debug_records_label_and_join_decisions(tmp_path):
    expression = tmp_path / "expr.tsv"
    metadata = tmp_path / "metadata.tsv"
    probe_map = tmp_path / "probe_map.tsv"
    axis = tmp_path / "axis.tsv"
    output = tmp_path / "axis_scores.tsv"
    coverage = tmp_path / "coverage.tsv"
    labels = tmp_path / "labels.tsv"
    join = tmp_path / "join.tsv"
    debug = tmp_path / "debug.tsv"
    log = tmp_path / "axis.log"

    expression.write_text("ID_REF\tGSM1\tGSM2\tGSM3\nprobe_snca\t1\t3\t9\nprobe_nefl\t2\t4\t9\n", encoding="utf-8")
    metadata.write_text(
        "sample_id\tgeo_accession\tsample_title\tdisease_state\tendpoint_status\tlabel__pd_vs_control\n"
        "GSM1\tGSM1\tControl sample\tignored\tunambiguous\t0\n"
        "GSM2\tGSM2\tPD sample\tignored\tunambiguous\t1\n"
        "GSM3\tGSM3\tAmbiguous sample\tignored\tambiguous\t1\n",
        encoding="utf-8",
    )
    probe_map.write_text("platform_id\tprobe_id\tgene_symbol\nGPL570\tprobe_snca\tSNCA\nGPL570\tprobe_nefl\tNEFL\n", encoding="utf-8")
    axis.write_text("axis_id\tgene_members\nneuronal_vulnerability_axis\tSNCA;NEFL\n", encoding="utf-8")

    subprocess.run(
        [
            sys.executable,
            str(BUILDER),
            "--expression-file",
            str(expression),
            "--sample-metadata",
            str(metadata),
            "--probe-map",
            str(probe_map),
            "--cohort-id",
            "mock",
            "--axis-registry",
            str(axis),
            "--output",
            str(output),
            "--coverage-output",
            str(coverage),
            "--label-summary-output",
            str(labels),
            "--join-output",
            str(join),
            "--metadata-debug-output",
            str(debug),
            "--log-file",
            str(log),
        ],
        cwd=ROOT,
        check=True,
    )

    debug_row = next(csv.DictReader(debug.open("r", encoding="utf-8"), delimiter="\t"))
    assert debug_row["selected_label_column"] == "label__pd_vs_control"
    assert debug_row["selected_endpoint_status_column"] == "endpoint_status"
    assert debug_row["selected_join_key"] == "geo_accession"
    assert debug_row["metadata_rows_before_filtering"] == "3"
    assert debug_row["metadata_rows_after_endpoint_status_filtering"] == "2"
    assert debug_row["expression_sample_count"] == "3"
    assert debug_row["matched_sample_count"] == "2"
    assert debug_row["final_labeled_matched_sample_count"] == "2"
    assert "GSM3" in debug_row["unmatched_expression_sample_examples"]
    assert "endpoint_status=ambiguous" in debug_row["filtered_out_metadata_examples"]
