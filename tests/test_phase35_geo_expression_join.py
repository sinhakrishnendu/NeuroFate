import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts/165_build_pd_axis_scores_from_geo_expression.py"


def test_geo_accession_join_and_sample_title_fallback(tmp_path):
    expression = tmp_path / "expr.tsv"
    metadata = tmp_path / "metadata.tsv"
    probe_map = tmp_path / "probe_map.tsv"
    axis = tmp_path / "axis.tsv"
    output = tmp_path / "axis_scores.tsv"
    coverage = tmp_path / "coverage.tsv"
    labels = tmp_path / "labels.tsv"
    join = tmp_path / "join.tsv"
    log = tmp_path / "axis.log"
    expression.write_text("ID_REF\tGSM1\tGSM2\nprobe_snca\t1\t3\nprobe_nefl\t2\t4\n", encoding="utf-8")
    metadata.write_text(
        "cohort_id\tsample_id\tgeo_accession\tsample_title\tdisease_state\ttissue_or_region\tlabel__pd_vs_control\nmock\tS1\tGSM1\tcontrol sample\tControl\tSN\t0\nmock\tS2\tGSM2\tpd sample\tParkinson's disease\tSN\t1\n",
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
            "--log-file",
            str(log),
        ],
        cwd=ROOT,
        check=True,
    )
    join_rows = list(csv.DictReader(join.open("r", encoding="utf-8"), delimiter="\t"))
    assert {row["selected_join_key"] for row in join_rows if row["join_status"] == "matched"} == {"geo_accession"}

    expression.write_text("ID_REF\tControl Sample\tPD Sample\nprobe_snca\t1\t3\nprobe_nefl\t2\t4\n", encoding="utf-8")
    metadata.write_text(
        "cohort_id\tsample_id\tgeo_accession\tsample_title\tdisease_state\ttissue_or_region\tlabel__pd_vs_control\nmock\tS1\t\tcontrol sample\tControl\tSN\t0\nmock\tS2\t\tpd sample\tParkinson's disease\tSN\t1\n",
        encoding="utf-8",
    )
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
            "--log-file",
            str(log),
        ],
        cwd=ROOT,
        check=True,
    )
    join_rows = list(csv.DictReader(join.open("r", encoding="utf-8"), delimiter="\t"))
    assert {row["selected_join_key"] for row in join_rows if row["join_status"] == "matched"} == {"normalized_sample_title"}
