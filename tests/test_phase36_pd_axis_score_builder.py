import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts/165_build_pd_axis_scores_from_geo_expression.py"


def test_phase36_builder_writes_sample_level_axis_scores_only(tmp_path):
    expression = tmp_path / "expr.tsv"
    metadata = tmp_path / "metadata.tsv"
    probe_map = tmp_path / "probe_map.tsv"
    axis = tmp_path / "axis.tsv"
    output = tmp_path / "phase36_gse20141_pd_snpc_lcm_axis_scores.tsv"
    coverage = tmp_path / "phase36_gse20141_pd_snpc_lcm_axis_feature_coverage.tsv"
    labels = tmp_path / "phase36_gse20141_pd_snpc_lcm_axis_label_summary.tsv"
    join = tmp_path / "phase36_gse20141_expression_metadata_join.tsv"
    debug = tmp_path / "phase36_gse20141_builder_metadata_debug.tsv"
    log = tmp_path / "axis.log"
    gsm = [f"GSM5039{50 + idx:02d}" for idx in range(18)]
    expression.write_text(
        "ID_REF\t" + "\t".join(gsm) + "\n"
        "probe_snca\t" + "\t".join(str(idx + 1) for idx in range(18)) + "\n"
        "probe_nefl\t" + "\t".join(str(idx + 2) for idx in range(18)) + "\n"
        "probe_not_neurofate\t" + "\t".join(["999"] * 18) + "\n",
        encoding="utf-8",
    )
    rows = [
        "cohort_id\tsample_id\tgeo_accession\tsample_title\tdisease_state\ttissue_or_region\tendpoint_status\tlabel__pd_vs_control"
    ]
    for idx, sample in enumerate(gsm):
        label = "0" if idx < 8 else "1"
        rows.append(f"gse20141_pd_snpc_lcm\t{sample}\t{sample}\tSample {idx}\t\tSNpc\tunambiguous\t{label}")
    metadata.write_text("\n".join(rows) + "\n", encoding="utf-8")
    probe_map.write_text("platform_id\tprobe_id\tgene_symbol\nGPL570\tprobe_snca\tSNCA\nGPL570\tprobe_nefl\tNEFL\n", encoding="utf-8")
    axis.write_text(
        "axis_id\tgene_members\n"
        "neuronal_vulnerability_axis\tSNCA;NEFL\n"
        "synuclein_mitochondrial_axis\tSNCA\n",
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
            "gse20141_pd_snpc_lcm",
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
            "--expected-sample-count",
            "18",
            "--log-file",
            str(log),
        ],
        cwd=ROOT,
        check=True,
    )

    axis_rows = list(csv.DictReader(output.open("r", encoding="utf-8"), delimiter="\t"))
    assert len(axis_rows) == 18
    assert set(axis_rows[0]) == {
        "cohort_id",
        "sample_id",
        "expression_sample_id",
        "metadata_join_key",
        "label__pd_vs_control",
        "disease_state",
        "tissue_or_region",
        "axis__neuronal_vulnerability_axis",
        "axis__synuclein_mitochondrial_axis",
    }
    assert "probe_not_neurofate" not in output.read_text(encoding="utf-8")
    label_counts = {row["label__pd_vs_control"]: row["count"] for row in csv.DictReader(labels.open("r", encoding="utf-8"), delimiter="\t")}
    assert label_counts == {"0": "8", "1": "10"}
    assert {row["selected_join_key"] for row in csv.DictReader(join.open("r", encoding="utf-8"), delimiter="\t") if row["join_status"] == "matched"} == {"geo_accession"}


def test_phase36_auto_orientation_skips_geo_series_preamble(tmp_path):
    expression = tmp_path / "series_matrix.txt.gz"
    metadata = tmp_path / "metadata.tsv"
    probe_map = tmp_path / "probe_map.tsv"
    axis = tmp_path / "axis.tsv"
    output = tmp_path / "axis_scores.tsv"
    coverage = tmp_path / "coverage.tsv"
    labels = tmp_path / "labels.tsv"
    join = tmp_path / "join.tsv"
    log = tmp_path / "axis.log"
    import gzip

    with gzip.open(expression, "wt", encoding="utf-8", newline="") as handle:
        handle.write(
            "!Series_title\t\"mock\"\n"
            "!Sample_geo_accession\t\"GSM1\"\t\"GSM2\"\n"
            "!series_matrix_table_begin\n"
            "\"ID_REF\"\t\"GSM1\"\t\"GSM2\"\n"
            "\"probe_snca\"\t1\t3\n"
            "\"probe_nefl\"\t2\t4\n"
            "!series_matrix_table_end\n"
        )
    metadata.write_text(
        "sample_id\tgeo_accession\tsample_title\tendpoint_status\tlabel__pd_vs_control\n"
        "GSM1\tGSM1\tControl\tunambiguous\t0\n"
        "GSM2\tGSM2\tPD\tunambiguous\t1\n",
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

    rows = list(csv.DictReader(output.open("r", encoding="utf-8"), delimiter="\t"))
    assert [row["expression_sample_id"] for row in rows] == ["GSM1", "GSM2"]
