import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts/165_build_pd_axis_scores_from_geo_expression.py"


def test_phase35_builder_keeps_18_gse20141_samples_and_axis_only_output(tmp_path):
    gsm = [f"GSM5039{50 + idx:02d}" for idx in range(18)]
    expression = tmp_path / "expr.tsv"
    metadata = tmp_path / "metadata.tsv"
    probe_map = tmp_path / "probe_map.tsv"
    axis = tmp_path / "axis.tsv"
    output = tmp_path / "axis_scores.tsv"
    coverage = tmp_path / "coverage.tsv"
    labels = tmp_path / "labels.tsv"
    join = tmp_path / "join.tsv"
    log = tmp_path / "axis.log"
    expression.write_text(
        "ID_REF\t" + "\t".join(gsm) + "\n"
        "probe_snca\t" + "\t".join(str(idx + 1) for idx in range(18)) + "\n"
        "probe_nefl\t" + "\t".join(str(idx + 2) for idx in range(18)) + "\n"
        "probe_random\t" + "\t".join(["100"] * 18) + "\n",
        encoding="utf-8",
    )
    metadata_rows = ["cohort_id\tsample_id\tgeo_accession\tsample_title\tdisease_state\ttissue_or_region\tlabel__pd_vs_control"]
    for idx, sample in enumerate(gsm):
        label = "0" if idx < 8 else "1"
        disease = "Control" if label == "0" else "Parkinson's disease"
        metadata_rows.append(f"gse20141_pd_snpc_lcm\t{sample}\t{sample}\tSample {idx}\t{disease}\tSNpc\t{label}")
    metadata.write_text("\n".join(metadata_rows) + "\n", encoding="utf-8")
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
            "--expected-sample-count",
            "18",
            "--log-file",
            str(log),
        ],
        cwd=ROOT,
        check=True,
    )
    rows = list(csv.DictReader(output.open("r", encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 18
    summary = {row["label__pd_vs_control"]: row["count"] for row in csv.DictReader(labels.open("r", encoding="utf-8"), delimiter="\t")}
    assert summary == {"0": "8", "1": "10"}
    assert "probe_random" not in output.read_text(encoding="utf-8")
    assert "axis__neuronal_vulnerability_axis" in rows[0]


def test_phase35_builder_excludes_ambiguous_metadata_labels(tmp_path):
    expression = tmp_path / "expr.tsv"
    metadata = tmp_path / "metadata.tsv"
    probe_map = tmp_path / "probe_map.tsv"
    axis = tmp_path / "axis.tsv"
    output = tmp_path / "axis_scores.tsv"
    coverage = tmp_path / "coverage.tsv"
    labels = tmp_path / "labels.tsv"
    join = tmp_path / "join.tsv"
    log = tmp_path / "axis.log"
    expression.write_text("ID_REF\tGSM1\tGSM2\tGSM3\nprobe_snca\t1\t3\t9\nprobe_nefl\t2\t4\t9\n", encoding="utf-8")
    metadata.write_text(
        "cohort_id\tsample_id\tgeo_accession\tsample_title\tdisease_state\ttissue_or_region\tlabel__pd_vs_control\nmock\tGSM1\tGSM1\tC\tControl\tSN\t0\nmock\tGSM2\tGSM2\tPD\tParkinson's disease\tSN\t1\nmock\tGSM3\tGSM3\tUnknown\tUnknown\tSN\t\n",
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
