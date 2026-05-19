import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "scripts/171_audit_gse7621_axis_direction_and_probe_mapping.py"


def test_phase38_direction_probe_audit_flags_synuclein_opposite_and_neuronal_direction_only(tmp_path):
    axis_scores = tmp_path / "scores.tsv"
    probe_map = tmp_path / "probe_map.tsv"
    coverage = tmp_path / "coverage.tsv"
    stats = tmp_path / "stats.tsv"
    output = tmp_path / "audit.tsv"
    preview = tmp_path / "audit.md"
    log = tmp_path / "audit.log"
    axis_scores.write_text(
        "sample_id\tlabel__pd_vs_control\taxis__synuclein_mitochondrial_axis\taxis__neuronal_vulnerability_axis\n"
        "C1\t0\t1.0\t1.0\n"
        "C2\t0\t0.8\t0.8\n"
        "P1\t1\t-1.0\t0.4\n"
        "P2\t1\t-0.8\t0.2\n",
        encoding="utf-8",
    )
    probe_map.write_text(
        "platform_id\tprobe_id\tgene_symbol\n"
        "GPL570\tp1\tSNCA\nGPL570\tp2\tSNCA\nGPL570\tp3\tMAPT\nGPL570\tp4\tNEFL\n",
        encoding="utf-8",
    )
    coverage.write_text(
        "cohort_id\taxis_id\tgenes_requested\tgenes_found\tgenes_missing\tfound_gene_members\tmissing_gene_members\tstatus\n"
        "mock\tsynuclein_mitochondrial_axis\t6\t5\t1\tAPOE;LRRK2;MAPT;PINK1;SNCA\tPRKN\tok\n"
        "mock\tneuronal_vulnerability_axis\t6\t1\t0\tNEFL\t\tok\n",
        encoding="utf-8",
    )
    stats.write_text(
        "cohort_id\taxis_id\teffect_size\tpvalue\tfdr\tevidence_label\n"
        "mock\tsynuclein_mitochondrial_axis\t-0.8\t0.001\t0.01\topposite_direction\n"
        "mock\tneuronal_vulnerability_axis\t-0.2\t0.4\t0.7\tdirectionally_consistent_but_not_significant\n",
        encoding="utf-8",
    )
    subprocess.run(
        [
            sys.executable,
            str(AUDIT),
            "--axis-scores",
            str(axis_scores),
            "--probe-map",
            str(probe_map),
            "--axis-feature-coverage",
            str(coverage),
            "--replication-stats",
            str(stats),
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
    rows = {row["axis_id"]: row for row in csv.DictReader(output.open("r", encoding="utf-8"), delimiter="\t")}
    assert rows["synuclein_mitochondrial_axis"]["phase38_direction_flag"] == "statistically_significant_opposite_direction"
    assert rows["synuclein_mitochondrial_axis"]["prkn_missing"] == "true"
    assert "SNCA:2" in rows["synuclein_mitochondrial_axis"]["focus_axis_probe_counts"]
    assert rows["neuronal_vulnerability_axis"]["phase38_direction_flag"] == "directionally_consistent_not_significant"
    assert "not shared AD/PD replication" in rows["synuclein_mitochondrial_axis"]["safe_interpretation"]
