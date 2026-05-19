import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts/165_build_pd_axis_scores_from_geo_expression.py"


def write_common_inputs(tmp_path: Path, labels: list[str], disease_states: list[str] | None = None) -> dict[str, Path]:
    gsm = [f"GSM5039{50 + idx:02d}" for idx in range(18)]
    paths = {
        "expression": tmp_path / "expr.tsv",
        "metadata": tmp_path / "metadata.tsv",
        "probe_map": tmp_path / "probe_map.tsv",
        "axis": tmp_path / "axis.tsv",
        "output": tmp_path / "axis_scores.tsv",
        "coverage": tmp_path / "coverage.tsv",
        "labels": tmp_path / "labels.tsv",
        "join": tmp_path / "join.tsv",
        "debug": tmp_path / "debug.tsv",
        "log": tmp_path / "axis.log",
    }
    paths["expression"].write_text(
        "ID_REF\t" + "\t".join(gsm) + "\n"
        "probe_snca\t" + "\t".join(str(idx + 1) for idx in range(18)) + "\n"
        "probe_nefl\t" + "\t".join(str(idx + 2) for idx in range(18)) + "\n"
        "probe_random\t" + "\t".join(["100"] * 18) + "\n",
        encoding="utf-8",
    )
    disease_states = disease_states or ["not_used"] * 18
    rows = [
        "cohort_id\tsample_id\tgeo_accession\tsample_title\tdisease_state\ttissue_or_region\tendpoint_status\tlabel__pd_vs_control"
    ]
    for idx, sample in enumerate(gsm):
        rows.append(
            f"gse20141_pd_snpc_lcm\t{sample}\t{sample}\tSample {idx}\t{disease_states[idx]}\tSNpc\tunambiguous\t{labels[idx]}"
        )
    paths["metadata"].write_text("\n".join(rows) + "\n", encoding="utf-8")
    paths["probe_map"].write_text(
        "platform_id\tprobe_id\tgene_symbol\nGPL570\tprobe_snca\tSNCA\nGPL570\tprobe_nefl\tNEFL\n",
        encoding="utf-8",
    )
    paths["axis"].write_text("axis_id\tgene_members\nneuronal_vulnerability_axis\tSNCA;NEFL\n", encoding="utf-8")
    return paths


def run_builder(paths: dict[str, Path]) -> None:
    subprocess.run(
        [
            sys.executable,
            str(BUILDER),
            "--expression-file",
            str(paths["expression"]),
            "--sample-metadata",
            str(paths["metadata"]),
            "--probe-map",
            str(paths["probe_map"]),
            "--cohort-id",
            "gse20141_pd_snpc_lcm",
            "--axis-registry",
            str(paths["axis"]),
            "--output",
            str(paths["output"]),
            "--coverage-output",
            str(paths["coverage"]),
            "--label-summary-output",
            str(paths["labels"]),
            "--join-output",
            str(paths["join"]),
            "--metadata-debug-output",
            str(paths["debug"]),
            "--expected-sample-count",
            "18",
            "--log-file",
            str(paths["log"]),
        ],
        cwd=ROOT,
        check=True,
    )


def test_numeric_label_column_is_accepted_without_disease_state_parsing(tmp_path):
    paths = write_common_inputs(tmp_path, ["0"] * 8 + ["1"] * 10)
    run_builder(paths)

    rows = list(csv.DictReader(paths["output"].open("r", encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 18
    assert {row["metadata_join_key"] for row in rows} == {"geo_accession"}
    summary = {row["label__pd_vs_control"]: row["count"] for row in csv.DictReader(paths["labels"].open("r", encoding="utf-8"), delimiter="\t")}
    assert summary == {"0": "8", "1": "10"}


def test_string_pd_control_label_variants_are_canonicalized(tmp_path):
    labels = ["Control", "control", "0", "0", "Control", "0", "control", "0"] + [
        "PD",
        "Parkinson's disease",
        "Parkinson disease",
        "1",
        "PD",
        "1",
        "Parkinson's disease",
        "Parkinson disease",
        "PD",
        "1",
    ]
    paths = write_common_inputs(tmp_path, labels)
    run_builder(paths)
    summary = {row["label__pd_vs_control"]: row["count"] for row in csv.DictReader(paths["labels"].open("r", encoding="utf-8"), delimiter="\t")}
    assert summary == {"0": "8", "1": "10"}
