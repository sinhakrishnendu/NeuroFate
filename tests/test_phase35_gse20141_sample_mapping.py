import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "scripts/168_audit_gse20141_sample_mapping.py"
PARSER = ROOT / "scripts/163_parse_pd_geo_series_matrix_metadata.py"


def write_mock_gse20141_series(path: Path, include_expression: bool = True) -> list[str]:
    gsm = [f"GSM5039{50 + idx:02d}" for idx in range(18)]
    titles = [f"Control {idx + 1}" for idx in range(8)] + [f"PD {idx + 1}" for idx in range(10)]
    labels = ["Control"] * 8 + ["Parkinson's disease"] * 10
    lines = [
        "!Sample_title\t" + "\t".join(titles),
        "!Sample_geo_accession\t" + "\t".join(gsm),
        "!Sample_source_name_ch1\t" + "\t".join(["substantia nigra pars compacta"] * 18),
        "!Sample_platform_id\t" + "\t".join(["GPL570"] * 18),
        "!Sample_characteristics_ch1\t" + "\t".join(f'"disease state: {label}"' for label in labels),
        "!Sample_characteristics_ch1\t" + "\t".join(['"brain region: substantia nigra pars compacta"'] * 18),
    ]
    if include_expression:
        lines.extend(["!series_matrix_table_begin", "ID_REF\t" + "\t".join(gsm), "probe1\t" + "\t".join(["1"] * 18), "!series_matrix_table_end"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return gsm


def test_gse20141_sample_mapping_audit_selects_geo_accession(tmp_path):
    series = tmp_path / "GSE20141_series_matrix.txt"
    metadata = tmp_path / "metadata.tsv"
    labels = tmp_path / "labels.tsv"
    platforms = tmp_path / "platforms.tsv"
    parsed_log = tmp_path / "parse.log"
    audit_out = tmp_path / "audit.tsv"
    preview = tmp_path / "preview.txt"
    audit_log = tmp_path / "audit.log"
    write_mock_gse20141_series(series)
    subprocess.run(
        [
            sys.executable,
            str(PARSER),
            "--series-matrix",
            str(series),
            "--cohort-id",
            "gse20141_pd_snpc_lcm",
            "--output",
            str(metadata),
            "--label-summary-output",
            str(labels),
            "--platform-output",
            str(platforms),
            "--log-file",
            str(parsed_log),
        ],
        cwd=ROOT,
        check=True,
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
            str(audit_out),
            "--preview-output",
            str(preview),
            "--log-file",
            str(audit_log),
        ],
        cwd=ROOT,
        check=True,
    )
    rows = list(csv.DictReader(audit_out.open("r", encoding="utf-8"), delimiter="\t"))
    best = next(row for row in rows if row["is_best_join_key"] == "true")
    assert best["best_join_key"] == "geo_accession"
    assert best["matched_sample_count"] == "18"
    assert best["unmatched_expression_samples"] == "0"
    assert best["pd_count"] == "10"
    assert best["control_count"] == "8"


def test_phase34_metadata_parser_writes_sample_id_as_geo_accession(tmp_path):
    series = tmp_path / "GSE20141_series_matrix.txt"
    metadata = tmp_path / "metadata.tsv"
    labels = tmp_path / "labels.tsv"
    platforms = tmp_path / "platforms.tsv"
    log = tmp_path / "parse.log"
    gsm = write_mock_gse20141_series(series, include_expression=False)
    subprocess.run(
        [
            sys.executable,
            str(PARSER),
            "--series-matrix",
            str(series),
            "--cohort-id",
            "gse20141_pd_snpc_lcm",
            "--output",
            str(metadata),
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
    rows = list(csv.DictReader(metadata.open("r", encoding="utf-8"), delimiter="\t"))
    assert rows[0]["sample_id"] == gsm[0]
    summary = {row["label__pd_vs_control"]: row["count"] for row in csv.DictReader(labels.open("r", encoding="utf-8"), delimiter="\t")}
    assert summary == {"0": "8", "1": "10"}
