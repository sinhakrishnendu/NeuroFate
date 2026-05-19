import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
READINESS = ROOT / "scripts/123_build_pnas_readiness_matrix.py"


def test_phase38_readiness_marks_mixed_pd_evidence_and_divergence(tmp_path):
    tables = tmp_path / "results" / "tables"
    reports = tmp_path / "results" / "reports"
    tables.mkdir(parents=True)
    reports.mkdir(parents=True)
    (tables / "phase37_gse7621_pd_sn_bulk_pd_axis_replication_statistics.tsv").write_text(
        "axis_id\teffect_size\tpvalue\tfdr\tevidence_label\n"
        "synuclein_mitochondrial_axis\t-0.76\t0.0018\t0.018\topposite_direction\n",
        encoding="utf-8",
    )
    (tables / "phase32_crosscohort_axis_evidence_summary.tsv").write_text(
        "axis_id\tcrosscohort_evidence_class\tgse174367_fdr\n"
        "neuronal_vulnerability_axis\tstrong_ad_axis_with_nominal_external_replication\t0.24\n",
        encoding="utf-8",
    )
    subprocess.run(
        [
            sys.executable,
            str(READINESS),
            "--output",
            str(reports / "phase23_pnas_readiness_matrix.tsv"),
            "--gap-report",
            str(reports / "phase23_pnas_gap_report.md"),
        ],
        cwd=tmp_path,
        check=True,
    )
    rows = {row["criterion"]: row for row in csv.DictReader((reports / "phase23_pnas_readiness_matrix.tsv").open("r", encoding="utf-8"), delimiter="\t")}
    assert rows["independent_pd_replication"]["status"] == "mixed_pd_evidence"
    assert rows["pd_divergent_axis_candidate"]["status"] == "present"
    assert rows["shared_ad_pd_axis_claim"]["status"] == "not_ready"
    assert rows["pnas_biological_claim"]["status"] == "promising_but_requires_pd_resolution"
