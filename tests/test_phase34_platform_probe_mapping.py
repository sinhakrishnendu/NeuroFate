import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MAPPER = ROOT / "scripts/164_prepare_phase34_platform_axis_probe_mapping.py"
BUILDER = ROOT / "scripts/165_build_pd_axis_scores_from_geo_expression.py"


def test_phase34_platform_mapper_keeps_only_axis_probes(tmp_path):
    platform = tmp_path / "GPL_mock.tsv"
    axis = tmp_path / "axis.tsv"
    aliases = tmp_path / "aliases.tsv"
    output = tmp_path / "probe_map.tsv"
    audit = tmp_path / "audit.tsv"
    log = tmp_path / "probe_map.log"
    platform.write_text("ID\tGene Symbol\nprobe_snca\tSNCA\nprobe_nefl\tNEFL\nprobe_random\tRANDOM\n", encoding="utf-8")
    axis.write_text("axis_id\tgene_members\nneuronal_vulnerability_axis\tSNCA;NEFL\n", encoding="utf-8")
    aliases.write_text("gene_symbol\tensembl_gene_id\talias_type\tsource_note\nSNCA\tENSG00000145335\tensembl_gene_id\ttest\n", encoding="utf-8")
    subprocess.run(
        [
            sys.executable,
            str(MAPPER),
            "--platform-file",
            str(platform),
            "--platform-id",
            "GPLMOCK",
            "--axis-registry",
            str(axis),
            "--alias-table",
            str(aliases),
            "--output",
            str(output),
            "--audit-output",
            str(audit),
            "--log-file",
            str(log),
        ],
        cwd=ROOT,
        check=True,
    )
    rows = list(csv.DictReader(output.open("r", encoding="utf-8"), delimiter="\t"))
    assert {row["probe_id"] for row in rows} == {"probe_snca", "probe_nefl"}
    assert "probe_random" not in output.read_text(encoding="utf-8")
    audit_rows = list(csv.DictReader(audit.open("r", encoding="utf-8"), delimiter="\t"))
    assert audit_rows[0]["probe_column"] == "ID"
    assert audit_rows[0]["mapped_gene_count"] == "2"


def test_phase34_platform_mapper_detects_geo_annotation_preamble(tmp_path):
    platform = tmp_path / "GPL570.annot"
    axis = tmp_path / "axis.tsv"
    aliases = tmp_path / "aliases.tsv"
    output = tmp_path / "probe_map.tsv"
    audit = tmp_path / "audit.tsv"
    log = tmp_path / "probe_map.log"
    platform.write_text(
        "\n".join(
            [
                "^Annotation",
                "!Annotation_platform = GPL570",
                "#ID = ID from Platform data table",
                "#Gene symbol = Entrez Gene symbol",
                "!platform_table_begin",
                "ID\tGene title\tGene symbol\tGene ID\tEnsembl\tGB_ACC",
                "probe1\talpha synuclein\tSNCA///RANDOM\t6622\tENSG00000145335\tNM_000345",
                "probe2\tneurofilament light\tNEFL // OTHER\t4747\t\tNM_006158",
                "probe3\tbackground\tNOTAXIS\t0\t\t",
                "!platform_table_end",
                "",
            ]
        ),
        encoding="utf-8",
    )
    axis.write_text("axis_id\tgene_members\nneuronal_vulnerability_axis\tSNCA;NEFL\n", encoding="utf-8")
    aliases.write_text("gene_symbol\tensembl_gene_id\talias_type\tsource_note\nSNCA\tENSG00000145335\tensembl_gene_id\ttest\nNEFL\tENSG00000104725\tensembl_gene_id\ttest\n", encoding="utf-8")
    subprocess.run(
        [
            sys.executable,
            str(MAPPER),
            "--platform-file",
            str(platform),
            "--platform-id",
            "GPL570",
            "--axis-registry",
            str(axis),
            "--alias-table",
            str(aliases),
            "--output",
            str(output),
            "--audit-output",
            str(audit),
            "--log-file",
            str(log),
        ],
        cwd=ROOT,
        check=True,
    )
    rows = list(csv.DictReader(output.open("r", encoding="utf-8"), delimiter="\t"))
    assert {row["probe_id"] for row in rows} == {"probe1", "probe2"}
    assert {row["gene_symbol"] for row in rows} == {"SNCA", "NEFL"}
    audit_row = list(csv.DictReader(audit.open("r", encoding="utf-8"), delimiter="\t"))[0]
    assert audit_row["parser_mode"] == "geo_platform_table"
    assert audit_row["header_line"] == "6"
    assert audit_row["probe_column"] == "ID"
    assert "Gene symbol" in audit_row["gene_columns"]


def test_phase34_expression_builder_outputs_axis_scores_not_genomewide_matrix(tmp_path):
    expression = tmp_path / "expr.tsv"
    metadata = tmp_path / "metadata.tsv"
    probe_map = tmp_path / "probe_map.tsv"
    axis = tmp_path / "axis.tsv"
    output = tmp_path / "axis_scores.tsv"
    coverage = tmp_path / "coverage.tsv"
    labels = tmp_path / "labels.tsv"
    log = tmp_path / "axis.log"
    expression.write_text("ID_REF\tGSM1\tGSM2\nprobe_snca\t1\t4\nprobe_nefl\t2\t5\nprobe_random\t99\t100\n", encoding="utf-8")
    metadata.write_text(
        "cohort_id\tsample_title\tgeo_accession\tdisease_state\ttissue_or_region\tlabel__pd_vs_control\nmock\tControl\tGSM1\tControl\tSN\t0\nmock\tPD\tGSM2\tParkinson's disease\tSN\t1\n",
        encoding="utf-8",
    )
    probe_map.write_text("platform_id\tprobe_id\tgene_symbol\tmulti_probe_gene\nGPL\tprobe_snca\tSNCA\tfalse\nGPL\tprobe_nefl\tNEFL\tfalse\n", encoding="utf-8")
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
            "mock_pd",
            "--axis-registry",
            str(axis),
            "--output",
            str(output),
            "--coverage-output",
            str(coverage),
            "--label-summary-output",
            str(labels),
            "--log-file",
            str(log),
        ],
        cwd=ROOT,
        check=True,
    )
    rows = list(csv.DictReader(output.open("r", encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 2
    assert "axis__neuronal_vulnerability_axis" in rows[0]
    assert "probe_random" not in output.read_text(encoding="utf-8")
