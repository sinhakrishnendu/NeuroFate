from pathlib import Path
import csv
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/164_prepare_geo_platform_gene_mapping.py"


def test_platform_mapper_keeps_only_neurofate_axis_probes(tmp_path):
    platform = tmp_path / "GPL_mock.tsv"
    axis_registry = tmp_path / "axis.tsv"
    aliases = tmp_path / "aliases.tsv"
    output = tmp_path / "mapping.tsv"
    log_file = tmp_path / "mapping.log"
    platform.write_text("ID\tGene Symbol\nprobe_snca\tSNCA\nprobe_random\tRANDOM\nprobe_nefl\tNEFL /// OTHER\n", encoding="utf-8")
    axis_registry.write_text("axis_id\tgene_members\nneuronal_vulnerability_axis\tSNCA;NEFL\n", encoding="utf-8")
    aliases.write_text("gene_symbol\tensembl_gene_id\talias_type\tsource_note\nSNCA\tENSG00000145335\tensembl_gene_id\ttest\n", encoding="utf-8")
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--platform-file",
            str(platform),
            "--axis-registry",
            str(axis_registry),
            "--alias-table",
            str(aliases),
            "--output",
            str(output),
            "--log-file",
            str(log_file),
        ],
        cwd=ROOT,
        check=True,
    )
    rows = list(csv.DictReader(output.open("r", encoding="utf-8"), delimiter="\t"))
    assert {row["probe_id"] for row in rows} == {"probe_snca", "probe_nefl"}
    assert {row["gene_symbol"] for row in rows} == {"SNCA", "NEFL"}
    assert "probe_random" not in output.read_text(encoding="utf-8")
