from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parents[1]


def read_registry():
    with (ROOT / "metadata/neurofate_axis_registry.tsv").open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_axis_registry_has_required_axes():
    rows = read_registry()
    axis_ids = {row["axis_id"] for row in rows}
    required = {
        "inflammatory_microglial_axis",
        "astrocyte_stress_axis",
        "myelin_oligodendrocyte_axis",
        "neuronal_vulnerability_axis",
        "synuclein_mitochondrial_axis",
        "amyloid_tau_axis",
        "immune_antigen_presentation_axis",
        "vascular_barrier_axis",
        "proteostasis_autophagy_axis",
        "global_neurodegeneration_axis",
    }
    assert required <= axis_ids


def test_axis_registry_schema_and_genes():
    rows = read_registry()
    assert rows
    expected_columns = {
        "axis_id",
        "axis_name",
        "biological_theme",
        "gene_members",
        "celltype_context",
        "expected_direction_in_ad",
        "expected_direction_in_pd",
        "primary_evidence_source",
        "interpretability_note",
        "overclaiming_risk",
    }
    assert expected_columns <= set(rows[0])
    all_genes = ";".join(row["gene_members"] for row in rows)
    for gene in ["SNCA", "MAPT", "APP", "TREM2", "GFAP", "MBP", "HLA-DRA", "PINK1", "PRKN"]:
        assert gene in all_genes
