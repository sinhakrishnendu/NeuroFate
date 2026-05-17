from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/57_audit_feature_leakage.py"


def load_script():
    spec = importlib.util.spec_from_file_location("leakage_audit", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_leakage_audit_declares_required_rules() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for token in [
        "LABEL_PREFIX",
        "donor_id",
        "cohort_id",
        "label__",
        "diagnosis",
        "pathology",
        "apoe_genotype",
        "cognitive",
        "feature_leakage_audit.tsv",
    ]:
        assert token in text


def test_leakage_audit_classifies_columns() -> None:
    module = load_script()
    assert module.classify_column("label__Cognitive_Status")["leakage_risk"] == "high"
    assert module.classify_column("donor_id")["leakage_risk"] == "high"
    assert module.classify_column("cohort_id")["leakage_risk"] == "medium"
    assert module.classify_column("gene_mean__APOE")["column_role"] == "predictor"


def test_phase12_scripts_avoid_single_cell_engines() -> None:
    for path in sorted((ROOT / "scripts").glob("5[7-9]_*.py")) + sorted(
        (ROOT / "scripts").glob("6[0-3]_*.py")
    ):
        text = path.read_text(encoding="utf-8").lower()
        for forbidden in ["scanpy", "read_h5ad", "anndata", "scvi", "scvelo", "cellrank", "import h5py", "torch"]:
            assert forbidden not in text, path
