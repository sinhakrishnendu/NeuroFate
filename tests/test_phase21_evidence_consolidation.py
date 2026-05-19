from __future__ import annotations

import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLAIMS = ROOT / "scripts/65_build_claim_strength_table.py"
AUDIT = ROOT / "scripts/54_no_overclaiming_audit.py"
METRICS = ROOT / "results/tables/phase20_gse243639_celltype_validation_metrics.tsv"
README = ROOT / "README.md"
DOC = ROOT / "docs/external_validation_expansion.md"
INTERPRETATION = ROOT / "RESULTS_INTERPRETATION.md"
MANUSCRIPT = ROOT / "manuscript/neurofate_landmark_manuscript.tex"


def load_claim_module():
    spec = importlib.util.spec_from_file_location("phase21_claims", CLAIMS)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_phase20_metrics_are_preliminary_not_moderate_or_clinical() -> None:
    with METRICS.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    row = rows[0]
    assert row["annotation_match_rate"] == "1.0"
    assert row["feature_count"] == "1590"
    assert row["auroc"] == "0.72777778"
    assert row["auprc"] == "0.79044444"
    assert row["balanced_accuracy"] == "0.64666667"
    assert row["empirical_permutation_pvalue"] == "0.10891089"
    assert row["reliability_flag"] == "preliminary_pd_internal_signal"
    assert "clinical" not in row["notes"].lower()


def test_claim_builder_prioritizes_phase20_over_phase17_phase18() -> None:
    text = CLAIMS.read_text(encoding="utf-8")
    assert "phase20_gse243639_celltype_validation_metrics.tsv" in text
    assert "phase20_claim_strength_delta.tsv" in text
    assert "superseded_by_phase20" in text
    assert "preliminary sample-level cell-type-aware PD internal signal" in text
    assert "Do not claim clinical PD prediction" in text
    module = load_claim_module()
    status = module.phase20_pd_status(module.read_tsv(METRICS), module.read_tsv(ROOT / "results/tables/phase20_gse243639_feature_group_counts.tsv"))
    assert status[0] == "preliminary_pd_internal_signal"
    assert status[1] == 29
    assert status[2] == "1590"
    assert status[3] == "1.0"


def test_phase21_docs_and_manuscript_use_conservative_phase20_language() -> None:
    combined = "\n".join(
        [
            README.read_text(encoding="utf-8"),
            DOC.read_text(encoding="utf-8"),
            INTERPRETATION.read_text(encoding="utf-8"),
            MANUSCRIPT.read_text(encoding="utf-8"),
        ]
    ).lower()
    assert "preliminary_pd_internal_signal" in combined
    assert "0.72777778" in combined
    assert "0.10891089" in combined
    assert "improves over phase 16" in combined
    assert "clinical validation" in combined
    assert "diagnostic classifier" in combined
    assert "causal" not in combined or "not" in combined


def test_no_overclaiming_audit_has_phase20_allowed_context() -> None:
    text = AUDIT.read_text(encoding="utf-8")
    assert "preliminary sample-level cell-type-aware pd internal signal" in text
    assert "phase 20 gse243639 pd evidence is preliminary" in text
    assert "clinical pd validation" in text
