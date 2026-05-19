import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLAIMS = ROOT / "results/reports/phase39_publication_claim_table.tsv"


def test_phase39_claim_table_uses_conservative_evidence_levels():
    rows = list(csv.DictReader(CLAIMS.open("r", encoding="utf-8"), delimiter="\t"))
    assert rows
    levels = {row["evidence_level"] for row in rows}
    assert {"strong", "nominal", "preliminary", "divergent", "insufficient"} <= levels
    neuronal = next(row for row in rows if row["claim_id"] == "claim_ad_neuronal_vulnerability")
    assert neuronal["evidence_level"] == "nominal"
    assert "AD-replicated neuronal vulnerability" in neuronal["allowed_language"]
    pd_div = next(row for row in rows if row["claim_id"] == "claim_pd_divergent_synuclein")
    assert pd_div["evidence_level"] == "divergent"
    assert "candidate PD-divergent" in pd_div["allowed_language"]


def test_phase39_claim_text_does_not_turn_divergence_into_shared_replication():
    rows = list(csv.DictReader(CLAIMS.open("r", encoding="utf-8"), delimiter="\t"))
    claim_text = "\n".join(row["claim_text"] + "\n" + row["allowed_language"] for row in rows).lower()
    assert "shared ad/pd replication" not in claim_text
    assert "validated shared" not in claim_text
    assert "clinical diagnostic" not in claim_text
    assert "causal mechanism" not in claim_text
