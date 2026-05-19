# ROSMAP / AD Knowledge Portal Manual Access Template

Do not run from Codex. Manual user execution only.

Expected local path:

```text
data/raw/external/rosmap_ad/
```

ROSMAP and related AD Knowledge Portal harmonized RNA-seq resources may require Synapse login, data-use certification, project approval, or controlled-access agreements.

Manual checklist:

1. Confirm the intended ROSMAP or AD Knowledge Portal dataset accession.
2. Confirm data-use terms and whether redistribution is prohibited.
3. Download only after approval, using official portal instructions.
4. Place approved files under `data/raw/external/rosmap_ad/`.
5. Record source, date accessed, license/terms, filename, size, and checksum in NeuroFate provenance tables.

Example placeholders, not executable commands:

```bash
# MANUAL_HEAVY only after approval:
# synapse get SYNAPSE_ID --downloadLocation data/raw/external/rosmap_ad/
# shasum -a 256 data/raw/external/rosmap_ad/OFFICIAL_FILE
```
