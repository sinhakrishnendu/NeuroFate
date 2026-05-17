# Reviewer Guidance

This guide explains how to inspect NeuroFate as a computational platform without requiring raw SEA-AD or Mathys data.

## Inspect Results

Start with:

- `RESULTS_INTERPRETATION.md`
- `results/reports/claim_strength_table.tsv`
- `results/reports/best_supported_claims.tsv`
- `results/reports/reviewer_audit_report.md`
- `results/reports/no_overclaiming_audit.tsv`

## Verify No Raw Data Are Bundled

Source release packages are built with:

```bash
python scripts/66_create_source_release_package.py
```

Results-review packages are built with:

```bash
python scripts/67_create_results_review_package.py
```

Both package builders exclude `data/`, raw H5AD files, model binaries, and large external files by rule.

## Run The Demo

```bash
python -m neurofate run-demo
```

The demo uses bundled synthetic data only and writes to `results/demo/`.

## Interpret Evidence Categories

- `strong_internal`: repeated donor-level benchmarks, low instability, permutation support, and no leakage.
- `moderate_internal`: useful internal evidence, but not enough for stronger claims.
- `exploratory_internal`: signal exists but needs replication.
- `preliminary_external_feasibility`: external pipeline transfer exists, but validation is underpowered.
- `insufficient_external_validation`: external evidence is missing or too weak.
- `failed_or_unstable`: weak or unstable benchmark support.

## Why Mathys Is Preliminary

Mathys 2019 is harmonized into six sample-level units in the current workflow. This is useful for feasibility and schema-transfer testing, but too small for definitive external validation.

## Full Cross-Cohort Validation Needs

- More external donor/sample-level units.
- Locked task labels.
- Repeated validation across cohorts.
- Permutation and ablation support.
- No high-severity overclaiming or unresolved leakage.
