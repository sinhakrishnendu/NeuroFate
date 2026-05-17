# Benchmarking And Validation

NeuroFate Phase 12 adds a leakage-aware benchmarking layer for donor-level analysis. It is designed for reproducibility, uncertainty reporting, and conservative claim strength.

## Donor-Level Validation

Benchmarks use `results/tables/phase5_donor_feature_table.tsv`. They do not read H5AD files, run Scanpy, or load single-cell matrices. Task labels are derived from `label__` columns, and predictors are restricted to approved feature prefixes such as `gene_mean__`, `gene_detection__`, `index__`, `cell_fraction__`, and `celltype_index__`.

## Leakage Prevention

Run:

```bash
python scripts/57_audit_feature_leakage.py --input results/tables/phase5_donor_feature_table.tsv --output results/reports/feature_leakage_audit.tsv
```

The audit flags label columns, donor/sample identifiers, cohort identifiers, and suspicious non-feature columns before modeling.

## Repeated Seeds

Run repeated classical baselines with:

```bash
python scripts/58_run_repeated_baseline_benchmarks.py --features results/tables/phase5_donor_feature_table.tsv --config configs/benchmark_config.yaml
```

The benchmark computes mean, standard deviation, and approximate confidence intervals for AUROC, AUPRC, balanced accuracy, and Brier score across configured seeds.

## Permutation Controls

Permutation controls shuffle labels within valid task rows and estimate empirical p-values:

```bash
python scripts/59_run_label_permutation_controls.py --features results/tables/phase5_donor_feature_table.tsv --config configs/benchmark_config.yaml
```

If real performance is not better than the null distribution, claims should be weakened.

## Feature Ablation

Feature ablation removes groups such as gene-level features, cell fractions, cell-type indices, inflammatory signatures, astrocyte signatures, neuronal signatures, and mitochondrial/neurodegeneration signatures:

```bash
python scripts/60_run_feature_ablation.py --features results/tables/phase5_donor_feature_table.tsv --config configs/benchmark_config.yaml
```

Feature groups should only be interpreted as robust when ablation effects are consistent across tasks and seeds.

## Evidence Categories

`scripts/63_classify_evidence_strength.py` assigns:

- `strong_internal`
- `moderate_internal`
- `preliminary_external`
- `insufficient`
- `failed_or_unstable`

Rules use sample size, repeated AUROC stability, permutation p-values, ablation consistency, external validation availability, and no-overclaiming status.

## Mathys Validation

Mathys 2019 should be described as preliminary external feasibility evidence unless sample-level granularity, feature overlap, and repeated validation support stronger claims.
