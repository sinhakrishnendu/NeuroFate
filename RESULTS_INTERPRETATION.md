# NeuroFate Results Interpretation

## 1. Executive Summary

NeuroFate is strongest as a reproducible, memory-safe, donor-level neurodegeneration modeling workflow. It integrates metadata-safe single-nucleus processing, targeted sparse gene-panel extraction, donor-level statistical analysis, interpretable machine learning, Apple Silicon MPS support, and reviewer-facing audits.

The current evidence supports internal donor-level modeling claims more strongly than external-validation claims. Mathys 2019 provides preliminary external feasibility after harmonization, but the harmonized Mathys table has six sample-level units, which is not enough for definitive cross-cohort validation.

## 2. What NeuroFate Currently Does Well

- Preserves raw data boundaries and avoids full H5AD matrix loading in metadata and reporting phases.
- Uses donor-level aggregation to reduce single-cell pseudo-replication.
- Provides leakage audits, no-overclaiming audits, reproducibility manifests, source packaging, and results-review packaging.
- Runs a tiny synthetic demo suitable for installation and release smoke testing.

## 3. Strongest Supported Biological/Computational Findings

The strongest current computational finding is internal donor-level stratification for dementia/reference-related tasks when Phase 12 repeated benchmarks support it. These findings should be described as internal donor-level evidence, not clinical prediction and not definitive external validation.

Microglial, astrocyte, neuronal, and neurodegeneration feature groups are useful as structured biological layers, but individual gene or feature claims should remain conservative unless supported by repeated benchmarks, permutation controls, and ablation consistency.

## 4. Internal Validation Status

Internal validation is based on donor-level tables and repeated classical-model benchmarking. Claim strength should be read from `results/reports/claim_strength_table.tsv` and `results/reports/best_supported_claims.tsv`, not from a single best metric.

## 5. External Validation Status

Mathys 2019 is currently external feasibility evidence only. The sample-level feature table contains six harmonized units, which is useful for testing schema transfer and pipeline behavior but insufficient for strong validation language.

## 6. Why Mathys Is Preliminary Feasibility, Not Definitive Validation

The Mathys CSV files require harmonization from cell-level covariates and count matrices into donor/sample-level NeuroFate features. After harmonization, there are six sample-level units. This is too small for reliable AUROC uncertainty, calibration, or stable subgroup interpretation. Use “preliminary external feasibility” rather than stronger validation language.

## 7. Leakage Audit Interpretation

Leakage audit flags for `label__` columns, donor identifiers, and cohort identifiers should be interpreted as detected and excluded when those columns are not used as predictors. A leakage audit is only fatal if label or identifier columns enter the predictor matrix.

## 8. Apple Silicon/MPS Neural-Model Interpretation

The MPS model is a real Apple Silicon implementation for small donor-level neural modeling. It is not a foundation model, not a deep single-cell model, and not evidence of state-of-the-art neural performance. It should be presented as an Apple-Silicon-ready implementation layer.

## 9. What Cannot Be Claimed

- NeuroFate is not clinical-grade.
- NeuroFate is not a diagnostic tool.
- NeuroFate does not infer causality.
- NeuroFate is not a foundation model.
- Current external evidence is insufficient for language claiming definitive validation across cohorts.
- Current outputs do not establish patient-level diagnosis.

## 10. Reviewer-Risk Table

| Topic | Reviewer Risk | Conservative Response |
| --- | --- | --- |
| External validation | High | Describe Mathys as preliminary external feasibility. |
| Clinical translation | High | State research software only. |
| Causality | High | Use association language. |
| MPS neural model | Medium | Present as small donor-level implementation. |
| Internal prediction | Medium | Report repeated-seed uncertainty and permutation controls. |
| Software release | Low | Provide PyPI metadata, CI, demo mode, and source package. |

## 11. Next Validation Needed Before Nature Computational Science Submission

- Add larger independent external cohorts with adequate donor/sample-level units.
- Lock task definitions and leakage-audit rules before final model comparison.
- Run full Phase 12 repeated benchmarks, permutation controls, and feature ablations.
- Confirm source and results packages exclude raw data and large private artifacts.
- Update the manuscript claim language from claim-strength tables.

## 12. Recommended Manuscript Claim Language

Recommended: “NeuroFate provides a reproducible, memory-safe, donor-level framework for neurodegeneration systems biology and internally benchmarked predictive modeling, with preliminary external feasibility demonstrated on harmonized Mathys 2019 data.”

Avoid stronger language unless future validation supports it.
