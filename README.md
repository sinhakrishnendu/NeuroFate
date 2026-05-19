# NeuroFate

NeuroFate is an Apple-Silicon-ready computational platform for neurodegeneration systems biology. It connects metadata-safe single-nucleus workflows, targeted sparse gene-panel extraction, donor-level statistical biology, interpretable machine learning, Apple Silicon MPS neural modeling, external cohort feasibility checks, and reviewer-aware reporting.

## What NeuroFate Does

- Builds reproducible neurodegeneration feature layers from metadata, sparse gene panels, and donor-level summaries.
- Supports SEA-AD and Mathys 2019 workflows while preserving room for Parkinson disease, microbiome/metabolite, protein interaction, and evolutionary modules.
- Produces manuscript-oriented tables, figures, model summaries, reports, reproducibility manifests, and no-overclaiming audits.
- Provides a guarded CLI so common workflows are discoverable without accidentally running heavy analyses.

## Who Should Use It

NeuroFate is intended for computational biologists, neurogenomics researchers, translational data scientists, and Apple Silicon users who need a reproducible, memory-conscious workflow for neurodegeneration cohort analysis. It is a research platform, not a clinical product.

## Accepted Data

- SEA-AD H5AD files for metadata-only inspection and carefully bounded sparse target-gene extraction.
- Mathys 2019 GEO CSV count and covariate files for external feasibility analysis.
- Donor-level feature tables for classical machine learning and small MPS neural models.
- TSV registries for datasets, provenance, feature definitions, manuscript modules, gene panels, and validation plans.

## Quickstart

Recommended Python: 3.11 or 3.12.

After the first PyPI release:

```bash
python -m pip install neurofate
neurofate check-system
neurofate run-demo
```

From the GitHub repository checkout:

```bash
conda env create -f environment.yml
conda activate neurofate
python -m pip install -e .
neurofate check-system
neurofate doctor
```

For a venv workflow:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
neurofate check-system
```

Run the bundled tiny demo without external data:

```bash
neurofate run-demo
```

If the editable install has not created the console entry point yet, use:

```bash
python -m neurofate run-demo
```

The demo uses `examples/tiny_demo/`, writes to `results/demo/`, and is meant for installation and CLI smoke testing only.

Generate end-user platform reports from existing outputs:

```bash
python scripts/51_generate_end_user_report.py --tables-dir results/tables --reports-dir results/reports
python scripts/52_generate_reproducibility_manifest.py --output results/reports/reproducibility_manifest.json
python scripts/53_validate_neurofate_outputs.py --output results/reports/output_validation_report.tsv
python scripts/54_no_overclaiming_audit.py --output results/reports/no_overclaiming_audit.tsv
```

## Apple Silicon Notes

NeuroFate keeps full single-cell matrices out of memory and moves modeling to donor-level tables. Optional PyTorch MPS workflows use small neural models only, with CPU fallback when MPS is unavailable. PyTorch is optional and should be installed deliberately:

```bash
python -m pip install -e ".[torch]"
```

## Safety And Memory Design

- No datasets are downloaded automatically.
- Heavy commands are written as manual templates.
- H5AD expression matrix arrays are never touched by metadata-only scripts.
- Sparse expression extraction is target-gene bounded, chunked, and protected against dense conversion.
- Model scripts operate on donor-level tables rather than full single-cell matrices.
- Generated reports summarize existing outputs and do not recompute analyses.

## Current Validation Status

SEA-AD supports the internal NeuroFate analysis path. Mathys 2019 currently serves as preliminary external feasibility evidence unless larger sample-level harmonization is added. Phase 15 adds a multi-cohort external-validation planning layer for larger AD and PD cohorts, but no claim should be upgraded until donor/sample-level external validation is reliable.

## Limitations

NeuroFate is not a clinical system, does not establish cause-and-effect biology by itself, and does not replace cohort-specific quality control. External validation remains limited by available harmonized cohorts, access constraints, sample counts, and feature overlap.

## Full Workflow Overview

1. Validate registries, safety flags, provenance templates, and manuscript-module alignment.
2. Manually acquire allowed external datasets from official sources.
3. Extract metadata-only SEA-AD tables, decode categories, and prepare Table 1.
4. Plan and run bounded sparse target-gene extraction.
5. Generate Phase 3/4 biological summaries and statistical findings.
6. Build donor-level Phase 5 classical ML and optional Phase 6 MPS neural models.
7. Use Mathys 2019 as external feasibility evidence while avoiding overclaiming.
8. Generate reports, manifests, output inventories, and no-overclaiming audits.

## What Is Bundled

Bundled: source code, configuration templates, registries, documentation, tests, CI, and a tiny synthetic demo dataset.

Not bundled: SEA-AD, Mathys GEO files, ROSMAP, Synapse-controlled data, STRING, HMDB, KEGG, trained real-data models, or large generated result artifacts.

## How To Cite

Use `CITATION.cff` for the NeuroFate software citation and cite each external dataset according to its source-specific instructions. NeuroFate outputs should acknowledge both the software and the underlying data resources.

## Release Status

Current release metadata target: `0.1.0`. This is a research software release candidate focused on reproducibility, safety, and reviewer-facing transparency.

## PyPI And GitHub Release

NeuroFate is prepared for PyPI distribution under the package name `neurofate`. The console entry point is:

```bash
neurofate
```

Before uploading to PyPI or pushing a public GitHub release, review:

- `PYPI_RELEASE_CHECKLIST.md`
- `RELEASE_CHECKLIST.md`
- `docs/pypi_release.md`
- `LICENSE`
- `CITATION.cff`
- `codemeta.json`

Manual release build commands:

```bash
python -m pip install -e ".[dev]"
python -m build
python -m twine check dist/*
```

## First Manual Check

Run this from inside the project root:

```bash
python3 scripts/00_check_system.py --dry-run --config configs/project_config.yaml --log-file results/logs/00_check_system.log
```

Please paste back the terminal output, `results/logs/00_check_system.log`, and any Python or architecture warnings.

## Phase 10: Platform Hardening And End-User Readiness

Phase 10 turns NeuroFate into a reusable platform layer: installable package metadata, a guarded CLI, user-facing documentation, report generation, reproducibility manifests, output validation, and no-overclaiming review.

### LIGHTWEIGHT: inspect the platform

```bash
neurofate check-system
neurofate doctor
```

### LIGHTWEIGHT: generate end-user artifacts from existing outputs

```bash
python scripts/51_generate_end_user_report.py --tables-dir results/tables --reports-dir results/reports
python scripts/52_generate_reproducibility_manifest.py --output results/reports/reproducibility_manifest.json
python scripts/53_validate_neurofate_outputs.py --output results/reports/output_validation_report.tsv
python scripts/54_no_overclaiming_audit.py --output results/reports/no_overclaiming_audit.tsv
```

These commands summarize existing files only. They do not download data, read full H5AD matrices, run Scanpy/scVI/scVelo, train models, or recompute expression analyses.

## Phase 11: Release Engineering And Tiny Demo

Phase 11 adds CI, citation metadata, release notes, contributor guidance, data-use documentation, runtime benchmark placeholders, an output inventory command, and the bundled tiny synthetic demo.

### LIGHTWEIGHT: run the tiny demo

```bash
neurofate run-demo
```

Expected outputs:

- `results/demo/demo_donor_feature_table.tsv`
- `results/demo/demo_model_metrics.tsv`
- `results/demo/demo_report.md`

### LIGHTWEIGHT: inventory generated outputs

```bash
python scripts/56_inventory_outputs.py --output results/reports/output_inventory.tsv
```

## Phase 12: Benchmarking, Leakage Audit, and Robustness Validation

Phase 12 adds donor-level leakage auditing, repeated classical-model benchmarks, label-permutation controls, feature-ablation studies, uncertainty reporting, benchmark figures, and evidence-strength classification. These scripts use donor-level tables only and do not touch H5AD files or single-cell matrices.

### LIGHTWEIGHT: audit feature leakage

```bash
python scripts/57_audit_feature_leakage.py --input results/tables/phase5_donor_feature_table.tsv --output results/reports/feature_leakage_audit.tsv
```

### MANUAL BENCHMARK: run repeated baselines

```bash
python scripts/58_run_repeated_baseline_benchmarks.py --features results/tables/phase5_donor_feature_table.tsv --config configs/benchmark_config.yaml
```

### MANUAL BENCHMARK: run robustness controls

```bash
python scripts/59_run_label_permutation_controls.py --features results/tables/phase5_donor_feature_table.tsv --config configs/benchmark_config.yaml
python scripts/60_run_feature_ablation.py --features results/tables/phase5_donor_feature_table.tsv --config configs/benchmark_config.yaml
python scripts/63_classify_evidence_strength.py
python scripts/61_generate_benchmark_uncertainty_report.py
python scripts/62_generate_phase12_benchmark_figures.py
```

## Phase 13: Evidence Cleanup, Claim Hardening, and Release Packaging

Phase 13 fixes Mathys reporting, builds conservative claim-strength tables, improves no-overclaiming checks, creates reviewer-facing interpretation files, and prepares clean source/results package builders.

### LIGHTWEIGHT: regenerate Phase 13 audit artifacts

```bash
python scripts/64_audit_mathys_gene_extraction.py
python scripts/45_generate_phase9_results_text.py --tables-dir results/tables --output results/tables/phase9_results_summary.txt
python scripts/65_build_claim_strength_table.py
python scripts/54_no_overclaiming_audit.py --output results/reports/no_overclaiming_audit.tsv
python scripts/66_create_source_release_package.py
python scripts/67_create_results_review_package.py
python scripts/68_generate_reviewer_audit_report.py
```

Phase 13 outputs include:

- `results/tables/phase13_mathys_gene_extraction_audit.tsv`
- `results/reports/claim_strength_table.tsv`
- `results/reports/best_supported_claims.tsv`
- `results/reports/claim_language_matrix.tsv`
- `results/reports/reviewer_audit_report.md`
- `release_artifacts/neurofate_source_release_<timestamp>.zip`
- `release_artifacts/neurofate_results_review_<timestamp>.zip`

`dist/` is reserved for PyPI artifacts generated by `python -m build`, such as `.whl` and `.tar.gz` files. Reviewer/source ZIP packages are written under `release_artifacts/`.

The Mathys wording is intentionally conservative: six harmonized sample-level units support preliminary external feasibility, not definitive external validation.

## Phase 15: Real External Validation Expansion

Phase 15 adds a generic external-cohort onboarding and validation framework for independent AD, PD, controlled-access, and bulk RNA-seq cohorts. It is planning-first: candidates are registered, local files can be inventoried after manual acquisition, metadata can be mapped safely, sparse extraction templates can be generated, and claim strength can be updated only when reliable external validation exists.

Registered candidates:

- `gse243639_pd_snpc`: priority A PD substantia nigra pars compacta snRNA-seq validation target.
- `gse174367_ad_multiomics`: priority B AD multi-omics validation target.
- `gse147528_ad_snrna`: priority C AD progression snRNA-seq validation target.
- `rosmap_ad_bulk_rnaseq`: priority D controlled/registered bulk validation target.
- `ad_knowledge_portal_harmonized_optional`: optional controlled AD Knowledge Portal harmonized cohort.

### LIGHTWEIGHT: triage candidate cohorts

```bash
python scripts/69_triage_external_validation_candidates.py
```

Expected outputs:

- `results/reports/phase15_external_dataset_triage.tsv`
- `results/reports/phase15_external_dataset_priority_summary.md`

### LIGHTWEIGHT: inspect manually acquired external files

```bash
python scripts/70_inspect_external_dataset_files.py \
  --dataset-id gse243639_pd_snpc \
  --input-dir data/raw/external/gse243639_pd_snpc \
  --output-summary results/reports/phase15_gse243639_pd_snpc_file_inventory.tsv \
  --format-output results/reports/phase15_gse243639_pd_snpc_format_recommendation.tsv \
  --log-file results/logs/70_inspect_external_dataset_files.log
```

### LIGHTWEIGHT: inspect external metadata safely

```bash
python scripts/71_inspect_external_metadata_safe.py \
  --dataset-id gse243639_pd_snpc \
  --metadata-file data/raw/external/gse243639_pd_snpc/METADATA_FILE \
  --format auto \
  --output results/reports/phase15_gse243639_pd_snpc_metadata_field_audit.tsv \
  --mapping-output results/reports/phase15_gse243639_pd_snpc_canonical_mapping_suggestions.tsv \
  --log-file results/logs/71_inspect_external_metadata_safe.log
```

### MANUAL TEMPLATE ONLY: plan external sparse extraction

```bash
python scripts/73_prepare_external_sparse_extraction_plan.py \
  --dataset-id gse243639_pd_snpc \
  --format h5ad_csr \
  --input-matrix data/raw/external/gse243639_pd_snpc/COUNTS_OR_CONTAINER_FILE \
  --metadata-file data/raw/external/gse243639_pd_snpc/METADATA_FILE \
  --feature-file data/raw/external/gse243639_pd_snpc/FEATURE_FILE \
  --panel metadata/target_gene_panel_v1.tsv \
  --output-plan results/tables/phase15_gse243639_pd_snpc_sparse_extraction_plan.tsv \
  --manual-script-output results/logs/manual_phase15_gse243639_pd_snpc_extraction_template.sh \
  --log-file results/logs/73_prepare_external_sparse_extraction_plan.log
```

### LIGHTWEIGHT: generate Phase 15 report

```bash
python scripts/77_generate_phase15_external_validation_report.py \
  --output results/reports/phase15_external_validation_report.md
```

CLI equivalents:

```bash
neurofate external-triage
neurofate inspect-external
neurofate plan-external-extraction
neurofate validate-multi-external
neurofate external-report
```

These commands do not download datasets, do not read full H5AD matrices, do not run Scanpy, and do not train models. Manual acquisition templates live in `scripts/manual_downloads/`.

## Phase 16: GSE243639 PD External Validation

Phase 16 onboards GSE243639 as the first serious Parkinson disease external cohort extension. The workflow parses the semicolon-delimited clinical file with header line 6, streams the genes-as-rows count CSV to extract only NeuroFate target genes, aggregates features to sample level, and runs conservative PD/control validation only when the user explicitly launches it.

This is independent PD cohort validation/extension. It is not clinical validation, not cause-and-effect inference, and not direct AD-to-PD disease-label transfer.

### MANUAL: extract GSE243639 target genes

```bash
python scripts/78_extract_gse243639_target_gene_panel.py \
  --counts data/raw/external/gse243639_pd_snpc/GSE243639_Filtered_count_table.csv.gz \
  --panel metadata/target_gene_panel_v1.tsv \
  --output data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz \
  --audit-output results/tables/phase16_gse243639_gene_extraction_audit.tsv \
  --log-file results/logs/78_extract_gse243639_target_gene_panel.log
```

### MANUAL: build GSE243639 sample-level features

```bash
python scripts/79_build_gse243639_feature_table.py \
  --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz \
  --clinical data/raw/external/gse243639_pd_snpc/GSE243639_Clinical_data.csv.gz \
  --phase5-schema results/tables/phase5_donor_feature_table.tsv \
  --output results/tables/phase16_gse243639_feature_table.tsv \
  --schema-output results/tables/phase16_gse243639_feature_schema_alignment.tsv \
  --label-summary-output results/tables/phase16_gse243639_label_summary.tsv \
  --log-file results/logs/79_build_gse243639_feature_table.log
```

### MANUAL: run conservative PD validation

```bash
python scripts/80_run_gse243639_pd_external_validation.py \
  --sea-ad-features results/tables/phase5_donor_feature_table.tsv \
  --pd-features results/tables/phase16_gse243639_feature_table.tsv \
  --metrics-output results/tables/phase16_gse243639_external_validation_metrics.tsv \
  --predictions-output results/tables/phase16_gse243639_external_predictions.tsv \
  --log-file results/logs/80_run_gse243639_pd_external_validation.log
```

## Phase 17: GSE243639 Cell-Type-Aware PD Refinement

Phase 17 uses the existing GSE243639 UMAP/annotation workbook only as an annotation table. NeuroFate does not recompute UMAP, does not cluster cells, and keeps all derived features at sample level before PD/control validation. The goal is to test whether cell-type-aware aggregation improves the exploratory PD signal from Phase 16.

### MANUAL: inspect the annotation workbook

```bash
python scripts/83_inspect_gse243639_umap_annotations.py \
  --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx \
  --output results/reports/phase17_gse243639_umap_annotation_audit.tsv \
  --preview-output results/reports/phase17_gse243639_umap_annotation_preview.tsv \
  --log-file results/logs/83_inspect_gse243639_umap_annotations.log
```

### MANUAL: build the cell annotation map

```bash
python scripts/84_build_gse243639_cell_annotation_map.py \
  --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx \
  --cell-sample-map data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv \
  --output data/interim/external/gse243639_pd_snpc/gse243639_cell_annotation_map.tsv \
  --summary-output results/tables/phase17_gse243639_cell_annotation_summary.tsv \
  --log-file results/logs/84_build_gse243639_cell_annotation_map.log
```

### MANUAL: build cell-type-aware sample features

```bash
python scripts/85_build_gse243639_celltype_feature_table.py \
  --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz \
  --annotations data/interim/external/gse243639_pd_snpc/gse243639_cell_annotation_map.tsv \
  --clinical data/raw/external/gse243639_pd_snpc/GSE243639_Clinical_data.csv.gz \
  --phase5-schema results/tables/phase5_donor_feature_table.tsv \
  --output results/tables/phase17_gse243639_celltype_feature_table.tsv \
  --schema-output results/tables/phase17_gse243639_celltype_schema_alignment.tsv \
  --label-summary-output results/tables/phase17_gse243639_celltype_label_summary.tsv \
  --log-file results/logs/85_build_gse243639_celltype_feature_table.log
```

### MANUAL: run robust cell-type-aware PD validation

```bash
python scripts/86_run_gse243639_celltype_pd_validation.py \
  --features results/tables/phase17_gse243639_celltype_feature_table.tsv \
  --metrics-output results/tables/phase17_gse243639_celltype_validation_metrics.tsv \
  --predictions-output results/tables/phase17_gse243639_celltype_predictions.tsv \
  --importance-output results/tables/phase17_gse243639_celltype_feature_importance.tsv \
  --log-file results/logs/86_run_gse243639_celltype_pd_validation.log
```

Phase 17 evidence remains exploratory unless robustness controls support a stronger category. It is not medical validation, not cause-and-effect inference, and not direct AD-to-PD disease-label transfer.

## Phase 18: GSE243639 Cell-ID Repair and Repaired Cell-Type-Aware PD Validation

Phase 18 repairs the Phase 17 annotation-join failure. Phase 17 should be treated as a technical audit result, not a biological conclusion, when annotation matching collapses the feature table. Phase 18 audits cell-ID formats, rebuilds the annotation map against expression cell IDs, preserves Phase 16 global gene features, and adds repaired cell-type/cluster-aware features only after successful matching.

Final PD interpretation should use Phase 18 rather than Phase 17 if the repair succeeds.

```bash
python scripts/90_audit_gse243639_cell_id_matching.py --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz --cell-sample-map data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx --output results/tables/phase18_gse243639_cell_id_matching_audit.tsv --preview-output results/reports/phase18_gse243639_cell_id_matching_preview.tsv --log-file results/logs/90_audit_gse243639_cell_id_matching.log
python scripts/84_build_gse243639_cell_annotation_map.py --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx --cell-sample-map data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv --output data/interim/external/gse243639_pd_snpc/gse243639_cell_annotation_map.tsv --summary-output results/tables/phase18_gse243639_annotation_match_summary.tsv --candidate-output results/reports/phase18_gse243639_annotation_column_candidates.tsv --log-file results/logs/84_build_gse243639_cell_annotation_map.log
python scripts/85_build_gse243639_celltype_feature_table.py --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz --annotations data/interim/external/gse243639_pd_snpc/gse243639_cell_annotation_map.tsv --clinical data/raw/external/gse243639_pd_snpc/GSE243639_Clinical_data.csv.gz --phase5-schema results/tables/phase5_donor_feature_table.tsv --output results/tables/phase18_gse243639_celltype_feature_table.tsv --schema-output results/tables/phase18_gse243639_celltype_schema_alignment.tsv --label-summary-output results/tables/phase18_gse243639_celltype_label_summary.tsv --feature-group-output results/tables/phase18_gse243639_feature_group_counts.tsv --log-file results/logs/85_build_gse243639_celltype_feature_table.log
python scripts/91_run_gse243639_repaired_celltype_pd_validation.py --features results/tables/phase18_gse243639_celltype_feature_table.tsv --metrics-output results/tables/phase18_gse243639_celltype_validation_metrics.tsv --predictions-output results/tables/phase18_gse243639_celltype_predictions.tsv --importance-output results/tables/phase18_gse243639_celltype_feature_importance.tsv --log-file results/logs/91_run_gse243639_repaired_celltype_pd_validation.log
```

## Phase 19: GSE243639 Annotation-Linkage Forensics

Phase 19 decides whether the `GSE243639_UMAP_coordinates.xlsx` workbook can be safely linked to expression cells at all. Phase 16 remains the valid global sample-level PD extension. Phase 17 and Phase 18 cell-type-aware outputs should not be interpreted biologically when annotation linkage is unsafe. The correct scientific decision is to retire the workbook annotation route rather than force a mapping.

Cell-type-aware PD validation requires one of these audited outcomes:

- direct ID linkage with high overlap,
- normalized ID linkage with high overlap,
- cautious row-order linkage only after independent row-count and sample-grouping checks.

If none passes, Phase 19 writes an unsafe/inconclusive decision and blocks creation of a safe annotation map.

```bash
python scripts/95_forensic_gse243639_workbook_cell_ids.py --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz --cell-sample-map data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx --output results/reports/phase19_gse243639_cell_id_forensic_preview.tsv --log-file results/logs/95_forensic_gse243639_workbook_cell_ids.log
python scripts/96_deep_gse243639_cell_id_normalization_audit.py --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz --cell-sample-map data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx --output results/tables/phase19_gse243639_deep_cell_id_overlap.tsv --best-rule-output results/reports/phase19_gse243639_best_normalization_rule.md --log-file results/logs/96_deep_gse243639_cell_id_normalization_audit.log
python scripts/97_audit_gse243639_row_order_annotation_link.py --cell-sample-map data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx --clinical data/raw/external/gse243639_pd_snpc/GSE243639_Clinical_data.csv.gz --output results/tables/phase19_gse243639_row_order_link_audit.tsv --preview-output results/reports/phase19_gse243639_row_order_link_preview.tsv --log-file results/logs/97_audit_gse243639_row_order_annotation_link.log
python scripts/98_decide_gse243639_annotation_linkage.py --normalization-audit results/tables/phase19_gse243639_deep_cell_id_overlap.tsv --row-order-audit results/tables/phase19_gse243639_row_order_link_audit.tsv --phase18-match-summary results/tables/phase18_gse243639_annotation_match_summary.tsv --output-md results/reports/phase19_gse243639_annotation_linkage_decision.md --output-tsv results/tables/phase19_gse243639_annotation_linkage_decision.tsv
python scripts/99_build_gse243639_safe_annotation_map_if_valid.py --decision results/tables/phase19_gse243639_annotation_linkage_decision.tsv --existing-annotation-map data/interim/external/gse243639_pd_snpc/gse243639_cell_annotation_map.tsv --output data/interim/external/gse243639_pd_snpc/gse243639_safe_cell_annotation_map.tsv --blocked-output results/reports/phase19_annotation_linkage_blocked.md
```

## Phase 20: Safe Annotation Map Consumption and Corrected Cell-Type-Aware PD Validation

Phase 19 found safe normalized ID linkage for GSE243639. Phase 20 repairs the feature-builder interface so the safe annotation map is consumed through `cell_id_expression`, with unique-cell annotation match rates reported correctly. Only Phase 20 should be used for final cell-type-aware PD interpretation if annotation match rate is high; Phase 16 remains the valid global sample-level PD extension.

Current Phase 20 evidence: annotation match rate `1.0`, feature count `1590`, logistic repeated AUROC `0.72777778`, AUPRC `0.79044444`, balanced accuracy `0.64666667`, empirical permutation p-value `0.10891089`, and reliability `preliminary_pd_internal_signal`. This improves over Phase 16 global features but remains preliminary and must not be described as clinical validation, a diagnostic classifier, causal evidence, or cross-disease validation.

## Phase 21: PNAS-Oriented NeuroFate Biological Discovery Framework

Phase 21 reframes NeuroFate as a systems neurobiology discovery framework centered on donor/sample-level neurodegeneration axes rather than clinical prediction. The central hypothesis is that AD and PD share candidate glial-inflammatory, myelin, and neuronal vulnerability axes while diverging in amyloid/tau- and synuclein-associated axis structure.

Phase 21 adds `metadata/neurofate_axis_registry.tsv`, `PNAS_DISCOVERY_STRATEGY.md`, and `docs/pnas_validation_strategy.md`. Axis scripts use only existing donor/sample-level tables and do not access raw H5AD files, run Scanpy, recompute UMAP, cluster cells, or train deep models.

Manual Phase 21 command sequence:

```bash
python scripts/106_build_neurofate_axis_scores.py && \
python scripts/107_test_neurofate_axis_associations.py && \
python scripts/108_compare_ad_pd_axis_patterns.py && \
python scripts/109_run_axis_randomization_controls.py --n-random 500 && \
python scripts/110_generate_phase21_pnas_axis_figures.py && \
python scripts/111_generate_phase21_pnas_biology_report.py
```

Phase 21 claim language must remain conservative: candidate shared axis, preliminary disease-specific axis, donor-level association, or exploratory cross-disease convergence. It must not be described as clinical validation, causal mechanism proof, definitive shared mechanism, or validated across diseases.

## Phase 22: Endpoint-Locked NeuroFate-Axis Validation

Phase 22 supersedes Phase 21 for PNAS-facing biological claims. Phase 21 remains useful for exploratory axis discovery, but its largest-effect-across-label comparison can mix non-equivalent endpoints. Phase 22 locks primary endpoints before testing and uses matched random-axis controls based on the same endpoint statistic.

Primary endpoints:

- SEA-AD AD endpoint: `label__Cognitive_Status`, Dementia versus No dementia.
- GSE243639 PD endpoint: `diagnosis`, Parkinson's versus Control.

Manual Phase 22 command sequence:

```bash
python scripts/112_test_axis_associations_endpoint_locked.py && \
python scripts/113_compare_ad_pd_axes_endpoint_locked.py && \
python scripts/114_run_endpoint_locked_random_axis_controls.py --n-random 1000 && \
python scripts/115_build_endpoint_locked_axis_evidence_table.py && \
python scripts/116_generate_phase22_endpoint_locked_figures.py && \
python scripts/117_generate_phase22_endpoint_locked_pnas_report.py
```

Use `results/tables/phase22_endpoint_locked_axis_evidence_table.tsv` and `results/reports/phase22_endpoint_locked_axis_claims.md` for PNAS-oriented axis claims. Do not claim clinical utility, diagnostic use, causal mechanisms, definitive shared mechanisms, or validation across diseases.

```bash
python scripts/100_audit_phase19_safe_annotation_map_schema.py --safe-map data/interim/external/gse243639_pd_snpc/gse243639_safe_cell_annotation_map.tsv --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz --output results/tables/phase20_safe_annotation_map_schema_audit.tsv --preview-output results/reports/phase20_safe_annotation_map_preview.tsv --log-file results/logs/100_audit_phase19_safe_annotation_map_schema.log
python scripts/101_build_gse243639_phase20_celltype_features.py --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz --annotations data/interim/external/gse243639_pd_snpc/gse243639_safe_cell_annotation_map.tsv --clinical data/raw/external/gse243639_pd_snpc/GSE243639_Clinical_data.csv.gz --phase5-schema results/tables/phase5_donor_feature_table.tsv --output results/tables/phase20_gse243639_celltype_feature_table.tsv --schema-output results/tables/phase20_gse243639_celltype_schema_alignment.tsv --label-summary-output results/tables/phase20_gse243639_celltype_label_summary.tsv --feature-group-output results/tables/phase20_gse243639_feature_group_counts.tsv --log-file results/logs/101_build_gse243639_phase20_celltype_features.log
python scripts/102_run_gse243639_phase20_celltype_pd_validation.py --features results/tables/phase20_gse243639_celltype_feature_table.tsv --metrics-output results/tables/phase20_gse243639_celltype_validation_metrics.tsv --predictions-output results/tables/phase20_gse243639_celltype_predictions.tsv --importance-output results/tables/phase20_gse243639_celltype_feature_importance.tsv --log-file results/logs/102_run_gse243639_phase20_celltype_pd_validation.log
python scripts/103_compare_phase16_17_18_20_pd_validation.py --output results/tables/phase20_pd_validation_comparison.tsv --summary-output results/reports/phase20_pd_repair_summary.md
python scripts/104_generate_phase20_gse243639_figures.py
python scripts/105_generate_phase20_gse243639_report.py
```

## Phase 1B: Safety and Metadata Registry Validation

Phase 1B registry design is manuscript-driven. The current scientific backbone is `manuscript/neurofate_landmark_manuscript.tex`, and the metadata registries mirror its major planned modules: single-cell transcriptomics, Alzheimer disease, Parkinson disease, gut-brain axis, microbiome/metabolite signatures, protein interaction/network biology, evolutionary conservation, positive selection, multimodal fate prediction, and interpretability/reporting.

### LIGHTWEIGHT: validate registries manually

```bash
python scripts/04_validate_registries.py --dry-run --log-file results/logs/04_validate_registries.log
```

### LIGHTWEIGHT: run Phase 1B registry tests manually

```bash
python -m pytest tests/test_dataset_registry.py tests/test_feature_registry.py tests/test_safety_flags.py tests/test_manuscript_module_map.py
```

Please paste back:

- terminal output,
- `results/logs/04_validate_registries.log`,
- any failed test names if pytest reports a failure.

## Phase 1C: Dataset Intake and Provenance Tracking

Phase 1C adds a manuscript-driven intake and provenance layer so future datasets can be traced from source to local file, biological module, feature layer, and planned result table/figure. These files are still templates: they do not grant access, fetch remote resources, verify real files, compute checksums, or process biological matrices.

### LIGHTWEIGHT: prepare or validate the dataset intake checklist

```bash
python scripts/05_prepare_dataset_intake_sheet.py --dry-run --log-file results/logs/05_prepare_dataset_intake_sheet.log
```

### LIGHTWEIGHT: validate provenance placeholders

```bash
python scripts/06_validate_provenance.py --dry-run --log-file results/logs/06_validate_provenance.log
```

### LIGHTWEIGHT: run Phase 1C tests manually

```bash
python -m pytest tests/test_intake_config.py tests/test_dataset_intake_checklist.py tests/test_provenance_template.py tests/test_provenance_validator.py
```

Please paste back:

- terminal output from both scripts,
- `results/logs/05_prepare_dataset_intake_sheet.log`,
- `results/logs/06_validate_provenance.log`,
- pytest output.

## Phase 1D: Dataset Selection Plan

Phase 1D maps each planned manuscript claim to the minimum dataset or dataset layer needed before any acquisition work begins. This is a planning layer only: it does not download datasets, access remote URLs, open biological files, or run analysis.

### LIGHTWEIGHT: validate the dataset selection plan

```bash
python scripts/07_validate_dataset_selection_plan.py --dry-run --log-file results/logs/07_validate_dataset_selection_plan.log
```

### LIGHTWEIGHT: run Phase 1D tests manually

```bash
python -m pytest tests/test_dataset_selection_plan.py
```

Please paste back:

- terminal output,
- `results/logs/07_validate_dataset_selection_plan.log`,
- pytest output.

## Phase 2A: Real Dataset Acquisition Plan

Phase 2A creates manual download templates and source metadata only. No download occurs automatically. Do not run these scripts from Codex.

SEA-AD processed single-nucleus RNA-seq from Allen/Brain Knowledge Platform/AWS Open Data is the first priority dataset. Mathys 2019 and ROSMAP are secondary AD resources and may require Synapse, controlled access, or data-use certification depending on the selected files. STRING human protein interaction links support the network layer. HMDB/KEGG/gut-brain metabolite signatures remain a manual curation placeholder.

Before any manual download, validate the local source plan:

```bash
python scripts/08_validate_real_dataset_sources.py --dry-run --log-file results/logs/08_validate_real_dataset_sources.log
```

Then review and edit the relevant guarded manual template:

```bash
less scripts/manual_downloads/download_sea_ad_manual.sh
```

Only after filling the official source URI from the official resource page should you run the SEA-AD template manually:

```bash
RUN_MANUAL_DOWNLOAD=YES SEA_AD_SOURCE_URI='s3://OFFICIAL_SEA_AD_BUCKET/OFFICIAL_PROCESSED_SNRNA_FILE' bash scripts/manual_downloads/download_sea_ad_manual.sh
```

The scripts intentionally keep the heavy command commented until you verify and uncomment exactly one official command.

Run Phase 2A tests:

```bash
python -m pytest tests/test_real_dataset_sources.py tests/test_manual_download_templates.py
```

Please paste back:

- terminal output from `scripts/08_validate_real_dataset_sources.py`,
- `results/logs/08_validate_real_dataset_sources.log`,
- pytest output,
- any official SEA-AD object URI you choose before attempting download.

## Phase 2C: SEA-AD Metadata-Only Extraction

Phase 2C prepares the first manuscript-ready SEA-AD metadata tables. This phase is metadata-only: it reads selected `obs` fields and `var` gene identifiers from the local H5AD file with explicit guards against `X` access. It does not run Scanpy, create an in-memory analysis object, read expression matrix arrays, normalize, cluster, run PCA/UMAP, or train models.

### METADATA-ONLY: extract SEA-AD obs/var metadata

```bash
python scripts/11_extract_sea_ad_metadata_only.py \
  --input data/raw/sea_ad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
  --outdir data/interim/sea_ad \
  --tables-dir results/tables \
  --log-file results/logs/11_extract_sea_ad_metadata_only.log
```

Expected outputs:

- `data/interim/sea_ad/sea_ad_obs_metadata_minimal.tsv`
- `data/interim/sea_ad/sea_ad_var_genes.tsv`
- `results/tables/sea_ad_metadata_summary.tsv`
- `results/tables/table1_sea_ad_cohort_cell_summary.tsv`

### METADATA-ONLY: summarize extracted SEA-AD metadata

```bash
python scripts/12_summarize_sea_ad_metadata.py \
  --metadata data/interim/sea_ad/sea_ad_obs_metadata_minimal.tsv \
  --tables-dir results/tables \
  --log-file results/logs/12_summarize_sea_ad_metadata.log
```

Expected outputs:

- `results/tables/sea_ad_donor_summary.tsv`
- `results/tables/sea_ad_celltype_by_ad_pathology.tsv`
- `results/tables/sea_ad_celltype_by_cognitive_status.tsv`

Please paste back:

- terminal output from both scripts,
- `results/logs/11_extract_sea_ad_metadata_only.log`,
- `results/logs/12_summarize_sea_ad_metadata.log`,
- the first 10 lines of each generated TSV.

## Phase 2D: Decode SEA-AD Categories and Table 1

Phase 2D decodes SEA-AD categorical metadata labels using `obs/__categories` from the local H5AD file and the metadata-only TSV from Phase 2C. This phase still does not read expression matrix arrays, run Scanpy, load an analysis object, normalize, cluster, run PCA/UMAP, or train models.

### METADATA-ONLY: decode SEA-AD categorical metadata

```bash
python scripts/13_decode_sea_ad_categories.py \
  --input data/raw/sea_ad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
  --metadata data/interim/sea_ad/sea_ad_obs_metadata_minimal.tsv \
  --outdir data/interim/sea_ad \
  --tables-dir results/tables \
  --log-file results/logs/13_decode_sea_ad_categories.log
```

Expected outputs:

- `data/interim/sea_ad/sea_ad_obs_metadata_decoded.tsv`
- `results/tables/table1_sea_ad_publication_ready.tsv`
- `results/tables/sea_ad_category_mapping.tsv`

Please paste back:

- terminal output,
- `results/logs/13_decode_sea_ad_categories.log`,
- the first 20 lines of `results/tables/table1_sea_ad_publication_ready.tsv`,
- the first 20 lines of `results/tables/sea_ad_category_mapping.tsv`.

## Phase 2E: Sparse Gene Extraction Planning

Phase 2E prepares a small targeted gene panel and guarded sparse extraction tooling. The planner reads only `sea_ad_var_genes.tsv` and does not open the H5AD file. The future extractor is manual-only, chunked, CSR-aware, and writes nonzero values for selected genes only. It does not create dense matrices, run Scanpy, normalize, cluster, run PCA/UMAP, or train models.

### LIGHTWEIGHT: plan target gene panel extraction

```bash
python scripts/14_plan_sparse_gene_extraction.py \
  --var data/interim/sea_ad/sea_ad_var_genes.tsv \
  --panel metadata/target_gene_panel_v1.tsv \
  --tables-dir results/tables \
  --log-file results/logs/14_sparse_gene_extraction_plan.log
```

Expected outputs:

- `results/tables/target_gene_panel_presence.tsv`
- `results/logs/manual_sparse_gene_extraction_template.sh`

### DRY RUN ONLY: validate future sparse extraction settings

```bash
python scripts/15_sparse_gene_extraction_safe.py \
  --input data/raw/sea_ad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
  --var data/interim/sea_ad/sea_ad_var_genes.tsv \
  --panel metadata/target_gene_panel_v1.tsv \
  --output data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz \
  --max-genes 64 \
  --m5-max-profile \
  --dry-run
```

### MANUAL EXPRESSION EXTRACTION: M5 Max high-memory profile

Run only after reviewing the planner output and dry run:

```bash
RUN_MANUAL_EXTRACTION=YES python scripts/15_sparse_gene_extraction_safe.py \
  --input data/raw/sea_ad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
  --var data/interim/sea_ad/sea_ad_var_genes.tsv \
  --panel metadata/target_gene_panel_v1.tsv \
  --output data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz \
  --max-genes 64 \
  --m5-max-profile \
  --execute
```

The M5 Max profile uses chunked sparse reads with up to 50,000 rows per chunk and a 32 GB memory cap. It still never densifies the matrix or reads the full expression matrix into memory.

## PHASE 3: Sparse Expression Disease Analysis

Phase 3 uses only the pre-extracted sparse gene panel table and decoded SEA-AD metadata. It performs chunked/tabular disease-associated summaries for the targeted 30-gene panel and generates first-pass manuscript figures from summary tables only. It does not read the H5AD file, load the full matrix, create dense matrices, run Scanpy, run PCA/UMAP/clustering, or train models.

### LIGHTWEIGHT/TABULAR: compute disease expression statistics

```bash
python scripts/16_compute_sparse_expression_statistics.py \
  --expression data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz \
  --metadata data/interim/sea_ad/sea_ad_obs_metadata_decoded.tsv \
  --panel metadata/target_gene_panel_v1.tsv \
  --tables-dir results/tables \
  --log-file results/logs/16_compute_sparse_expression_statistics.log
```

Expected outputs:

- `results/tables/gene_by_celltype_summary.tsv`
- `results/tables/gene_by_ad_pathology.tsv`
- `results/tables/gene_by_cognitive_status.tsv`
- `results/tables/microglial_activation_signature.tsv`
- `results/tables/astrocyte_stress_signature.tsv`
- `results/tables/neuronal_signature_summary.tsv`
- `results/tables/neurodegeneration_signature_summary.tsv`

### LIGHTWEIGHT/TABULAR: generate manuscript figures

```bash
python scripts/17_generate_phase3_figures.py \
  --tables-dir results/tables \
  --figures-dir results/figures
```

Expected outputs:

- `results/figures/figure1_celltype_composition.png`
- `results/figures/figure2_microglial_activation.png`
- `results/figures/figure3_neurodegeneration_signatures.png`
- `results/figures/figure4_ad_pathology_gene_trends.png`

### LIGHTWEIGHT/TABULAR: draft Phase 3 results text

```bash
python scripts/18_generate_phase3_results_text.py \
  --tables-dir results/tables \
  --output results/tables/phase3_results_summary.txt
```

Expected output:

- `results/tables/phase3_results_summary.txt`

Please paste back:

- terminal output from the three Phase 3 commands,
- `results/logs/16_compute_sparse_expression_statistics.log`,
- the first 20 lines of each Phase 3 TSV table,
- `results/tables/phase3_results_summary.txt`,
- the generated figure files if visual review is needed.

## PHASE 4: Statistical Neurodegeneration Biology

Phase 4 converts the sparse target-gene table and decoded SEA-AD metadata into statistically testable neurodegeneration findings. It uses donor-aware summaries where feasible, rank-based statistics, Benjamini-Hochberg correction, confidence intervals, APOE genotype comparisons, mixed-pathology comparisons, and composite biological indices. It still does not read the H5AD file, load the full expression matrix, create dense matrices, run Scanpy/scVI/scVelo, run PCA/UMAP/clustering, or train models.

### LIGHTWEIGHT/TABULAR: compute Phase 4 statistical tables

```bash
python scripts/19_compute_phase4_statistics.py \
  --expression data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz \
  --metadata data/interim/sea_ad/sea_ad_obs_metadata_decoded.tsv \
  --tables-dir results/tables \
  --log-file results/logs/19_compute_phase4_statistics.log
```

Expected outputs:

- `results/tables/phase4_gene_statistics.tsv`
- `results/tables/phase4_celltype_vulnerability.tsv`
- `results/tables/phase4_apoe_analysis.tsv`
- `results/tables/phase4_mixed_pathology.tsv`
- `results/tables/phase4_composite_indices.tsv`

### LIGHTWEIGHT/TABULAR: generate Phase 4 manuscript figures

```bash
python scripts/20_generate_phase4_figures.py \
  --tables-dir results/tables \
  --figures-dir results/figures
```

Expected outputs:

- `results/figures/figure5_braak_associations.png`
- `results/figures/figure6_apoe_microglia.png`
- `results/figures/figure7_celltype_vulnerability_heatmap.png`
- `results/figures/figure8_composite_indices.png`

### LIGHTWEIGHT/TABULAR: draft Phase 4 results text

```bash
python scripts/21_generate_phase4_results_text.py \
  --tables-dir results/tables \
  --output results/tables/phase4_results_summary.txt
```

Expected output:

- `results/tables/phase4_results_summary.txt`

Please paste back:

- terminal output from the three Phase 4 commands,
- `results/logs/19_compute_phase4_statistics.log`,
- the first 20 lines of each Phase 4 TSV table,
- `results/tables/phase4_results_summary.txt`,
- the generated figure files if visual review is needed.

## PHASE 5: Predictive Neurodegeneration Modeling

Phase 5 builds donor-level features and runs lightweight interpretable predictive models. The workflow uses only the sparse target-gene expression table, decoded metadata TSV, and Phase 3/4 output tables. Donor aggregation happens before modeling, labels are excluded from features, train/test splitting is donor-level, and random seeds are fixed for reproducibility. It does not read the H5AD file, load the full expression matrix, create dense single-cell matrices, run Scanpy/scVI/scVelo, run PCA/UMAP/clustering, or train deep-learning models.

### LIGHTWEIGHT/TABULAR: build donor-level feature table

```bash
python scripts/22_build_donor_feature_table.py \
  --expression data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz \
  --metadata data/interim/sea_ad/sea_ad_obs_metadata_decoded.tsv \
  --output results/tables/phase5_donor_feature_table.tsv \
  --log-file results/logs/22_build_donor_feature_table.log
```

Expected output:

- `results/tables/phase5_donor_feature_table.tsv`

### LIGHTWEIGHT/DONOR-LEVEL: run interpretable predictive models

```bash
python scripts/23_run_phase5_models.py \
  --features results/tables/phase5_donor_feature_table.tsv \
  --tables-dir results/tables \
  --log-file results/logs/23_run_phase5_models.log
```

Expected outputs:

- `results/tables/phase5_model_metrics.tsv`
- `results/tables/phase5_feature_importance.tsv`
- `results/tables/phase5_neurofate_scores.tsv`

The NeuroFate Neurodegeneration Risk Score is defined as the mean out-of-fold predicted risk probability across available donor-level Phase 5 tasks.

### LIGHTWEIGHT/TABULAR: generate Phase 5 manuscript figures

```bash
python scripts/24_generate_phase5_figures.py \
  --tables-dir results/tables \
  --figures-dir results/figures
```

Expected outputs:

- `results/figures/figure9_model_performance.png`
- `results/figures/figure10_feature_importance.png`
- `results/figures/figure11_neurofate_score_distribution.png`
- `results/figures/figure12_donor_risk_heatmap.png`

### LIGHTWEIGHT/TABULAR: draft Phase 5 results text

```bash
python scripts/25_generate_phase5_results_text.py \
  --tables-dir results/tables \
  --output results/tables/phase5_results_summary.txt
```

Expected output:

- `results/tables/phase5_results_summary.txt`

Please paste back:

- terminal output from the four Phase 5 commands,
- `results/logs/22_build_donor_feature_table.log`,
- `results/logs/23_run_phase5_models.log`,
- the first 20 lines of each Phase 5 TSV table,
- `results/tables/phase5_results_summary.txt`,
- the generated figure files if visual review is needed.

## PHASE 6: Apple Silicon Metal Neural NeuroFate Model

Phase 6 trains a small donor-level NeuroFate MLP with PyTorch and Apple Silicon Metal/MPS acceleration when available. It uses only `results/tables/phase5_donor_feature_table.tsv`; it does not read single-cell files, load expression matrices, run Scanpy/scVI/scVelo, run PCA/UMAP/clustering, use transformers, or train a large model.

### LIGHTWEIGHT: check PyTorch MPS device

```bash
python scripts/26_check_mps_device.py
```

Expected output:

- torch version,
- MPS build status,
- MPS availability,
- selected device,
- small tensor test result.

### DONOR-LEVEL NEURAL MODEL: train NeuroFate MPS MLP

```bash
python scripts/27_train_neurofate_mps_model.py \
  --config configs/neurofate_mps_model_config.yaml \
  --features results/tables/phase5_donor_feature_table.tsv \
  --all-tasks \
  --tables-dir results/tables \
  --models-dir results/models \
  --log-file results/logs/27_train_neurofate_mps_model.log
```

Expected outputs:

- `results/models/neurofate_mps_<task>.pt`
- `results/tables/phase6_mps_model_metrics.tsv`
- `results/tables/phase6_mps_training_log.tsv`
- `results/tables/phase6_mps_predictions.tsv`

### LIGHTWEIGHT/TABULAR: generate Phase 6 manuscript figures

```bash
python scripts/28_generate_phase6_figures.py \
  --tables-dir results/tables \
  --figures-dir results/figures
```

Expected outputs:

- `results/figures/figure13_mps_model_performance.png`
- `results/figures/figure14_mps_training_curves.png`
- `results/figures/figure15_mps_prediction_distribution.png`

### LIGHTWEIGHT/TABULAR: draft Phase 6 results text

```bash
python scripts/29_generate_phase6_results_text.py \
  --tables-dir results/tables \
  --output results/tables/phase6_results_summary.txt
```

Expected output:

- `results/tables/phase6_results_summary.txt`

Please paste back:

- terminal output from the MPS check,
- terminal output from the training command,
- `results/logs/27_train_neurofate_mps_model.log`,
- the first 20 lines of each Phase 6 TSV table,
- `results/tables/phase6_results_summary.txt`,
- the generated Phase 6 figures if visual review is needed.

## PHASE 7: Cross-Cohort External Validation

Phase 7 prepares NeuroFate for external validation across AD, PD, and optional reference cohorts. It introduces an external validation registry, gene/metadata overlap planning, guarded manual sparse extraction templates, donor-level feature harmonization, and cross-cohort validation modes. External downloads, extraction, and modeling remain manual. The workflow does not run Scanpy/scVI/scVelo, does not load full H5AD matrices, does not create dense single-cell matrices, and keeps validation donor-level whenever feasible.

### LIGHTWEIGHT: prepare external validation plan

```bash
python scripts/30_prepare_external_validation_plan.py \
  --registry metadata/external_validation_registry.tsv \
  --panel metadata/target_gene_panel_v1.tsv \
  --tables-dir results/tables \
  --log-file results/logs/30_prepare_external_validation_plan.log
```

Expected outputs:

- `results/tables/external_validation_gene_overlap.tsv`
- `results/tables/external_validation_metadata_overlap.tsv`

### MANUAL_HEAVY TEMPLATE: external sparse target-gene extraction

Run only after an external cohort has been downloaded manually and a lightweight external var metadata TSV is available:

```bash
RUN_MANUAL_EXTERNAL_EXTRACTION=YES python scripts/31_sparse_external_gene_extraction.py \
  --dataset-id mathys_2019_ad \
  --input data/raw/external/mathys_2019/example_external.h5ad \
  --var data/interim/external/mathys_2019/var_genes.tsv \
  --panel metadata/target_gene_panel_v1.tsv \
  --output data/interim/external/mathys_2019/sparse_gene_panel_expression.tsv.gz \
  --log-file results/logs/31_sparse_external_gene_extraction_mathys.log \
  --execute
```

### LIGHTWEIGHT/DONOR-LEVEL: build harmonized cross-cohort feature tables

```bash
python scripts/32_build_crosscohort_feature_tables.py \
  --feature-table sea_ad=results/tables/phase5_donor_feature_table.tsv \
  --feature-table mathys_2019_ad=results/tables/mathys_2019_phase5_donor_feature_table.tsv \
  --output results/tables/crosscohort_donor_feature_table.tsv \
  --overlap-output results/tables/crosscohort_feature_overlap.tsv \
  --log-file results/logs/32_build_crosscohort_feature_tables.log
```

Expected outputs:

- `results/tables/crosscohort_donor_feature_table.tsv`
- `results/tables/crosscohort_feature_overlap.tsv`

### DONOR-LEVEL MODEL VALIDATION: cross-cohort generalization

```bash
python scripts/33_run_crosscohort_validation.py \
  --features results/tables/crosscohort_donor_feature_table.tsv \
  --tables-dir results/tables \
  --log-file results/logs/33_run_crosscohort_validation.log
```

Validation modes:

- train SEA-AD, test external cohort,
- leave-one-cohort-out,
- pooled multi-cohort training.

Expected outputs:

- `results/tables/phase7_crosscohort_metrics.tsv`
- `results/tables/phase7_generalization_summary.tsv`

### LIGHTWEIGHT/TABULAR: generate Phase 7 manuscript figures

```bash
python scripts/34_generate_phase7_figures.py \
  --tables-dir results/tables \
  --figures-dir results/figures
```

Expected outputs:

- `results/figures/figure16_crosscohort_generalization.png`
- `results/figures/figure17_cohort_transfer_performance.png`
- `results/figures/figure18_feature_stability.png`
- `results/figures/figure19_multicohort_neurofate_scores.png`

### LIGHTWEIGHT/TABULAR: draft Phase 7 results text

```bash
python scripts/35_generate_phase7_results_text.py \
  --tables-dir results/tables \
  --output results/tables/phase7_results_summary.txt
```

Expected output:

- `results/tables/phase7_results_summary.txt`

Please paste back:

- terminal output from the Phase 7 commands you run,
- `results/logs/30_prepare_external_validation_plan.log`,
- `results/logs/32_build_crosscohort_feature_tables.log`,
- `results/logs/33_run_crosscohort_validation.log`,
- the first 20 lines of each Phase 7 TSV table,
- `results/tables/phase7_results_summary.txt`,
- generated Phase 7 figures if visual review is needed.

## PHASE 8: Mathys 2019 External Validation Cohort

Phase 8 onboards Mathys et al. 2019 AD snRNA-seq (`GSE138852`) as the first real external NeuroFate validation cohort. The workflow prepares manual GEO acquisition, metadata-only inspection, target-gene overlap, sparse extraction planning, and future donor-level feature harmonization. It does not download automatically, does not process full matrices automatically, does not run Scanpy/scVI/scVelo, and does not create dense matrices.

### MANUAL_HEAVY: review Mathys GEO acquisition template

```bash
RUN_MANUAL_DOWNLOAD=YES bash scripts/manual_downloads/download_mathys2019_geo_manual.sh
```

GEO references:

- `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138852`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/`

Manual download templates to adapt after reviewing the GEO supplementary file names:

```bash
curl -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/FILE_NAME" \
  -o data/raw/external/mathys_2019/FILE_NAME

wget -O data/raw/external/mathys_2019/FILE_NAME \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138852/suppl/FILE_NAME"

md5 data/raw/external/mathys_2019/FILE_NAME > metadata/checksums/mathys2019_FILE_NAME.md5
shasum -a 256 data/raw/external/mathys_2019/FILE_NAME > metadata/checksums/mathys2019_FILE_NAME.sha256
```

Expected local structure:

- `data/raw/external/mathys_2019/`
- `data/interim/external/mathys_2019/`
- `results/tables/mathys2019_metadata_summary.tsv`
- `data/interim/external/mathys_2019/mathys_var_genes.tsv`

### METADATA-ONLY: inspect Mathys H5AD metadata

```bash
python scripts/36_inspect_external_h5ad_metadata.py \
  --input data/raw/external/mathys_2019/mathys2019_external.h5ad \
  --output results/tables/mathys2019_metadata_summary.tsv \
  --var-output data/interim/external/mathys_2019/mathys_var_genes.tsv \
  --log-file results/logs/36_inspect_mathys2019.log
```

Expected outputs:

- `results/tables/mathys2019_metadata_summary.tsv`
- `data/interim/external/mathys_2019/mathys_var_genes.tsv`
- `results/logs/36_inspect_mathys2019.log`

### LIGHTWEIGHT: compare Mathys genes with NeuroFate target panel

```bash
python scripts/37_prepare_mathys_gene_panel_overlap.py \
  --var data/interim/external/mathys_2019/mathys_var_genes.tsv \
  --panel metadata/target_gene_panel_v1.tsv \
  --overlap-output results/tables/mathys_gene_overlap.tsv \
  --missing-output results/tables/mathys_missing_target_genes.tsv \
  --log-file results/logs/37_prepare_mathys_gene_panel_overlap.log
```

Expected outputs:

- `results/tables/mathys_gene_overlap.tsv`
- `results/tables/mathys_missing_target_genes.tsv`

### LIGHTWEIGHT: prepare Mathys sparse extraction plan

```bash
python scripts/38_prepare_mathys_sparse_extraction.py \
  --overlap results/tables/mathys_gene_overlap.tsv \
  --input data/raw/external/mathys_2019/mathys2019_external.h5ad \
  --var data/interim/external/mathys_2019/mathys_var_genes.tsv \
  --output data/interim/external/mathys_2019/sparse_gene_panel_expression.tsv.gz \
  --plan-output results/tables/mathys_sparse_extraction_plan.tsv \
  --manual-script-output results/logs/manual_mathys_sparse_extraction_template.sh \
  --log-file results/logs/38_prepare_mathys_sparse_extraction.log
```

Expected outputs:

- `results/tables/mathys_sparse_extraction_plan.tsv`
- `results/logs/manual_mathys_sparse_extraction_template.sh`

### FUTURE MANUAL_HEAVY: run Mathys sparse extraction

Run only after reviewing the plan and local paths:

```bash
RUN_MANUAL_EXTERNAL_EXTRACTION=YES bash results/logs/manual_mathys_sparse_extraction_template.sh
```

### DONOR-LEVEL: build future Mathys NeuroFate feature table

```bash
python scripts/39_build_mathys_donor_feature_table.py \
  --metadata data/interim/external/mathys_2019/mathys_obs_metadata_decoded.tsv \
  --expression data/interim/external/mathys_2019/sparse_gene_panel_expression.tsv.gz \
  --donor-field donor_id \
  --celltype-field cell_type \
  --output results/tables/mathys_2019_phase5_donor_feature_table.tsv \
  --schema-output results/tables/mathys_2019_feature_schema_alignment.tsv \
  --phase5-schema results/tables/phase5_donor_feature_table.tsv \
  --log-file results/logs/39_build_mathys_donor_feature_table.log
```

Expected outputs:

- `results/tables/mathys_2019_phase5_donor_feature_table.tsv`
- `results/tables/mathys_2019_feature_schema_alignment.tsv`

Please paste back:

- terminal output from the metadata inspection and gene-overlap commands,
- `results/logs/36_inspect_mathys2019.log`,
- first 40 lines of `results/tables/mathys2019_metadata_summary.tsv`,
- first 40 lines of `results/tables/mathys_gene_overlap.tsv`,
- first 40 lines of `results/tables/mathys_missing_target_genes.tsv`.

## PHASE 9: Mathys 2019 CSV External Validation

Phase 9 adapts NeuroFate to the real Mathys `GSE138852` GEO files:

- `data/raw/external/mathys_2019/GSE138852_counts.csv.gz`
- `data/raw/external/mathys_2019/GSE138852_covariates.csv.gz`

The workflow is CSV-native and does not use Scanpy, AnnData, H5AD conversion, PCA/UMAP/clustering, or dense single-cell matrix expansion. The count table is inspected for orientation, the NeuroFate 30-gene panel is extracted into a sparse-like long TSV, Mathys covariates are mapped to sample-level labels, and the resulting donor/sample feature table is aligned to the SEA-AD Phase 5 schema for external validation.

### LIGHTWEIGHT/CSV: inspect Mathys count and covariate structure

```bash
python scripts/40_inspect_mathys_csv_structure.py \
  --counts data/raw/external/mathys_2019/GSE138852_counts.csv.gz \
  --covariates data/raw/external/mathys_2019/GSE138852_covariates.csv.gz \
  --panel metadata/target_gene_panel_v1.tsv \
  --summary-output results/tables/mathys_csv_structure_summary.tsv \
  --preview-output data/interim/external/mathys_2019/mathys_covariates_preview.tsv \
  --log-file results/logs/40_inspect_mathys_csv_structure.log
```

Expected outputs:

- `results/tables/mathys_csv_structure_summary.tsv`
- `data/interim/external/mathys_2019/mathys_covariates_preview.tsv`

### CSV/SPARSE-LIKE: extract Mathys NeuroFate target genes

```bash
python scripts/41_extract_mathys_target_gene_panel.py \
  --counts data/raw/external/mathys_2019/GSE138852_counts.csv.gz \
  --panel metadata/target_gene_panel_v1.tsv \
  --output data/interim/external/mathys_2019/mathys_sparse_gene_panel_expression.tsv.gz \
  --orientation auto \
  --log-file results/logs/41_extract_mathys_target_gene_panel.log
```

The extractor supports both orientations:

- genes as rows and cells as columns,
- cells as rows and genes as columns.

Expected output:

- `data/interim/external/mathys_2019/mathys_sparse_gene_panel_expression.tsv.gz`

### DONOR/SAMPLE-LEVEL: build Mathys feature table

```bash
python scripts/42_build_mathys_feature_table.py \
  --expression data/interim/external/mathys_2019/mathys_sparse_gene_panel_expression.tsv.gz \
  --covariates data/raw/external/mathys_2019/GSE138852_covariates.csv.gz \
  --phase5-schema results/tables/phase5_donor_feature_table.tsv \
  --output results/tables/mathys_2019_phase5_donor_feature_table.tsv \
  --schema-output results/tables/mathys_2019_feature_schema_alignment.tsv \
  --label-summary-output results/tables/mathys_2019_label_summary.tsv \
  --log-file results/logs/42_build_mathys_feature_table.log
```

Covariate mapping:

- `cell_id`: first covariate column,
- diagnosis: `oupSample.batchCond` or `oupSample.subclustCond`,
- cell type: `oupSample.cellType`,
- subcluster: `oupSample.subclustID`,
- sample ID: inferred from cell barcode suffix where possible, otherwise a clearly logged pseudo-donor fallback.

Expected outputs:

- `results/tables/mathys_2019_phase5_donor_feature_table.tsv`
- `results/tables/mathys_2019_feature_schema_alignment.tsv`
- `results/tables/mathys_2019_label_summary.tsv`

### DONOR/SAMPLE-LEVEL: run Mathys external validation

```bash
python scripts/43_run_mathys_external_validation.py \
  --sea-ad-features results/tables/phase5_donor_feature_table.tsv \
  --mathys-features results/tables/mathys_2019_phase5_donor_feature_table.tsv \
  --metrics-output results/tables/phase9_mathys_external_validation_metrics.tsv \
  --predictions-output results/tables/phase9_mathys_external_predictions.tsv \
  --log-file results/logs/43_run_mathys_external_validation.log
```

Validation modes:

- train SEA-AD baseline, test Mathys,
- Mathys internal train/test diagnostic if enough sample units exist.

Expected outputs:

- `results/tables/phase9_mathys_external_validation_metrics.tsv`
- `results/tables/phase9_mathys_external_predictions.tsv`

### LIGHTWEIGHT/TABULAR: generate Phase 9 figures

```bash
python scripts/44_generate_phase9_figures.py \
  --tables-dir results/tables \
  --figures-dir results/figures
```

Expected outputs:

- `results/figures/figure20_mathys_gene_overlap.png`
- `results/figures/figure21_mathys_external_validation.png`
- `results/figures/figure22_mathys_celltype_composition.png`

### LIGHTWEIGHT/TABULAR: draft Phase 9 results text

```bash
python scripts/45_generate_phase9_results_text.py \
  --tables-dir results/tables \
  --output results/tables/phase9_results_summary.txt
```

Expected output:

- `results/tables/phase9_results_summary.txt`

Please paste back:

- terminal output from all Phase 9 commands you run,
- `results/logs/40_inspect_mathys_csv_structure.log`,
- `results/logs/41_extract_mathys_target_gene_panel.log`,
- `results/logs/42_build_mathys_feature_table.log`,
- `results/logs/43_run_mathys_external_validation.log`,
- first 40 lines of each Phase 9 TSV table,
- `results/tables/phase9_results_summary.txt`.

## Manual Command Templates

### LIGHTWEIGHT: prepare dataset registry preview

```bash
python scripts/01_prepare_dataset_registry.py --dry-run --datasets configs/datasets.yaml --output metadata/dataset_registry.tsv --log-file results/logs/01_prepare_dataset_registry.log
```

### LIGHTWEIGHT: write dataset registry from config

```bash
python scripts/01_prepare_dataset_registry.py --write --datasets configs/datasets.yaml --output metadata/dataset_registry.tsv --log-file results/logs/01_prepare_dataset_registry.log
```

### LIGHTWEIGHT: validate declared local input paths

```bash
python scripts/02_validate_input_files.py --dry-run --registry metadata/dataset_registry.tsv --log-file results/logs/02_validate_input_files.log
```

### LIGHTWEIGHT: preview a gene universe from small text inputs

```bash
python scripts/03_make_gene_universe.py --dry-run --input metadata/feature_registry.tsv --output data/interim/gene_universe.tsv --log-file results/logs/03_make_gene_universe.log
```

### HEAVY: future single-cell processing placeholder

Do not run until a later project step explicitly defines inputs, memory limits, and expected outputs.

```bash
# HEAVY - TEMPLATE ONLY
# python -m neurofate.singlecell --input data/raw/example.h5ad --output data/processed/example_processed.h5ad
```

### HEAVY: future model training placeholder

Do not run until model design, dataset splits, and compute plan are approved.

```bash
# HEAVY - TEMPLATE ONLY
# python -m neurofate.models --train --config configs/model_config.yaml
```

## Repository Layout

- `configs/`: project, dataset, model, and intake configuration.
- `data/`: raw, interim, processed, and external data locations. Large files are ignored by git.
- `metadata/`: lightweight TSV registries, including the manuscript module map, provenance template, and intake checklist.
- `scripts/`: safe setup and validation scripts.
- `neurofate/`: package modules with placeholder APIs.
- `notebooks/`: project overview notebook with no executed analysis.
- `results/`: logs, figures, tables, and model output directories.
- `manuscript/`: LaTeX manuscript stub and references.
- `tests/`: lightweight import/config tests.

## Current Status

No analysis has been executed. The repository is ready for the first manual system check.

## Phase 23: Independent Replication Cohort Onboarding

Phase 23 prepares the replication layer needed before PNAS-level biological claims. Phase 22 provides endpoint-locked candidate axes, but candidate axes are not enough for strong claims without independent AD and PD replication.

Replication priorities:

- PD priority 1: `GSE184950`, human substantia nigra single-cell/single-nucleus transcriptomics in Parkinson disease.
- AD priority 1: `GSE174367`, AD multi-omics with bulk/sample-level data prioritized first if usable.
- AD priority 2: `GSE147528`, AD progression snRNA-seq across caudal entorhinal cortex and superior frontal gyrus.

Manual planning commands:

```bash
python scripts/118_triage_replication_cohort_files.py
python scripts/123_build_pnas_readiness_matrix.py
```

Downloads remain manual only through guarded templates in `scripts/manual_downloads/`. Phase 23 does not run SRA tools, does not process raw SRA, does not load H5AD, and does not upgrade claims until replication statistics exist.

## Phase 24: GSE184950 PD Replication Onboarding

Phase 24 focuses on GSE184950 as the first independent PD replication cohort. The local `GSE184950_add2.xlsx` workbook is metadata only; it confirms sample fields and processed per-sample file names but does not contain expression values. Axis replication requires manual acquisition and inspection of `GSE184950_RAW.tar`.

Phase 24 workflow:

```bash
python scripts/124_parse_gse184950_geo_metadata_workbook.py
# Manual download only:
# RUN_MANUAL_DOWNLOAD=YES bash scripts/manual_downloads/download_gse184950_raw_manual.sh
python scripts/125_list_gse184950_raw_archive.py
python scripts/126_plan_gse184950_processed_matrix_extraction.py
```

Prefer processed 10x matrices if present in the RAW archive. FASTQ/SRA processing is avoided unless absolutely necessary and is not run by NeuroFate/Codex.

## Phase 25: GSE184950 Series-Matrix Metadata and Replication Planning

Phase 25 replaces the incomplete `GSE184950_add2.xlsx` workbook metadata with the GEO series matrix as the primary GSE184950 metadata source. The series matrix provides all 34 sample labels, including 10 Unaffected Control, 6 Parkinson's Disease, and 18 Parkinson's Disease Dementia samples.

Primary replication endpoint:

- `label__pd_pdd_vs_control`: Parkinson's Disease and Parkinson's Disease Dementia are positive; Unaffected Control is negative.

Safe planning commands:

```bash
python scripts/131_parse_gse184950_series_matrix.py
python scripts/132_reconcile_gse184950_archive_with_series_metadata.py
python scripts/133_plan_gse184950_selective_tar_extraction.py
```

Phase 25 still does not extract the RAW archive, process FASTQ/SRA files, run Scanpy, create H5AD/AnnData objects, recompute UMAP, cluster cells, train models, or load dense matrices. It prepares a selective processed-matrix route for future sample-level NeuroFate-Axis replication.

## Phase 26: GSE184950 PD Replication Execution

Phase 26 adds the guarded execution route for GSE184950. It first inspects nested per-sample archives inside `GSE184950_RAW.tar` without extracting files, then allows a user-confirmed selective extraction of processed 10x matrix files only if they are present.

Safe commands:

```bash
python scripts/135_inspect_gse184950_nested_archives.py
RUN_MANUAL_GSE184950_EXTRACTION=YES python scripts/136_extract_gse184950_processed_matrices_selective.py --execute
python scripts/127_extract_gse184950_axis_genes_from_10x.py --run-manual-extraction YES --audit-output results/tables/phase26_gse184950_axis_gene_extraction_audit.tsv
python scripts/128_build_gse184950_axis_scores.py
python scripts/129_test_gse184950_axis_replication.py
```

The replication endpoint remains sample-level PD/PDD versus Unaffected Control. FASTQ processing, SRA tools, Scanpy, H5AD/AnnData creation, clustering, UMAP, dense full-matrix loading, deep learning, and clinical/diagnostic claims remain out of scope.

## Phase 27: GSE184950 Clean Replication and Conservative Evidence Integration

Phase 27 repairs the GSE184950 replication outputs by removing non-sample rows such as `processed_matrices` from sample-level axis-score tables. Clean GSE184950 outputs must contain exactly 34 biological samples, with 24 PD/PDD positives and 10 Unaffected Control negatives.

Phase 27 also makes evidence integration more conservative. Directionally consistent axis effects are not treated as replication unless they have statistical support (`p < 0.05` or `FDR < 0.1`) with adequate sample and class counts. If GSE184950 effects remain weak by FDR, the cohort supports independent PD replication feasibility and infrastructure, not validated AD/PD shared biology.

## Phase 28: Independent AD Replication Onboarding

Phase 28 addresses the next major PNAS bottleneck: independent AD replication. SEA-AD provides strong internal AD axis evidence, but an additional AD cohort is required before making strong biological replication claims.

Priority cohorts:

- `GSE174367`: first-priority AD multi-omics cohort, with bulk/sample-level matrices preferred first.
- `GSE147528`: secondary AD progression snRNA-seq cohort.
- `GSE157827`: optional backup AD single-nucleus cohort.

Safe onboarding commands:

```bash
python scripts/144_triage_ad_replication_files.py
# After manual review and acquisition of a GEO series matrix:
python scripts/143_parse_geo_series_matrix_generic.py --series-matrix <series_matrix.txt.gz> --cohort-id gse174367_ad_multiomics
```

Downloads remain guarded manual templates only. Phase 28 does not process SRA/FASTQ, run Scanpy, create H5AD/AnnData, run UMAP/clustering, load dense matrices, or train models.

## Phase 29: GSE174367 Bulk RNA Independent AD Replication

Phase 29 adds the first executable independent AD replication lane using `GSE174367_bulkRNA_processed.rda.gz`. The parsed GSE174367 series matrix contains 230 samples, including 118 AD and 112 Control records. The bulk RNA file is used before the larger single-nucleus resources because NeuroFate-Axis is donor/sample-level and the immediate PNAS bottleneck is independent AD replication.

Safe Phase 29 sequence:

```bash
python scripts/149_inspect_gse174367_bulk_rda.py
python scripts/150_convert_gse174367_bulk_rda_to_axis_matrix.py
python scripts/151_build_gse174367_bulk_axis_scores.py
python scripts/152_test_gse174367_bulk_ad_axis_replication.py
```

The conversion step writes only NeuroFate axis genes, not a genome-wide converted matrix. If the RDA structure or sample mapping is ambiguous, it stops for manual review. Phase 29 does not read the single-nucleus H5 matrix, process single-nucleus expression, run Scanpy, create AnnData/H5AD, run UMAP/clustering, train models, or make unsupported clinical-use, diagnosis-use, cause-effect, or journal-readiness claims.

## Phase 32: Cross-Cohort Evidence Consolidation

Phase 32 consolidates endpoint-locked NeuroFate-Axis evidence after successful GSE174367 bulk-RNA Ensembl mapping and AD replication. The strongest consolidated finding is the `neuronal_vulnerability_axis`, which shows directionally consistent SEA-AD discovery support and nominal independent AD replication in GSE174367.

Current interpretation:

- SEA-AD remains the strong AD discovery anchor.
- GSE174367 adds nominal independent AD replication, strongest for neuronal vulnerability.
- GSE243639 remains a preliminary PD extension.
- GSE184950 provides clean PD replication infrastructure, but current axis-level replication is weak or direction-only.

PNAS readiness is closer, but not complete. The major remaining biological bottleneck is stronger independent PD axis replication.

## PNAS Submission Package Status

The current manuscript route is now framed for a conservative PNAS-style biological paper rather than an overextended shared-mechanism claim. The defensible center is:

- endpoint-locked SEA-AD AD axis discovery,
- nominal independent AD replication of the `neuronal_vulnerability_axis` in GSE174367 bulk RNA,
- preliminary GSE243639 PD convergence,
- clean but weak GSE184950 PD/PDD replication infrastructure,
- explicit claim boundaries blocking clinical, diagnostic, causal, and definitive shared AD/PD language.

Use:

```bash
python scripts/159_build_crosscohort_axis_evidence_summary.py
python scripts/123_build_pnas_readiness_matrix.py
python scripts/160_generate_phase32_pnas_decision_report.py
python scripts/162_build_pnas_submission_package_report.py
```

The official PNAS-template manuscript is `manuscript/Research_report.tex`. Stronger shared AD/PD biology still requires statistically supported independent PD axis replication.

## Phase 33: PD Replication Expansion

Phase 33 addresses the current PNAS bottleneck directly: AD replication is nominally supported, but PD axis replication remains weak. The new lane prioritizes public donor/sample-level PD expression cohorts that can be analyzed without FASTQ/SRA processing:

- `GSE20141`: first-priority laser-captured substantia nigra pars compacta neuron cohort for testing the `neuronal_vulnerability_axis`.
- `GSE20186`: PD expression superseries to triage for usable processed subcohorts.
- `GSE7621`, `GSE8397`, and `GSE20292`: backup substantia nigra bulk/microarray cohorts.
- `GSE157783`: optional small midbrain snRNA cohort only if processed matrices are available.

Safe first command:

```bash
python scripts/162_triage_pd_replication_geo_files.py
```

After manual acquisition of a series matrix and any needed platform annotation, use:

```bash
python scripts/163_build_axis_scores_from_geo_series_matrix.py --series-matrix <GSE_series_matrix.txt.gz> --cohort-id <cohort_id>
python scripts/165_test_pd_axis_replication_microarray.py --axis-scores results/tables/phase33_<cohort_id>_axis_scores.tsv --cohort-id <cohort_id> --output results/tables/phase33_<cohort_id>_pd_axis_replication_statistics.tsv --fdr-output results/tables/phase33_<cohort_id>_pd_axis_replication_fdr.tsv --log-file results/logs/165_<cohort_id>_pd_axis_replication.log
```

Phase 33 extracts only NeuroFate axis genes and keeps outputs sample-level. It does not download data automatically, process FASTQ/SRA, use Scanpy, create H5AD/AnnData, run UMAP/clustering, write dense genome-wide outputs, train models, or make clinical/diagnostic/causal/shared-mechanism claims.

## Phase 34: PD Microarray/Bulk Replication Expansion

Phase 34 makes the PD replication route more explicit for small public substantia nigra, SNpc, LCM, and microarray/bulk cohorts. It prioritizes:

- `GSE20141`: first manual target; 10 PD and 8 controls; laser-dissected SNpc neuron-focused expression.
- `GSE7621`: second target; reported 16 PD and 9 control substantia nigra replicates.
- `GSE8397`: region-aware PD/control brain expression, used only if substantia nigra labels can be separated cleanly.
- `GSE20186`: superseries/subseries route if metadata and platform mapping are clear.

Safe first metadata acquisition command for GSE20141:

```bash
RUN_MANUAL_DOWNLOAD=YES scripts/manual_downloads/download_gse20141_manual.sh
```

Then triage locally:

```bash
python scripts/162_triage_phase34_pd_geo_files.py
```

Phase 34 uses GEO series matrix, SOFT/MINiML metadata, platform annotations, or processed expression tables. It keeps only NeuroFate axis probes/genes, keeps all outputs sample-level, and does not run FASTQ/SRA processing or single-cell workflows.

## Phase 37: GSE7621 PD Replication

Phase 37 adds `GSE7621` as the next independent PD replication attempt after GSE20141 produced direction-only, non-significant support. GSE7621 is treated as a donor/sample-level substantia nigra bulk/microarray cohort.

The Phase 37 route:

- manually acquires the GEO series matrix and any required platform annotation;
- parses sample metadata into unambiguous PD/control labels where possible;
- audits expression-column to metadata joins before scoring;
- maps platform probes only to NeuroFate axis genes;
- builds sample-level axis scores without writing genome-wide expression output;
- tests endpoint-locked PD versus Control axis replication conservatively.

No clinical, diagnostic, causal, or validated shared-mechanism claims are allowed. Direction-only PD evidence remains preliminary, and the shared AD/PD axis claim remains blocked unless statistically supported PD replication aligns with an AD-supported axis.
