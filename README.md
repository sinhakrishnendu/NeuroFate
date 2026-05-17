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

SEA-AD supports the internal NeuroFate analysis path. Mathys 2019 currently serves as preliminary external feasibility evidence unless larger sample-level harmonization is added. Reports and manuscript text should keep this distinction visible.

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
