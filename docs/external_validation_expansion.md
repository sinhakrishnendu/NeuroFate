# External Validation Expansion

NeuroFate Phase 15 expands the platform from Mathys 2019 feasibility testing toward a reusable multi-cohort validation layer. The goal is independent donor/sample-level validation without weakening the project safety rules: NeuroFate does not download external data automatically, does not run Scanpy-based analysis, and does not read full expression matrices during planning.

## Why External Validation Matters

Internal SEA-AD performance is useful for model development, but reviewer-grade claims require independent cohorts, transparent provenance, donor/sample-level harmonization, and conservative reliability labels. Phase 15 separates acquisition, file inspection, metadata mapping, sparse extraction planning, feature-table construction, validation, and claim-strength updating.

## Dataset Priority Table

| Priority | Dataset | Disease | Role |
| --- | --- | --- | --- |
| A | GSE243639 | Parkinson disease | First serious PD external validation target using human substantia nigra pars compacta snRNA-seq. |
| B | GSE174367 | Alzheimer disease | AD multi-omics external validation candidate if processed count and metadata files are usable. |
| C | GSE147528 | Alzheimer disease | AD progression snRNA-seq candidate pending file-format and donor metadata triage. |
| D | ROSMAP / AD Knowledge Portal | Alzheimer disease | Controlled or registered bulk donor-level validation option. |

## AD Cohorts

GSE174367 and GSE147528 are treated as manually acquired external AD candidates. They must pass local file inventory, metadata-field mapping, target-gene overlap, and donor/sample-level feature-table checks before any validation result is interpreted.

## PD Cohorts

GSE243639 is the highest-priority PD candidate. It is intended to test whether NeuroFate target-gene and donor-level signatures can transfer to Parkinson disease without pretending that SEA-AD Alzheimer labels are equivalent to PD labels.

## Phase 16: GSE243639 PD External Validation

Phase 16 onboards GSE243639 as an independent Parkinson disease cohort extension. The clinical file is parsed with a semicolon delimiter and a 1-based header line of 6 because the first five lines contain prose and a blank line. The count table is treated as a genes-as-rows CSV gzip file, with cell IDs such as `s.0096_AAACCCAAGTACGAGC.1`; sample IDs are derived from the prefix before the underscore.

The Phase 16 workflow is intentionally conservative:

- target-gene extraction is streaming and NeuroFate-panel-only,
- expression is aggregated to sample level before validation,
- Parkinson's/control labels are evaluated within GSE243639 as an independent PD cohort internal validation,
- SEA-AD to GSE243639 output is described as feature-space transfer feasibility, not direct AD-to-PD disease-label transfer,
- the workflow is not medical validation and does not establish cause-and-effect biology.

Manual commands:

```bash
python scripts/78_extract_gse243639_target_gene_panel.py \
  --counts data/raw/external/gse243639_pd_snpc/GSE243639_Filtered_count_table.csv.gz \
  --panel metadata/target_gene_panel_v1.tsv \
  --output data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz \
  --audit-output results/tables/phase16_gse243639_gene_extraction_audit.tsv \
  --log-file results/logs/78_extract_gse243639_target_gene_panel.log

python scripts/79_build_gse243639_feature_table.py \
  --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz \
  --clinical data/raw/external/gse243639_pd_snpc/GSE243639_Clinical_data.csv.gz \
  --phase5-schema results/tables/phase5_donor_feature_table.tsv \
  --output results/tables/phase16_gse243639_feature_table.tsv \
  --schema-output results/tables/phase16_gse243639_feature_schema_alignment.tsv \
  --label-summary-output results/tables/phase16_gse243639_label_summary.tsv \
  --log-file results/logs/79_build_gse243639_feature_table.log

python scripts/80_run_gse243639_pd_external_validation.py \
  --sea-ad-features results/tables/phase5_donor_feature_table.tsv \
  --pd-features results/tables/phase16_gse243639_feature_table.tsv \
  --metrics-output results/tables/phase16_gse243639_external_validation_metrics.tsv \
  --predictions-output results/tables/phase16_gse243639_external_predictions.tsv \
  --log-file results/logs/80_run_gse243639_pd_external_validation.log
```

## Phase 17: GSE243639 Cell-Type-Aware PD Refinement

Phase 17 uses `GSE243639_UMAP_coordinates.xlsx` only as an existing annotation and coordinate table. NeuroFate does not recompute UMAP, does not cluster cells, and does not create H5AD/AnnData objects. The aim is to test whether sample-level cell-type-aware aggregation improves the exploratory PD signal relative to Phase 16 global sample features.

Manual commands:

```bash
python scripts/83_inspect_gse243639_umap_annotations.py \
  --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx \
  --output results/reports/phase17_gse243639_umap_annotation_audit.tsv \
  --preview-output results/reports/phase17_gse243639_umap_annotation_preview.tsv \
  --log-file results/logs/83_inspect_gse243639_umap_annotations.log

python scripts/84_build_gse243639_cell_annotation_map.py \
  --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx \
  --cell-sample-map data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv \
  --output data/interim/external/gse243639_pd_snpc/gse243639_cell_annotation_map.tsv \
  --summary-output results/tables/phase17_gse243639_cell_annotation_summary.tsv \
  --log-file results/logs/84_build_gse243639_cell_annotation_map.log

python scripts/85_build_gse243639_celltype_feature_table.py \
  --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz \
  --annotations data/interim/external/gse243639_pd_snpc/gse243639_cell_annotation_map.tsv \
  --clinical data/raw/external/gse243639_pd_snpc/GSE243639_Clinical_data.csv.gz \
  --phase5-schema results/tables/phase5_donor_feature_table.tsv \
  --output results/tables/phase17_gse243639_celltype_feature_table.tsv \
  --schema-output results/tables/phase17_gse243639_celltype_schema_alignment.tsv \
  --label-summary-output results/tables/phase17_gse243639_celltype_label_summary.tsv \
  --log-file results/logs/85_build_gse243639_celltype_feature_table.log

python scripts/86_run_gse243639_celltype_pd_validation.py \
  --features results/tables/phase17_gse243639_celltype_feature_table.tsv \
  --metrics-output results/tables/phase17_gse243639_celltype_validation_metrics.tsv \
  --predictions-output results/tables/phase17_gse243639_celltype_predictions.tsv \
  --importance-output results/tables/phase17_gse243639_celltype_feature_importance.tsv \
  --log-file results/logs/86_run_gse243639_celltype_pd_validation.log
```

Reliability remains conservative: Phase 17 can be `moderate_pd_internal_validation`, `preliminary_pd_internal_signal`, or `weak_pd_signal`. The result is not medical validation, not cause-and-effect inference, and not direct AD-to-PD disease-label transfer.

## Phase 18: GSE243639 Cell-ID Repair and Repaired Cell-Type-Aware PD Validation

Phase 18 treats the weak Phase 17 result as a technical audit outcome because annotation matching collapsed the cell-type-aware feature space. The repair starts by comparing expression cell IDs, the cell-sample map, and workbook cell IDs under multiple normalization rules. It then rebuilds the annotation map from expression/cell-map IDs so that every matched annotation remains traceable to an expression cell ID.

Final PD interpretation should use Phase 18 rather than Phase 17 if the repaired annotation join succeeds.

Manual commands:

```bash
python scripts/90_audit_gse243639_cell_id_matching.py \
  --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz \
  --cell-sample-map data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv \
  --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx \
  --output results/tables/phase18_gse243639_cell_id_matching_audit.tsv \
  --preview-output results/reports/phase18_gse243639_cell_id_matching_preview.tsv \
  --log-file results/logs/90_audit_gse243639_cell_id_matching.log

python scripts/84_build_gse243639_cell_annotation_map.py \
  --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx \
  --cell-sample-map data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv \
  --output data/interim/external/gse243639_pd_snpc/gse243639_cell_annotation_map.tsv \
  --summary-output results/tables/phase18_gse243639_annotation_match_summary.tsv \
  --candidate-output results/reports/phase18_gse243639_annotation_column_candidates.tsv \
  --log-file results/logs/84_build_gse243639_cell_annotation_map.log

python scripts/85_build_gse243639_celltype_feature_table.py \
  --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz \
  --annotations data/interim/external/gse243639_pd_snpc/gse243639_cell_annotation_map.tsv \
  --clinical data/raw/external/gse243639_pd_snpc/GSE243639_Clinical_data.csv.gz \
  --phase5-schema results/tables/phase5_donor_feature_table.tsv \
  --output results/tables/phase18_gse243639_celltype_feature_table.tsv \
  --schema-output results/tables/phase18_gse243639_celltype_schema_alignment.tsv \
  --label-summary-output results/tables/phase18_gse243639_celltype_label_summary.tsv \
  --feature-group-output results/tables/phase18_gse243639_feature_group_counts.tsv \
  --log-file results/logs/85_build_gse243639_celltype_feature_table.log

python scripts/91_run_gse243639_repaired_celltype_pd_validation.py \
  --features results/tables/phase18_gse243639_celltype_feature_table.tsv \
  --metrics-output results/tables/phase18_gse243639_celltype_validation_metrics.tsv \
  --predictions-output results/tables/phase18_gse243639_celltype_predictions.tsv \
  --importance-output results/tables/phase18_gse243639_celltype_feature_importance.tsv \
  --log-file results/logs/91_run_gse243639_repaired_celltype_pd_validation.log
```

If annotation match rate remains poor or repaired feature count remains too small, Phase 18 reports `technical_failure_annotation_join` and no biological PD conclusion should be drawn from the cell-type-aware layer.

## Phase 19: GSE243639 Annotation-Linkage Forensics

Phase 19 addresses the stricter question raised by Phase 18: whether the workbook annotation IDs can be safely linked to expression cells. Current Phase 18 evidence indicates that expression extraction and the cell-sample map are internally consistent, while workbook `CELL_ID` values have zero overlap under direct and normalized comparisons. Phase 16 therefore remains the valid global sample-level PD extension unless Phase 19 proves a safe annotation linkage rule.

Phase 17 and Phase 18 should not be interpreted biologically when annotation linkage is unsafe. Cell-type-aware PD validation requires safe annotation linkage; if linkage is unsafe, the correct scientific decision is not to force it.

Phase 19 evaluates:

- forensic previews of expression, cell-sample-map, and workbook cell IDs,
- deep normalization rules including Seurat-style conversions,
- row-order linkage only as an audited hypothesis,
- a final annotation-linkage decision that can block annotation-map creation.

Manual commands:

```bash
python scripts/95_forensic_gse243639_workbook_cell_ids.py \
  --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz \
  --cell-sample-map data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv \
  --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx \
  --output results/reports/phase19_gse243639_cell_id_forensic_preview.tsv \
  --log-file results/logs/95_forensic_gse243639_workbook_cell_ids.log

python scripts/96_deep_gse243639_cell_id_normalization_audit.py \
  --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz \
  --cell-sample-map data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv \
  --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx \
  --output results/tables/phase19_gse243639_deep_cell_id_overlap.tsv \
  --best-rule-output results/reports/phase19_gse243639_best_normalization_rule.md \
  --log-file results/logs/96_deep_gse243639_cell_id_normalization_audit.log

python scripts/97_audit_gse243639_row_order_annotation_link.py \
  --cell-sample-map data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv \
  --xlsx data/raw/external/gse243639_pd_snpc/GSE243639_UMAP_coordinates.xlsx \
  --clinical data/raw/external/gse243639_pd_snpc/GSE243639_Clinical_data.csv.gz \
  --output results/tables/phase19_gse243639_row_order_link_audit.tsv \
  --preview-output results/reports/phase19_gse243639_row_order_link_preview.tsv \
  --log-file results/logs/97_audit_gse243639_row_order_annotation_link.log

python scripts/98_decide_gse243639_annotation_linkage.py \
  --normalization-audit results/tables/phase19_gse243639_deep_cell_id_overlap.tsv \
  --row-order-audit results/tables/phase19_gse243639_row_order_link_audit.tsv \
  --phase18-match-summary results/tables/phase18_gse243639_annotation_match_summary.tsv \
  --output-md results/reports/phase19_gse243639_annotation_linkage_decision.md \
  --output-tsv results/tables/phase19_gse243639_annotation_linkage_decision.tsv

python scripts/99_build_gse243639_safe_annotation_map_if_valid.py \
  --decision results/tables/phase19_gse243639_annotation_linkage_decision.tsv \
  --existing-annotation-map data/interim/external/gse243639_pd_snpc/gse243639_cell_annotation_map.tsv \
  --output data/interim/external/gse243639_pd_snpc/gse243639_safe_cell_annotation_map.tsv \
  --blocked-output results/reports/phase19_annotation_linkage_blocked.md
```

## Phase 20: Safe Annotation Map Consumption and Corrected Cell-Type-Aware PD Validation

Phase 19 can establish a safe normalized ID linkage, but the feature builder still has to consume the safe map correctly. Phase 20 repairs that interface by auditing the safe-map schema, joining expression `cell_id` values directly to `safe_map.cell_id_expression`, retaining only matched safe-map rows, and computing annotation match rate from unique expression cells rather than sparse expression rows.

Only Phase 20 should be used for final cell-type-aware PD interpretation if the annotation match rate is high. Phase 16 remains the valid global sample-level PD extension, and Phase 17/18 are retained as technical diagnostics.

The current successful Phase 20 run has annotation match rate `1.0`, feature count `1590`, logistic repeated AUROC `0.72777778`, AUPRC `0.79044444`, balanced accuracy `0.64666667`, empirical permutation p-value `0.10891089`, and reliability `preliminary_pd_internal_signal`. The correct interpretation is that safe-map cell-type-aware features improve the GSE243639 PD signal relative to Phase 16 global features, but the evidence remains preliminary because permutation support is not significant at 0.05 and there is not yet independent PD-cohort replication.

Manual commands:

```bash
python scripts/100_audit_phase19_safe_annotation_map_schema.py \
  --safe-map data/interim/external/gse243639_pd_snpc/gse243639_safe_cell_annotation_map.tsv \
  --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz \
  --output results/tables/phase20_safe_annotation_map_schema_audit.tsv \
  --preview-output results/reports/phase20_safe_annotation_map_preview.tsv \
  --log-file results/logs/100_audit_phase19_safe_annotation_map_schema.log

python scripts/101_build_gse243639_phase20_celltype_features.py \
  --expression data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz \
  --annotations data/interim/external/gse243639_pd_snpc/gse243639_safe_cell_annotation_map.tsv \
  --clinical data/raw/external/gse243639_pd_snpc/GSE243639_Clinical_data.csv.gz \
  --phase5-schema results/tables/phase5_donor_feature_table.tsv \
  --output results/tables/phase20_gse243639_celltype_feature_table.tsv \
  --schema-output results/tables/phase20_gse243639_celltype_schema_alignment.tsv \
  --label-summary-output results/tables/phase20_gse243639_celltype_label_summary.tsv \
  --feature-group-output results/tables/phase20_gse243639_feature_group_counts.tsv \
  --log-file results/logs/101_build_gse243639_phase20_celltype_features.log

python scripts/102_run_gse243639_phase20_celltype_pd_validation.py \
  --features results/tables/phase20_gse243639_celltype_feature_table.tsv \
  --metrics-output results/tables/phase20_gse243639_celltype_validation_metrics.tsv \
  --predictions-output results/tables/phase20_gse243639_celltype_predictions.tsv \
  --importance-output results/tables/phase20_gse243639_celltype_feature_importance.tsv \
  --log-file results/logs/102_run_gse243639_phase20_celltype_pd_validation.log
```

## Controlled-Access Cohorts

ROSMAP and AD Knowledge Portal resources may require registration, Synapse approvals, or data-use certifications. NeuroFate records these as controlled-access candidates and does not bundle or retrieve them.

## Accepted Formats

The Phase 15 planning scripts support local files in H5AD/H5, MTX, CSV/TSV, text, RDS, and LOOM containers at the inventory level. Extraction planning currently creates manual templates for H5AD sparse CSR, CSV genes-as-rows, CSV cells-as-rows, MTX plus features/barcodes, and bulk RNA-seq matrices.

## Metadata Harmonization

External metadata are mapped into canonical fields where possible: donor ID, sample ID, cell ID, diagnosis, disease status, pathology, brain region, age, sex, cell type, batch, APOE genotype, and sequencing platform. Ambiguous fields remain low-confidence suggestions until manually reviewed.

## Donor/Sample-Level Aggregation

External validation must use donor-level, sample-level, pseudo-donor, or bulk-sample units. Cell-level rows are never used directly as independent validation observations. Scripts write warnings when sample counts are too small for strong claims.

## Reliability Categories

- `reliable_external_validation`: sufficient sample size, both classes present, and adequate feature overlap.
- `preliminary_external_feasibility`: an external test exists but is too small for a validation claim.
- `insufficient_sample_size`: too few donor/sample units.
- `insufficient_feature_overlap`: target features do not align well enough.
- `failed_label_mapping`: labels cannot support a case-control or disease-status evaluation.

## How Claims Are Upgraded

Phase 15 can upgrade claim strength only when a reliable external validation result exists. Mathys n=6 or any similarly small external check remains preliminary feasibility and cannot justify wording such as "validated across cohorts."

## Why Mathys Remains Preliminary

The current Mathys harmonization has six sample-level units. This is useful for testing file handling, feature alignment, and external transfer mechanics, but it is not enough for robust cross-cohort validation.

## How To Add A New Cohort

1. Add a row to `metadata/phase15_external_validation_candidates.tsv`.
2. Create or adapt a guarded manual acquisition template.
3. Run `scripts/70_inspect_external_dataset_files.py` on the local directory.
4. Run `scripts/71_inspect_external_metadata_safe.py` on metadata or container files.
5. Run `scripts/72_plan_external_target_gene_overlap.py` on feature/gene files.
6. Generate a manual extraction plan with `scripts/73_prepare_external_sparse_extraction_plan.py`.
7. Build a donor/sample-level table with `scripts/74_build_external_feature_table_generic.py`.
8. Run multi-external validation only after feature tables and labels have been reviewed.

## Phase 21: PNAS-Oriented Biological Discovery Framework

Phase 21 shifts NeuroFate from validation engineering toward a biological discovery question: whether AD and PD share conserved donor/sample-level neurodegeneration axes while retaining disease-specific amyloid/tau and synuclein/mitochondrial structure.

The Phase 21 axis framework uses only existing donor/sample-level tables:

- SEA-AD donor-level features from the AD anchor workflow.
- GSE243639 Phase 20 sample-level cell/cluster-aware PD features.
- Mathys feasibility outputs only as context, not as strong validation.

The new `metadata/neurofate_axis_registry.tsv` defines curated axes for inflammatory microglial activation, astrocyte stress, myelin/oligodendrocyte biology, neuronal vulnerability, synuclein/mitochondrial stress, amyloid/tau biology, antigen presentation, vascular/barrier biology, proteostasis/autophagy, and global neurodegeneration burden.

Phase 21 adds axis-level claim categories:

- `shared_ad_pd_candidate`
- `ad_enriched_axis`
- `pd_enriched_axis`
- `inconclusive_axis`
- `insufficient_coverage`

These are biological association categories, not clinical validation labels. Acceptable language includes candidate shared axis, preliminary disease-specific axis, donor-level association, and exploratory cross-disease convergence. NeuroFate should not claim causal axes, proven disease mechanisms, clinical biomarkers, definitive shared mechanisms, or validation across diseases without independent replication.

## Phase 23: Independent Replication Cohort Onboarding

Phase 23 extends the endpoint-locked framework toward PNAS readiness by preparing independent replication cohorts. GSE184950 is the priority PD replication target. GSE174367 is the priority AD replication target, with bulk/sample-level expression planned first if available. GSE147528 is a secondary AD progression-validation target.

The Phase 23 scripts are planning and donor/sample-level scaffolds only:

- local file inventory and next-action triage,
- guarded manual acquisition templates,
- generic sample-matrix axis scoring,
- snRNA-seq extraction planning templates,
- endpoint-locked binary replication association tests,
- replication evidence integration,
- PNAS readiness matrix.

No Phase 23 script downloads data, runs SRA tools, processes raw SRA, loads H5AD expression matrices, recomputes UMAP, clusters cells, or trains models. Replication status cannot upgrade a NeuroFate-Axis claim unless directionally consistent independent cohort statistics exist.

## Phase 24: GSE184950 RAW Archive Inspection and Extraction Planning

Phase 24 prepares the GSE184950 PD replication cohort after the metadata workbook has been downloaded. The workbook alone is insufficient for replication because it does not include expression data. It provides sample metadata and names of processed per-sample files such as `A22.tar.gz` and `D14.tar.gz`.

The safe Phase 24 plan is:

- parse workbook metadata from the `METADATA TEMPLATE` sheet,
- manually download `GSE184950_RAW.tar`,
- list archive members with Python `tarfile` without extraction,
- identify processed 10x matrices, barcodes, features/genes, and FASTQ files,
- generate a manual selective extraction template for processed matrices,
- avoid FASTQ/SRA processing by default,
- compute sample-level axis scores only after processed matrices are manually prepared.

## Phase 25: GSE184950 Series-Matrix Metadata and Replication Route

Phase 25 makes GSE184950 usable for real endpoint-locked replication planning by replacing incomplete workbook metadata with the GEO series matrix. The series matrix contains all 34 sample records and per-sample supplementary tar links, while the workbook alone is insufficient for replication.

The primary endpoint is:

- PD/PDD positive: Parkinson's Disease and Parkinson's Disease Dementia.
- Control negative: Unaffected Control.

Phase 25 scripts parse the series matrix, reconcile expected tar files with the listed RAW archive inventory, and create a guarded selective extraction plan for processed matrices. They do not download data, extract archives, process FASTQ/SRA files, use Scanpy, create H5AD/AnnData, recompute UMAP, cluster cells, load dense matrices, or train models.

## Phase 26: GSE184950 Nested Archive Inspection and Axis Replication

Phase 26 adds a safe execution path for the GSE184950 replication cohort. The outer RAW archive is opened for listing only, and each nested per-sample processed archive is inspected through a tar stream without writing files. Selective extraction is available only for processed matrix files and only with the manual guard `RUN_MANUAL_GSE184950_EXTRACTION=YES`.

Allowed extracted file names are limited to processed 10x matrix components:

- `matrix.mtx.gz`
- `features.tsv.gz`
- `genes.tsv.gz`
- `barcodes.tsv.gz`

The downstream NeuroFate-Axis extraction reads only axis genes and writes sample-level summaries before endpoint-locked PD/PDD-vs-Control replication statistics are calculated. The result is an independent PD replication test, not a clinical or diagnostic validation.

## Phase 27: GSE184950 Clean Sample QC and Conservative Evidence Integration

Phase 27 audits GSE184950 sample IDs against the series-matrix metadata and removes technical pseudo-samples from clean axis-score tables. Clean replication statistics must use exactly 34 biological samples.

Replication evidence is integrated conservatively. Directional consistency without p/FDR support is labelled as a preliminary signal, not replication. Claim upgrades remain blocked unless clean GSE184950 statistics show statistically supported independent replication.

## Phase 28: Independent AD Replication Onboarding

Phase 28 adds AD replication candidates to close the largest remaining PNAS gap. GSE174367 is prioritized because it may provide bulk/sample-level AD expression suitable for donor/sample-level NeuroFate-Axis replication. GSE147528 and GSE157827 are secondary or optional cohorts.

Phase 28 scripts parse GEO series metadata, triage local files, build axis scores from sample-level matrices, and run endpoint-locked AD replication tests. They do not download data, extract archives, process SRA/FASTQ, use Scanpy, create H5AD/AnnData, run UMAP or clustering, load dense matrices, or train models.

## Phase 29: GSE174367 Bulk RNA Independent AD Replication

GSE174367 now has a bulk RNA replication route using the processed `GSE174367_bulkRNA_processed.rda.gz` file and the parsed 230-sample series matrix metadata. The parsed endpoint contains 118 AD and 112 Control records.

The Phase 29 scripts:

- inspect RDA structure before conversion,
- extract only genes present in the NeuroFate axis registry,
- map expression samples back to series-matrix metadata,
- build sample-level axis scores,
- run endpoint-locked AD-versus-Control replication tests,
- keep single-nucleus H5 resources optional for later work.

If object selection or sample mapping is ambiguous, Phase 29 stops and reports the ambiguity instead of guessing. It does not make clinical, diagnostic, causal, or definitive PNAS-ready claims.

## Phase 32: Cross-Cohort Evidence Consolidation

Phase 32 consolidates existing outputs rather than processing new data. It uses:

- Phase 22 SEA-AD endpoint-locked axis evidence,
- Phase 31 GSE174367 bulk-RNA AD replication statistics,
- Phase 27 GSE184950 clean PD/PDD replication statistics,
- Phase 20 GSE243639 PD extension metrics.

The current strongest axis is `neuronal_vulnerability_axis`, supported as a nominal independent AD replication candidate. PD evidence remains preliminary, so stronger PD axis replication is the next external-validation priority.

## Phase 33: PD Replication Expansion

Phase 33 expands PD replication beyond GSE243639 and GSE184950 using public donor/sample-level expression cohorts. The purpose is to test whether NeuroFate axes, especially the neuronal vulnerability axis, reproduce in independent PD substantia nigra or laser-captured neuron datasets.

Registered targets:

- `GSE20141` PD substantia nigra pars compacta laser-captured neurons.
- `GSE20186` PD expression superseries.
- `GSE7621`, `GSE8397`, and `GSE20292` PD substantia nigra bulk/microarray backups.
- `GSE157783` optional small midbrain snRNA cohort if processed matrices are accessible.

Phase 33 acquisition remains manual and guarded. Scripts prefer GEO series matrix files and processed expression tables, use platform annotation only for conservative probe-to-gene mapping, extract only NeuroFate axis genes, and keep all analyses at donor/sample level. It does not download data automatically, process FASTQ/SRA, use Scanpy, create H5AD/AnnData, run UMAP or clustering, write dense genome-wide matrices, or train models.

## Phase 34: PD Microarray/Bulk Replication Expansion

Phase 34 adds an explicit PD microarray/bulk replication lane for the remaining PNAS bottleneck. The registry focuses on public sample-level cohorts:

- `GSE20141` SNpc laser-captured neuron-focused cohort.
- `GSE7621` substantia nigra bulk/microarray cohort.
- `GSE8397` region-aware PD/control brain expression cohort.
- `GSE20186` superseries/subseries route.

The workflow is:

1. Manually acquire series matrix, SOFT/MINiML, and platform annotation files.
2. Parse sample metadata and confirm unambiguous PD/control labels.
3. Map platform probes only to NeuroFate axis genes.
4. Build sample-level axis scores without writing genome-wide converted expression.
5. Run endpoint-locked PD-versus-Control axis replication statistics.

Claim integration remains conservative. Direction-only PD evidence remains preliminary, while statistically supported PD replication requires direction consistency and `p < 0.05` or `FDR < 0.1`.

## Phase 37: GSE7621 PD Replication

Phase 37 extends the PD replication expansion to GSE7621, a substantia nigra bulk/microarray cohort. It is intended to test whether the current direction-only PD evidence strengthens in another independent sample-level cohort.

The route uses only GEO series-matrix expression values and platform probe mappings. Raw CEL/CHP files, FASTQ/SRA files, Scanpy, H5AD/AnnData, UMAP, clustering, dense genome-wide converted outputs, and model training remain out of scope.

Claim language remains conservative: GSE7621 can support an independent PD replication candidate only if an axis is directionally consistent and has p/FDR support. It cannot by itself establish clinical utility, diagnostic performance, causality, or a validated shared AD/PD mechanism.

## Phase 38: GSE7621 Direction-Aware PD Replication Interpretation

Phase 38 audits the successful GSE7621 result. The dataset is technically usable: GEO accession sample joins are complete, PD/control labels are available, and NeuroFate axis coverage is strong. The biological result is mixed rather than confirmatory.

The neuronal vulnerability axis is directionally consistent but not significant, so it does not provide statistically supported PD replication of the AD-replicated axis. The synuclein-mitochondrial axis is statistically strong but opposite in direction, which supports only a candidate PD-divergent interpretation after probe and direction audit.

This improves the biological story by identifying a possible PD-specific divergence, but it does not validate a shared AD/PD mechanism.

## Phase 39: Manuscript-Readiness Consolidation

Phase 39 consolidates the external-validation evidence for manuscript writing. The feasible target is eLife, with NAR Genomics and Bioinformatics as backup. The manuscript should present NeuroFate-Axis as an endpoint-locked donor/sample-level framework, not as a clinical or diagnostic system.

The external-validation story is:

- GSE174367 provides nominal independent AD replication of neuronal vulnerability.
- GSE243639 provides preliminary PD signal.
- GSE184950 and GSE20141 show weak or direction-only PD replication.
- GSE7621 shows candidate PD-divergent synuclein--mitochondrial remodeling.

The correct conclusion is conservative: AD evidence is strongest, PD evidence is biologically informative but unresolved, and shared AD/PD mechanism claims remain unsupported.
