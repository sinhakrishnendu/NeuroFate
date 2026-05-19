# PNAS Validation Strategy

## Why Software Evidence Is Not Enough

NeuroFate is increasingly usable as a software platform, but a PNAS-oriented biological discovery paper needs a clear biological question, interpretable biological axes, negative controls, replicated cohort support, and conservative claim language. Installation, CLI quality, and release packaging support reproducibility; they do not by themselves establish a neurodegeneration discovery.

## NeuroFate-Axis Biological Reframing

The Phase 21 reframing treats NeuroFate outputs as donor/sample-level biological axes rather than only prediction features. Each axis summarizes a curated neurodegeneration layer such as microglial inflammation, astrocyte stress, myelin/oligodendrocyte biology, neuronal vulnerability, synuclein/mitochondrial stress, amyloid/tau biology, antigen presentation, vascular/barrier features, proteostasis/autophagy, and global neurodegeneration burden.

## Current Cohort Hierarchy

- SEA-AD: internal AD anchor cohort.
- Mathys 2019: AD external feasibility only because harmonized sample count is small.
- GSE243639: independent PD cohort with Phase 20 corrected cell/cluster-aware sample-level features and preliminary PD internal signal.
- Future AD cohorts: GSE174367, GSE147528, and ROSMAP/AD Knowledge Portal when access and formats are resolved.
- Future PD cohorts: at least one additional donor/sample-level PD brain cohort is needed for replication.

## Validation Ladder

1. Axis feature coverage in each cohort.
2. Donor/sample-level axis association with disease or pathology labels.
3. FDR-corrected axis statistics.
4. Random-axis negative controls.
5. AD/PD same-direction or disease-enriched comparison.
6. Independent replication in additional AD and PD cohorts.
7. Manuscript claims linked to the claim-strength table.

## Claim-Strength Rules

- `shared_ad_pd_candidate`: same-direction effect in AD and PD with adequate sample-level coverage, still requiring replication.
- `disease_specific_candidate`: stronger association in one disease cohort and weak or inconclusive evidence in the other.
- `axis_level_preliminary_evidence`: association-level support only.
- `axis_level_insufficient_validation`: no biological claim should be made.

## Forbidden Claim Shapes

Do not describe NeuroFate axes as causal axes, proven disease mechanisms, clinical biomarkers, definitive shared mechanisms, or validated across diseases. Acceptable language includes candidate shared axis, preliminary disease-specific axis, donor-level association, and exploratory cross-disease convergence.

## Future Cohorts Needed

The highest priority is one additional PD cohort to replicate the GSE243639 Phase 20 signal and one larger AD external cohort beyond Mathys feasibility. Controlled-access ROSMAP or AD Knowledge Portal data can strengthen donor-level external validation if provenance and data-use requirements are satisfied.

## Phase 22 Endpoint-Locked Validation

Phase 22 repairs the main statistical weakness of the initial Phase 21 axis comparison. Phase 21 selected the strongest axis effect across available labels, which is useful for exploration but not safe for PNAS-facing AD/PD biological claims.

Phase 22 locks endpoints before testing:

- `sea_ad_cognitive_dementia`: primary AD endpoint, Dementia versus No dementia.
- `sea_ad_ad_pathology_ordinal`: secondary AD pathology endpoint.
- `gse243639_pd_diagnosis`: primary PD endpoint, Parkinson's versus Control.
- Secondary GSE243639 Lewy/CERAD endpoints remain within-cohort pathology analyses.

Only `phase22_endpoint_locked_axis_evidence_table.tsv` and `phase22_endpoint_locked_axis_claims.md` should drive PNAS-facing axis claims. Phase 21 outputs remain exploratory and should not be used to claim shared AD/PD biology.

## Phase 23 Replication Cohort Onboarding

Phase 23 adds the missing replication scaffold. The registry `metadata/phase23_replication_cohort_registry.tsv` tracks:

- `gse184950_pd_sn`: first PD replication priority.
- `gse174367_ad_multiomics`: first AD replication priority, with bulk/sample-level data preferred initially.
- `gse147528_ad_progression`: secondary AD progression cohort.

Replication evidence is integrated only after donor/sample-level axis association statistics are available. Directionally consistent replication can produce `preliminary_replicated_candidate` status, but it still does not justify clinical, diagnostic, causal, or validated-across-diseases language. Stronger claims require at least one independent replication cohort and ideally two independent cohorts per disease.

## Phase 24 GSE184950 Replication Workflow

GSE184950 is the first PD replication cohort to onboard. The workbook `GSE184950_add2.xlsx` provides metadata fields such as disease state, donor ID, age, gender, postmortem interval, Braak stage, processed data file, and raw file. It does not provide expression data.

The RAW archive workflow is:

1. Parse the workbook metadata.
2. Manually download `GSE184950_RAW.tar`.
3. List archive contents without extracting members.
4. Prefer processed 10x matrices if the archive contains per-sample tarballs or matrix/barcode/feature files.
5. Generate selective extraction templates only.
6. Build sample-level axis scores and endpoint-locked PD/PDD versus Control statistics after processed matrices are manually prepared.

## Phase 25 Series-Matrix Repair For GSE184950

Phase 25 supersedes workbook-derived GSE184950 metadata for replication. The add2 workbook is useful as a schema hint but exposes only representative rows. The GEO series matrix is the primary metadata source because it provides all 34 samples, disease labels, donor-related characteristics, and per-sample supplementary tar links.

The Phase 25 endpoint is locked before analysis:

- `label__pd_pdd_vs_control = 1`: Parkinson's Disease or Parkinson's Disease Dementia.
- `label__pd_pdd_vs_control = 0`: Unaffected Control.

Archive reconciliation compares the listed RAW archive members to the series-matrix supplementary tar manifest without extracting files. Selective extraction templates target processed matrix files only. FASTQ/SRA processing remains out of scope unless a separate manual route is approved.

## Phase 26 GSE184950 Replication Execution

Phase 26 inspects nested per-sample archives and permits guarded extraction of processed 10x matrix files only. It does not process raw sequence files. The intended route is:

1. Inspect nested archives without extraction.
2. Selectively extract processed matrix, feature/gene, and barcode files only after `RUN_MANUAL_GSE184950_EXTRACTION=YES`.
3. Extract NeuroFate axis genes and aggregate directly to sample level.
4. Build axis scores.
5. Run endpoint-locked PD/PDD-vs-Control replication statistics.

GSE184950 strengthens the independent PD replication layer only if processed matrices are available and endpoint-locked axis effects are directionally consistent. Claims remain preliminary unless replication support is strong.

## Phase 27 Clean Replication Criteria

Phase 27 removes non-sample rows from GSE184950 axis-score tables and regenerates endpoint-locked statistics from the clean 34-sample table. The clean endpoint contains 24 PD/PDD positive samples and 10 Unaffected Control samples.

Evidence integration uses the following conservative rule:

- statistically supported replication: directionally consistent and `p < 0.05` or `FDR < 0.1`;
- directionally consistent preliminary signal: directionally consistent but statistically weak;
- weak or no replication: weak, missing, or unsupported effect;
- opposite direction: directionally inconsistent with the Phase 22 PD signal.

Direction-only evidence is not enough for PNAS-level replication claims.

## Phase 28 Independent AD Replication

Phase 28 adds the independent AD replication onboarding layer. GSE174367 is the first-priority cohort, with bulk/sample-level RNA-seq preferred before single-nucleus routes because the endpoint-locked NeuroFate-Axis framework operates at donor/sample level.

The AD replication workflow is:

1. Register candidate AD cohorts.
2. Use guarded manual acquisition templates.
3. Parse GEO series-matrix metadata when present.
4. Triage local files without opening large matrices.
5. Build axis scores from sample-level matrices when available.
6. Run endpoint-locked AD-vs-Control or pathology-defined replication tests.

Independent AD replication is marked `statistically_supported` only when at least one axis is directionally consistent with SEA-AD Phase 22 and passes p/FDR support. Otherwise it remains `available_but_preliminary` or `missing`.

## Phase 29: GSE174367 Bulk RNA AD Replication

Phase 29 turns the first-priority AD replication target into a guarded execution lane. The bulk processed RDA is inspected before conversion, and conversion is limited to NeuroFate axis genes. This avoids creating a genome-wide converted expression file and keeps the analysis sample-level.

The endpoint is AD versus Control (`label__ad_vs_control`: AD = 1, Control = 0). Replication statistics are endpoint-locked and compared with the SEA-AD Phase 22 dementia direction. The single-nucleus H5 matrix is not read in Phase 29; it remains a backup route only if bulk RNA is insufficient.

PNAS-facing interpretation remains conservative: statistically supported independent AD replication can strengthen candidate axis evidence, while direction-only support remains preliminary.

## Phase 32: Cross-Cohort Evidence Consolidation

Phase 32 combines the Phase 22 SEA-AD endpoint-locked axis table, Phase 31 GSE174367 AD replication, Phase 27 GSE184950 PD replication, and Phase 20 GSE243639 PD extension. The consolidated evidence table separates nominal AD replication from stronger FDR-robust replication.

Allowed Phase 32 language:

- nominal independent AD replication,
- directionally consistent AD replication candidate,
- preliminary PD extension,
- not yet ready for shared AD/PD mechanism claims.

The PNAS readiness matrix now marks independent AD replication as nominally supported when Phase 32 detects GSE174367 p-value support without FDR robustness. Shared AD/PD axis claims remain not ready until stronger independent PD axis replication is available.

## Phase 33: PD Replication Expansion

Phase 33 adds a public bulk/LCM/microarray replication route for PD. This route is designed to remove the main remaining PNAS bottleneck without relying on FASTQ/SRA processing or single-cell reanalysis.

Priority order:

1. `GSE20141` laser-dissected substantia nigra pars compacta neurons.
2. `GSE20186` PD expression superseries and usable subseries.
3. `GSE7621`, `GSE8397`, and `GSE20292` substantia nigra bulk/microarray cohorts.
4. `GSE157783` only if processed sample-level or matrix files are available.

The Phase 33 extractor reads GEO series-matrix expression sections only when present, retains only NeuroFate axis genes, and stops if probe identifiers require a platform annotation that has not been supplied. The platform mapper keeps only probes mapping to NeuroFate axis genes. Replication claims require endpoint-locked PD-versus-Control statistics; direction-only evidence cannot upgrade claim strength.

## Phase 34: PD Microarray/Bulk Replication Expansion

Phase 34 operationalizes the PD replication push with small public cohorts that are compatible with donor/sample-level NeuroFate-Axis testing. It adds a registry, guarded acquisition templates, metadata parsing, platform probe mapping, sample-level axis scoring, and conservative PD-versus-Control statistics.

The first recommended dataset is `GSE20141`, because it is laser-dissected SNpc neuron-focused and directly tests the neuronal vulnerability axis. `GSE7621` is the next bulk substantia nigra target. `GSE8397` requires region-aware filtering, and `GSE20186` requires clean subseries triage.

Phase 34 does not process FASTQ/SRA files and does not use single-cell tooling. It extracts or maps only NeuroFate axis genes/probes and does not write genome-wide converted matrices.

## Phase 37: GSE7621 PD Replication

Phase 37 adds GSE7621 as the next independent PD replication target after GSE20141 showed directionally consistent but non-significant support for the neuronal vulnerability axis. The workflow remains endpoint locked and sample level:

1. Manually acquire the GSE7621 series matrix.
2. Parse sample metadata and platform identifiers.
3. Audit expression sample IDs against parsed metadata without reading full expression rows.
4. Reuse an existing platform probe map when possible or build a compact NeuroFate-only probe map.
5. Build sample-level axis scores.
6. Test PD versus Control replication with conservative p/FDR rules.

The PNAS bottleneck improves only if GSE7621 or another independent PD cohort shows directionally consistent statistical support. Direction-only support remains preliminary.

## Phase 38: GSE7621 Direction-Aware Interpretation

Phase 38 separates technical success from shared-axis evidence. GSE7621 passes sample-matching and axis-coverage checks, but its statistically strongest signal is an opposite-direction synuclein-mitochondrial axis. That result should be treated as a candidate PD-divergent axis, not a shared AD/PD replication result.

The neuronal vulnerability axis remains directionally consistent but statistically weak in GSE7621. Therefore the PNAS-facing interpretation remains conservative: stronger PD replication is still needed, and the synuclein-mitochondrial result requires independent confirmation and platform/probe review.

## Phase 39: Manuscript-Readiness Consolidation

Phase 39 shifts the near-term writing target from PNAS to eLife. The PNAS ladder remains useful for discipline, but the current evidence supports a realistic systems-biology manuscript rather than a definitive shared-mechanism claim.

The eLife-ready claim is: NeuroFate-Axis identifies endpoint-locked neurodegeneration axes, with strongest support for an AD-replicated neuronal vulnerability axis and additional PD evidence that is preliminary or divergent. PNAS remains a future goal after stronger PD replication.
