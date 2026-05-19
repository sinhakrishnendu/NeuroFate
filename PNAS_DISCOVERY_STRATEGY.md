# NeuroFate PNAS Discovery Strategy

## 1. Central Biological Question

Do Alzheimer disease and Parkinson disease share donor-level neurodegeneration fate axes that capture conserved glial-inflammatory, myelin, and neuronal vulnerability programs, while diverging in amyloid/tau- and synuclein-associated axis structure?

## 2. Why AD and PD Can Be Compared

AD and PD are distinct diseases, but both involve age-associated human brain degeneration, glial activation, neuronal vulnerability, proteostasis stress, mitochondrial dysfunction, and inflammatory remodeling. NeuroFate compares them at donor/sample level rather than cell level, reducing pseudo-replication and allowing disease-cohort differences to be interpreted as systems-level associations.

## 3. NeuroFate Axis Concept

A NeuroFate axis is a curated, interpretable donor-level score built from available gene, cell-type, and signature features. Axes are not clinical biomarkers. They are structured biological summaries designed to test whether disease-associated variance converges on shared or disease-specific neurodegeneration programs.

## 4. Shared Neurodegeneration Axes

Candidate shared axes include inflammatory microglial activation, antigen presentation, astrocyte stress, myelin/oligodendrocyte disruption, neuronal vulnerability, mitochondrial stress, and global neurodegeneration burden. A shared axis requires same-direction donor/sample-level association in both AD and PD cohorts with sufficient coverage and conservative statistical support.

## 5. Disease-Specific Axes

Candidate disease-specific axes include amyloid/tau-associated structure in AD and synuclein/mitochondrial structure in PD. A disease-specific axis requires stronger association in one disease cohort and weak, inconsistent, or insufficient evidence in the other. NeuroFate should describe such axes as candidates until replicated.

## 6. Required Validation Cohorts

Current cohort hierarchy:

- SEA-AD: internal AD anchor cohort.
- Mathys 2019: preliminary AD external feasibility due small harmonized sample count.
- GSE243639: independent PD cohort with Phase 20 corrected cell/cluster-aware donor-level features.
- Future AD cohorts: GSE174367, GSE147528, ROSMAP/AD Knowledge Portal if access allows.
- Future PD cohorts: at least one additional donor/sample-level PD brain cohort for replication.

## 7. What Can Be Claimed Now

NeuroFate can claim that it provides a donor/sample-level framework for defining and testing candidate neurodegeneration axes across AD and PD cohorts. GSE243639 Phase 20 provides preliminary PD internal signal after safe annotation linkage, with annotation match rate 1.0, 1590 features, AUROC approximately 0.728, AUPRC approximately 0.790, balanced accuracy approximately 0.647, empirical permutation p approximately 0.109, and reliability `preliminary_pd_internal_signal`.

## 8. What Cannot Be Claimed Yet

NeuroFate cannot yet claim clinical diagnosis, clinical-grade validation, causal disease mechanisms, deployed biomarkers, definitive shared mechanisms, or validation across diseases. Axis language must remain candidate, preliminary, donor-level, and association-based until independent replication is complete.

## 9. PNAS Readiness Criteria

PNAS readiness requires:

- axis-level AD and PD association tables,
- random-axis negative controls,
- conservative multiple-testing correction,
- clear separation of shared versus disease-specific candidate axes,
- at least one robust AD anchor cohort and one robust PD cohort,
- ideally independent replication in an additional AD and PD cohort,
- no leakage or overclaiming flags,
- manuscript claims aligned to evidence-strength tables.

## 10. Next Validation Priorities

1. Run Phase 21 axis scores on existing SEA-AD and GSE243639 donor/sample-level tables.
2. Add an additional PD cohort to replicate Phase 20 axis-level findings.
3. Add one larger AD external cohort beyond Mathys n=6 feasibility.
4. Run random-axis controls and FDR correction before elevating any shared-axis claim.
5. Update the manuscript only from `phase21_axis_claim_strength.tsv` and no-overclaiming outputs.

## Phase 22 Endpoint-Locked Repair

Phase 21 is now treated as an exploratory axis-discovery layer because it compared the strongest absolute axis effect across heterogeneous labels. PNAS-facing biological claims must instead use Phase 22 endpoint-locked tables.

Primary Phase 22 endpoints are fixed before testing:

- SEA-AD primary AD endpoint: `sea_ad_cognitive_dementia`, comparing Dementia versus No dementia from `label__Cognitive_Status`.
- GSE243639 primary PD endpoint: `gse243639_pd_diagnosis`, comparing Parkinson's versus Control.

Secondary pathology endpoints remain within-cohort support unless explicitly labelled. Matched random-axis controls must use the same endpoint and the same association statistic as the curated axis. Phase 22 replaces Phase 21 for biological claim language.

## Phase 23 Independent Replication

Phase 23 defines the next validation step required for PNAS-level claims. The priority PD replication cohort is GSE184950, a human substantia nigra Parkinson disease single-cell transcriptomics dataset. The priority AD replication cohort is GSE174367, with bulk/sample-level RNA-seq planned first if available because donor/sample-level axis replication is the immediate need. GSE147528 is retained as a secondary AD progression-validation cohort.

Phase 23 does not download or process data automatically. It creates guarded acquisition templates, local file triage, sample-matrix axis scoring support, snRNA extraction planning templates, endpoint-locked replication association tests, replication integration, and a PNAS readiness matrix. Axis claims should remain candidate or preliminary until independent replication is directionally consistent and adequately powered.

## Phase 24 GSE184950 RAW Archive Planning

Phase 24 specializes the PD replication plan for GSE184950. The `GSE184950_add2.xlsx` workbook is useful for metadata parsing and processed-file names, but it is insufficient for axis replication because it does not contain expression data. The next required step is manual download of `GSE184950_RAW.tar`, safe archive listing without extraction, and selective planning for processed 10x matrix members. FASTQ/SRA processing is deliberately avoided unless no processed matrices exist and a separate manual plan is approved.

## Phase 25 GSE184950 Series-Matrix Metadata

Phase 25 corrects the metadata source for GSE184950. The add2 workbook is incomplete for replication, while the GEO series matrix provides the full 34-sample cohort and per-sample supplementary tar links. NeuroFate now treats the series matrix as the primary GSE184950 metadata source.

The endpoint-locked replication contrast is PD/PDD versus Unaffected Control:

- positive: Parkinson's Disease and Parkinson's Disease Dementia,
- negative: Unaffected Control.

This remains a biological replication endpoint, not a clinical diagnostic endpoint. FASTQ/SRA processing remains out of scope unless processed matrices are unavailable and a separate manual preprocessing plan is approved.

## Phase 26 GSE184950 Nested Archive Replication Route

Phase 26 turns GSE184950 from a candidate into an executable replication route while preserving safety. The nested per-sample archives inside `GSE184950_RAW.tar` are inspected without extraction. If processed 10x matrices are present, the user may manually enable selective extraction of only `matrix.mtx.gz`, `features.tsv.gz` or `genes.tsv.gz`, and `barcodes.tsv.gz`.

The resulting axis-gene extraction aggregates directly to sample-level mean expression and detection rate before endpoint-locked PD/PDD-vs-control statistics are computed. Any replication interpretation remains candidate or preliminary unless effects are directionally consistent, adequately covered, and statistically stable.

## Phase 27 GSE184950 Clean Replication Integration

Phase 27 cleans the GSE184950 sample-level outputs so only the 34 biological samples from the series matrix are retained. The spurious `processed_matrices` pseudo-sample is treated as a technical reporting artifact and excluded.

For PNAS-facing interpretation, direction-only agreement is not enough. GSE184950 axes are considered statistically supported replication only when direction is consistent and p-value or FDR support is adequate. Weak FDR support should be described as independent PD replication feasibility, not validated shared AD/PD mechanism.

## Phase 28 Independent AD Replication Onboarding

Phase 28 starts the missing AD replication layer. The priority target is GSE174367 because NeuroFate-Axis needs an independent donor/sample-level AD cohort to test whether SEA-AD dementia-associated axes replicate outside the discovery cohort. GSE147528 and GSE157827 are secondary or optional AD cohorts.

The preferred first route is bulk/sample-level expression, because endpoint-locked NeuroFate-Axis claims are donor/sample-level. Single-nucleus resources are planned only through processed-matrix templates, with no raw sequence processing. PNAS readiness remains blocked until at least one independent AD cohort shows statistically supported endpoint-locked replication.

## Phase 29 GSE174367 Bulk AD Replication Lane

Phase 29 implements the preferred GSE174367 bulk RNA route. The series matrix supplies the sample-level endpoint context, with 118 AD and 112 Control records, while `GSE174367_bulkRNA_processed.rda.gz` is inspected and converted only to a NeuroFate axis-gene matrix.

The Phase 29 endpoint is locked as AD versus Control through `label__ad_vs_control`. Axis effects are compared with the SEA-AD Phase 22 dementia endpoint direction, and direction-only agreement remains preliminary. A claim upgrade requires direction consistency plus statistical support. The single-nucleus H5 resource remains optional for a later processed-matrix route if bulk RNA is insufficient.

## Phase 32 Cross-Cohort Evidence Consolidation

Phase 32 reframes the current evidence after GSE174367 bulk RNA replication. NeuroFate-Axis now has independent AD replication evidence, but it is nominal rather than FDR-robust. The neuronal vulnerability axis is the leading AD replicated candidate because it is directionally consistent between SEA-AD and GSE174367.

The PD side remains the bottleneck. GSE243639 provides a preliminary PD extension, and GSE184950 provides clean independent PD/PDD replication infrastructure, but current axis-level support is not strong enough for a shared AD/PD biological claim. The next PNAS-critical step is stronger PD axis replication or an additional PD cohort with better endpoint/pathology support.

## Conservative Submission Position

The project can now be written as a PNAS-style systems-biology resource and discovery paper if the central claim is narrowed to the evidence:

- NeuroFate-Axis provides an endpoint-locked, donor/sample-level framework for neurodegeneration axis discovery.
- SEA-AD identifies robust AD-associated axes under a locked dementia endpoint.
- GSE174367 provides nominal independent AD replication of the neuronal vulnerability axis.
- GSE243639 and GSE184950 provide preliminary PD convergence and replication infrastructure, not a confirmed shared AD/PD mechanism.

The stronger statement that NeuroFate-Axis has established a replicated shared AD/PD mechanism remains blocked. The correct next biological move is stronger independent PD axis replication, while the correct manuscript move is to submit the conservative AD-replication-centered story and do not claim clinical utility, diagnostic utility, causality, or definitive shared-mechanism evidence.

## Phase 33 PD Replication Expansion

Phase 33 targets the remaining biological bottleneck for a stronger PNAS claim: independent PD axis replication. The strategy deliberately shifts from additional single-nucleus discovery toward public sample-level PD expression cohorts that can test NeuroFate axes without raw sequence processing.

The first recommended dataset is `GSE20141`, because laser-captured substantia nigra pars compacta neurons provide a direct donor/sample-level test of the neuronal vulnerability axis. `GSE20186` and related substantia nigra bulk/microarray cohorts are retained as follow-up replication sources if platform annotation can safely map probes to NeuroFate axis genes.

Evidence rules remain conservative:

- a PD axis is considered statistically supported only with direction consistency and `p < 0.05` or `FDR < 0.1`;
- direction-only PD evidence is labelled preliminary;
- shared AD/PD language remains blocked until AD and PD evidence converge on the same axis with adequate independent support;
- clinical, diagnostic, causal, or definitive shared-mechanism claims remain out of scope.

## Phase 34 PD Microarray/Bulk Expansion

Phase 34 narrows the next experimental move to public sample-level PD datasets that can be analyzed rapidly and conservatively. The practical priority is `GSE20141`, followed by `GSE7621`, `GSE8397`, and clear `GSE20186` subseries.

The key scientific question is whether the AD-replicated neuronal vulnerability axis shows statistically supported directionally consistent effects in an independent PD substantia nigra/SNpc cohort. If Phase 34 produces only direction-only support, the manuscript should remain in the current conservative position. If Phase 34 produces p/FDR-supported PD replication aligned with the AD-replicated axis, the claim can move to a replicated candidate cross-disease axis while still avoiding definitive mechanism language.

## Phase 37 GSE7621 PD Replication

Phase 37 moves from the technically successful but statistically weak GSE20141 result to the next PD cohort, `GSE7621`. The cohort is used only through sample-level GEO series-matrix expression values and conservative platform probe mappings.

The decision rule does not change: an independent PD replication result can improve the PNAS bottleneck only when an axis is directionally consistent and has `p < 0.05` or `FDR < 0.1`. Direction-only support remains preliminary. A shared AD/PD axis remains a candidate at most unless the PD-supported axis aligns with an AD-supported or AD-replicated axis.
