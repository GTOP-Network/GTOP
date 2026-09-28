# GTOP Analysis Pipeline — Detailed Documentation

## 1. Scope

This document provides a detailed map of the computational workflows in `analysis_pipeline/`, including the updated transcript-discovery, ASE, QTL, frequency-difference, portability, colocalization, and SMR workflows.

## 2. Workflow map

| Stage | Directory | Major purpose |
|---|---|---|
| Transcript discovery | `00_transcript_detection/01_transcript_discovery/` | Seven complementary long-read transcript discovery approaches |
| Transcript integration | `00_transcript_detection/02_merge_filter_transcript/` | FLNC support filtering, TAMA merging, SQANTI3 annotation, enhanced GTF construction |
| Transcript quantification | `00_transcript_detection/03_quantification/` | FLAIR transcript quantification |
| Peptide validation | `00_transcript_detection/04_peptide_validation/` | DIA-NN-based proteomic support |
| RNA phenotype preparation | `01_data_preparation/` | Gene, splice-junction, and transcript-level phenotypes |
| Variant calling | `02_Variant_calling/` | LRS/SRS small variants, population structure, and VEP annotation |
| ASE/ASTS/ASJ | `03_ASE_ASTS/` | Short-read ASE, long-read ASE/ASTS, allele-specific splicing, and ASJ |
| QTL mapping | `04_QTL_mapping/` | eQTL, sQTL, TR-xQTL, and MAJIQTL |
| Fine-mapping | `05_finemap/` | SuSiE and cross-ancestry SuSiEx |
| Frequency differences | `06_frequency_differences/` | Frequency-differentiated QTLs |
| eQTL portability | `07_eQTL_portability/` | GTOP–GTEx portability and mashR |
| Enrichment | `08_QTL_enrichment/` | torus and S-LDSC |
| Colocalization | `09_colocalization/` | GWAS/QTL fine-mapping and colocalization |
| SMR | `10_SMR/` | Summary-data-based Mendelian randomization |

## 3. Transcript discovery and validation

### 3.1 Transcript discovery

`00_transcript_detection/01_transcript_discovery/` contains seven discovery workflows:

| Tool | Directory |
|---|---|
| Bambu | `bambu/` |
| FLAIR | `flair/` |
| FLAMES | `flames/` |
| IsoQuant | `isoquant/` |
| Iso-Seq | `isoseq/` |
| IsoTools | `isotools/` |
| TALON | `talon/` |

Each workflow provides pipeline/run scripts where applicable. The discovery outputs are subsequently integrated rather than treating a single caller as the sole transcript reference.

### 3.2 Transcript integration and filtering

`02_merge_filter_transcript/` contains workflows for FLNC support filtering, raw transcript QC, chromosome-wise processing, TAMA-based merging, cross-tool transcript merging, SQANTI3 structural annotation, FLAIR quantification, alternative first-exon analysis, short-read junction analysis, and enhanced GTF construction.

### 3.3 Quantification and proteomic validation

`03_quantification/` contains the FLAIR-based transcript quantification workflow. `04_peptide_validation/` prepares transcript-derived protein FASTA files and runs DIA-NN-related workflows to evaluate peptide support.

## 4. Molecular phenotype preparation

### Short-read gene expression

`01_data_preparation/SRS_RNA_quantification/` implements FastQC → STAR alignment → RNA-SeQC gene-level quantification.

### Splicing

`01_data_preparation/splicing_quantification/` implements alignment, intron-usage extraction, phenotype preprocessing, and LeafCutter/TensorQTL preparation.

### Long-read transcript quantification

`01_data_preparation/transcript_quantification/` contains Salmon/RSEM workflows and utilities for transcript-level molecular phenotypes.

## 5. Variant calling

`02_Variant_calling/` provides separate LRS and SRS small-variant workflows, followed by PCA/ADMIXTURE population analysis and VEP functional annotation.

## 6. Allele-specific analyses

`03_ASE_ASTS/` provides complementary short-read and long-read approaches.

### lorals

`run_LRS_ase_lorals.sh` performs haplotype-aware long-read mapping and quantifies ASE and ASTS.

### isoLASER

The four `isolaser/run.step*.sh` scripts prepare transcript references, annotate BAM files, create sample metadata, and run isoLASER jointly for allele-specific splicing.

### Short-read ASE

`run_SRS_ase.sh` performs variant-aware short-read ASE analysis with allele-count aggregation and quality control.

### longcallR

`run_longcallR.sh` provides an additional long-read workflow:

```text
RNA reads/BAM → longcallR SNP calling/phasing → phased BAM + RNA VCF
                                            ├→ high-confidence ASE
                                            └→ high-confidence ASJ
```

Matched DNA variants are used for high-confidence ASE/ASJ calls.

## 7. QTL mapping

`04_QTL_mapping/01_data_preparation/` prepares phenotypes, covariates, and genotypes. Association workflows are:

```text
run_eQTL_mapping.sh
run_sQTL_mapping.sh
run_TR_xQTL_mapping.sh
run_MAJIQTL.sh
```

The TR workflows additionally provide TensorQTL implementations for TR-eQTL and TR-sQTL mapping. `run_MAJIQTL.sh` provides a MAJIQTL-based splicing-QTL workflow using SNP, SV, and TR VCF inputs and tissue-specific splicing phenotypes/covariates.

## 8. Fine-mapping

### SuSiE

`05_finemap/SuSiE_Fine_mapping/` contains gene-level genotype/phenotype preparation, standard SuSiE fine-mapping, joint fine-mapping, result summarization, and SV/TR annotation. Main wrappers are `run.sh` and `run_joint.sh`.

### Cross-ancestry SuSiEx

`05_finemap/Cross_ancestries_Fine_mapping_SuSiEX/` harmonizes GTOP and external ancestry-specific QTL datasets, identifies shared eGenes, updates SNP metadata, and summarizes SuSiEx results.

## 9. Frequency-differentiated QTLs

`06_frequency_differences/` follows four sequential stages:

1. `a1.gtop_gnomad_af_comparison.R`: allele-frequency comparison;
2. `a2.qtl_fine_mapping_summary.R`: fine-mapping and LD-expanded QTL summary;
3. `a3.qtl_mege_to_loci.R`: cross-tissue credible-set/locus aggregation;
4. `a4.fdQTL_table.R`: fdQTL identification and table generation.

## 10. eQTL portability and mashR

`07_eQTL_portability/` evaluates GTOP–GTEx overlap, mash-adjusted cross-tissue effects, portability using multiple metrics, and downstream mash analyses. The `MashR/` directory provides workflows for strong pairs, random pairs, and fine-mapping pairs.

## 11. QTL enrichment

`08_QTL_enrichment/` implements torus-based genomic annotation enrichment and S-LDSC-based heritability enrichment. The torus workflow converts TensorQTL nominal results to torus input; the S-LDSC workflow constructs annotations from significant and fine-mapped QTL variants, calculates LD scores, and estimates enrichment.

## 12. GWAS-QTL colocalization

`09_colocalization/` follows:

```text
GWAS fine-mapping → gene/locus preparation → GTOP QTL fine-mapping → SuSiE-coloc → result merging
```

The principal scripts are `a1.finemapping_for_GWAS.R`, `a2.prepare_genes.R`, `a3.finemapping_for_QTL_GTOP_LD.R`, `a4.coloc.R`, and `a5.merge_result.R`, coordinated by `run.sh`.

## 13. SMR

`10_SMR/run_SMR.sh` provides the SMR workflow for eQTL, sQTL, and transcript-usage QTL summary statistics. Study-specific paths and SMR installations must be adapted before execution.

## 14. Practical reproducibility notes

Complete reproduction requires the corresponding input datasets and reference resources. Before execution, inspect the module README and entry-point script, replace absolute paths, activate the required environment, confirm reference versions, verify consistent sample IDs, and retain intermediate files required downstream.
