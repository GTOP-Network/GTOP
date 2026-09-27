# GTOP Analysis Pipeline — Detailed Documentation

## 1. Scope

This document provides a detailed map of the computational workflows in `analysis_pipeline/`. It is intended to help users identify the appropriate entry point, understand dependencies between modules, and reproduce individual analyses without requiring the entire pipeline to be rerun.

## 2. Workflow map

| Stage | Directory | Major outputs / purpose |
|---|---|---|
| Transcript discovery | `00_transcript_detection/01_transcript_discovery/` | Transcript models from Bambu, FLAIR and Iso-Seq |
| Transcript integration | `00_transcript_detection/02_merge_filter_transcript/` | Integrated and filtered transcript annotation |
| Transcript quantification | `00_transcript_detection/03_quantification/` | Transcript-level abundance |
| Peptide validation | `00_transcript_detection/04_peptide_validation/` | Proteomic support for translated transcripts |
| RNA phenotype preparation | `01_data_preparation/` | Gene, splice-junction and transcript phenotypes |
| Variant calling | `02_Variant_calling/` | Small variants, SV/TR resources, population structure and annotation |
| ASE/ASTS | `03_ASE_ASTS/` | Allele-specific expression and transcript-structure signals |
| QTL mapping | `04_QTL_mapping/` | eQTL, sQTL, TR-eQTL and TR-sQTL associations |
| Fine-mapping | `05_finemap/` | Credible sets and PIP estimates |
| Population-specific QTL | `06_population-specific_eQTL/` | Frequency differentiation, portability and heterogeneity |
| Enrichment | `07_QTL_enrichment/` | Genomic annotation and heritability enrichment |
| Colocalization | `08_colocalization/` | GWAS-QTL colocalization |
| SMR | `09_SMR/` | Summary-data-based Mendelian randomization |

## 3. Transcript detection

### 3.1 Transcript discovery

`01_transcript_discovery/` contains three complementary discovery workflows:

- `bambu/`: Bambu-based transcript discovery and preprocessing.
- `flair/`: FLAIR read processing, splice correction and transcript reconstruction.
- `isoseq/`: Python implementation of the Iso-Seq processing workflow.

The three approaches generate complementary transcript models that are subsequently integrated.

### 3.2 Transcript integration and filtering

`02_merge_filter_transcript/` contains custom Python workflows for:

- merging transcript models by intron-chain identity;
- assigning transcripts to genes;
- integrating outputs from different discovery tools;
- FLAIR-based transcript quantification;
- SQANTI3-based structural annotation;
- generation of enhanced GTOP transcript references.

### 3.3 Quantification

`03_quantification/` contains FLAIR-based transcript quantification workflows used after construction of the transcript reference.

### 3.4 Peptide validation

`04_peptide_validation/` prepares transcript-derived FASTA files and runs the proteomics workflow used to assess peptide support.

## 4. Molecular phenotype preparation

### 4.1 Short-read RNA-seq

`01_data_preparation/SRS_RNA_quantification/` contains the short-read gene-expression quantification entry point.

### 4.2 Splicing

`01_data_preparation/splicing_quantification/` contains sequential scripts for:

1. read alignment;
2. intron-usage extraction;
3. phenotype mapping/preparation;
4. preparation of LeafCutter-derived phenotypes for TensorQTL.

### 4.3 Long-read transcript expression

`01_data_preparation/transcript_quantification/` contains Salmon/RSEM-related Python workflows and utilities for generating transcript-level molecular phenotypes.

## 5. Variant calling

`02_Variant_calling/` contains four major entry points:

- `01_LRS_WGS_Variant_Calling.sh`
- `02_SRS_WGS_Variant_Calling.sh`
- `03_PCA_ADMIXTURE.sh`
- `04_vep_annotation.sh`

The resulting variant resources support population analysis, allele-specific analyses, QTL mapping and fine-mapping.

## 6. Allele-specific analyses

`03_ASE_ASTS/` contains separate short-read and long-read workflows.

The long-read branch includes haplotype preparation, allele-aware alignment, allelic coverage calculation and donor-level aggregation. The `isolaser/` subdirectory provides the multi-step isoLASER workflow for allele-specific splicing.

## 7. QTL mapping

Before association testing, `04_QTL_mapping/01_data_preparation/` prepares:

- molecular phenotypes;
- covariates;
- tissue-specific expression matrices;
- genotype inputs;
- residualized phenotypes where applicable.

The top-level scripts then perform:

```text
run_eQTL_mapping.sh
run_sQTL_mapping.sh
run_TR_xQTL_mapping.sh
```

The two TensorQTL Python scripts provide tandem-repeat-specific association workflows.

## 8. Fine-mapping

### 8.1 SuSiE

`05_finemap/SuSiE_Fine_mapping/` contains the main fine-mapping implementation.

Important workflow components include:

- `prepare_genotype_by_gene*.R`: gene-centered genotype preparation;
- `prepare_pheno_by_gene*.R`: phenotype preparation;
- `generate_joint_GT_by_gene.R`: joint genotype construction;
- `finemapping.R`: fine-mapping;
- `summarize_susie_results_by_tissue*.R`: result summarization;
- `anno_SV_in_finemap_res.R`: SV annotation;
- TR-specific preparation and normalization scripts.

`run.sh` and `run_joint.sh` serve as the principal execution entry points for the corresponding workflows.

### 8.2 Cross-ancestry fine-mapping

`Cross_ancestries_Fine_mapping_SuSiEX/` prepares GTOP and external ancestry-specific QTL datasets for SuSiEx-based cross-ancestry fine-mapping. Supporting scripts harmonize SNP metadata, identify shared eGenes, update RSID versions, and summarize results.

## 9. Population-specific eQTL

The `a1`–`a4` scripts implement successive analyses of:

- ancestry-specific allele frequencies;
- frequency annotation of fine-mapped variants;
- replication in GTEx;
- comparison of GTOP and GTEx effect sizes;
- heterogeneity QTL preparation and calculation;
- mashR-based multivariate QTL modeling and fine-mapping.

## 10. Enrichment, colocalization and SMR

The final modules connect GTOP regulatory variation to functional genomic annotations and complex traits.

`07_QTL_enrichment/` provides torus and S-LDSC workflows.

`08_colocalization/` separates GWAS locus preparation, GWAS fine-mapping, QTL fine-mapping, colocalization, and result merging.

`09_SMR/run_SMR.sh` provides the SMR integration workflow.

## 11. Practical reproducibility notes

The repository contains analysis logic and executable workflows, but complete reproduction requires the corresponding input datasets and reference resources. Users should replace absolute paths and environment-specific settings in the scripts before execution.

For large-scale analyses, it is recommended to execute one module at a time and retain the intermediate files required by downstream modules.

