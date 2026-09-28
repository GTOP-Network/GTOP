# GTOP Analysis Pipeline

This directory contains the computational workflows used to process, analyze, and interpret the **GTOP Phase-I long-read multi-omics resource**.

The pipeline is organized into modular workflows covering long-read transcript discovery and validation, molecular phenotype preparation, genome-wide variant calling, allele-specific analyses, molecular QTL mapping, fine-mapping, frequency-differentiated QTLs, eQTL portability, enrichment analyses, GWAS-QTL colocalization, and SMR.

## Analysis workflow

```text
Long-read RNA-seq
        │
        ├── Transcript discovery
        │     ├── Bambu
        │     ├── FLAIR
        │     ├── FLAMES
        │     ├── IsoQuant
        │     ├── Iso-Seq
        │     ├── IsoTools
        │     └── TALON
        │
        ├── Transcript integration and filtering
        ├── Transcript quantification
        └── Proteomic validation
        │
        ▼
Genome-wide sequencing
        ├── long-read WGS variant calling
        ├── short-read WGS variant calling
        ├── population structure
        └── VEP annotation
        │
        ▼
Allele-specific analyses
        ├── short-read ASE
        ├── long-read ASE / ASTS with lorals
        ├── allele-specific splicing with isoLASER
        └── long-read ASE / ASJ with longcallR
        │
        ▼
Molecular QTL mapping
        ├── eQTL
        ├── sQTL
        ├── TR-eQTL
        ├── TR-sQTL
        └── MAJIQTL-based splicing QTL analysis
        │
        ▼
Fine-mapping
        ├── SuSiE
        └── cross-ancestry SuSiEx
        │
        ▼
Downstream genetic analyses
        ├── frequency-differentiated QTLs
        ├── eQTL portability and mashR
        ├── QTL enrichment
        ├── GWAS-QTL colocalization
        └── SMR
```

## Directory structure

| Module | Main purpose |
|---|---|
| `00_transcript_detection/` | Long-read transcript discovery, integration, filtering, quantification, and proteomic validation |
| `01_data_preparation/` | Short-read gene expression, splicing, and long-read transcript-level phenotype preparation |
| `02_Variant_calling/` | LRS/SRS variant calling, population structure analysis, and VEP annotation |
| `03_ASE_ASTS/` | ASE, ASTS, allele-specific splicing, and allele-specific junction analysis |
| `04_QTL_mapping/` | eQTL, sQTL, TR-eQTL, TR-sQTL, and MAJIQTL-based QTL analyses |
| `05_finemap/` | SuSiE and cross-ancestry SuSiEx fine-mapping |
| `06_frequency_differences/` | Frequency-differentiated QTL and ancestry-aware fine-mapping analyses |
| `07_eQTL_portability/` | GTOP–GTEx eQTL portability, cross-tissue sharing, and mashR analyses |
| `08_QTL_enrichment/` | Genomic annotation enrichment with torus and heritability enrichment with S-LDSC |
| `09_colocalization/` | GWAS and cis-QTL fine-mapping and colocalization |
| `10_SMR/` | Summary-data-based Mendelian randomization |

## Module overview

### 00. Transcript detection

`00_transcript_detection/` reconstructs and evaluates transcript models from PacBio long-read RNA-seq using **seven complementary approaches**: Bambu, FLAIR, FLAMES, IsoQuant, Iso-Seq, IsoTools, and TALON. The resulting models are merged, filtered, structurally annotated, quantified, and evaluated using proteomic evidence.

### 01. Data preparation

`01_data_preparation/` generates molecular phenotypes for downstream QTL analyses, including short-read gene expression, LeafCutter-based splicing phenotypes, and long-read transcript-level quantification using Salmon/RSEM-related workflows.

### 02. Variant calling

`02_Variant_calling/` processes both long- and short-read WGS data. The main workflows are `01_LRS_WGS_Variant_Calling.sh`, `02_SRS_WGS_Variant_Calling.sh`, `03_PCA_ADMIXTURE.sh`, and `04_vep_annotation.sh`.

### 03. ASE and ASTS

`03_ASE_ASTS/` contains short-read ASE, long-read ASE/ASTS using lorals, allele-specific splicing using isoLASER, and long-read ASE/ASJ using longcallR. The longcallR workflow calls and phases RNA SNPs and uses matched DNA variants for high-confidence ASE/ASJ analysis.

### 04. QTL mapping

`04_QTL_mapping/` supports cis-eQTL, cis-sQTL, TR-eQTL, TR-sQTL, and MAJIQTL-based splicing QTL analysis. Main entry points include `run_eQTL_mapping.sh`, `run_sQTL_mapping.sh`, `run_TR_xQTL_mapping.sh`, and `run_MAJIQTL.sh`.

### 05. Fine-mapping

`05_finemap/` contains SuSiE fine-mapping and cross-ancestry SuSiEx workflows, including gene-centered genotype/phenotype preparation, credible-set/PIP summarization, and SV/TR annotation.

### 06. Frequency-differentiated QTLs

`06_frequency_differences/` compares allele frequencies across ancestral populations, summarizes fine-mapped QTLs and LD-expanded variants, merges credible sets across tissues, and identifies frequency-differentiated QTLs through the sequential `a1`–`a4` workflow.

### 07. eQTL portability

`07_eQTL_portability/` evaluates GTOP–GTEx eQTL sharing and portability and uses mashR to model cross-tissue effect patterns. The `MashR/` subdirectory contains strong-pair, random-pair, and fine-mapping-pair workflows.

### 08. QTL enrichment

`08_QTL_enrichment/` implements torus-based genomic annotation enrichment and S-LDSC-based heritability enrichment.

### 09. GWAS-QTL colocalization

`09_colocalization/` performs GWAS fine-mapping, GWAS–gene/locus pairing, GTOP QTL fine-mapping, SuSiE-coloc, and result integration.

### 10. SMR

`10_SMR/` contains the workflow for integrating eQTL, sQTL, and transcript-usage QTL summary statistics with GWAS summary statistics using SMR.

## Reproducibility

The workflows use Bash, Python, and R and are designed primarily for HPC environments. Many scripts contain study-specific paths, software environments, and scheduler settings that should be adapted before execution. Most modules require external reference resources such as genome assemblies, annotations, genotype data, LD resources, phenotype matrices, and GWAS/QTL summary statistics.

## Recommended execution order

```text
00 → 01 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → 10
```

Individual downstream modules can also be run independently when their required upstream inputs are available.

## Key software

The workflows incorporate, among others: Bambu, FLAIR, FLAMES, IsoQuant, Iso-Seq, IsoTools, TALON, TAMA, SQANTI3, Salmon, RSEM, LeafCutter, TensorQTL, MAJIQTL, lorals, isoLASER, longcallR, SuSiE, SuSiEx, mashR, torus, LDSC/S-LDSC, SuSiE-coloc, and SMR.

## Citation

If you use the GTOP analysis workflows or derived results, please cite the GTOP study and the original software packages used in the corresponding analyses.
