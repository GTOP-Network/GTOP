# GTOP Analysis Pipeline

This directory contains the computational workflows used to process, analyze, and interpret the **GTOP Phase-I long-read multi-omics resource**.

The pipeline is organized as a sequence of modular workflows covering long-read transcript discovery, transcript annotation and quantification, genome-wide variant calling, allele-specific analyses, molecular QTL mapping, fine-mapping, ancestry-aware QTL analyses, enrichment analyses, GWAS-QTL colocalization, and SMR.

## Analysis workflow

```text
Long-read RNA-seq
        │
        ├── Transcript discovery
        │     ├── Bambu
        │     ├── FLAIR
        │     └── Iso-Seq
        │
        ├── Transcript integration and filtering
        │     ├── transcript merging
        │     ├── FLAIR quantification
        │     └── SQANTI3 annotation
        │
        ├── Transcript quantification
        └── Peptide validation
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
        ├── long-read ASE
        ├── ASTS
        └── allele-specific splicing
        │
        ▼
Molecular QTL mapping
        ├── eQTL
        ├── sQTL
        ├── TR-eQTL
        └── TR-sQTL
        │
        ▼
Fine-mapping
        ├── SuSiE
        └── cross-ancestry SuSiEx
        │
        ▼
Downstream analyses
        ├── population-specific QTLs
        ├── QTL enrichment
        ├── GWAS-QTL colocalization
        └── SMR
```

## Directory structure

| Module | Main purpose |
|---|---|
| `00_transcript_detection/` | Long-read transcript discovery, integration, filtering, quantification, and peptide validation |
| `01_data_preparation/` | Preparation of gene, splicing, and transcript-level molecular phenotypes |
| `02_Variant_calling/` | Variant discovery, population structure analysis, and functional annotation |
| `03_ASE_ASTS/` | Allele-specific expression, allele-specific transcript structure, and allele-specific splicing |
| `04_QTL_mapping/` | eQTL, sQTL, TR-eQTL, and TR-sQTL mapping |
| `05_finemap/` | SuSiE and cross-ancestry fine-mapping |
| `06_population-specific_eQTL/` | Ancestry-dependent allele frequencies, QTL portability, heterogeneity QTLs, and mashR analyses |
| `07_QTL_enrichment/` | Genomic annotation and heritability enrichment |
| `08_colocalization/` | GWAS fine-mapping and GWAS-QTL colocalization |
| `09_SMR/` | Summary-data-based Mendelian randomization |
| `Pipeline.md` | Detailed overview of the analysis pipeline |

## Module overview

### 00. Transcript detection

`00_transcript_detection/` reconstructs and evaluates transcript models from long-read RNA-seq.

The workflow integrates:

- **Bambu** for reference-guided transcript discovery.
- **FLAIR** for reference-guided isoform reconstruction and quantification.
- **Iso-Seq** for transcript reconstruction.
- **SQANTI3** for transcript structural annotation and quality assessment.
- Custom scripts for transcript merging, filtering, gene assignment, and construction of enhanced transcript references.
- **DIANN**-based workflows for peptide support/validation.

The main submodules are:

```text
00_transcript_detection/
├── 01_transcript_discovery/
├── 02_merge_filter_transcript/
├── 03_quantification/
└── 04_peptide_validation/
```

### 01. Data preparation

`01_data_preparation/` generates molecular phenotypes for downstream QTL analyses.

It includes:

- Short-read RNA-seq gene-expression quantification.
- Splicing quantification using LeafCutter-based workflows.
- Long-read transcript quantification using Salmon/RSEM-related workflows.
- Construction of molecular phenotype matrices.

### 02. Variant calling

`02_Variant_calling/` processes long- and short-read WGS data and generates variant resources for downstream analyses.

Main workflows:

```text
01_LRS_WGS_Variant_Calling.sh
02_SRS_WGS_Variant_Calling.sh
03_PCA_ADMIXTURE.sh
04_vep_annotation.sh
```

These workflows cover genome-wide variant discovery, population structure analysis using PCA/ADMIXTURE, and functional annotation with VEP.

### 03. ASE and ASTS

`03_ASE_ASTS/` identifies allele-specific regulatory effects.

The workflow includes:

- Short-read ASE.
- Long-read ASE.
- Allele-specific transcript structure (ASTS).
- Allele-specific splicing using `isoLASER`.
- Haplotype-aware alignment and donor-level aggregation.

The main entry points are:

```text
run_SRS_ase.sh
run_LRS_ase_lorals.sh
```

### 04. QTL mapping

`04_QTL_mapping/` performs molecular QTL mapping.

Supported analyses include:

- cis-eQTL
- cis-sQTL
- TR-eQTL
- TR-sQTL

The workflow first prepares phenotypes, covariates, and genotypes, followed by QTL association testing. TensorQTL-based Python workflows are provided for tandem-repeat QTL analyses.

### 05. Fine-mapping

`05_finemap/` performs variant-level fine-mapping.

Two major workflows are provided:

1. **SuSiE fine-mapping**
   - Gene/tissue-level genotype and phenotype preparation.
   - Variant-class-specific and joint fine-mapping.
   - Credible-set and PIP summarization.
   - SV/TR annotation and downstream summaries.

2. **Cross-ancestry fine-mapping with SuSiEx**
   - Preparation of GTOP and external ancestry-specific QTL inputs.
   - Shared eGene identification.
   - Cross-ancestry fine-mapping.
   - Summary and harmonization of fine-mapping results.

### 06. Population-specific eQTL

`06_population-specific_eQTL/` investigates ancestry-dependent regulatory effects.

The workflow includes:

- Allele-frequency comparisons across ancestry groups.
- Fine-mapped QTL frequency annotation.
- Replication and comparison with GTEx.
- Identification of heterogeneity QTLs.
- mashR-based multivariate effect-size modeling and fine-mapping.

Scripts are organized into analytical stages `a1`–`a4`, with additional workflows under `MashR/`.

### 07. QTL enrichment

`07_QTL_enrichment/` evaluates the functional enrichment of QTLs.

Two complementary approaches are implemented:

- **torus** for enrichment of QTL associations in genomic annotations.
- **S-LDSC** for enrichment of complex-trait heritability in QTL-derived annotations.

### 08. GWAS-QTL colocalization

`08_colocalization/` integrates GTOP molecular QTLs with GWAS signals.

The workflow includes:

1. GWAS locus definition and annotation.
2. GWAS fine-mapping.
3. QTL fine-mapping preparation.
4. GWAS-QTL colocalization.
5. Result integration and summarization.

### 09. SMR

`09_SMR/` contains the SMR workflow for integrating molecular QTL summary statistics with GWAS summary statistics.

## Reproducibility

The workflows were developed for high-performance computing environments and use a combination of **Bash, Python, and R**.

Most workflows require external reference data, including genome assemblies, gene annotations, genotype files, phenotype matrices, GWAS summary statistics, and/or external QTL resources. Paths in shell and R scripts should therefore be adapted to the local environment before execution.

For a complete description of a specific workflow, consult:

- the module-level README, where available;
- the corresponding shell/Python/R entry-point script;
- `Pipeline.md` for the overall analysis order.

## Recommended execution order

For reproducing the complete GTOP analysis, the modules are generally intended to be followed in this order:

```text
00  Transcript detection
 ↓
01  Molecular phenotype preparation
 ↓
02  Variant calling and annotation
 ↓
03  ASE / ASTS
 ↓
04  QTL mapping
 ↓
05  Fine-mapping
 ↓
06  Population-specific QTL analyses
 ↓
07  QTL enrichment
 ↓
08  GWAS-QTL colocalization
 ↓
09  SMR
```

Individual downstream modules can also be run independently when the required upstream inputs are available.

## Software environment

The repository combines R, Python, and command-line workflows and relies on third-party software packages used by individual modules. Because the exact software requirements can differ between modules, users should inspect the corresponding scripts before execution and install the relevant tools in their HPC environment.

Key software/tools represented in the workflows include Bambu, FLAIR, Iso-Seq, SQANTI3, Salmon, RSEM, LeafCutter, TensorQTL, SuSiE, SuSiEx, mashR, torus, S-LDSC, locus-level colocalization workflows, and SMR.

## Notes

The scripts in this repository document the computational workflows used for the GTOP study. They are provided primarily to facilitate transparency and reproducibility of the published analyses. Large intermediate files, reference resources, and computationally intensive outputs are not necessarily distributed with the repository.

## Citation

If you use the GTOP analysis workflows or derived results, please cite the GTOP study and the original software packages used in the corresponding analyses.

> Zou X, Li X, Zhang T, et al. *A long-read multi-omics atlas broadens discovery of regulatory variation across 33 human tissues.*
