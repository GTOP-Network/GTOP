# GTOP Figure Generation

This directory contains the R scripts and processed input files used to generate the main, extended, and supplementary figures for the GTOP Phase-I manuscript.

## Directory structure

```text
figure_generation/
├── Figure_1/
├── Figure_2/
├── Figure_3/
├── Figure_4/
├── Figure_5/
├── Figure_6/
├── Extended_Figures/
└── Supp_Figures/
```

## Main figures

| Figure | Script | Main content |
|---|---|---|
| Figure 1 | `Figure_1/Figure1.R` | Population/sample structure and genome-wide genetic variation |
| Figure 2 | `Figure_2/Figure2.R` | Long-read transcriptome characterization, peptide support, and ASE/ASTS |
| Figure 3 | `Figure_3/Figure3.R` | Molecular QTL landscape across variant classes |
| Figure 4 | `Figure_4/Figure4.R` | Population-specific QTLs, fine-mapping, frequency-differentiated QTLs, and eQTL portability |
| Figure 5 | `Figure_5/Figure5.R` | SV/TR tagging and joint fine-mapping |
| Figure 6 | `Figure_6/Figure6.R` | Disease association, enrichment, colocalization, and SV/TR disease-related signals |

### Figure 1

Covers genotype PCA, RNA-seq sample structure, tissue clustering, variant number/length summaries, per-genome variant counts, SV comparisons, and LRS/SRS TR comparisons.

### Figure 2

Covers long-/short-read sample structure, novel transcript discovery, SQANTI3 annotation/coding status, alternative splicing, transcript tissue breadth, peptide support, WGCNA/GO analysis, and ASE/ASTS summaries.

### Figure 3

Covers eQTL/sQTL summaries, fine-mapped QTL distance to TSS, GTOP–GTEx effect-size comparison, eGene sharing across variant classes, SV effect size versus length, pathogenic TR-QTLs, representative junction signals, and MASH-based analyses.

### Figure 4

Covers credible-set size/PIP comparisons, frequency-differentiated QTLs, representative fdQTL loci, and eQTL portability.

### Figure 5

Evaluates SV/TR tagging by small variants, poorly tagged signals, and the composition/enrichment of joint fine-mapping credible sets.

### Figure 6

Integrates S-LDSC enrichment, GWAS/QTL colocalization, comparisons with GTEx/JCTF/MAGE, GTOP-specific signals, GWAS enrichment, rare-variant analyses, disease-associated SV/TR signals, and representative SV-associated loci.

## Extended figures

| Figure | Main content |
|---|---|
| Extended Figure 1 | Genetic variation and pathogenic TR-related analyses |
| Extended Figure 2 | Additional long-read/short-read transcriptome characterization |
| Extended Figure 3 | Long-read ASE and ASTS analyses |
| Extended Figure 4 | Molecular QTL characteristics and functional enrichment |
| Extended Figure 5 | QTL tissue sharing and tissue specificity |
| Extended Figure 6 | Cross-ancestry fine-mapping and MPRA-related analyses |
| Extended Figure 7 | Frequency-differentiated sQTLs and mash-based eQTL portability |
| Extended Figure 8 | SV-sQTL and TR-sQTL signals poorly tagged by small variants |
| Extended Figure 9 | SMR and GWAS/QTL colocalization analyses |

## Supplementary figures

`Supp_Figures/` contains the plotting scripts currently distributed for Supplementary Figures 2–42. The repository does not necessarily contain a script for every integer in that range; the table below reflects the scripts actually included.

| Figure | Description |
|---|---|
| S2 | Whole-genome sequencing statistics |
| S3 | Short-read RNA-seq quality control |
| S4 | Long-read RNA-seq quality |
| S5 | ADMIXTURE ancestry components |
| S6 | LRS/SRS variant detection comparison |
| S7 | Variant discovery and genomic variation |
| S8 | VEP genomic annotation |
| S11 | Consensus transcript discovery using seven tools |
| S12 | Transcript discovery and filtering |
| S14 | Peptide support and predicted protein structures |
| S15 | Cross-platform gene/transcript expression consistency |
| S16 | Iso-Seq sampling-depth evaluation |
| S17 | Transcript-level co-expression modules |
| S18 | Genotype phasing comparison |
| S19 | ASE and ASTS across tissues |
| S20 | Allele-specific gene expression and splicing |
| S21 | Allele-specific splicing |
| S22 | PCA-based gene-expression covariate selection |
| S23 | PCA-based splicing covariate selection |
| S24 | PCA-based transcript-usage covariate selection |
| S25 | Inferred factors and measured covariates |
| S26 | Cross-tissue sQTL consistency |
| S27 | SV/TR-associated QTL representation and tagging |
| S28 | sQTL enrichment in genomic features |
| S29 | eQTL effect-size correlation across pancreas sites |
| S30 | eQTL and ASE effect-size correlation |
| S31 | GTOP–GTEx eQTL effect-size correlation |
| S32 | Representative SV associated with IRGM expression |
| S34 | GTOP–GTEx eQTL sharing and portability |
| S35 | Effects of allele frequency and sample size on portability |
| S36 | Fine-mapped sQTL credible sets containing SV/TR lead variants |
| S37 | Joint fine-mapping after coverage filtering |
| S38 | Robustness of SV/TR prioritization |
| S39 | Representative SV prioritized by joint fine-mapping |
| S40 | Heritability enrichment and GWAS colocalization |
| S41 | GWAS loci with and without cis-eQTL colocalization |
| S42 | Frequency-differentiated QTLs at disease-associated loci |

## Inputs and reproducibility

Processed figure inputs are stored in the corresponding `input/` directories. They include text tables, compressed tables, `.RDS` objects, and `.RData` objects.

```text
analysis_pipeline/
       ↓
processed analysis results
       ↓
figure-specific input/
       ↓
R plotting script
       ↓
manuscript figure
```

Most scripts can be run with:

```bash
Rscript Figure3.R
```

after adapting local paths and installing the packages required by the target script.

Common dependencies include `ggplot2`, `data.table`, `dplyr`, `tidyverse`, `ggpubr`, `patchwork`, `cowplot`, `ComplexHeatmap`, `locuszoomr`, `EnsDb.Hsapiens.v86`, `ggvenn`, `ggradar`, and `ggrastr`. Individual scripts may require additional packages.
