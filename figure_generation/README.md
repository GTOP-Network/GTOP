# GTOP Figure Generation

This directory contains the R scripts and input tables used to generate the main, extended, and supplementary figures for the GTOP Phase-I manuscript.

The figure-generation code is organized by manuscript figure and is intended to make the published visualizations traceable to the underlying analysis outputs.

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
├── Supp_Figures/
└── shared plotting utilities
```

Each figure directory contains one or more R scripts and an `input/` directory containing the tables or serialized R objects required by the plotting code. Several directories also contain a figure-specific README.

## Main figures

| Figure | Main content |
|---|---|
| Figure 1 | GTOP cohort, tissue transcriptomic landscape, and genome-wide genetic variation |
| Figure 2 | Long-read transcriptome characterization and allele-specific analyses |
| Figure 3 | Molecular QTL landscape across variant classes |
| Figure 4 | Population-specific QTLs, fine-mapping, and QTL portability |
| Figure 5 | SV/TR tagging and joint fine-mapping |
| Figure 6 | Disease association, GWAS enrichment, colocalization and SMR |

### Figure 1

`Figure_1/Figure1.R`

Covers genotype population structure, RNA-seq sample structure, tissue clustering, genome-wide variant counts and lengths, variant counts per sample, SV concordance, and TR comparison between long- and short-read sequencing.

### Figure 2

`Figure_2/Figure2.R`

Covers transcript-level sample structure, transcript discovery relative to reference annotations, SQANTI3 transcript classification and coding potential, alternative splicing, tissue distribution of transcripts, peptide support, cardiac co-expression modules, ASE/ASTS summary statistics, and an example locus.

### Figure 3

`Figure_3/Figure3.R`

Covers eQTL/sQTL discovery by variant class, variant-to-TSS distance, GTOP–GTEx eQTL effect-size concordance, overlap among eGenes, SV-eQTL effect sizes by SV length, disease-associated tandem repeats, a representative splice-junction example, and tissue sharing of molecular QTLs.

### Figure 4

`Figure_4/Figure4.R`

Covers cross-ancestry fine-mapping, credible-set size and PIP, an example fine-mapped locus, geographic allele-frequency differences, GWAS/QTL colocalization examples, GTOP–GTEx allele-frequency comparisons, and eQTL portability.

### Figure 5

`Figure_5/Figure5.R`

Covers tagging of SV/TR eQTLs by nearby SNVs, the composition and size of joint fine-mapping credible sets, the distribution of credible sets across tissues, and enrichment of SV/TR variants among fine-mapped signals.

### Figure 6

`Figure_6/Figure6.R`

Covers GWAS heritability enrichment, colocalized GWAS loci, comparisons with GTEx/JCTF/MAGE, representative locus-level colocalization, GWAS enrichment, overlap between GWAS risk variants and QTLs, fine-mapped SV/TR colocalization, and representative sQTL loci.

## Extended figures

`Extended_Figures/` contains scripts for Extended Figures 1–9:

| Figure | Description |
|---|---|
| Extended Figure 1 | Comparison of genetic variants detected by long- and short-read sequencing |
| Extended Figure 2 | Transcript discovery, quantification and peptide validation |
| Extended Figure 3 | Long-read ASE and ASTS across tissues |
| Extended Figure 4 | Molecular QTL characteristics across variant classes |
| Extended Figure 5 | Tissue sharing and context-dependent QTL patterns |
| Extended Figure 6 | Cross-ancestry fine-mapping of eQTLs |
| Extended Figure 7 | Shared and frequency-differentiated regulatory effects |
| Extended Figure 8 | SV-sQTL and TR-sQTL signals poorly tagged by small variants |
| Extended Figure 9 | Characteristics of risk-variant-associated genes |

The directory also contains shared plotting utilities, including `multi-circle_function.R`.

## Supplementary figures

`Supp_Figures/` contains scripts for the supplementary figures used in the manuscript. The current repository includes scripts corresponding to Supplementary Figures 2–42, with figure numbers following the manuscript numbering.

Major categories include:

- WGS and RNA-seq quality control;
- ancestry inference;
- long-/short-read variant comparison;
- VEP annotation;
- transcript discovery and filtering;
- peptide support and predicted protein structures;
- cross-platform expression consistency;
- transcript sampling-depth evaluation;
- co-expression modules;
- genotype phasing;
- ASE, ASTS and allele-specific splicing;
- QTL covariate selection;
- sQTL consistency;
- SV/TR QTL representation and tagging;
- genomic enrichment;
- eQTL replication and portability;
- joint fine-mapping;
- robustness analyses;
- GWAS enrichment and colocalization;
- frequency-differentiated QTLs.

## Inputs

Figure-specific input files are stored under each figure directory, for example:

```text
Figure_3/input/
Figure_4/input/
Extended_Figures/input/
Supp_Figures/input/
```

Inputs include plain-text tables, compressed tables, and serialized R objects (`.RData`, `.RDS`). These files generally represent processed results generated by the upstream analysis pipeline rather than raw sequencing data.

## Running a figure script

Most figure scripts can be run directly after adapting the working directory and input paths:

```r
setwd("/path/to/GTOP/figure_generation/Figure_3")
source("Figure3.R")
```

or from the shell:

```bash
Rscript Figure3.R
```

The exact command and input requirements may differ between figures; consult the figure-specific README and the corresponding R script.

## Shared dependencies

The figure-generation scripts use a range of R packages for statistical graphics, data manipulation, genomic visualization, and multi-panel composition. Examples include:

- `ggplot2`
- `data.table`
- `dplyr`
- `tidyverse`
- `ggpubr`
- `patchwork`
- `cowplot`
- `ComplexHeatmap`
- `locuszoomr`
- `EnsDb.Hsapiens.v86`
- `ggvenn`
- `ggradar`
- `ggrastr`

Additional packages are loaded by individual figure scripts as needed.

## Reproducibility

The figure scripts are designed to regenerate the manuscript figures from the corresponding processed input files. They are not intended to reproduce upstream statistical analyses; those workflows are documented under `analysis_pipeline/`.

For reproducibility:

1. Prepare the required upstream analysis outputs.
2. Place the corresponding processed inputs in each figure's `input/` directory.
3. Update environment-specific paths in the R scripts.
4. Install the packages required by the target figure.
5. Run the corresponding figure script.

## Citation

If you reuse the GTOP figure-generation code or derived results, please cite the GTOP study and the original software packages used in the analyses.
