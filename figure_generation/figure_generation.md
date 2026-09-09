# GTOP Figure Generation


This directory contains the R scripts and supporting input files used to generate the main, extended, and supplementary figures for the GTOP Phase-I manuscript.


The figure-generation workflow is organized according to the manuscript figures. Each figure directory generally contains:


* An R script for figure generation

* An `input/` directory containing required intermediate data

* A figure-specific README describing the panels and analyses


The scripts are intended to reproduce the major figures and panels presented in the GTOP manuscript.


---


## Directory Structure


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


---


# Main Figures


## Figure 1 — GTOP cohort, genomic variation, and long-read sequencing


```text

Figure_1/

├── input/

├── Figure1.R

└── README-FIG.1.md

```


Figure 1 summarizes the GTOP cohort and genomic variation identified using long- and short-read sequencing.


The panels include:


* **Figure 1B:** PCA of GTOP, GTEx, and 1000 Genomes genotypes to characterize continental ancestry structure.

* **Figure 1C:** MDS of GTOP RNA-seq samples based on gene-expression profiles.

* **Figure 1D:** Hierarchical clustering of tissues based on gene expression.

* **Figure 1E:** Genome-wide cumulative variant length and variant counts.

* **Figure 1F:** Numbers of small variants, structural variants, and tandem repeats detected from long-read and matched short-read sequencing.

* **Figure 1G–I:** Comparison of small variants, structural variants, and tandem repeats identified from long-read and short-read sequencing.


The corresponding implementation is available in `Figure_1/Figure1.R`.


---


## Figure 2 — Transcriptome diversity revealed by long-read RNA sequencing


```text

Figure_2/

├── input/

├── Figure2.R

└── README-FIG.2.md

```


Figure 2 characterizes transcript discovery and transcript-level regulatory diversity enabled by GTOP long-read RNA sequencing.


The panels include:


* **Figure 2A:** MDS of full-length RNA-seq samples based on transcript-level expression.

* **Figure 2B:** Comparison of GTOP transcript discovery with GENCODE v47.

* **Figure 2C:** Transcript classification and coding-potential assessment using SQANTI3.

* **Figure 2D:** Alternative splicing patterns of GTOP novel transcripts relative to GENCODE annotation.

* **Figure 2E:** Tissue distribution of expressed novel and annotated transcripts.

* **Figure 2F:** Peptide support for novel transcripts across tissues.

* **Figure 2G:** WGCNA transcript co-expression modules across cardiac tissues.

* **Figure 2H:** Numbers of genes analyzed for ASE and ASTS and the corresponding significant events.

* **Figure 2I:** Distribution of ASE- and ASTS-associated genes across the 33 tissues.


The corresponding implementation is available in `Figure_2/Figure2.R`.


---


## Figure 3 — Molecular QTL landscape across variant classes


```text

Figure_3/

├── input/

├── Figure3.R

└── README-FIG.3.md

```


Figure 3 describes the landscape of molecular QTLs identified across GTOP tissues and variant classes.


The panels include:


* **Figure 3A:** Numbers of cis-eQTLs and cis-sQTLs by tissue and variant class.

* **Figure 3B:** Genomic distance between eQTL variants and the transcription start sites of their target genes.

* **Figure 3C:** Correlation of eQTL effect sizes between GTOP and matched GTEx tissues.

* **Figure 3D:** Overlap of eGenes identified by SNV-, SV-, and TR-eQTLs.

* **Figure 3E:** Disease-associated tandem repeats associated with molecular phenotypes across GTOP tissues.

* **Figure 3F:** Example of a pathogenic tandem repeat associated with MUC1 expression.

* **Figure 3G:** Tissue sharing of eQTL signals.


The corresponding implementation is available in `Figure_3/Figure3.R`.


---


## Figure 4 — Population-specific QTLs and cross-ancestry fine-mapping


```text

Figure_4/

├── input/

├── Figure4.R

└── README-FIG.4.md

```


Figure 4 focuses on ancestry-dependent QTL architecture and the benefits of cross-ancestry fine-mapping.


The panels include:


* **Figure 4A:** Geographic allele-frequency distributions of GTOP independent eQTLs.

* **Figure 4B:** Comparison of allele frequencies between GTOP and European populations.

* **Figure 4C:** Classification of fd-QTLs according to allele-frequency differences between East Asian and European populations.

* **Figure 4D:** Number of genes regulated by fd-QTLs across tissues.

* **Figure 4E:** Example locus showing colocalization between a disease-associated signal and molecular QTLs.

* **Figure 4F:** Number of heterogeneity QTLs across GTOP tissues.

* **Figure 4G:** Comparison of he-QTL effect-size statistics between GTOP and GTEx.

* **Figure 4H:** Size of 95% credible sets and comparison between cross-ancestry and single-population fine-mapping.

* **Figure 4I:** Proportion of credible sets containing variants with high maximum PIP.

* **Figure 4J:** PIP distribution of variants in fine-mapped credible sets for the ABO eGene in liver.


The corresponding implementation is available in `Figure_4/Figure4.R`.


---


## Figure 5 — Fine-mapping of SVs and tandem repeats


```text

Figure_5/

├── input/

├── Figure5.R

└── README-FIG.5.md

```


Figure 5 investigates the contribution of structural variants and tandem repeats to molecular QTL signals and fine-mapped credible sets.


The panels include:


* **Figure 5A:** LD/tagging correlation between SV/TR eQTLs and their best-tagging SNVs.

* **Figure 5B:** Number of fine-mapped variants and credible sets across tissues.

* **Figure 5C:** Number of credible sets and the proportion led by SVs or TRs.

* **Figure 5D:** Composition of fine-mapped eQTL credible sets containing SVs, TRs, or only SNVs.

* **Figure 5E:** Enrichment of SVs and TRs among variants prioritized at different PIP thresholds.


The corresponding implementation is available in `Figure_5/Figure5.R`.


---


## Figure 6 — Linking molecular QTLs to complex traits and diseases


```text

Figure_6/

├── input/

├── Figure6.R

└── README-FIG.6.md

```


Figure 6 integrates GTOP molecular QTLs with GWAS and disease-associated genetic signals.


The panels include:


* **Figure 6A:** GWAS enrichment of cis-eQTLs and cis-sQTLs.

* **Figure 6B:** GWAS loci colocalized with eQTLs, sQTLs, or both, together with the overlap between colocalization and SMR-identified genes.

* **Figure 6C:** Comparison of colocalized genes identified by GTOP, GTEx, JCTF, and MAGE.

* **Figure 6D:** Example of GTOP-specific colocalization between ADAP1 eQTLs and Crohn's disease.

* **Figure 6E:** Proportion of GWAS variants overlapping molecular QTLs.

* **Figure 6F:** GWAS/QTL overlap stratified by variant frequency across ancestry groups.

* **Figure 6G:** Fine-mapped SV/TR-containing colocalized signals.

* **Figure 6H:** Example hyperlipidemia locus showing the relationship between GWAS, SNV-eQTL, and SV-eQTL signals.


The corresponding implementation is available in `Figure_6/Figure6.R`.


---


# Extended Figures


```text

Extended_Figures/

├── input/

├── Extended_Figure1.R

├── Extended_Figure2.R

├── ...

├── Extended_Figure9.R

└── README.md

```


The Extended Figures provide additional analyses supporting the major conclusions of the GTOP manuscript.


### Extended Figure 1


Comparison of genetic variants detected using long-read and short-read sequencing.


### Extended Figure 2


Transcript discovery, transcript quantification, and peptide validation.


### Extended Figure 3


Long-read identification of allele-specific expression (ASE) and allele-specific transcript structure (ASTS).


### Extended Figure 4


Characteristics of molecular QTLs across different variant classes.


### Extended Figure 5


Tissue sharing and context-dependent patterns of molecular QTLs.


### Extended Figure 6


Characterization of fd-QTLs and he-QTLs.


### Extended Figure 7


Cross-ancestry fine-mapping of eQTLs using GTOP and GTEx.


### Extended Figure 8


Evaluation of SV-eQTLs and TR-eQTLs that are poorly tagged by small-variant eQTLs.


### Extended Figure 9


Characterization of genes associated with disease-risk variants.


The corresponding R scripts are provided as `Extended_Figure1.R` through `Extended_Figure9.R`.


---


# Supplementary Figures


```text

Supp_Figures/

├── input/

├── Supp_Figure_2.R

├── Supp_Figure_3.R

├── ...

├── Supp_Figure_29.R

└── README.md

```


The supplementary figure scripts provide additional quality-control, validation, methodological, and supporting analyses.


These include:


* **Supp. Figure 2:** Whole-genome sequencing statistics

* **Supp. Figure 3:** Short-read RNA-seq quality control

* **Supp. Figure 4:** Long-read RNA-seq quality assessment

* **Supp. Figure 5:** ADMIXTURE-based ancestry inference

* **Supp. Figure 6:** Variant detection statistics and mutation patterns

* **Supp. Figure 7:** VEP-based genomic annotation of variants

* **Supp. Figure 9:** Transcript discovery, annotation, and filtering

* **Supp. Figure 10:** DIA-MS-based peptide validation

* **Supp. Figure 11:** Consistency of gene/transcript quantification across sequencing platforms and analytical methods

* **Supp. Figure 12:** Transcript-level co-expression modules

* **Supp. Figure 13:** Genotype phasing statistics

* **Supp. Figure 14:** ASE and ASTS analyses

* **Supp. Figure 15:** Allele-specific splicing

* **Supp. Figures 16–18:** PCA-based selection of molecular phenotype covariates

* **Supp. Figure 19:** Selection of genotype PCs for QTL mapping

* **Supp. Figure 20:** Genomic-feature enrichment of sQTLs

* **Supp. Figure 21:** Concordance of eQTL effects across pancreas sampling sites

* **Supp. Figure 22:** Concordance between eQTL and ASE effects


Additional supplementary figure scripts are included in the same directory.


---


# Running the Figure Generation Scripts


The figure-generation scripts are designed to operate on processed intermediate results generated by the GTOP analysis pipelines.


A typical workflow is:


```text

GTOP raw / processed data

          │

          ▼

analysis_pipeline/

          │

          ├── Transcript processing

          ├── Variant calling

          ├── ASE / ASTS

          ├── QTL mapping

          ├── Fine-mapping

          ├── Enrichment

          └── Colocalization

          │

          ▼

Intermediate analysis results

          │

          ▼

figure_generation/

          │

          ├── Figure_1/

          ├── Figure_2/

          ├── ...

          ├── Figure_6/

          ├── Extended_Figures/

          └── Supp_Figures/

          │

          ▼

Manuscript figures

```


Each figure directory contains its own `input/` directory and figure-specific R script. Users should first ensure that the required intermediate data are available in the corresponding `input/` directory before executing the R script.


---


# Reproducibility


The figure-generation scripts reproduce the major visualizations reported in the GTOP manuscript. The scripts rely on intermediate results generated by the corresponding upstream analysis modules.


For reproducible execution:


1. Prepare the required GTOP processed data and analysis results.

2. Place or link the required inputs in the corresponding `input/` directory.

3. Install the R packages required by the individual script.

4. Run the corresponding figure-generation script.

5. Inspect the generated figures against the manuscript panels.


Because the figure scripts depend on intermediate datasets and external annotations, users should refer to the corresponding figure-specific README and script for the exact input files and analysis settings.


---


## Citation


If you use these figure-generation scripts or reproduce GTOP analyses, please cite the GTOP manuscript:


> Zou X, Li X, Zhang T, et al. *A long-read multi-omics atlas broadens discovery of regulatory variation across 33 human tissues.*


Please also cite the original software packages and external datasets used in the corresponding analyses.


---


## License


This project is released under the MIT License.


