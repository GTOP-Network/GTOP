# GTOP Analysis Pipeline


This directory contains the computational pipelines used for processing, analyzing, and interpreting the GTOP Phase-I multi-omics resource. The workflows cover long-read transcript discovery and quantification, genome-wide variant identification, allele-specific analyses, molecular QTL mapping, fine-mapping, ancestry-specific QTL analyses, QTL enrichment, GWAS-QTL colocalization, and SMR analyses.


The pipeline is organized into sequential modules, from primary molecular data processing to downstream functional genomics analyses.


---


## Overview


The GTOP analysis framework consists of the following major stages:


```text

Long-read RNA-seq

      │

      ├── Transcript discovery

      │      ├── Bambu

      │      ├── FLAIR

      │      └── Iso-Seq

      │

      ├── Transcript integration and filtering

      │      ├── Transcript merging

      │      ├── FLAIR quantification

      │      ├── SQANTI3 annotation

      │      └── High-confidence transcript reference

      │

      └── Transcript quantification / peptide validation

      │

      ▼

Genome-wide sequencing

      │

      ├── Small variant calling

      ├── Structural variant calling

      ├── Tandem repeat genotyping

      ├── Population structure

      └── Variant annotation

      │

      ▼

Allele-specific analyses

      │

      ├── Short-read ASE

      ├── Long-read ASE

      ├── ASTS

      └── Allele-specific splicing

      │

      ▼

Molecular QTL mapping

      │

      ├── eQTL

      ├── sQTL

      └── TR-eQTL / TR-sQTL

      │

      ▼

Fine-mapping

      │

      ├── SuSiE

      └── Cross-ancestry fine-mapping

      │

      ▼

Downstream functional analyses

      │

      ├── Population-specific QTLs

      ├── QTL enrichment

      ├── GWAS-QTL colocalization

      └── SMR

```


---


## Directory Structure


```text

analysis_pipeline/

├── 00_transcript_detection/

├── 01_data_preparation/

├── 02_Variant_calling/

├── 03_ASE_ASTS/

├── 04_QTL_mapping/

├── 05_finemap/

├── 06_population-specific_eQTL/

├── 07_QTL_enrichment/

├── 08_colocalization/

└── 09_SMR/

```


Each module contains the scripts required for the corresponding stage of the GTOP analysis. Where applicable, module-specific README files provide additional details on individual workflows and execution commands.


---


# 00. Transcript Detection


`00_transcript_detection/` contains the long-read RNA-seq transcript discovery and annotation workflow.


The pipeline integrates three complementary transcript discovery approaches:


* **Bambu** — reference-guided transcript discovery

* **FLAIR** — reference-guided isoform reconstruction with splice-site correction

* **Iso-Seq** — reference-free transcript reconstruction


These complementary approaches are used to maximize transcript discovery while improving structural confidence and reproducibility.


### 00.1 Transcript discovery


```text

00_transcript_detection/

└── 01_transcript_discovery/

    ├── bambu/

    ├── flair/

    ├── isoseq/

    └── README.md

```


The transcript discovery workflow includes:


1. Preparation of genome annotation objects

2. Read-class generation

3. Reference-guided transcript discovery with Bambu

4. Read alignment and splice-site correction with FLAIR

5. Isoform collapsing with FLAIR

6. Reference-free transcript reconstruction using Iso-Seq

7. Chromosome-level and genome-wide transcript annotation merging


The three approaches generate complementary transcript models that are subsequently integrated.


### 00.2 Transcript merging and filtering


```text

00_transcript_detection/

└── 02_merge_filter_transcript/

```


Transcript models identified by Bambu, FLAIR, and Iso-Seq are integrated according to intron-chain identity.


The workflow includes:


1. Integration of transcript models from the three discovery approaches

2. Prioritization of transcripts supported by multiple methods

3. Removal of low-confidence and mono-exonic transcript models

4. Transcript-level quantification using FLAIR

5. Transcript annotation using SQANTI3

6. Structural and expression-based filtering

7. Construction of GTOP and enhanced transcript references


The resulting high-confidence transcript annotation is used for downstream expression, transcript usage, ASE/ASTS, and QTL analyses.


### 00.3 Transcript quantification


```text

00_transcript_detection/

└── 03_quantification/

```


This module contains scripts for transcript-level expression quantification based on the GTOP transcript reference.


### 00.4 Peptide validation


```text

00_transcript_detection/

└── 04_peptide_validation/

```


This module contains scripts used to evaluate transcript-level peptide support using proteomics data.


---


# 01. Data Preparation


`01_data_preparation/` contains preprocessing workflows for molecular phenotypes used in downstream analyses.


```text

01_data_preparation/

├── SRS_RNA_quantification/

├── splicing_quantification/

└── transcript_quantification/

```


The module prepares:


* Short-read RNA-seq gene expression

* Splicing phenotypes

* Long-read transcript expression

* Molecular phenotype matrices for downstream QTL analyses


The resulting phenotype matrices are used as inputs for the QTL mapping workflow.


---


# 02. Variant Calling


`02_Variant_calling/` contains workflows for identifying and annotating genetic variants from long- and short-read sequencing data.


```text

02_Variant_calling/

├── 01_LRS_WGS_Variant_Calling.sh

├── 02_SRS_WGS_Variant_Calling.sh

├── 03_PCA_ADMIXTURE.sh

└── 04_vep_annotation.sh

```


The workflow includes:


### Long-read WGS variant calling


Variant discovery from long-read WGS data.


### Short-read WGS variant calling


Variant discovery from matched short-read WGS data.


### Population structure


PCA and ADMIXTURE analyses are performed to characterize sample ancestry and population structure.


### Variant annotation


Variants are functionally annotated using VEP.


These variant sets form the basis for allele-specific analyses, QTL mapping, fine-mapping, and downstream functional interpretation.


---


# 03. ASE and ASTS


`03_ASE_ASTS/` contains workflows for allele-specific expression (ASE), allele-specific transcript structure (ASTS), and allele-specific splicing.


```text

03_ASE_ASTS/

├── isolaser/

├── scripts/

├── run_LRS_ase_lorals.sh

├── run_SRS_ase.sh

└── README.md

```


Three complementary analyses are implemented:


### Short-read ASE


Variant-aware read mapping is performed using STAR together with WASP filtering to reduce reference-mapping bias.


### Long-read ASE and ASTS


Long-read allele-specific expression and transcript-structure analyses are performed using `lorals`.


The long-read workflow includes:


1. Processing donor-specific VCFs

2. Construction of haplotype-specific references

3. Haplotype-aware read mapping

4. Transcriptome alignment

5. Allelic coverage calculation

6. Transcript-level allele-specific quantification

7. Filtering and statistical testing


### Allele-specific splicing


`isoLASER` is used to identify allele-specific transcript splicing events from long-read RNA-seq data.


---


# 04. QTL Mapping


`04_QTL_mapping/` contains the molecular QTL mapping workflows used in GTOP.


```text

04_QTL_mapping/

├── 01_data_preparation/

├── run_eQTL_mapping.sh

├── run_sQTL_mapping.sh

├── run_TR_xQTL_mapping.sh

├── run_tensorqtl_TR_eQTL.py

└── run_tensorqtl_TR_sQTL.py

```


The main QTL analyses include:


* **cis-eQTL mapping**

* **cis-sQTL mapping**

* **TR-eQTL mapping**

* **TR-sQTL mapping**


The workflow includes phenotype and genotype preparation followed by QTL association testing. TensorQTL-based scripts are provided for tandem-repeat QTL analyses.


---


# 05. Fine-mapping


`05_finemap/` contains fine-mapping workflows for identifying variants with high posterior probability of driving molecular QTL associations.


```text

05_finemap/

├── SuSiE_Fine_mapping/

└── Cross_ancestries_Fine_mapping_SuSiEX/

```


The analyses include:


### SuSiE fine-mapping


SuSiE is used to identify credible sets and estimate posterior inclusion probabilities (PIP) for candidate variants.


### Cross-ancestry fine-mapping


Cross-ancestry fine-mapping is performed using GTOP and external ancestry-specific datasets to improve resolution of causal QTL signals.


---


# 06. Population-specific eQTL


`06_population-specific_eQTL/` contains analyses investigating ancestry-dependent properties of QTLs.


```text

06_population-specific_eQTL/

├── MashR/

├── a1.different_ancestries_freq.R

├── a1.fm_SNV_QTL_add_freq.R

├── a2.QTL_replication_in_GTEx.R

├── a3.compare_with_GTEx.R

├── a4.calculate_he_QTLs.R

├── a4.data_preparation.R

├── a4.finemapping_data.R

├── a4.mash_finemapping_QTLs.R

└── he-QTL_mash.sh

```


The workflow covers:


* Allele-frequency comparison across ancestry groups

* Identification and characterization of population-specific QTLs

* Fine-mapping with ancestry-aware allele-frequency information

* Replication in GTEx

* Comparison of GTOP and GTEx QTL effect sizes

* Estimation and fine-mapping of heterogeneity QTLs (he-QTLs)

* Multivariate effect-size analysis using mashR


The scripts are organized into analytical stages (`a1`–`a4`) corresponding to variant-frequency analysis, QTL characterization, GTEx comparison, and heterogeneity/fine-mapping analyses.


---


# 07. QTL Enrichment


`07_QTL_enrichment/` evaluates whether GTOP QTLs are enriched in genomic annotations and disease-relevant heritability.


```text

07_QTL_enrichment/

├── README.md

├── run_torus.sh

└── run_sldsc.sh

```


Two complementary approaches are implemented:


### Genomic annotation enrichment


Bayesian enrichment analysis is performed using **torus**.


The workflow prepares tensorQTL nominal association results and estimates enrichment parameters for different genomic annotation categories.


### Heritability enrichment


Stratified LD score regression (**S-LDSC**) is used to quantify enrichment of trait heritability in QTL-derived annotations.


Variant annotations can include:


* Significant QTL variants

* Fine-mapped credible-set variants

* Credible-set variants weighted by maximum PIP


These analyses connect molecular QTLs with genomic regulatory annotations and complex-trait genetic architecture.


---


# 08. GWAS-QTL Colocalization


`08_colocalization/` contains integrative analyses combining GWAS and molecular QTL signals.


```text

08_colocalization/

├── scripts/

├── colocalization.sh

└── README.md

```


The workflow consists of:


1. Identification and annotation of GWAS loci

2. GWAS fine-mapping

3. Preparation of gene-based regulatory annotations

4. QTL fine-mapping

5. GWAS-QTL colocalization

6. Integration of fine-mapping and colocalization results


The final integrated results are used to identify genes and regulatory variants that may mediate disease-associated genetic effects.


---


# 09. SMR


`09_SMR/` contains scripts for Summary-data-based Mendelian Randomization (SMR) analyses.


```text

09_SMR/

└── run_SMR.sh

```


SMR analyses integrate molecular QTL signals with GWAS summary statistics to prioritize genes whose molecular phenotypes may contribute to complex traits and diseases.


---


# Reproducibility


The scripts in this directory were developed for high-performance computing environments and include shell, Python, and R workflows.


Because individual modules depend on external reference genomes, annotations, genotype files, phenotype matrices, and software packages, users should consult the README and scripts within each module before execution.


In general, the recommended execution strategy is:


1. Prepare reference genomes and annotations.

2. Process raw sequencing data and generate molecular phenotypes.

3. Generate and annotate genetic variants.

4. Perform ASE/ASTS analyses.

5. Prepare genotype and phenotype matrices.

6. Perform QTL mapping.

7. Perform fine-mapping.

8. Conduct population-specific and functional enrichment analyses.

9. Integrate QTLs with GWAS using colocalization and SMR.


The workflows are modular and can be executed independently when the corresponding upstream inputs are available.


---


## Citation


If you use the GTOP analysis pipelines or derived results in your work, please cite the GTOP manuscript:


> Zou X, Li X, Zhang T, et al. *A long-read multi-omics atlas broadens discovery of regulatory variation across 33 human tissues.*


Please also cite the original software packages used in individual analyses.


---


## License


This project is released under the MIT License.


