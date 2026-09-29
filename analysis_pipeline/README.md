# GTOP Analysis Pipeline

This directory contains the code for GTOP data processing and downstream analyses, from long-read transcript discovery and RNA phenotype preparation to variant calling, allele-specific analyses, QTL mapping, fine-mapping, cross-ancestry analyses, enrichment, colocalization, and SMR.

Detailed instructions for running individual workflows are provided in the README files within the corresponding subdirectories.


## Pipeline stages

###Transcript discovery 

`00_transcript_detection/01_transcript_discovery/` | Seven complementary long-read transcript discovery approaches 

**Input**
PacBio full-length non-chimeric (FLNC) reads in FASTQ/BAM format and their alignments to hg38
the hg38 reference genome FASTA; GENCODE v47 annotation GTF 

**Output**
Per-sample transcript annotations in GTF format generated independently by Bambu, FLAIR, FLAMES, IsoQuant, Iso-Seq, IsoTools, and TALON, together with caller-specific read-support information 

| Transcript integration | `00_transcript_detection/02_merge_filter_transcript/` | FLNC support filtering, TAMA merging, SQANTI3 annotation, and enhanced GTF construction | Transcript annotations and read-support information from the seven discovery tools; PacBio FLNC reads and genome alignments; short-read splice-junction evidence from STAR; hg38 genome and GENCODE v47 annotation; CAGE/TSS annotations and poly(A) motif list | Merged, quality-filtered GTOP transcript catalogue and an enhanced reference combining GTOP novel transcripts with GENCODE v47; transcript GTF/FASTA, annotation and QC tables, and predicted protein sequences |
| Transcript quantification | `00_transcript_detection/03_quantification/` | FLAIR transcript quantification | Per-sample PacBio FLNC FASTQ files and the final GTOP transcript reference, including transcript sequences and annotations | Transcript-level read-count and TPM matrices across samples, plus gene-level count and TPM matrices aggregated from transcripts belonging to the same gene |
| Peptide validation | `00_transcript_detection/04_peptide_validation/` | DIA-NN-based proteomic support | Tissue-matched DIA-MS data in mzML format; GENCODE v47 and predicted GTOP protein sequences; transcript annotations; long-read transcript TPM estimates; tissue metadata | Tissue-specific protein databases; DIA-NN peptide-evidence reports; per-tissue protein abundance matrices |
| RNA phenotype preparation | `01_data_preparation/` | Gene-, splice-junction-, and transcript-level phenotypes | RNA-seq FASTQ files; donor-specific DNA VCFs; genome/reference annotations; enhanced GTOP transcript reference; tissue and sample metadata; exon/gene annotations | Gene-level counts/TPM; WASP-filtered alignments and splicing phenotypes; transcript-level count/TPM matrices and filtered, imputed, normalized transcript-usage phenotypes |
| Variant calling | `02_Variant_calling/` | LRS/SRS small variants, population structure, and VEP annotation | LRS HiFi WGS BAM files and SRS WGS FASTQ files | Per-sample genotype VCFs from LRS and SRS WGS, plus a cohort-level merged and site-filtered VCF |
| ASE/ASTS/ASJ | `03_ASE_ASTS/` | Short-read ASE, long-read ASE/ASTS, allele-specific splicing, and ASJ | WASP-filtered RNA BAMs, donor-specific heterozygous VCFs, reference/annotation files; long-read RNA/genomic data and phased genotypes; inputs for lorals, isoLASER, and longcallR | Allelic counts and ASE estimates; lorals ASE/ASTS results; isoLASER summaries and merged gVCF; longcallR phased RNA VCF/BAM and DNA-supported ASE/ASJ results |
| QTL mapping | `04_QTL_mapping/` | eQTL, sQTL, TR-xQTL, and MAJIQTL | Phenotypes in BED/Parquet format; covariates; genotype data preferably in PLINK2 PGEN/PVAR/PSAM format | Nominal cis-QTL summary statistics for variant–phenotype pairs and permutation-mode phenotype-level statistics with empirical P-values for genome-wide FDR |
| Fine-mapping | `05_finemap/` | SuSiE and cross-ancestry SuSiEx fine-mapping | Individual-level normalized phenotypes for each eGene locus and genotype matrices for variants within 1 Mb of the target eGene TSS | Tissue-level fine-mapping tables containing locus ID, variant ID, PIP, credible-set ID, credible-set size, credible-set purity, and tissue |
| Frequency differences | `06_frequency_differences/` | Frequency-differentiated QTLs | gnomAD allele frequencies across ancestry groups; SuSiE fine-mapping results containing gene ID, variant ID, PIP, and credible-set ID | fd-QTL list, including independent loci and frequency-difference classifications across ancestry groups |
| eQTL portability | `07_eQTL_portability/` | GTOP–GTEx portability and mashR analyses | GTOP and GTEx eQTL summary statistics | Six portability metrics and estimated proportions of portable eQTLs |
| Enrichment | `08_QTL_enrichment/` | torus and S-LDSC enrichment analyses | torus: tensorQTL nominal QTL statistics and variant annotations; S-LDSC: significant QTL pairs, fine-mapping tables, reference genotypes, baselineLD resources, weights, and GWAS summary statistics | torus annotation-enrichment estimates with confidence intervals; S-LDSC heritability contributions, enrichment statistics, and annotation coefficients |
| Colocalization | `09_colocalization/` | GWAS/QTL fine-mapping and colocalization | GWAS summary statistics and fine-mapping results; QTL summary statistics and fine-mapping results | Colocalization results from SuSiE-coloc, supplemented by coloc.abf analyses |
| SMR | `10_SMR/` | Summary-data-based Mendelian randomization | GWAS and QTL summary statistics | SMR results including gene ID, SMR P-value, and HEIDI P-value |

## Directory structure

```text
analysis_pipeline/
├── 00_transcript_detection/
│   ├── 01_transcript_discovery/
│   ├── 02_merge_filter_transcript/
│   ├── 03_quantification/
│   └── 04_peptide_validation/
├── 01_data_preparation/
├── 02_Variant_calling/
├── 03_ASE_ASTS/
├── 04_QTL_mapping/
├── 05_finemap/
├── 06_frequency_differences/
├── 07_eQTL_portability/
├── 08_QTL_enrichment/
├── 09_colocalization/
└── 10_SMR/
```

## Main analysis components

### Transcript discovery

GTOP long-read transcript discovery uses seven complementary approaches: **Bambu, FLAIR, FLAMES, IsoQuant, Iso-Seq, IsoTools, and TALON**. Independent transcript annotations are subsequently integrated, filtered, and annotated to generate the GTOP transcript reference.

### Data preparation and variant calling

Short-read RNA-seq data are processed to generate gene-level, splice-junction, and transcript-level phenotypes. LRS and SRS whole-genome sequencing data are processed for small-variant discovery, population structure analysis, and functional annotation.

### Allele-specific analyses

The ASE/ASTS module contains workflows based on short-read ASE, lorals, isoLASER, and **longcallR**, supporting allele-specific expression, allele-specific transcript usage/splicing, and allele-specific junction analyses.

### QTL mapping

The QTL module supports eQTL, sQTL, TR-xQTL, and **MAJIQTL** analyses. Phenotypes, covariates, and genotype data are prepared in formats compatible with the corresponding QTL workflows.

### Fine-mapping and downstream analyses

Fine-mapping is performed using SuSiE and cross-ancestry SuSiEx. Downstream modules characterize frequency-differentiated QTLs, eQTL portability, QTL enrichment, GWAS/QTL colocalization, and SMR associations.

## Reproducibility

The pipelines are organized as modular workflows so that individual stages can be rerun independently when their required inputs are available. Please consult the README files in each subdirectory for software requirements, configuration, input preparation, and execution commands.
