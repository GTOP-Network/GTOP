# GTOP Analysis Pipeline: Detailed Documentation

This document provides a detailed description of the GTOP data-processing and downstream-analysis workflows, including the principal inputs, analytical procedures, and expected outputs for each stage.

## 1. Transcript detection

### Transcript discovery

**Input:** PacBio full-length non-chimeric (FLNC) reads in FASTQ or BAM format, their alignments to hg38, the hg38 reference genome FASTA, and the GENCODE v47 annotation GTF.

**Workflow:** This workflow, located in `00_transcript_detection/01_transcript_discovery/`, applies seven complementary long-read transcript discovery approaches: Bambu, FLAIR, FLAMES, IsoQuant, Iso-Seq, IsoTools, and TALON.

**Expected output:** Per-sample transcript annotations in GTF format generated independently by the seven tools, together with caller-specific read-support information for subsequent transcript filtering and integration.

### Transcript integration

**Input:** Transcript annotations and read-support information from the seven discovery tools, PacBio FLNC reads and genome alignments, short-read splice-junction evidence from STAR, the hg38 genome, the GENCODE v47 annotation, CAGE/TSS annotations, and a poly(A) motif list.

**Workflow:** This workflow, located in `00_transcript_detection/02_merge_filter_transcript/`, performs FLNC support filtering, TAMA-based transcript merging, SQANTI3 annotation, and construction of the enhanced GTOP transcript reference.

**Expected output:** A merged and quality-filtered GTOP transcript catalogue and an enhanced reference combining GTOP novel transcripts with GENCODE v47 annotations, together with transcript GTF and FASTA files, annotation and quality-control tables, and predicted protein sequences.

### Transcript quantification

**Input:** Per-sample PacBio FLNC FASTQ files and the final GTOP transcript reference, including transcript sequences and associated annotations.

**Workflow:** This workflow, located in `00_transcript_detection/03_quantification/`, uses FLAIR to quantify GTOP transcripts.

**Expected output:** Transcript-level read counts quantified using FLAIR and combined across samples into count and TPM matrices, together with gene-level count and TPM matrices obtained by aggregating transcripts belonging to the same gene.

### Peptide validation

**Input:** Tissue-matched DIA mass-spectrometry data in mzML format, GENCODE v47 protein sequences, predicted GTOP protein sequences, transcript annotations, long-read transcript TPM estimates, and tissue metadata.

**Workflow:** This workflow, located in `00_transcript_detection/04_peptide_validation/`, uses DIA-NN-based proteomic evidence to evaluate protein support for GTOP transcripts.

**Expected output:** Tissue-specific protein databases containing coding transcripts with TPM > 5 in at least one corresponding long-read sample, DIA-NN reports providing peptide evidence for the searched proteins, and per-tissue protein abundance matrices.

## 2. RNA phenotype preparation

### Gene-level quantification

**Input:** Raw paired-end RNA-seq FASTQ data for each sample, together with the associated sample and reference information required by the workflow.

**Workflow:** This workflow, implemented in `01_data_preparation/`, generates gene-expression phenotypes from short-read RNA-seq data.

**Expected output:** Gene-level read-count and transcripts-per-million (TPM) matrices.

### Splicing quantification

**Input:** Paired-end RNA-seq FASTQ files, donor-specific DNA VCFs, a STAR genome index, the reference genome FASTA, a gene annotation GTF, tissue metadata, sample-to-participant lookup tables, exon annotations containing `chr`, `start`, `end`, `strand`, and `gene_id`, and gene metadata containing `gene_name`, `gene_type`, `gene_id`, `chr`, `start`, `end`, and `strand`.

**Workflow:** This workflow, implemented in `01_data_preparation/`, generates splice-junction and intron-cluster phenotypes for downstream QTL analyses.

**Expected output:** WASP-filtered alignments, junction counts, intron-cluster ratios, and filtered and normalized splicing phenotypes.

### Transcript quantification

**Input:** Paired-end short-read RNA-seq FASTQ files, the hg38 reference genome FASTA, the enhanced GTOP novel–GENCODE v47 transcript reference in GTF and FASTA format with associated gene–transcript annotations, and tissue metadata.

**Workflow:** This workflow, implemented in `01_data_preparation/`, generates transcript-level expression and transcript-usage phenotypes from short-read RNA-seq data using Salmon and RSEM.

**Expected output:** Transcript-level count and TPM matrices, tissue-specific transcript TPM matrices obtained by integrating the two quantification methods after concordance filtering, and filtered, imputed, and normalized transcript-usage phenotypes for downstream QTL analysis.

## 3. Variant calling

### LRS/SRS variant calling

**Input:** Long-read PacBio HiFi WGS BAM files and short-read WGS FASTQ files for each sample.

**Workflow:** This workflow, located in `02_Variant_calling/`, processes long-read and short-read whole-genome sequencing data to generate cohort-level variant calls.

**Expected output:** Per-sample genotype VCF files generated from LRS and SRS WGS data, together with a cohort-level merged and site-filtered VCF containing variants and genotypes across all samples.

### Population structure analysis

**Input:** Cohort-level genotype data derived from whole-genome sequencing.

**Workflow:** This workflow, implemented within `02_Variant_calling/`, performs principal component analysis (PCA) and admixture analysis to characterize population structure.

**Expected output:** PCA coordinates, admixture proportions, and related ancestry-inference results for downstream analyses.

### Variant annotation

**Input:** The cohort-level variant set generated by the variant-calling workflow.

**Workflow:** This workflow, implemented within `02_Variant_calling/`, applies the Variant Effect Predictor (VEP) to annotate genetic variants.

**Expected output:** A VEP-annotated variant set for downstream variant interpretation and functional analyses.

The principal entry points for this module are `01_LRS_WGS_Variant_Calling.sh`, `02_SRS_WGS_Variant_Calling.sh`, `03_PCA_ADMIXTURE.sh`, and `04_vep_annotation.sh`.

## 4. ASE, ASTS, and ASJ

### Short-read ASE

**Input:** WASP-filtered RNA BAM files, donor-specific heterozygous DNA VCFs, the reference FASTA, gene BED annotations containing chromosome, start, end, and gene ID, and an ASE manifest containing `sample_id`, `tissue_site_detail`, and `ase_readcount_file`.

**Workflow:** This workflow is implemented in `03_ASE_ASTS/run_SRS_ase.sh` and identifies allele-specific expression events from short-read RNA-seq data.

**Expected output:** Per-site allelic counts, donor-level LAMP estimates, and `<individual_id>.ase_table.tsv.gz` files.

### lorals

**Input:** Long-read RNA BAM and FASTQ files, phased donor genotype VCFs, a one-column donor list, genome and transcriptome FASTA files, gene BED annotations, and a gene-to-transcript mapping TSV.

**Workflow:** This workflow is implemented in `03_ASE_ASTS/run_LRS_ase_lorals.sh` and identifies allele-specific expression and allele-specific transcript-structure events from long-read RNA-seq data.

**Expected output:** Processed allele-specific expression (ASE) and allele-specific transcript-structure (ASTS) results.

### isoLASER

**Input:** Long-read genomic BAM files, the genome FASTA, a GTF annotation, a transcriptome reference, and a sample list.

**Workflow:** This workflow is located in `03_ASE_ASTS/isolaser/` and uses isoLASER to identify allele-specific splicing events.

**Expected output:** Per-sample `.mi_summary.tab` and `.mi_summary.filtered.tab` files, a merged genotyped gVCF, and joint splicing-linkage summaries.

### longcallR

**Input:** Long-read RNA FASTQ files, the reference FASTA, the GTF annotation, and matched DNA VCFs.

**Workflow:** This workflow is implemented in `03_ASE_ASTS/run_longcallR.sh` and performs high-confidence allele-specific expression and allele-specific splicing/junction analyses.

**Expected output:** Phased RNA VCF and BAM files, DNA-supported `<sample>.ase.tsv` files, and DNA-supported `<sample>.asj.tsv` files.

## 5. QTL mapping

### eQTL mapping

**Input:** Molecular phenotypes in BED, BED-compressed, or compatible Parquet format; covariates in tab-delimited text or dataframe format; and genotype data, preferably in PLINK2 PGEN/PVAR/PSAM format. Phenotype BED files should contain `chr`, `start`, `end`, and `phenotype_id` as the first four columns, followed by sample columns matching the genotype input. A cis-window can be defined around the TSS or around the supplied start and end positions. A BED template can be generated from a GTF annotation using `pyqtl`'s `io.gtf_to_tss_bed` function.

**Workflow:** This workflow is implemented in `04_QTL_mapping/run_eQTL_mapping.sh` and performs cis-eQTL mapping using nominal and permutation-based association testing.

**Expected output:** Nominal-mode (`cis_nominal`) summary statistics for all variant–phenotype pairs in Parquet format and permutation-mode (`cis`) phenotype-level summary statistics with empirical P-values for genome-wide FDR calculation.

### sQTL mapping

**Input:** Normalized splicing phenotypes, covariates, and matched genotype data prepared according to the QTL mapping input conventions described for eQTL mapping.

**Workflow:** This workflow is implemented in `04_QTL_mapping/run_sQTL_mapping.sh` and performs cis-sQTL mapping using the prepared splicing phenotypes.

**Expected output:** Cis-sQTL summary statistics and empirical significance estimates for downstream fine-mapping and functional analyses.

### TR-xQTL mapping

**Input:** Transcript-repeat phenotypes together with matched genotype and covariate data prepared according to the QTL mapping input conventions.

**Workflow:** This workflow is implemented in `04_QTL_mapping/run_TR_xQTL_mapping.sh` and maps genetic associations with transcript-repeat phenotypes.

**Expected output:** TR-xQTL association results for downstream fine-mapping and functional analyses.

### MAJIQTL mapping

**Input:** MAJIQ-derived splicing phenotypes and matched genotype and covariate data.

**Workflow:** This workflow is implemented in `04_QTL_mapping/run_MAJIQTL.sh` and maps genetic associations with MAJIQ-derived splicing phenotypes.

**Expected output:** MAJIQTL association results for downstream comparison and integration with other QTL modalities.

## 6. Fine-mapping

### SuSiE fine-mapping

**Input:** For each eGene identified from eQTL permutation results, a tab-delimited file containing individual-level normalized phenotypes with three columns corresponding to individual ID, individual ID, and the normalized phenotype value, together with a genotype matrix containing all variants within 1 Mb of the target eGene TSS.

**Workflow:** This workflow is located in `05_finemap/SuSiE_Fine_mapping/` and performs SuSiE fine-mapping to identify variants with high posterior inclusion probabilities within eQTL loci.

**Expected output:** A tissue-level tab-delimited table containing `locus_id`, `variant_id`, `pip`, `cs`, `cs_size`, `cs_purity`, and `Tissue` for all fine-mapped eGenes.

### Cross-ancestry SuSiEx fine-mapping

**Input:** Cross-ancestry genotype and phenotype summary data prepared for joint fine-mapping.

**Workflow:** This workflow is located in `05_finemap/Cross_ancestries_Fine_mapping_SuSiEX/` and applies SuSiEx to perform cross-ancestry joint fine-mapping.

**Expected output:** Cross-ancestry fine-mapping results and credible sets for comparison with ancestry-specific fine-mapping results and downstream analyses.

## 7. Downstream analyses

### Frequency differences

**Input:** Ancestry-group-specific allele frequencies from gnomAD and SuSiE fine-mapping results containing gene ID, variant ID, posterior inclusion probability (PIP), and credible-set (CS) ID.

**Workflow:** This workflow is located in `06_frequency_differences/` and identifies frequency-differentiated QTLs (fd-QTLs), compares allele frequencies, and consolidates signals into independent loci. The principal scripts are `a1.gtop_gnomad_af_comparison.R`, `a2.qtl_fine_mapping_summary.R`, `a3.qtl_mege_to_loci.R`, and `a4.fdQTL_table.R`.

**Expected output:** A list of fd-QTLs, including independent loci and their frequency-difference classifications across ancestry groups.

### eQTL portability

**Input:** GTOP eQTL summary statistics and GTEx eQTL summary statistics.

**Workflow:** This workflow is located in `07_eQTL_portability/` and evaluates eQTL portability using six portability metrics. The directory also contains mashR workflows for cross-tissue and cross-resource effect-size comparisons.

**Expected output:** Portability metrics, estimates of the proportion of portable eQTLs, and mashR-based cross-tissue and cross-resource effect-size comparison results.

### QTL enrichment

**Input:** For torus, gzip-compressed tab-delimited nominal QTL statistics generated by tensorQTL and a compatible variant annotation file. For S-LDSC, significant QTL pairs, fine-mapping tables containing variant ID, PIP, and credible-set assignments, chromosome-specific PLINK reference genotypes, baselineLD resources, weights, and munged GWAS `.sumstats.gz` files.

**Workflow:** The analyses are located in `08_QTL_enrichment/` and use torus and stratified LD score regression (S-LDSC) to evaluate the functional and disease relevance of QTLs.

**Expected output:** Torus annotation-enrichment parameter estimates with confidence intervals and S-LDSC trait-level regression results containing heritability contributions, enrichment statistics, and annotation coefficients.

### GWAS-QTL colocalization

**Input:** GWAS summary statistics and GWAS fine-mapping results, together with molecular QTL summary statistics and QTL fine-mapping results.

**Workflow:** This workflow is located in `09_colocalization/` and integrates GWAS and molecular QTL signals through fine-mapping and colocalization analyses. The principal scripts are `a1.finemapping_for_GWAS.R`, `a2.prepare_genes.R`, `a3.finemapping_for_QTL_GTOP_LD.R`, `a4.coloc.R`, `a5.merge_result.R`, and `run.sh`.

**Expected output:** GWAS-QTL colocalization results generated using SuSiE-coloc, supplemented with coloc.abf analyses where applicable.

### SMR

**Input:** GWAS summary statistics and QTL summary statistics.

**Workflow:** This workflow is located in `10_SMR/` and performs summary-data-based Mendelian randomization (SMR) and HEIDI analyses to evaluate genetically predicted molecular-trait associations with complex traits. The principal entry point is `run_SMR.sh`.

**Expected output:** SMR results containing gene ID, SMR P-value, and HEIDI P-value.

## Recommended execution order

The recommended order for running the GTOP analysis pipeline is:

1. `00_transcript_detection/`
2. `01_data_preparation/`
3. `02_Variant_calling/`
4. `03_ASE_ASTS/`
5. `04_QTL_mapping/`
6. `05_finemap/`
7. `06_frequency_differences/`
8. `07_eQTL_portability/`
9. `08_QTL_enrichment/`
10. `09_colocalization/`
11. `10_SMR/`

Individual downstream modules can also be run independently when their required upstream inputs are available.
