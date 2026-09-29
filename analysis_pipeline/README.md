# GTOP Analysis Pipeline: Detailed Documentation

This document provides a detailed description of the GTOP data-processing and downstream-analysis workflows, including the principal inputs, analytical procedures, and expected outputs for each stage.

## 1. Transcript detection

### Transcript discovery

The transcript discovery workflow is located in `00_transcript_detection/01_transcript_discovery/` and applies seven complementary long-read transcript discovery approaches: Bambu, FLAIR, FLAMES, IsoQuant, Iso-Seq, IsoTools, and TALON. The inputs include PacBio full-length non-chimeric (FLNC) reads in FASTQ or BAM format, their alignments to hg38, the hg38 reference genome FASTA, and the GENCODE v47 annotation GTF. The expected outputs are per-sample transcript annotations in GTF format generated independently by each of the seven tools, together with caller-specific read-support information for subsequent transcript filtering and integration.

### Transcript integration

The transcript integration workflow is located in `00_transcript_detection/02_merge_filter_transcript/` and performs FLNC support filtering, TAMA-based transcript merging, SQANTI3 annotation, and construction of the enhanced GTOP transcript reference. The inputs include transcript annotations and read-support information from the seven discovery tools, PacBio FLNC reads and their genome alignments, short-read splice-junction evidence from STAR, the hg38 genome, and the GENCODE v47 annotation. CAGE/TSS annotations and a poly(A) motif list are also used for transcript quality assessment. The expected outputs are a merged and quality-filtered GTOP transcript catalogue and an enhanced reference combining GTOP novel transcripts with GENCODE v47 annotations. Additional outputs include transcript annotations in GTF format, transcript sequences in FASTA format, associated annotation and quality-control tables, and predicted protein sequences.

### Transcript quantification

The transcript quantification workflow is located in `00_transcript_detection/03_quantification/` and uses FLAIR to quantify GTOP transcripts. The inputs include per-sample PacBio FLNC FASTQ files and the final GTOP transcript reference, including transcript sequences and associated annotations. The expected outputs are transcript-level read counts quantified using FLAIR and combined across samples into count and TPM matrices, together with gene-level count and TPM matrices obtained by aggregating transcripts belonging to the same gene.

### Peptide validation

The peptide validation workflow is located in `00_transcript_detection/04_peptide_validation/` and uses DIA-NN-based proteomic evidence to evaluate protein support for GTOP transcripts. The inputs include tissue-matched DIA mass-spectrometry data in mzML format, GENCODE v47 protein sequences, predicted GTOP protein sequences, transcript annotations, long-read transcript TPM estimates, and tissue metadata. The expected outputs are tissue-specific protein databases containing coding transcripts with TPM > 5 in at least one corresponding long-read sample, DIA-NN reports providing peptide evidence for the searched proteins, and per-tissue protein abundance matrices.

## 2. RNA phenotype preparation

### Gene-level quantification

The gene-level quantification workflow is implemented in `01_data_preparation/` and generates gene-expression phenotypes from short-read RNA-seq data. The input is raw paired-end RNA-seq FASTQ data for each sample, together with the associated sample and reference information required by the workflow. The expected outputs are gene-level read-count and transcripts-per-million (TPM) matrices.

### Splicing quantification

The splicing quantification workflow is implemented in `01_data_preparation/` and generates splice-junction and intron-cluster phenotypes for downstream QTL analyses. The inputs include paired-end RNA-seq FASTQ files, donor-specific DNA VCFs, a STAR genome index, the reference genome FASTA, and a gene annotation GTF. Tissue metadata are provided as a headerless two-column file containing tissue name and tissue code, and the sample-to-participant lookup is provided as a headerless TSV containing sample ID and participant ID. Exon annotations are supplied as a headered TSV containing `chr`, `start`, `end`, `strand`, and `gene_id`, while gene metadata are supplied as a headerless TSV containing `gene_name`, `gene_type`, `gene_id`, `chr`, `start`, `end`, and `strand`. The expected outputs are WASP-filtered alignments, junction counts, intron-cluster ratios, and filtered and normalized splicing phenotypes.

### Transcript quantification

The transcript quantification workflow is implemented in `01_data_preparation/` and generates transcript-level expression and transcript-usage phenotypes from short-read RNA-seq data. The inputs include paired-end short-read RNA-seq FASTQ files, the hg38 reference genome FASTA, the enhanced GTOP novel–GENCODE v47 transcript reference in GTF and FASTA format with associated gene–transcript annotations, and tissue metadata. The expected outputs are transcript-level count and TPM matrices generated using Salmon and RSEM, tissue-specific transcript TPM matrices obtained by averaging the two methods' estimates and retaining transcripts that pass concordance filters, and filtered, imputed, and normalized transcript-usage phenotypes for downstream QTL analysis.

## 3. Variant calling

### LRS/SRS variant calling

The variant-calling workflow is located in `02_Variant_calling/` and processes both long-read and short-read whole-genome sequencing data. The inputs are raw whole-genome sequencing data for each sample, including long-read PacBio HiFi WGS BAM files and short-read WGS FASTQ files. The expected outputs are per-sample genotype VCF files generated from LRS and SRS WGS data, together with a cohort-level merged and site-filtered VCF containing variants and genotypes across all samples.

### Population structure analysis

Population-structure analysis is performed within `02_Variant_calling/` using genotype data derived from whole-genome sequencing. The workflow generates PCA and admixture results for characterizing population structure and for downstream analyses requiring ancestry information.

### Variant annotation

Variant annotation is performed within `02_Variant_calling/` using the cohort-level variant set as input. The workflow applies VEP to generate an annotated variant set for downstream interpretation.

The principal entry points for this module are `01_LRS_WGS_Variant_Calling.sh`, `02_SRS_WGS_Variant_Calling.sh`, `03_PCA_ADMIXTURE.sh`, and `04_vep_annotation.sh`.

## 4. ASE, ASTS, and ASJ

### Short-read ASE

The short-read allele-specific expression workflow is implemented in `03_ASE_ASTS/run_SRS_ase.sh`. The inputs include WASP-filtered RNA BAM files, donor-specific heterozygous DNA VCFs, the reference FASTA, gene BED annotations containing chromosome, start, end, and gene ID, and an ASE manifest containing `sample_id`, `tissue_site_detail`, and `ase_readcount_file`. The expected outputs are per-site allelic counts, donor-level LAMP estimates, and `<individual_id>.ase_table.tsv.gz` files.

### lorals

The long-read allele-specific expression and allele-specific transcript-structure workflow is implemented in `03_ASE_ASTS/run_LRS_ase_lorals.sh`. The inputs include long-read RNA BAM and FASTQ files, phased donor genotype VCFs, a one-column donor list, genome and transcriptome FASTA files, gene BED annotations, and a gene-to-transcript mapping TSV. The expected outputs are processed allele-specific expression (ASE) and allele-specific transcript-structure (ASTS) results.

### isoLASER

The allele-specific splicing workflow using isoLASER is located in `03_ASE_ASTS/isolaser/`. The inputs include long-read genomic BAM files, the genome FASTA, a GTF annotation, a transcriptome reference, and a sample list. The expected outputs are per-sample `.mi_summary.tab` and `.mi_summary.filtered.tab` files, a merged genotyped gVCF, and joint splicing-linkage summaries.

### longcallR

The longcallR workflow is implemented in `03_ASE_ASTS/run_longcallR.sh` and is used for high-confidence allele-specific expression and allele-specific splicing/junction analyses. The inputs include long-read RNA FASTQ files, the reference FASTA, the GTF annotation, and matched DNA VCFs. The expected outputs are phased RNA VCF and BAM files, DNA-supported `<sample>.ase.tsv` files, and DNA-supported `<sample>.asj.tsv` files.

## 5. QTL mapping

### eQTL mapping

The eQTL mapping workflow is implemented in `04_QTL_mapping/run_eQTL_mapping.sh`. Phenotypes must be provided in BED or BED-compressed format, with a single header line beginning with `#`; the first four columns must correspond to `chr`, `start`, `end`, and `phenotype_id`, followed by sample columns whose identifiers match those in the genotype input. BED input in Parquet format is also supported. The BED file can specify the center of the cis-window, usually the TSS, using `start == end - 1`, or can provide start and end positions, in which case the cis-window is defined as `[start-window, end+window]`. A BED template can be generated from a GTF annotation using `pyqtl`'s `io.gtf_to_tss_bed` function. Covariates can be supplied as a tab-delimited text file or dataframe, with row and column headers. Genotypes are preferably provided in PLINK2 PGEN/PVAR/PSAM format. The expected outputs include nominal-mode (`cis_nominal`) summary statistics for all variant–phenotype pairs in Parquet format and permutation-mode (`cis`) phenotype-level summary statistics with empirical P-values, enabling calculation of genome-wide FDR.

### sQTL mapping

The sQTL mapping workflow is implemented in `04_QTL_mapping/run_sQTL_mapping.sh` and uses normalized splicing phenotypes, covariates, and matched genotype data. The same phenotype, covariate, and genotype input conventions described for eQTL mapping apply. The expected outputs are cis-sQTL summary statistics and empirical significance estimates for downstream analyses.

### TR-xQTL mapping

The transcript-repeat QTL workflow is implemented in `04_QTL_mapping/run_TR_xQTL_mapping.sh`. The inputs include transcript-repeat phenotypes together with matched genotype and covariate data prepared according to the QTL mapping input conventions. The expected outputs are TR-xQTL mapping results for downstream fine-mapping and functional analyses.

### MAJIQTL mapping

The MAJIQTL workflow is implemented in `04_QTL_mapping/run_MAJIQTL.sh` and maps genetic associations with MAJIQ-derived splicing phenotypes. The inputs include MAJIQ-derived splicing phenotypes and matched genotype and covariate data. The expected output is a set of MAJIQTL association results for downstream comparison and integration with other QTL modalities.

## 6. Fine-mapping

### SuSiE fine-mapping

The SuSiE fine-mapping workflow is located in `05_finemap/SuSiE_Fine_mapping/` and identifies variants with high posterior inclusion probabilities within eQTL loci. For each eGene identified from eQTL permutation results, the first input is a tab-delimited file containing individual-level normalized phenotypes, with three columns corresponding to individual ID, individual ID, and the normalized phenotype value (for example, TMM-normalized expression). The second input is a genotype matrix containing genotypes for all variants within 1 Mb of the TSS of the target eGene. The expected output is a tissue-level tab-delimited table summarizing all eGenes, with columns for `locus_id`, `variant_id`, `pip`, `cs`, `cs_size`, `cs_purity`, and `Tissue`.

### Cross-ancestry SuSiEx fine-mapping

The cross-ancestry fine-mapping workflow is located in `05_finemap/Cross_ancestries_Fine_mapping_SuSiEX/` and applies SuSiEx to cross-ancestry genotype and phenotype summary data prepared for joint analysis. The expected outputs are cross-ancestry fine-mapping results and credible sets that can be compared with ancestry-specific fine-mapping results and used in downstream analyses.

## 7. Downstream analyses

### Frequency differences

The frequency-difference analysis is located in `06_frequency_differences/` and identifies frequency-differentiated QTLs (fd-QTLs). The inputs include ancestry-group-specific allele frequencies from gnomAD and SuSiE fine-mapping results containing gene ID, variant ID, posterior inclusion probability (PIP), and credible-set (CS) ID. The expected outputs are a list of fd-QTLs, including independent loci and their frequency-difference classifications across ancestry groups.

The principal scripts are `a1.gtop_gnomad_af_comparison.R`, `a2.qtl_fine_mapping_summary.R`, `a3.qtl_mege_to_loci.R`, and `a4.fdQTL_table.R`.

### eQTL portability

The eQTL portability workflow is located in `07_eQTL_portability/` and evaluates the portability of GTOP eQTL effects relative to GTEx. The inputs are GTOP eQTL summary statistics and GTEx eQTL summary statistics. The expected outputs include evaluation results for six portability metrics and estimates of the proportion of portable eQTLs. The directory also contains mashR workflows for cross-tissue and cross-resource effect-size comparisons.

### QTL enrichment

The enrichment analyses are located in `08_QTL_enrichment/` and use torus and stratified LD score regression (S-LDSC) to assess the functional and disease relevance of QTLs. For torus, the inputs are gzip-compressed, tab-delimited nominal QTL statistics generated by tensorQTL and a compatible variant annotation file. For S-LDSC, the inputs include significant QTL pairs with variant IDs in column 2, a headered fine-mapping table with variant ID, PIP, and credible-set assignment in columns 2, 3, and 4, respectively, chromosome-specific PLINK reference genotypes, baselineLD resources, weights, and munged GWAS `.sumstats.gz` files. The expected torus outputs are annotation-enrichment parameter estimates with confidence intervals, whereas S-LDSC produces trait-level regression results containing heritability contributions, enrichment statistics, and annotation coefficients.

### GWAS-QTL colocalization

The colocalization workflow is located in `09_colocalization/` and integrates GWAS and molecular QTL signals through fine-mapping and colocalization analyses. The inputs include GWAS summary statistics and GWAS fine-mapping results, together with QTL summary statistics and QTL fine-mapping results. The expected outputs are colocalization results generated using SuSiE-coloc, supplemented with coloc.abf analyses where applicable.

The principal scripts are `a1.finemapping_for_GWAS.R`, `a2.prepare_genes.R`, `a3.finemapping_for_QTL_GTOP_LD.R`, `a4.coloc.R`, `a5.merge_result.R`, and `run.sh`.

### SMR

The summary-data-based Mendelian randomization workflow is located in `10_SMR/`. The inputs are GWAS summary statistics and QTL summary statistics. The workflow performs SMR and HEIDI analyses to evaluate genetically predicted molecular-trait associations with complex traits. The expected outputs include SMR results containing gene ID, SMR P-value, and HEIDI P-value.

The principal entry point is `run_SMR.sh`.

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
