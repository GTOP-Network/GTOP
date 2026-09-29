# GTOP Analysis Pipeline: Detailed Documentation

This document provides a detailed map of the GTOP data-processing and downstream-analysis workflows, including the principal inputs and expected outputs for each stage.

## 1. Transcript detection

| Module | Directory | Input | Expected output |
|---|---|---|---|
| Transcript discovery | `00_transcript_detection/01_transcript_discovery/` | PacBio FLNC reads (FASTQ/BAM), hg38-aligned reads, hg38 reference FASTA, GENCODE v47 GTF | Independent per-sample transcript GTFs and read-support information from Bambu, FLAIR, FLAMES, IsoQuant, Iso-Seq, IsoTools, and TALON |
| Transcript integration | `00_transcript_detection/02_merge_filter_transcript/` | Seven-tool transcript annotations/support; FLNC reads/alignments; STAR splice-junction evidence; hg38; GENCODE v47; CAGE/TSS annotations; poly(A) motif list | Merged and quality-filtered GTOP transcript catalogue; enhanced GTF/FASTA reference; annotation/QC tables; predicted proteins |
| Transcript quantification | `00_transcript_detection/03_quantification/` | PacBio FLNC FASTQ files and final GTOP transcript reference | Transcript-level counts and TPM; gene-level counts and TPM |
| Peptide validation | `00_transcript_detection/04_peptide_validation/` | Tissue-matched DIA mzML files; GENCODE/GTOP protein sequences; transcript annotations; LR transcript TPM; tissue metadata | Tissue-specific protein databases, DIA-NN peptide evidence, and protein abundance matrices |

## 2. RNA phenotype preparation

| Module | Directory | Input | Expected output |
|---|---|---|---|
| Gene-level quantification | `01_data_preparation/` | Paired-end RNA-seq FASTQ files and associated sample/reference information | Gene-level read counts and TPM |
| Splicing quantification | `01_data_preparation/` | Paired-end RNA-seq FASTQ files, donor-specific VCFs, STAR genome index, reference FASTA/GTF, tissue/sample metadata, exon/gene annotations | WASP-filtered alignments, junction counts, intron-cluster ratios, and filtered/normalized splicing phenotypes |
| Transcript quantification | `01_data_preparation/` | Paired-end RNA-seq FASTQ files, hg38 FASTA, enhanced GTOP–GENCODE v47 transcript reference, tissue metadata | Salmon/RSEM transcript count and TPM matrices; tissue-specific transcript TPM; filtered, imputed, normalized transcript-usage phenotypes |

## 3. Variant calling

| Module | Directory | Input | Expected output |
|---|---|---|---|
| LRS/SRS variant calling | `02_Variant_calling/` | LRS HiFi WGS BAMs and SRS WGS FASTQs | Per-sample genotype VCFs and cohort-level merged/site-filtered VCF |
| Population structure | `02_Variant_calling/` | Genotype data derived from WGS | PCA/admixture results used for population-structure characterization |
| Variant annotation | `02_Variant_calling/` | Cohort variant set | VEP-annotated variant set |

Principal entry points include `01_LRS_WGS_Variant_Calling.sh`, `02_SRS_WGS_Variant_Calling.sh`, `03_PCA_ADMIXTURE.sh`, and `04_vep_annotation.sh`.

## 4. ASE, ASTS, and ASJ

| Workflow | Directory | Input | Expected output |
|---|---|---|---|
| Short-read ASE | `03_ASE_ASTS/run_SRS_ase.sh` | WASP-filtered RNA BAMs, donor-specific heterozygous VCFs, reference FASTA, gene BED, ASE manifest | Per-site allelic counts, donor-level LAMP estimates, and `<individual_id>.ase_table.tsv.gz` |
| lorals | `03_ASE_ASTS/run_LRS_ase_lorals.sh` | Long-read RNA BAM/FASTQ, phased donor genotype VCFs, donor list, genome/transcriptome FASTA, gene BED, gene–transcript mapping | Processed allele-specific expression and allele-specific transcript usage results |
| isoLASER | `03_ASE_ASTS/isolaser/` | Long-read genomic BAMs, genome FASTA, GTF, transcriptome reference, sample list | Per-sample `.mi_summary.tab` and `.mi_summary.filtered.tab`, merged genotyped gVCF, and joint splicing-linkage summaries |
| longcallR | `03_ASE_ASTS/run_longcallR.sh` | Long-read RNA FASTQ, reference FASTA, GTF, matched DNA VCFs | Phased RNA VCF/BAM, DNA-supported `<sample>.ase.tsv`, and DNA-supported `<sample>.asj.tsv` |

## 5. QTL mapping

| Module | Directory | Input | Expected output |
|---|---|---|---|
| eQTL | `04_QTL_mapping/run_eQTL_mapping.sh` | Phenotypes in BED/Parquet, covariates, PLINK2 genotype data | Nominal cis-eQTL summary statistics and permutation-mode statistics with empirical P-values |
| sQTL | `04_QTL_mapping/run_sQTL_mapping.sh` | Splicing phenotypes, covariates, genotype data | cis-sQTL summary statistics and empirical significance estimates |
| TR-xQTL | `04_QTL_mapping/run_TR_xQTL_mapping.sh` | Transcript-repeat phenotypes and matched genotype/covariate data | TR-xQTL mapping results |
| MAJIQTL | `04_QTL_mapping/run_MAJIQTL.sh` | MAJIQ-derived splicing phenotypes and matched genotype/covariate data | MAJIQTL QTL results |

Phenotype BED files use a header beginning with `#`; the first four columns are `chr`, `start`, `end`, and `phenotype_id`, followed by sample columns whose identifiers match the genotype input. Covariates can be supplied as tab-delimited text or a dataframe. PLINK2 PGEN/PVAR/PSAM is the preferred genotype format.

## 6. Fine-mapping

| Module | Directory | Input | Expected output |
|---|---|---|---|
| SuSiE | `05_finemap/SuSiE_Fine_mapping/` | Per-locus normalized individual-level phenotypes and genotype matrices for variants within 1 Mb of the target eGene TSS | Tissue-level tables containing locus ID, variant ID, PIP, credible-set ID, credible-set size, credible-set purity, and tissue |
| Cross-ancestry SuSiEx | `05_finemap/Cross_ancestries_Fine_mapping_SuSiEX/` | Cross-ancestry genotype/phenotype summary data prepared for SuSiEx | Cross-ancestry fine-mapping results and credible sets |

## 7. Downstream analyses

| Module | Directory | Input | Expected output |
|---|---|---|---|
| Frequency differences | `06_frequency_differences/` | gnomAD ancestry-specific allele frequencies and SuSiE fine-mapping results | Frequency-differentiated QTLs and ancestry-group classifications |
| eQTL portability | `07_eQTL_portability/` | GTOP and GTEx eQTL summary statistics | Portability metrics and estimated proportion of portable eQTLs |
| Enrichment | `08_QTL_enrichment/` | QTL statistics/annotations for torus; QTL pairs, fine-mapping, LD/reference resources, and GWAS summary statistics for S-LDSC | torus enrichment estimates; S-LDSC heritability/enrichment and annotation coefficients |
| Colocalization | `09_colocalization/` | GWAS and QTL summary statistics and fine-mapping results | SuSiE-coloc and coloc.abf colocalization results |
| SMR | `10_SMR/` | GWAS and QTL summary statistics | SMR and HEIDI results |

## Principal scripts

### `06_frequency_differences/`

- `a1.gtop_gnomad_af_comparison.R`
- `a2.qtl_fine_mapping_summary.R`
- `a3.qtl_mege_to_loci.R`
- `a4.fdQTL_table.R`

### `07_eQTL_portability/`

Contains GTOP–GTEx overlap, portability, and mashR workflows.

### `08_QTL_enrichment/`

Contains torus and S-LDSC analyses for assessing QTL enrichment and heritability-related signals.

### `09_colocalization/`

- `a1.finemapping_for_GWAS.R`
- `a2.prepare_genes.R`
- `a3.finemapping_for_QTL_GTOP_LD.R`
- `a4.coloc.R`
- `a5.merge_result.R`
- `run.sh`

### `10_SMR/`

- `run_SMR.sh`

## Execution order

The recommended order is:

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

Individual downstream modules can be run independently when their required upstream inputs are available.
