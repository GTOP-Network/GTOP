# GTOP ASE and ASTS analysis

This directory contains code for the ASE and ASTS analyses performed as part of the GTOP publication. The workflow is divided into four main parts:

1. **LRS ASE/ASTS analysis with lorals**: Long-read ASE and ASTS quantification using [lorals](https://github.com/LappalainenLab/lorals)
2. **LRS splicing linkage analysis with isoLASER**: Long-read allele-specific transcript splicing using [isoLASER](https://github.com/gxiaolab/isoLASER)
3. **SRS ASE analysis**: Variant-aware read mapping with STAR and WASP filtering
4. **LRS ASE/ASJ analysis with longcallR**: RNA SNP calling and phasing, followed by high-confidence ASE and ASJ analysis using matched DNA variants

---

## Part 1: LRS ASE/ASTS analysis with lorals

The `run_LRS_ASE_lorals.sh` script processes long-read sequencing data to quantify ASE and ASTS using [lorals](https://github.com/LappalainenLab/lorals). The pipeline consists of six main steps, each implemented as a function within the `run_LRS_ASE_lorals.sh` master script.

|Function|Description|
|---|---|
|`process_vcf`|Processes the VCF file containing variants from all donors and generates a reference genome per haplotype for each donor.|
|`hap_map`|Performs haplotype-aware mapping of long reads to the reference genome per haplotype of each donor.|
|`hap_map_trans`|Aligns reads to the transcriptome.|
|`ase_cal_chr`|Calculates and annotates the allelic coverage of each variant.|
|`asts_cal_quant`|Calculates the number of reads containing the reference or alternate allele assigned to each transcript.|
|`process_asts`|Aggregates and processes the ASTS quantification results, performing filtering and statistical tests.|

---

## Part 2: LRS splicing linkage analysis with isoLASER

The `isolaser` directory contains scripts that perform splicing linkage analysis using [isoLASER](https://github.com/gxiaolab/isoLASER).

|Script|Description|
|---|---|
|`run.step1.sh`|Uses the GTF file to generate a transcriptome reference for alignment and extract exonic parts.|
|`run.step2.sh`|Annotates the BAM file and runs isoLASER.|
|`run.step3.sh`|Creates fofn.txt containing information about the individual samples.|
|`run.step4.sh`|Runs isoLASER in joint mode.|

---

## Part 3: SRS ASE analysis with STAR and WASP

The `run_SRS_align_ASE.sh` script completes ASE analysis for Short-Read Sequencing, adapted from the GTEx project's ASE workflow. The pipeline consists of three main steps, each implemented as a function within the `run_SRS_align_ASE.sh` master script.

|Function|Description|
|---|---|
|`ase_snp_level`|Uses GATK to calculate reference and alternate allele counts at heterozygous sites for each sample.|
|`ase_snp_cal_lamp`|Calculates the global foreign allele frequency (lamp value) per individual, based on GTEx's `ase_calculate_lamp.py`.|
|`ase_snp_sum`|Integrates ASE data for each donor, applies quality filters and performs statistical tests for each sample, inspired by GTEx's `ase_aggregate_by_individual.py`.|

---

## Part 4: LRS ASE/ASJ analysis with longcallR

The `run_longcallR.sh` script processes long-read sequencing data to identify allele-specific expression (ASE) and allele-specific junctions (ASJ) using [longcallR](https://github.com/huangnengCSU/longcallR). The pipeline consists of four main steps, each implemented as a function within the `run_longcallR.sh` master script. Both ASE and ASJ analyses use the high-confidence mode with a matched DNA VCF.

|Function|Description|
|---|---|
|`align_reads`|Aligns long reads to the reference genome using minimap2 and generates a sorted, indexed BAM file. This step is skipped when an aligned BAM is provided.|
|`snp_call`|Calls and phases RNA SNPs using longcallR, generating a phased RNA VCF and a haplotype-tagged BAM file.|
|`ase_cal`|Quantifies high-confidence ASE using the phased BAM, phased RNA VCF and matched DNA VCF (`--vcf1` and `--vcf3`).|
|`asj_cal`|Identifies high-confidence ASJ using the phased BAM, phased RNA VCF and matched DNA VCF (`--rna-vcf` and `--dna-vcf`).|


