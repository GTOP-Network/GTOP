# GTOP ASE and ASTS analysis

This directory contains code for the ASE and ASTS analyses performed as part of the GTOP publication. The workflow is divided into three main:

1. **SRS ASE analysis**: Variant-aware read mapping with STAR and WASP filtering
    
2. **LRS ASE/ASTS analysis with lorals**: Long-read ASE and ASTS quantification using [lorals](https://github.com/LappalainenLab/lorals)
    
3. **LRS splicing linkage analysis with isoLASER**: Long-read allele-specific transcript splicing using [isoLASER](https://github.com/gxiaolab/isoLASER)
    
  

---


## Part 1: LRS ASE/ASTS analysis with lorals

The `run_LRS_ASE_lorals.sh` script processes long-read sequencing data to quantify ASE and ASTS using the [lorals](https://github.com/LappalainenLab/lorals). The pipeline consists of six main steps, each implemented as a function within the `run_LRS_ASE_lorals.sh` master script.


|Function|Description|
|---|---|
|`process_vcf`|Processes the VCF file containing all donors variants and generates a reference genome per haplotype of each donor.|
|`hap_map`|Performs haplotype-aware mapping of long reads to the reference genome per haplotype of each donor.|
|`hap_map_trans`|Aligns reads to the transcriptome.|
|`ase_cal_chr`|Calculates the allelic coverage of each variant and annotated.|
|`asts_cal_quant`|Calculates the number of reads containing the ref or alt allele assigned to each transcript.|
|`process_asts`|Aggregates and processes the ASTS quantification results, performing filtering and statistical tests.|

---

## Part 2: LRS splicing linkage analysis with isoLASER

The isolaser directory contains script performs splicing linkage analysis using [isoLASER](https://github.com/gxiaolab/isoLASER). 

|script|Description|
|---|---|
|`run.step1.sh`|Use the GTF file to generate a transcriptome reference for alignment and extract exonic parts from GTF.|
|`run.step2.sh`|Annotate bam file and run IsoLASER..|
|`run.step3.sh`|Make fofn.txt which contains informaiton of the individual samples|
|`run.step4.sh`|Run IsoLASER joint mode.|


---
## Part 3: SRS ASE analysis with STAR and WASP

The `run_SRS_align_ASE.sh` script completes ASE analysis for Short-Read Sequencing, adapted from the GTEx project's ASE workflow. The pipeline consists of three main steps, each implemented as a function within the `run_SRS_align_ASE.sh` master script.

|Function|Description|
|---|---|
|`ase_snp_level`|Uses GATK to calculate reference and alternate allele counts at heterozygous sites for each sample|
|`ase_snp_cal_lamp`|Calculates the global foreign allele frequency (lamp value) per individual, based on GTEx's `ase_calculate_lamp.py`|
|`ase_snp_sum`|Integrates ASE data for each donor, applies quality filters and performs statistical tests for each sample inspired by GTEx's `ase_aggregate_by_individual.py`|


