#!/usr/bin/env bash
# Edit the settings below, then run: bash run_longcallR.sh
set -euo pipefail

# Sample and input files. Use either bam or reads; leave the other empty.
sample="sample01"
bam="/path/to/sample01.sorted.bam"
reads=""
ref="/path/to/genome.fa"
gtf="/path/to/annotation.gtf"
dna_vcf="/path/to/sample01.dna.vcf.gz"

# Analysis settings.
preset="hifi-isoseq"
threads=8
outdir="results"

# Output paths.
sample_dir="$outdir/$sample"
snp_prefix="$sample_dir/snp_call/$sample"
ase_prefix="$sample_dir/ase/$sample"
asj_prefix="$sample_dir/asj/$sample"
phased_bam="$snp_prefix.phased.bam"
rna_vcf="$snp_prefix.vcf"

# Optional alignment of RNA FASTA/FASTQ reads.
align_reads() {
    bam="$snp_prefix.aligned.bam"
    if [ "$preset" = "hifi-isoseq" ]; then
        minimap2 -t "$threads" -ax splice:hq "$ref" "$reads" |
            samtools sort -@ "$threads" -o "$bam" -
    elif [ "$preset" = "hifi-masseq" ]; then
        minimap2 -t "$threads" -ax splice:hq -uf "$ref" "$reads" |
            samtools sort -@ "$threads" -o "$bam" -
    elif [ "$preset" = "ont-drna" ]; then
        minimap2 -t "$threads" -ax splice -uf -k14 "$ref" "$reads" |
            samtools sort -@ "$threads" -o "$bam" -
    else
        minimap2 -t "$threads" -ax splice "$ref" "$reads" |
            samtools sort -@ "$threads" -o "$bam" -
    fi
    samtools index -@ "$threads" "$bam"
}

# Call RNA SNPs, phase variants, and assign reads to haplotypes.
snp_call() {
    longcallR -b "$bam" -f "$ref" -o "$snp_prefix" \
        -t "$threads" -p "$preset"
    samtools index -@ "$threads" "$phased_bam"
}

# High-confidence ASE using matched DNA variants.
ase_cal() {
    longcallR ase -b "$phased_bam" -a "$gtf" -t "$threads" \
        --vcf1 "$rna_vcf" --vcf3 "$dna_vcf" -o "$ase_prefix"
}

# High-confidence ASJ using matched DNA variants.
asj_cal() {
    longcallR asj -b "$phased_bam" -a "$gtf" -f "$ref" -t "$threads" \
        --rna-vcf "$rna_vcf" --dna-vcf "$dna_vcf" -o "$asj_prefix"
}

# Run all steps sequentially. Stop immediately if a command fails.
mkdir -p "$outdir"
mkdir "$sample_dir"
mkdir "$sample_dir/snp_call" "$sample_dir/ase" "$sample_dir/asj"
if [ -n "$reads" ]; then
    echo "Aligning RNA reads..."
    align_reads
fi
echo "Calling and phasing RNA SNPs..."
snp_call
echo "Running high-confidence ASE..."
ase_cal
echo "Running high-confidence ASJ..."
asj_cal
echo "Completed: $sample_dir"
