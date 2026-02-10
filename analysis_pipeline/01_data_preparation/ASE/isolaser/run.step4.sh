#!/bin/bash

#SBATCH --job-name=isol
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=7-00:00:00
#SBATCH --output=tmp/isol/%x_%A_%a_%N_%j.o
#SBATCH --error=tmp/isol/%x_%A_%a_%N_%j.e

set -euo pipefail

THREADS=$SLURM_CPUS_PER_TASK

workdir="/flashfs1/scratch.global/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/isolaser"
refg="${workdir}/ref/genome.fa"
outdir="${workdir}/output/isolaser/concat"
mkdir -p $outdir

# Run the combine vcf step
function combine_vcf() {
    echo "#### start combine_vcf at $(date "+%Y-%m-%d %H:%M:%S") ####"
    isolaser_combine_vcf -i ${workdir}/output/isolaser/fofn.txt \
    -f ${refg} -o ${outdir}/isolaser && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}
# Perform joint analysis
function exec_joint() {
    echo "#### start exec_joint at $(date "+%Y-%m-%d %H:%M:%S") ####"
    isolaser_joint -i ${workdir}/output/isolaser/fofn.txt \
    -t ${workdir}/output/db \
    -v ${outdir}/isolaser.merged.genotyped.gvcf.gz \
    -o ${outdir}/isolaser.merged.mi_summary.tab && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}


# main func
function main() {
    source /lustre/software/anaconda/anaconda3-2019.10-py37/bin/activate /lustre/home/cxue/.conda/envs/isoLASER
    module load samtools gatk
    combine_vcf
    exec_joint
}
main
