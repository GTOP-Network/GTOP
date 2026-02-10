#!/bin/bash

#SBATCH --job-name=ilp
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=7-00:00:00
#SBATCH --output=tmp/ilp/%x_%A_%a_%N_%j.o
#SBATCH --error=tmp/ilp/%x_%A_%a_%N_%j.e

set -euo pipefail

THREADS=$SLURM_CPUS_PER_TASK

workdir="/flashfs1/scratch.global/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/isolaser"
refg="${workdir}/ref/genome.fa"

gtf="${workdir}/input/enhanced_gtf/GTOP.sorted.gtf.gz"
fasta="${workdir}/input/enhanced_gtf/GTOP.fa"
outdir="${workdir}/output/db"
mkdir -p $outdir

# generate transcriptome reference from GTF
function convert_gtf_to_fasta() {
    echo "#### start convert_gtf_to_fasta at $(date "+%Y-%m-%d %H:%M:%S") ####"
    isolaser_convert_gtf_to_fasta -g ${gtf} -f ${refg} -o ${fasta} && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}
# Extract exonic parts from GTF
function extract_exonic_parts() {
    echo "#### start extract_exonic_parts at $(date "+%Y-%m-%d %H:%M:%S") ####"
    isolaser_extract_exon_parts -g ${gtf} -o ${outdir} && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}


# main func
function main() {
    source /lustre/software/anaconda/anaconda3-2019.10-py37/bin/activate /lustre/home/cxue/.conda/envs/isoLASER
    extract_exonic_parts
}
main
