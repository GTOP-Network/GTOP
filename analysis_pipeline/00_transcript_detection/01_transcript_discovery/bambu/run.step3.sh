#!/bin/bash

#SBATCH --job-name=bb
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=7-00:00:00
#SBATCH --output=tmp/bb/%x_%A_%a_%N_%j.o
#SBATCH --error=tmp/bb/%x_%A_%a_%N_%j.e

set -euo pipefail

THREADS=$SLURM_CPUS_PER_TASK

workdir="/flashfs1/scratch.global/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/bambu"
samples=`ls /lustre/home/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/bambu/output/bambu/2026-01-31/AGTEX-*/*.rds`
outdate=$(date "+%Y-%m-%d")
outdir="${workdir}/output/bambu/${outdate}"
mkdir -p $outdir

# bambu
function run_bambu() {
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] start bambu workflow ..."
    cp ${workdir}/bambu.R ${outdir}/bambu.R
    echo "/lustre/home/lhgong/anaconda3/envs/longrw_bambu3/bin/Rscript ${outdir}/bambu.R ${THREADS} ${outdir} ${samples}" > ${outdir}/cmd.txt
    source /lustre/software/proxy/set_proxy.sh
    /lustre/home/lhgong/anaconda3/envs/longrw_bambu3/bin/Rscript \
    ${outdir}/bambu.R ${THREADS} ${outdir} ${samples} > ${outdir}/run.o 2> ${outdir}/run.e && \
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] done! "
}


# main func
function main() {
    source /lustre/home/lhgong/anaconda3/bin/activate longrw_bambu3
    run_bambu
}
main
