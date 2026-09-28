#!/bin/bash

#SBATCH --job-name=dnn
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=7-00:00:00
#SBATCH --array=1-33
#SBATCH --output=tmp/dnn/%x_%A_%a_%N_%j.o
#SBATCH --error=tmp/dnn/%x_%A_%a_%N_%j.e

set -euo pipefail

THREADS=$SLURM_CPUS_PER_TASK

workdir='/path/to/myprojs/gtop/20251108-LR-RNAseq/diann'
image='/path/to/toolkits/bioapps/diann/v2.0/diann_docker.img'
diann='/mnt/v1.8.1/diann-1.8.1'
tissues="${workdir}/input/tissues.txt"
tissue=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$tissues")
tissue=$(basename "$tissue" .faa)
outdate=$(date "+%Y-%m-%d")
outdir="${workdir}/output/diann/${outdate}/${tissue}"
mkdir -p $outdir

# diann
function run_diann() {
    echo "#### start diann workflow at $(date "+%Y-%m-%d %H:%M:%S") ####"

    singularity exec --bind ${workdir}:/mnt $image bash -c "
    $diann \
    --predictor \
    --fasta /mnt/input/tissue_pr_faa.copy/${tissue}.faa \
    --fasta-search \
    --out-lib /mnt/output/diann/${outdate}/${tissue}/${tissue} \
    --cut K*,R*,!*P \
    --missed-cleavages 2 \
    --fixed-mod Carbamidomethylation,57.021464,C \
    --var-mod Oxidation,15.994915,M \
    --var-mod Acetylation,42.010565,*n \
    --var-mods 2 \
    --threads ${THREADS} && \
    $diann \
    --lib /mnt/output/diann/${outdate}/${tissue}/${tissue}.predicted.speclib \
    --fasta /mnt/input/tissue_pr_faa.copy/${tissue}.faa \
    --dir /mnt/input/rawfiles/${tissue} \
    --out /mnt/output/diann/${outdate}/${tissue}/${tissue} \
    --matrices \
    --qvalue 0.01 \
    --pg-level 0 \
    --matrix-qvalue 0.01 \
    --threads ${THREADS}"

    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}

# main func
function main() {
    module load singularity
    run_diann
}
main
