#!/bin/bash

#SBATCH --job-name=pbb
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:00:00
#SBATCH --array=1-176
#SBATCH --output=tmp/pbb/%x_%A_%a_%N_%j.o
#SBATCH --error=tmp/pbb/%x_%A_%a_%N_%j.e

set -euo pipefail

THREADS=$SLURM_CPUS_PER_TASK

workdir="/lustre/home/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/bambu"

sampleids="${workdir}/input/sampleids.txt"
sampleid=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$sampleids")
bam="/lustre/home/cxue/project/GMTiP-RNA/20251031/long_read/HPC/output/sample_based/${sampleid}/${sampleid}.flnc_mapped.bam"
outdate=$(date "+%Y-%m-%d")
outdir="${workdir}/output/bambu/${outdate}"
subdir="${outdir}/${sampleid}"
mkdir -p $subdir

# bambu
function run_bambu_preprocessing() {
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] start bambu-preprocessing workflow ..."
    cp ${workdir}/preprocessing.R ${subdir}/preprocessing.R
    echo "/lustre/home/lhgong/anaconda3/envs/longrw_bambu3/bin/Rscript ${subdir}/preprocessing.R ${THREADS} ${subdir} ${bam}" > ${subdir}/cmd.txt
    source /lustre/software/proxy/set_proxy.sh
    /lustre/home/lhgong/anaconda3/envs/longrw_bambu3/bin/Rscript \
    ${subdir}/preprocessing.R ${THREADS} ${subdir} ${bam} > ${subdir}/run.o 2> ${subdir}/run.e && \
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] done! "
}


# main func
function main() {
    source /lustre/home/lhgong/anaconda3/bin/activate longrw_bambu3
    run_bambu_preprocessing
}
main
