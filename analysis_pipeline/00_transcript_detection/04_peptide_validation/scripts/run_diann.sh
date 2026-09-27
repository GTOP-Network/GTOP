#!/usr/bin/env bash
#SBATCH --job-name=diann
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --time=7-00:00:00
set -euo pipefail
: "${GTOP_PROJECT_DIR:?Set GTOP_PROJECT_DIR}"
: "${DIANN_IMAGE:?Set the Singularity/Apptainer image containing DIA-NN 1.8.1}"
: "${DIANN_BIN:?Set the path of DIA-NN 1.8.1 inside the container}"
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    CURRENT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
    list="$GTOP_PROJECT_DIR/output/MS/diann/input/tissue_fa_list.txt"
    n=$(awk 'NF {n++} END {print n+0}' "$list")
    [[ "$n" -gt 0 ]] || { echo 'No tissue FASTA files' >&2; exit 1; }
    log_dir="$GTOP_PROJECT_DIR/run_log/diann"
    mkdir -p "$log_dir"
    sbatch --partition="${SLURM_PARTITION:-cpuPartition}" --array="1-$n" \
      --output="$log_dir/%x_%A_%a.out" --error="$log_dir/%x_%A_%a.err" "$CURRENT_DIR/run_diann.sh"
    exit 0
fi
THREADS="${SLURM_CPUS_PER_TASK:-30}"
workdir="$GTOP_PROJECT_DIR/output/MS/diann"
run_tag="${DIANN_RUN_TAG:-run}"
fa_file_name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$workdir/input/tissue_fa_list.txt")
[[ -n "$fa_file_name" ]] || { echo 'Invalid tissue array index' >&2; exit 1; }
tissue=$(basename "$fa_file_name" .faa)
prefix="/mnt/output/diann/tissues/$run_tag/$tissue/$tissue"
mkdir -p "$workdir/output/diann/tissues/$run_tag/$tissue"
container=("${CONTAINER_RUNTIME:-singularity}" exec --bind "$workdir:/mnt" "$DIANN_IMAGE" "$DIANN_BIN")
"${container[@]}" --predictor --fasta "/mnt/input/tissue_faa/$fa_file_name" --fasta-search \
  --out-lib "$prefix" --cut 'K*,R*,!*P' --missed-cleavages 2 \
  --fixed-mod Carbamidomethylation,57.021464,C --var-mod Oxidation,15.994915,M \
  --var-mod Acetylation,42.010565,*n --var-mods 2 --threads "$THREADS" --no-quant-files
mv "$workdir/output/diann/tissues/$run_tag/$tissue/$tissue.log.txt" \
   "$workdir/output/diann/tissues/$run_tag/$tissue/$tissue.lib.log.txt"
"${container[@]}" --lib "$prefix.predicted.speclib" --fasta "/mnt/input/tissue_faa/$fa_file_name" \
  --dir "/mnt/input/rawfiles/tissues/$tissue" --out "$prefix" --matrices --qvalue 0.01 \
  --pg-level 0 --matrix-qvalue 0.01 --threads "$THREADS" --no-quant-files
