#!/usr/bin/env bash
set -euo pipefail
: "${GTOP_PROJECT_DIR:?Set GTOP_PROJECT_DIR using config.example.sh}"
CURRENT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
OUT_DIR="$GTOP_PROJECT_DIR/output/assembly/LRS/tx_merge/merged_tools/sqanti3_filter_1"
SRS_JUNCTION_DIR="${SRS_JUNCTION_DIR:-$GTOP_PROJECT_DIR/input/srs_junctions}"
mkdir -p "$OUT_DIR/SRS_junction_filter"
cd "$OUT_DIR/SRS_junction_filter"
case "${1:-}" in
prepare)
    python "$CURRENT_DIR/junction_srs.py" build-index --gtf "$OUT_DIR/filtered.gtf" --out-dir junction_index
    ;;
submit)
    python "$CURRENT_DIR/junction_srs_run.py" --sj-dir "$SRS_JUNCTION_DIR" \
      --index-dir "$OUT_DIR/SRS_junction_filter/junction_index" \
      --root-dir "$OUT_DIR/SRS_junction_filter/srs_jobs" \
      --partition "${SLURM_PARTITION:-cpuPartition}" --cpu-per-node 32 --sample-processes 32 --samples-per-job 32
    ;;
merge)
    mkdir -p srs_jobs/merged
    python "$CURRENT_DIR/junction_srs.py" merge --transcript-junctions junction_index/transcript_junctions.tsv \
      --count-list srs_jobs/all_count_files.list --out-prefix srs_jobs/merged/srs
    python "$CURRENT_DIR/junction_srs_run.py" --merged-file srs_jobs/merged/srs.transcript_junction_summary.tsv \
      --qc-min-samples 3 --qc-out srs_jobs/merged/srs.transcript_qc.tsv
    ;;
*) echo 'Usage: bash 13_junction_srs_run.sh prepare|submit|merge' >&2; exit 2 ;;
esac
