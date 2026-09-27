#!/usr/bin/env bash
set -euo pipefail
: "${GTOP_PROJECT_DIR:?Set GTOP_PROJECT_DIR using config.example.sh}"
MERGE_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
OUT_DIR="$GTOP_PROJECT_DIR/output/assembly/LRS/tx_merge/merged_tools/sqanti3_filter_1"
REF_GTF="${HG38_GTF:-${GTOP_REF_DIR:-$GTOP_PROJECT_DIR/reference}/gencode.v47.annotation.gtf}"
mkdir -p "$OUT_DIR/AF_filter"
cd "$OUT_DIR/AF_filter"
case "${1:-}" in
prepare)
    python "$MERGE_DIR/alt_novel_first_exon_qc.py" --log-file prepare.log prepare \
      --ref-gtf "$REF_GTF" --lrs-gtf "$OUT_DIR/filtered.gtf" \
      --sqanti "$OUT_DIR/filtered.RulesFilter_result_classification.txt" --out-prefix novel_first_exon
    ;;
submit)
    python "$MERGE_DIR/alt_first_exon_run.py"
    ;;
merge)
    python "$MERGE_DIR/alt_novel_first_exon_qc.py" --log-file final.merge.log merge \
      --transcripts novel_first_exon.transcripts.tsv --regions novel_first_exon.regions.tsv \
      --scan-results scan_results/GTOP*/out.sample_region_mismatch.tsv --out-prefix final_novel_first_exon
    ;;
*) echo 'Usage: bash 12_alt_first_exon_run.sh prepare|submit|merge' >&2; exit 2 ;;
esac
