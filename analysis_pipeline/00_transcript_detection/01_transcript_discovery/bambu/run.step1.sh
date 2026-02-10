#!/bin/bash

workdir="/lustre/home/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/bambu"

# bambu
function run_bambu_preparing() {
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] start bambu-preparing workflow ..."
    source /lustre/software/proxy/set_proxy.sh
    /lustre/home/lhgong/anaconda3/envs/longrw_bambu3/bin/Rscript ${workdir}/preparing.R && \
    echo "[$(date "+%Y-%m-%d %H:%M:%S")] done! "
}


# main func
function main() {
    source /lustre/home/lhgong/anaconda3/bin/activate longrw_bambu3
    run_bambu_preparing
}
main
