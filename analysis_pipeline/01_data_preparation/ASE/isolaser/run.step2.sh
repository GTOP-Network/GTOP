#!/bin/bash

#SBATCH --job-name=isoL
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:00:00
#SBATCH --array=1-176
#SBATCH --output=tmp/isol/%x_%A_%a_%N_%j.o
#SBATCH --error=tmp/isol/%x_%A_%a_%N_%j.e

set -euo pipefail

THREADS=$SLURM_CPUS_PER_TASK

workdir="/flashfs1/scratch.global/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/isolaser"
refg="${workdir}/ref/genome.fa"

gtop_gtf="${workdir}/input/enhanced_gtf/GTOP.sorted.gtf.gz"
gtop_fa="${workdir}/input/enhanced_gtf/GTOP.fa"
gtop_db="${workdir}/output/db"

sampleids="${workdir}/input/sampleids.txt"
sampleid=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$sampleids")
pbmm2_bam="/lustre/home/cxue/project/GMTiP-RNA/20251031/long_read/HPC/output/sample_based/${sampleid}/${sampleid}.flnc_mapped.bam"
outdir="${workdir}/output/isolaser/${sampleid}"
chrdir="${outdir}/chrs"
mkdir -p $chrdir

# convert bam file to fastq for re-alignment
function exec_bam2fq() {
    echo "#### start flnc.bam-to-fastq at $(date "+%Y-%m-%d %H:%M:%S") ####"
    samtools fastq --threads ${THREADS} ${pbmm2_bam} > ${outdir}/${sampleid}.pbmm2.fastq && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}
# re-align against newly generated transcriptome reference
function exec_minimap2() {
    echo "#### start minimap2 at $(date "+%Y-%m-%d %H:%M:%S") ####"
    minimap2 -t ${THREADS} -ax map-hifi --MD ${gtop_fa} ${outdir}/${sampleid}.pbmm2.fastq > ${outdir}/${sampleid}.sam && \
    rm -rf ${outdir}/${sampleid}.pbmm2.fastq && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}
# annotating with transcript ids
function filter_and_annotate() {
    echo "#### start filter_and_annotate at $(date "+%Y-%m-%d %H:%M:%S") ####"
    isolaser_annotate -b ${pbmm2_bam} -t ${outdir}/${sampleid}.sam -g ${gtop_gtf} -o ${outdir}/${sampleid}.anno.bam && \
    rm -rf ${outdir}/${sampleid}.sam && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}
# sort & index bam
function bam_sort_index() {
    echo "#### start bam_sort_index at $(date "+%Y-%m-%d %H:%M:%S") ####"
    samtools sort -@ ${THREADS} ${outdir}/${sampleid}.anno.bam -o ${outdir}/${sampleid}.anno.sorted.bam && \
    samtools index -@ ${THREADS} ${outdir}/${sampleid}.anno.sorted.bam && \
    rm -rf ${outdir}/${sampleid}.anno.bam && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}
# Run IsoLASER
function exec_isolaser() {
    echo "#### start isolaser at $(date "+%Y-%m-%d %H:%M:%S") ####"
    isolaser -b ${outdir}/${sampleid}.anno.sorted.bam -o ${outdir}/${sampleid} -t ${gtop_db} -f ${refg} \
    --nodes=${THREADS} --platform=PacBio && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}
function exec_bgzip() {
    sed -n '1p' ${outdir}/${sampleid}.mi_summary.tab > ${outdir}/${sampleid}.mi_summary.tab.header && \
    sed -n '2,$p' ${outdir}/${sampleid}.mi_summary.tab > ${outdir}/${sampleid}.mi_summary.tab.content && \
    sort -k1,1V -k2,2n -k3,3n ${outdir}/${sampleid}.mi_summary.tab.content > ${outdir}/${sampleid}.mi_summary.tab.content.sorted && \
    cat ${outdir}/${sampleid}.mi_summary.tab.header ${outdir}/${sampleid}.mi_summary.tab.content.sorted > ${outdir}/${sampleid}.mi_summary.tab.sorted && \
    bgzip -c -f ${outdir}/${sampleid}.mi_summary.tab.sorted > ${outdir}/${sampleid}.mi_summary.tab.gz && \
    tabix -f -p bed ${outdir}/${sampleid}.mi_summary.tab.gz && \
    rm -rf ${outdir}/${sampleid}.mi_summary.tab.header \
    ${outdir}/${sampleid}.mi_summary.tab.content \
    ${outdir}/${sampleid}.mi_summary.tab.content.sorted \
    ${outdir}/${sampleid}.mi_summary.tab.sorted
}
# obtain the significant allele-specific events (cis-directed splicing events)
function exec_isolaser_filter() {
    echo "#### start isolaser_filter at $(date "+%Y-%m-%d %H:%M:%S") ####"
    isolaser_filter -m ${outdir}/${sampleid}.mi_summary.tab -o ${outdir}/${sampleid}.mi_summary.filtered.tab && \
    echo "#### end at $(date "+%Y-%m-%d %H:%M:%S") ####"
}


# main func
function main() {
    source /lustre/home/lhgong/anaconda3/bin/activate isoquant
    exec_bam2fq
    exec_minimap2
    source /lustre/software/anaconda/anaconda3-2019.10-py37/bin/activate /lustre/home/cxue/.conda/envs/isoLASER
    filter_and_annotate
    source /lustre/home/lhgong/anaconda3/bin/activate isoquant
    bam_sort_index
    source /lustre/software/anaconda/anaconda3-2019.10-py37/bin/activate /lustre/home/cxue/.conda/envs/isoLASER
    exec_isolaser
    exec_isolaser_filter
    source /lustre/home/lhgong/anaconda3/bin/activate base
    module load samtools
    exec_bgzip
}
main
