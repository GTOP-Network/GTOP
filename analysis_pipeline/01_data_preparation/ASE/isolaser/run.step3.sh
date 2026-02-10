#!/bin/bash

workdir="/flashfs1/scratch.global/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/isolaser"
outdir="${workdir}/output/isolaser"
sampleids="${workdir}/input/sampleids.txt"
mapfile -t sampleid_arr < $sampleids
echo -n > ${outdir}/fofn.txt
for sid in "${sampleid_arr[@]}"; do
    bamf="${outdir}/${sid}/${sid}.anno.sorted.bam"
    mif="${outdir}/${sid}/${sid}.mi_summary.tab.gz"
    vcf="${outdir}/${sid}/${sid}.gvcf"
    printf "${{bamf}}\t${{mif}}\t${{vcf}}\n" >> ${outdir}/fofn.txt
done
cat ${outdir}/fofn.txt
