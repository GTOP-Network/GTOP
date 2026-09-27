#!/bin/bash


main(){

pre_processing_alignment
TR_calling
SNV_calling

}

#Function to processing bam and align reads to human hg38 reference fasta
function pre_processing_alignment(){

    # 1. adapter_trim
    fastp -i  ${fq1} -I  ${fq2} -o  ${out_dir}/${sampleID}_1.fq.gz -O  ${out_dir}/${sampleID}_2.fq.gz --thread ${THREAD_PER_JOB} -j ${out_dir}/${sampleID}.json -h ${out_dir}/${sampleID}.html

    # 2. Mapping
    bwa mem -t 16 -R '@RG\tID:foo_lane\tPL:illumina\tLB:library\tSM:${sampleID}' ${ref_index} ${fq_dir}/${sampleID}_1.fq.gz ${fq_dir}/${sampleID}_2.fq.gz |samtools view -@ 16 -S -b - > ${out_dir}/${sampleID}.bam
    samtools sort -@ 16 -O bam -o ${out_dir}/${sampleID}.sorted.bam ${out_dir}/${sampleID}.bam -T ${out_dir}/tmp/${sampleID}

    # 3. Mark_Duplicates
    java -Djava.io.tmpdir=${out_dir}/tmp/${sampleID} -jar picard.jar MarkDuplicates I=${bam_dir}/${sampleID}.sorted.bam O=${out_dir}/${sampleID}.sorted.markdup.bam M=${out_dir}/${sampleID}.markdup_metrics.txt TMP_DIR=${out_dir}/tmp/${sampleID}
    samtools index -@ 32 ${out_dir}/${sampleID}.sorted.markdup.bam 

    # 4. BQSR
    gatk BaseRecalibrator \
            -R reference.fasta \
            -I ${bam_dir}/${sampleID}.sorted.markdup.bam \
            --known-sites ${REF_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz \
            --known-sites ${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            --known-sites ${REF_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz \
            --known-sites ${REF_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
            --known-sites ${REF_DIR}/dbsnp_156.hg38.vcf.gz \
            -O ${out_dir}/${sampleID}.recal_data.table

    gatk ApplyBQSR \
            -R reference.fasta \
            -I ${bam_dir}/${sampleID}.sorted.markdup.bam \
            --bqsr-recal-file ${out_dir}/${sampleID}.recal_data.table \
            --create-output-bam-index false \
            -O ${out_dir}/${sampleID}.sorted.markdup.realign.BQSR.bam

}

# function to genotype the tandem repeat
function TR_calling(){

############################################################################################################
#    1. each tool to genotype TRs and generate raw calls in VCF format, with a single VCF file per sample  #
############################################################################################################

    # 1.1 GangSTR 
    GangSTR --bam ${sampleID}.sorted.markdup.realign.BQSR.bam \
            --ref reference.fasta \
            --regions repeat_catalog_v1.hg38.1_to_1000bp_motifs.GangSTR.bed \
            --verbose --grid-threshold 250 \
            --out ${sampleID}.GangSTR
            
    # 1.2 HipSTR
    HipSTR --bams ${sampleID}.sorted.markdup.realign.BQSR.bam \
            --fasta reference.fasta \
            --regions repeat_catalog_v1.hg38.1_to_1000bp_motifs.HipSTR.bed \
            --min-reads 20 --max-reads 2000000  --output-filters \
            --str-vcf ${sampleID}.HipSTR.vcf.gz

    # 1.3 ExpansionHunter
    ExpansionHunter --reads ${sampleID}.sorted.markdup.realign.BQSR.bam \
            --reference reference.fasta \
            --variant-catalog repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.json \
            --analysis-mode seeking --threads 16 \
            --output-prefix ${sampleID}.ExpansionHunter


############################################################################################################
#    2. merge individual output VCF files and dump low-quality sites of each tools                         # 
############################################################################################################

    # 2.1 merge_GangSTR_and_Filter_vcf
    mergeSTR --vcfs $FILE --out ${CUR_DIR}/mergeSTR/GTOP_merge_GangSTR --vcftype gangstr 
    dumpSTR --gangstr-filter-spanbound-only \
            --gangstr-filter-badCI \
            --gangstr-max-call-DP 1000 \
            --gangstr-min-call-DP 20 \
            --min-locus-callrate 0.85 \
            --gangstr-min-call-Q 0.95 \
            --min-locus-hwep 0.000001 \
            --filter-regions UCSC_hg38_genomic_segmental_duplications.sorted.bed.gz \
            --filter-regions-names SEGDUP \
            --vcftype gangstr \
            --zip \
            --vcf GTOP_merge_GangSTR.vcf.gz \
            --out GTOP_merge_GangSTR_dump

    bcftools filter --threads 4 -i 'FILTER=="PASS"' --no-version -Oz -o GTOP_merge_GangSTR_dump_filted.vcf.gz GTOP_merge_GangSTR_dump.vcf.gz 

    # 2.2 merge_HipSTR_and_Filter_vcf
    mergeSTR --vcfs $FILE --out GTOP_merge_HipSTR --vcftype hipstr 

    dumpSTR --hipstr-max-call-DP 1000 \
            --hipstr-min-call-DP 20 \
            --hipstr-min-call-Q 0.95 \
            --hipstr-max-call-flank-indel 0.15 \
            --hipstr-max-call-stutter 0.15 \
            --min-locus-callrate 0.85 \
            --min-locus-hwep 0.000001 \
            --filter-regions UCSC_hg38_genomic_segmental_duplications.sorted.bed.gz \
            --filter-regions-names SEGDUP \
            --vcftype hipstr \
            --zip \
            --vcf GTOP_merge_HipSTR.vcf.gz \
            --out GTOP_merge_HipSTR_dump
    bcftools filter --threads 4 -i 'FILTER=="PASS"' --no-version -Oz -o ${CUR_DIR}/dumpSTR/GTOP_merge_HipSTR_dump_filted.vcf.gz ${CUR_DIR}/dumpSTR/GTOP_merge_HipSTR_dump.vcf.gz 

    # 2.3 merge_ExpansionHunter_and_Filter_vcf
    bcftools filter ${sampleID}.ExpansionHunter.vcf --threads 4 -i 'FILTER=="PASS"' --no-version -O z  -o ${sampleID}.ExpansionHunter_filter.vcf.gz
    mergeSTR --vcfs $FILE --out GTOP_merge_ExpansionHunter --vcftype eh 
    dumpSTR --eh-min-call-LC 20 \
            --min-locus-callrate 0.85 \
            --min-locus-hwep 0.000001 \
            --filter-regions UCSC_hg38_genomic_segmental_duplications.sorted.bed.gz \
            --filter-regions-names SEGDUP \
            --vcftype eh \
            --zip \
            --vcf GTOP_merge_ExpansionHunter.vcf.gz \
            --out GTOP_merge_ExpansionHunter_dump

############################################################################################################
#    3. Ensemble genotyping                                                                                #
#  EnsembleTR takes VCF files from multiple TR genotypers as input and outputs a merged consensus callset. # 
############################################################################################################
    ensembletr --out ${CUR_DIR}/GTOP_SRS_TR.raw.vcf --ref reference.fasta --vcfs GTOP_merge_ExpansionHunter_dump.vcf.gz,GTOP_merge_HipSTR_dump.vcf.gz,GTOP_merge_GangSTR_dump.vcf.gz --exclude-single TURE

}


# function to genotype the SNVs and small Indels(< 50bp)
function SNV_calling(){

    # 1. make gVCF
    gatk HaplotypeCaller -R reference.fasta \
        -I ${bam_dir}/${sampleID}.sorted.markdup.realign.BQSR.bam \
        -O ${gVCF_dir}/${sampleID}.g.vcf.gz \
        -ERC GVCF
    
    # 2. Merge multi-sample VCF
    gatk --java-options \"-Xmx10g -XX:ParallelGCThreads=2 -Dsamjdk.compression_level=5 \" \
        CombineGVCFs \
        -R ${ref_dir}/reference.fasta \
        ${line} \
        --intervals chr${chr} \
        --tmp-dir ${out_dir}/tmp/chr${chr} \
        -O ${gVCF_dir}/cohort.chr${chr}.g.vcf.gz
        
    ### 3. joint_genotype_chrs

    gatk GenotypeGVCFs \
        -R ${ref_dir}/reference.fasta \
        -V ${gVCF_dir}/cohort.chr${chr}.g.vcf.gz \
        -O ${gVCF_dir}/cohort.chr${chr}.vcf.gz \

    bcftools concat ${line} -Oz -o ${out_gVCF_dirdir}/cohort.chrs.vcf.gz && tabix -p vcf  ${gVCF_dir}/cohort.chrs.vcf.gz 

    ### 4. Variant filter
    gatk  --java-options \"-Djava.io.tmpdir=`pwd`/tmp -Xms4G -Xmx4G -XX:ParallelGCThreads=2\" \
        VariantRecalibrator \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 \
    -tranche 99.8 \
    -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
    -tranche 95.0 -tranche 94.0 \
    -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -R ${ref_dir}/reference.fasta \
    -V ${out_dir}/cohort.chrs.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${ref_dir}/hapmap_3.3.hg38.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 ${ref_dir}/1000G_omni2.5.hg38.vcf.gz\
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${ref_dir}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${ref_dir}/dbsnp_156.hg38.vcf.gz  \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
    -mode SNP \
    -O ${out_dir}/cohort.SNP.recal \
    --tranches-file ${out_dir}/cohort.SNP.tranches \
    --rscript-file ${out_dir}/cohort.SNP.plots.R 

    gatk ApplyVQSR \
    -R ${ref_dir}/reference.fasta \
    -V ${out_dir}/cohort.chrs.vcf.gz \
    -O ${out_dir}/cohort.VQSR.SNP.99.8.chrs.vcf.gz \
    -ts-filter-level 99.8 \
    --tranches-file ${out_dir}/cohort.SNP.tranches \
    --recal-file ${out_dir}/cohort.SNP.recal \
    -mode SNP \

    gatk --java-options \"-Djava.io.tmpdir=`pwd`/tmp  -Xms4G -Xmx4G -XX:ParallelGCThreads=2 \" \
    VariantRecalibrator \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 \
    -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
    -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
    -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -R   ${ref_dir}/reference.fasta    \
    -V   ${out_dir}/cohort.chrs.vcf.gz     \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 ${ref_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz  \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${ref_dir}/dbsnp_156.hg38.vcf.gz \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
    -mode INDEL \
    -O ${out_dir}/cohort.INDEL.recal \
    --tranches-file ${out_dir}/cohort.INDEL.tranches    \
    --rscript-file  ${out_dir}/cohort.INDEL.plots.R 

    gatk ApplyVQSR \
    -R ${ref_dir}/reference.fasta \
    -V ${out_dir}/cohort.VQSR.SNP.99.8.chrs.vcf.gz \
    -O ${out_dir}/GTOP_SRS.VQSR.SNP.99.8.INDEL.99.95.chrs.vcf.gz \
    -ts-filter-level 99.95 \
    --tranches-file ${out_dir}/cohort.INDEL.tranches \
    --recal-file ${out_dir}/cohort.INDEL.recal \
    -mode INDEL \

    bcftools annotate -a GCF_000001405.40.vcf.gz -c ID    ${out_dir}/GTOP_SRS.VQSR.SNP.99.8.INDEL.99.95.chrs.vcf.gz   -Oz -o ${out_dir}/GTOP_SRS.VQSR.SNP.99.8.INDEL.99.95.dbSNP156.chrs.vcf.gz && tabix -p vcf  ${out_dir}/GTOP_SRS.VQSR.SNP.99.8.INDEL.99.95.dbSNP156.chrs.vcf.gz
    bcftools view -f PASS  ${out_dir}/GTOP_SRS.VQSR.SNP.99.5.INDEL.99.dbSNP156.ID.chrs.vcf.gz   -Oz -o  ${out_dir}/GTOP_SRS.VQSR.SNP.99.5.INDEL.99.dbSNP156.ID.PASS.chrs.vcf.gz  && tabix -p vcf  ${out_dir}/GTOP_SRS.VQSR.SNP.99.5.INDEL.99.dbSNP156.ID.PASS.chrs.vcf.gz  

}

main