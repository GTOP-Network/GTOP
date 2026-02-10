# Adipose 脂肪组织
# Muscle 肌肉
# Skin 皮肤：人体是最大的器官
# Adrenal_Gland  肾上腺
# Gallbladder 胆囊
# Liver 肝
# Pancreas_Body 胰体
# Pancreas_Head 胰头
# Pancreas_Tail 胰尾
# Spleen 脾脏：体内最大的淋巴器官，主要由淋巴组织构成，富含血管和血窦，具有免疫应答、滤血、储血和造血等功能，属于免疫系统的一部分
# Whole_Blood 全血

main(){
    copy_data
}


copy_data(){

    TISSUE="Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood"

    ## create folders
    for xQTL_type in SNV_eQTL SNV_juQTL SNV_tuQTL SV_eQTL SV_juQTL SV_tuQTL TR_eQTL TR_juQTL TR_tuQTL
    do
        mkdir -p data/tensorqtl/${xQTL_type}/{clean_xGene,nominal,nominal_slim,permutation}
    done

    ## SNV-eQTL  ------------------------------------------------------
    ### HPC SZBL ZT
    ## /lustre/home/xdzou/2024-10-21-GTBMap/2025-02-13-expression-quant/output/QTL_mapping_131/covariates/PC_0_no_Batch
    CURRENT_VERSION=QTL_mapping_131
    SNV_eQTL_path=/lustre/home/xdzou/2024-10-21-GTBMap/2025-02-13-expression-quant/output/${CURRENT_VERSION}/covariates/PC_0_no_Batch
    SNVeQTLPATH=/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/SNV_eQTL

    rsync -auxvLP ${SNV_eQTL_path}/cis_QTL/text_format/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.cis_eQTL.all_pairs.txt.gz zhangt@10.6.109.182:${SNVeQTLPATH}/nominal
    rsync -auxvLP ${SNV_eQTL_path}/permutation/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.perm.cis_qtl.txt.gz zhangt@10.6.109.182:${SNVeQTLPATH}/permutation
    rsync -auxvLP ${SNV_eQTL_path}/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.covariates.txt zhangt@10.6.109.182:${SNVeQTLPATH}/covariate
    rsync -auxvLP /lustre/home/xdzou/2024-10-21-GTBMap/2025-02-13-expression-quant/output/phenotype/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.phenotype.bed zhangt@10.6.109.182:${SNVeQTLPATH}/phenotype
    rsync -auxvLP /lustre/home/xdzou/2024-10-21-GTBMap/2025-02-13-expression-quant/input/Donor_age_sex_20250828.csv zhangt@10.6.109.182:/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl

    ## SV-eQTL  ------------------------------------------------------
    ### HPC SZBL ZT
    SVeQTLPATH=/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/SV_eQTL

    rsync -auxvLP ${SNV_eQTL_path}/SV_eQTL/cis_QTL/text_format/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.cis_eQTL.all_pairs.txt.gz zhangt@10.6.109.182:${SVeQTLPATH}/nominal
    rsync -auxvLP ${SNV_eQTL_path}/SV_eQTL/permutation/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.perm.cis_qtl.txt.gz zhangt@10.6.109.182:${SVeQTLPATH}/permutation

    for i in ${TISSUE//,/ }
    do
    mv ${SNVeQTLPATH}/phenotype/${i}.phenotype.bed ${SNVeQTLPATH}/phenotype/${i}.bed && bgzip ${SNVeQTLPATH}/phenotype/${i}.bed
    mv ${SNVeQTLPATH}/covariate/${i}.covariates.txt ${SNVeQTLPATH}/covariate/${i}.txt
    mv ${SNVeQTLPATH}/nominal/${i}.cis_eQTL.all_pairs.txt.gz ${SNVeQTLPATH}/nominal/${i}.txt.gz
    mv ${SNVeQTLPATH}/permutation/${i}.perm.cis_qtl.txt.gz ${SNVeQTLPATH}/permutation/${i}.txt.gz
    mv ${SVeQTLPATH}/nominal/${i}.cis_eQTL.all_pairs.txt.gz ${SVeQTLPATH}/nominal/${i}.txt.gz
    mv ${SVeQTLPATH}/permutation/${i}.perm.cis_qtl.txt.gz ${SVeQTLPATH}/permutation/${i}.txt.gz
    done


    ## sQTL: SNV SV   -------------------------------------------------------
    for variant_type in SNV SV
    do
    for sQTL in juQTL tuQTL
    do
    if [[ $variant_type == "SNV" ]];then
        DATAPATH=/media/bora_A/wangyn/2024-11-15-GTBMap-sqtl/2025-10-07-all-sqtl/output_old/${sQTL,,}/QTL_mapping_SNV_Indel
    elif [[ $variant_type == "SV" ]];then
        DATAPATH=/media/bora_A/wangyn/2024-11-15-GTBMap-sqtl/2025-10-07-all-sqtl/output_old/${sQTL,,}/QTL_mapping_SV_only
    fi
    sQTLPATH=/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/${variant_type}_${sQTL}
    for tissue in ${TISSUE//,/ }
    do
    rm ${sQTLPATH}/permutation/${tissue}.txt.gz
    ln -s ${DATAPATH}/cis_QTL/text_format/${tissue}.cis_eQTL.all_pairs.txt.gz ${sQTLPATH}/nominal/${tissue}.txt.gz
    ln -s ${DATAPATH}/permutation/${tissue}.cis_qtl.txt.gz ${sQTLPATH}/permutation/${tissue}.txt.gz
    done
    done
    done

    tuQTLPATH=/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/SNV_tuQTL
    juQTLPATH=/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/SNV_juQTL

    for tissue in ${TISSUE//,/ }
    do
    ln -s /media/bora_A/wangyn/2024-11-15-GTBMap-sqtl/2025-10-07-all-sqtl/input/ori_bed/tuqtl/${tissue}.filter.bed.gz ${tuQTLPATH}/phenotype/${tissue}.bed.gz
    ln -s /media/bora_A/wangyn/2024-11-15-GTBMap-sqtl/2025-02-16-Splicing/leafcutter2/splittissue_0906filter/${tissue}_regtools/leafcutter2.filter.leafcutter.bed.gz ${juQTLPATH}/phenotype/${tissue}.bed.gz
    ln -s /media/bora_A/wangyn/2024-11-15-GTBMap-sqtl/2025-02-16-Splicing/leafcutter2_filter_covariates_0906/${tissue}_regtools.GPC0.covariates.txt ${juQTLPATH}/covariate/${tissue}.txt
    done

    rsync -auxvLP /flashfs1/scratch.global/ynwang/ynwang/2024-11-15-GTBMap/2025-10-07-all-sqtl/output/covariates/tuqtl/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.GPC0.covariates.txt zhangt@10.6.109.182:/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/SNV_tuQTL/covariate
    
    ls /media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/SNV_tuQTL/covariate/*.txt | while read line;do mv ${line} ${line/.GPC0.covariates/};done


    ## TR    -------------------------------------------------------
    DATAPATH=/media/london_B/lixing/2024-08-29-TR-AsianGTEX/2025-03-14-LRS-TR-QTL/output/cis_eQTL
    TReQTLPATH=/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/TR_eQTL
    ls ${DATAPATH}/text_format/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.cis_eQTL.all_pairs.txt.gz | while read line;do ln -s $line ${TReQTLPATH}/nominal/`basename ${line/.cis_eQTL.all_pairs/}`;done
    cp ${DATAPATH}/permutation/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.cis_qtl.txt ${TReQTLPATH}/permutation
    ls ${TReQTLPATH}/permutation/*.txt | while read line;do mv ${line} ${line/.cis_qtl/}; bgzip ${line/.cis_qtl/};done


    DATAPATH=/media/london_B/lixing/2024-08-29-TR-AsianGTEX/2025-03-14-LRS-TR-QTL/output/cis_sQTL
    TRjuQTLPATH=/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/TR_juQTL
    ls ${DATAPATH}/text_format/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.cis_sQTL.all_pairs.txt.gz | while read line;do ln -s $line ${TRjuQTLPATH}/nominal/`basename ${line/.cis_sQTL.all_pairs/}`;done
    cp ${DATAPATH}/permutation/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.cis_qtl.txt ${TRjuQTLPATH}/permutation
    ls ${TRjuQTLPATH}/permutation/*.txt | while read line;do mv ${line} ${line/.cis_qtl/}; bgzip ${line/.cis_qtl/};done

    DATAPATH=/lustre/home/xlma/lixing/2024-11-14-GTBMap/2025-03-31-TR-xQTL/output/cis_tuQTL
    TRtuQTLPATH=/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/TR_tuQTL
    rsync -auxvLP ${DATAPATH}/permutation/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.cis_qtl.txt zhangt@10.6.109.182:${TRtuQTLPATH}/permutation
    ls ${TRtuQTLPATH}/permutation/*.txt | while read line;do mv ${line} ${line/.cis_qtl/}; bgzip ${line/.cis_qtl/};done
}


sample_size_file(){
    mkdir -p data/tensorqtl

    rm data/tensorqtl/SNV_eQTL/selected_samples.txt
    ls data/tensorqtl/SNV_eQTL/phenotype | while read line
    do
        zcat data/tensorqtl/SNV_eQTL/phenotype/${line} | head -n 1  | echo -e "${line/.bed.gz/}\t"`awk '{print NF-4}'` >> data/tensorqtl/SNV_eQTL/selected_samples.txt
    done

    ls data/tensorqtl/SNV_juQTL/phenotype | while read line
    do
        zcat data/tensorqtl/SNV_juQTL/phenotype/${line} | head -n 1  | echo -e "${line/.bed.gz/}\t"`awk '{print NF-4}'` >> data/tensorqtl/SNV_juQTL/selected_samples.txt
    done
}


cis_xQTL_slim(){
    TISSUE="Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood"

    for xQTL_type in SNV_eQTL SNV_juQTL SNV_tuQTL
    do
    parallel -j4 Rscript bin/tensorqtl/cis_slim.R ${xQTL_type},{} ::: ${TISSUE//,/ }
    done

    for xQTL_type in SV_eQTL SV_juQTL SV_tuQTL
    do
    parallel -j4 Rscript bin/tensorqtl/cis_slim.R ${xQTL_type},{} ::: ${TISSUE//,/ }
    done

    for xQTL_type in TR_eQTL TR_juQTL TR_tuQTL
    do
    parallel -j4 Rscript bin/tensorqtl/cis_slim.R ${xQTL_type},{} ::: ${TISSUE//,/ }
    done
}


merge_variants(){
    bcftools query -l ${snp_file} > $CUR_DIR/output/GWAS_LD/SNP_samples.txt
    bcftools query -l ${TR_file} > $CUR_DIR/output/GWAS_LD/TR_samples.txt
    bcftools query -l ${SV_file} > $CUR_DIR/output/GWAS_LD/SV_samples.txt
    diff $CUR_DIR/output/GWAS_LD/SNP_samples.txt $CUR_DIR/output/GWAS_LD/SV_samples.txt
    diff $CUR_DIR/output/GWAS_LD/SNP_samples.txt $CUR_DIR/output/GWAS_LD/TR_samples.txt

    python /media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/bin/tensorqtl/all_variant_susie_genotype.py
}


fine_mapping(){
    TISSUE="Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood"

    for xQTL_type in SNV_juQTL SNV_tuQTL SV_eQTL SV_juQTL SV_tuQTL TR_eQTL TR_juQTL TR_tuQTL
    do
    VARTYPE=`echo $xQTL_type | cut -f1 -d"_"`
    QTLTYPE=`echo $xQTL_type | cut -f2 -d"_"`
    rm -r data/tensorqtl/${xQTL_type}/susie
    mkdir -p data/tensorqtl/${xQTL_type}/susie/phenotype

    echo ${TISSUE//,/ } | tr " " "\n" | while read line; do cat data/tensorqtl/${xQTL_type}/clean_xGene/${line}.txt | grep -v "lead_variant" | cut -f 1 | sort | uniq -c | awk '$1>=1 {print $2}' > data/tensorqtl/${xQTL_type}/susie/phenotype/${line}.txt;done
    echo ${TISSUE//,/ } | tr " " "\n" | xargs -I {} echo "zcat data/tensorqtl/SNV_${QTLTYPE}/phenotype/{}.bed.gz | awk 'NR==FNR{a[\$1]=\$1;next}{if(a[\$4]||\$4==\"phenotype_id\"||\$4==\"ID\"){print \$0}}' data/tensorqtl/${xQTL_type}/susie/phenotype/{}.txt - | awk 'BEGIN{OFS=\"\\t\"} NR==1; NR>1 {print | \"sort -k1,1 -k2,2n\"}' | sed \"s/^chr//\" | bgzip > data/tensorqtl/${xQTL_type}/susie/phenotype/{}.bed.gz" | parallel -j 4 {}
    done

    for xQTL_type in TR_eQTL TR_juQTL TR_tuQTL
    do
    VARTYPE=`echo $xQTL_type | cut -f1 -d"_"`
    QTLTYPE=`echo $xQTL_type | cut -f2 -d"_"`
    rm -r data/tensorqtl/${xQTL_type}/susie
    mkdir -p data/tensorqtl/${xQTL_type}/susie/phenotype
    echo ${TISSUE//,/ } | tr " " "\n" | while read line; do cat data/tensorqtl/${xQTL_type}/clean_xGene/${line}.txt | grep -v "lead_variant" | cut -f 1 | sort | uniq -c | awk '$1>=1 {print $2}' > data/tensorqtl/${xQTL_type}/susie/phenotype/${line}.txt;done
    echo ${TISSUE//,/ } | tr " " "\n" | xargs -I {} echo "zcat data/tensorqtl/SNV_${QTLTYPE}/phenotype/{}.bed.gz | awk 'NR==FNR{a[\$1]=\$1;next}{if(a[\$4]||\$4==\"phenotype_id\"||\$4==\"ID\"){print \$0}}' data/tensorqtl/${xQTL_type}/susie/phenotype/{}.txt - | awk 'BEGIN{OFS=\"\\t\"} NR==1; NR>1 {print | \"sort -k1,1 -k2,2n\"}' | bgzip > data/tensorqtl/${xQTL_type}/susie/phenotype/{}.bed.gz" | parallel -j 4 {}
    done

    for xQTL_type in SV_eQTL SV_juQTL SV_tuQTL
    do
    VARTYPE=`echo $xQTL_type | cut -f1 -d"_"`
    QTLTYPE=`echo $xQTL_type | cut -f2 -d"_"`
    echo ${TISSUE//,/ } | tr " " "\n" | parallel -j 3 /media/dubai/home/dingruofan/anaconda3/envs/work/bin/python bin/tensorqtl/susie_map.py -i data/tensorqtl/${xQTL_type} -i2 data/tensorqtl/SNV_${QTLTYPE} -g data/genotype/GMTiP_${VARTYPE} -t {}
    done

    for xQTL_type in TR_tuQTL
    do
    VARTYPE=`echo $xQTL_type | cut -f1 -d"_"`
    QTLTYPE=`echo $xQTL_type | cut -f2 -d"_"`
    echo ${TISSUE//,/ } | tr " " "\n" | parallel -j 3 /media/dubai/home/dingruofan/anaconda3/envs/work/bin/python bin/tensorqtl/TR_susie_map.py -i data/tensorqtl/${xQTL_type} -i2 data/tensorqtl/SNV_${QTLTYPE} -t {}
    done
}


copy_susie(){
mkdir -p data/tensorqtl/{SNV_eQTL,SNV_juQTL,SNV_tuQTL,SV_eQTL,SV_juQTL,SV_tuQTL,TR_eQTL,TR_juQTL,TR_tuQTL}/{susie/phenotype,covariate}
mkdir -p data/genotype

rsync -auxvLP zhangt@10.6.109.182:/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/bin/tensorqtl /lustre/home/tzhang/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/bin

rsync -auxvLP zhangt@10.6.109.182:/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/genotype/GMTiP_TR_* /lustre/home/tzhang/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/genotype

# SNV_eQTL SNV_juQTL SNV_tuQTL SV_eQTL SV_juQTL SV_tuQTL #
for xQTL_type in SNV_eQTL SNV_juQTL SNV_tuQTL
do
    rsync -auxvLP zhangt@10.6.109.182:/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/${xQTL_type}/covariate/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.txt /lustre/home/tzhang/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/${xQTL_type}/covariate
done

for xQTL_type in TR_eQTL TR_juQTL TR_tuQTL
do
    rsync -auxvLP zhangt@10.6.109.182:/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/${xQTL_type}/susie/phenotype/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.bed.gz /lustre/home/tzhang/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/${xQTL_type}/susie/phenotype
done

for xQTL_type in TR_eQTL TR_juQTL TR_tuQTL
do
    rsync -auxvLP /lustre/home/tzhang/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/${xQTL_type}/susie/{Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood}.txt zhangt@10.6.109.182:/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/tensorqtl/${xQTL_type}/susie
done

}


merge_fine_mapping(){
    for xQTL_type in SNV_eQTL
    do
    VARTYPE=`echo $xQTL_type | cut -f1 -d"_"`
    QTLTYPE=`echo $xQTL_type | cut -f2 -d"_"`
    rm -r data/tensorqtl/${xQTL_type}/combined_susie
    mkdir -p data/tensorqtl/${xQTL_type}/combined_susie/phenotype

    echo ${TISSUE//,/ } | tr " " "\n" | while read line; do cat data/tensorqtl/{SNV_${QTLTYPE},TR_${QTLTYPE},SV_${QTLTYPE}}/clean_xGene/${line}.txt | grep -v "lead_variant" | cut -f 1 | sort | uniq -c | awk '$1>=1 {print $2}' > data/tensorqtl/${xQTL_type}/combined_susie/phenotype/${line}.txt;done

    echo ${TISSUE//,/ } | tr " " "\n" | xargs -I {} echo "zcat data/tensorqtl/${xQTL_type}/phenotype/{}.bed.gz | awk 'NR==FNR{a[\$1]=\$1;next}{if(a[\$4]||\$4==\"phenotype_id\"||\$4==\"ID\"){print \$0}}' data/tensorqtl/${xQTL_type}/combined_susie/phenotype/{}.txt - | bgzip > data/tensorqtl/${xQTL_type}/combined_susie/phenotype/{}.bed.gz" | parallel -j 4 {}

    echo ${TISSUE//,/ } | tr " " "\n" | parallel -j 3 /media/dubai/home/dingruofan/anaconda3/envs/work/bin/python bin/tensorqtl/all_variant_susie_map.py -i data/tensorqtl/${xQTL_type} -t {}
    done
}


merge_SNV_SV_fine_mapping(){
    # /media/bora_A/zhangt/miniconda3/bin/python /media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/bin/tensorqtl/SNV_SV_susie_genotype.py

    TISSUE="Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood"

    for xQTL_type in SNV_eQTL SNV_juQTL SNV_tuQTL
    do
    VARTYPE=`echo $xQTL_type | cut -f1 -d"_"`
    QTLTYPE=`echo $xQTL_type | cut -f2 -d"_"`
    rm -r data/tensorqtl/${xQTL_type}/combined_SNV_SV_susie
    mkdir -p data/tensorqtl/${xQTL_type}/combined_SNV_SV_susie/phenotype

    echo ${TISSUE//,/ } | tr " " "\n" | while read line; do cat data/tensorqtl/{SNV_${QTLTYPE},SV_${QTLTYPE}}/clean_xGene/${line}.txt | grep -v "lead_variant" | cut -f 1 | sort | uniq -c | awk '$1>=1 {print $2}' > data/tensorqtl/${xQTL_type}/combined_SNV_SV_susie/phenotype/${line}.txt;done

    echo ${TISSUE//,/ } | tr " " "\n" | xargs -I {} echo "zcat data/tensorqtl/${xQTL_type}/phenotype/{}.bed.gz | awk 'NR==FNR{a[\$1]=\$1;next}{if(a[\$4]||\$4==\"phenotype_id\"||\$4==\"ID\"){print \$0}}' data/tensorqtl/${xQTL_type}/combined_SNV_SV_susie/phenotype/{}.txt - | awk 'BEGIN{OFS=\"\\t\"} NR==1; NR>1 {print | \"sort -k1,1 -k2,2n\"}' | bgzip > data/tensorqtl/${xQTL_type}/combined_SNV_SV_susie/phenotype/{}.bed.gz" | parallel -j 4 {}

    echo ${TISSUE//,/ } | tr " " "\n" | parallel -j 3 /media/dubai/home/dingruofan/anaconda3/envs/work/bin/python bin/tensorqtl/SNV_SV_susie_map.py -i data/tensorqtl/${xQTL_type} -t {}
    done

    # /media/bora_A/zhangt/miniconda3/bin/python /media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/bin/tensorqtl/SNV_TR_susie_genotype.py

    # TISSUE="Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood"

    for xQTL_type in SNV_eQTL SNV_juQTL SNV_tuQTL
    do
    VARTYPE=`echo $xQTL_type | cut -f1 -d"_"`
    QTLTYPE=`echo $xQTL_type | cut -f2 -d"_"`
    rm -r data/tensorqtl/${xQTL_type}/combined_SNV_TR_susie
    mkdir -p data/tensorqtl/${xQTL_type}/combined_SNV_TR_susie/phenotype

    echo ${TISSUE//,/ } | tr " " "\n" | while read line; do cat data/tensorqtl/{SNV_${QTLTYPE},TR_${QTLTYPE}}/clean_xGene/${line}.txt | grep -v "lead_variant" | cut -f 1 | sort | uniq -c | awk '$1>=1 {print $2}' > data/tensorqtl/${xQTL_type}/combined_SNV_TR_susie/phenotype/${line}.txt;done

    echo ${TISSUE//,/ } | tr " " "\n" | xargs -I {} echo "zcat data/tensorqtl/${xQTL_type}/phenotype/{}.bed.gz | awk 'NR==FNR{a[\$1]=\$1;next}{if(a[\$4]||\$4==\"phenotype_id\"||\$4==\"ID\"){print \$0}}' data/tensorqtl/${xQTL_type}/combined_SNV_TR_susie/phenotype/{}.txt - | awk 'BEGIN{OFS=\"\\t\"} NR==1; NR>1 {print | \"sort -k1,1 -k2,2n\"}' | bgzip > data/tensorqtl/${xQTL_type}/combined_SNV_TR_susie/phenotype/{}.bed.gz" | parallel -j 4 {}

    echo ${TISSUE//,/ } | tr " " "\n" | parallel -j 3 /media/dubai/home/dingruofan/anaconda3/envs/work/bin/python bin/tensorqtl/SNV_TR_susie_map.py -i data/tensorqtl/${xQTL_type} -t {}
    done
}


merge_SNV_modifyTR_fine_mapping(){
    /media/bora_A/zhangt/miniconda3/bin/python /media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/bin/tensorqtl/SNV_modifyTR_susie_genotype.py

    TISSUE="Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood"

    for xQTL_type in SNV_tuQTL
    do
    VARTYPE=`echo $xQTL_type | cut -f1 -d"_"`
    QTLTYPE=`echo $xQTL_type | cut -f2 -d"_"`
    rm -r data/tensorqtl/${xQTL_type}/combined_SNV_modifyTR_susie
    mkdir -p data/tensorqtl/${xQTL_type}/combined_SNV_modifyTR_susie/phenotype

    echo ${TISSUE//,/ } | tr " " "\n" | while read line; do cat data/tensorqtl/{SNV_${QTLTYPE},TR_${QTLTYPE}}/clean_xGene/${line}.txt | grep -v "lead_variant" | cut -f 1 | sort | uniq -c | awk '$1>=1 {print $2}' > data/tensorqtl/${xQTL_type}/combined_SNV_modifyTR_susie/phenotype/${line}.txt;done

    echo ${TISSUE//,/ } | tr " " "\n" | xargs -I {} echo "zcat data/tensorqtl/${xQTL_type}/phenotype/{}.bed.gz | awk 'NR==FNR{a[\$1]=\$1;next}{if(a[\$4]||\$4==\"phenotype_id\"||\$4==\"ID\"){print \$0}}' data/tensorqtl/${xQTL_type}/combined_SNV_modifyTR_susie/phenotype/{}.txt - | awk 'BEGIN{OFS=\"\\t\"} NR==1; NR>1 {print | \"sort -k1,1 -k2,2n\"}' | bgzip > data/tensorqtl/${xQTL_type}/combined_SNV_modifyTR_susie/phenotype/{}.bed.gz" | parallel -j 4 {}

    echo ${TISSUE//,/ } | tr " " "\n" | parallel -j 3 /media/dubai/home/dingruofan/anaconda3/envs/work/bin/python bin/tensorqtl/SNV_modifyTR_susie_map.py -i data/tensorqtl/${xQTL_type} -t {}
    done
}


pQTL_run(){
    /media/dubai/home/dingruofan/anaconda3/envs/work/bin/python bin/tensorqtl/cis_map.py -i data/tensorqtl/SNV_pQTL -o data/tensorqtl/SNV_pQTL -g data/genotype/GMTiP -w 1000000 -t Whole_Blood
    /media/dubai/home/dingruofan/anaconda3/envs/work/bin/python bin/tensorqtl/calculate_qvalues.py -i data/tensorqtl/SNV_pQTL/tmp_permutation -o data/tensorqtl/SNV_pQTL/permutation -t Whole_Blood

    /media/dubai/home/dingruofan/anaconda3/envs/work/bin/python bin/tensorqtl/cis_map.py -i data/tensorqtl/SNV_pQTL -o data/tensorqtl/SV_pQTL -g data/genotype/GMTiP_SV -w 1000000 -t Whole_Blood
    /media/dubai/home/dingruofan/anaconda3/envs/work/bin/python bin/tensorqtl/calculate_qvalues.py -i data/tensorqtl/SV_pQTL/tmp_permutation -o data/tensorqtl/SV_pQTL/permutation -t Whole_Blood
}


#### Identify EAS specific SNP, TR, and SV
AF_SNV(){
    ## SNP --------------------------
    mkdir -p input/GMTiP/SNP
    CURRENTSNP=GMTiP_LRS131.snv_indel.maf05.reheader.vcf.gz
    SNPfile=/media/london_B/zouxudong/2024-10-21-aGTEx-main/Genotype/${CURRENTSNP}
    GTOPSNP=input/GMTiP/SNP/GMTiP.maf05.reheader.vcf.gz

    ln -s ${SNPfile} ${GTOPSNP} && bcftools index -t ${GTOPSNP}

    /media/bora_A/zhangt/src/bin/plink2 --max-alleles 2 --update-sex input/GMTiP/sex_info.txt --vcf-half-call missing --split-par 'hg38' --make-bed --vcf input/GMTiP/SNP/GMTiP.maf05.reheader.vcf.gz --out input/GMTiP/SNP/GMTiP
    plink2 --max-alleles 2 --update-sex input/GMTiP/sex_info.txt --vcf-half-call missing --make-pgen --vcf input/GMTiP/SNP/GMTiP.maf05.reheader.vcf.gz --out input/GMTiP/SNP/GMTiP

    ### PCA
    /media/bora_A/zhangt/src/QTLtools_1.2_Ubuntu16.04_x86_64/QTLtools_1.2_Ubuntu16.04_x86_64 pca --vcf input/GMTiP/SNP/GMTiP.maf05.reheader.vcf.gz --scale --center --distance 50000 --out input/GMTiP/SNP/GMTiP_geno_PCA

    zcat input/GMTiP/SNP/GMTiP.maf05.reheader.vcf.gz | awk '$0!~/^#/{print $1"_"$2"_"$4"_"$5"\t"$3}' | sort -u -k1,1 -k2,2 > input/GMTiP/SNP/GMTiP_SNP_chrpos_rsid.txt

    ## GMTiP EAS population
    /media/bora_A/zhangt/src/bin/plink2 --split-par 'hg38' --max-alleles 2 --update-sex input/GMTiP/sex_info.txt --freq --vcf-half-call missing --out input/GMTiP/SNP/SNP_EAS --vcf ${GTOPSNP}

    ## GTEx EUR population
    # GTEx_vcf=/media/pacific/share/Datasets/GTEx_data/GTEx_Genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz
    # awk '$3 == "European" {print $1}' input/GTEx/GTEx_V8.individual_race.txt > input/GTEx/GTEx_EUR_inds.txt
    # plink2 --vcf ${GTEx_vcf} \
    #         --keep input/GTEx/GTEx_EUR_inds.txt \
    #          --freq --out input/GTEx/GTEx_EUR
    # plink2 --freq --out 1000G_${POPNAME} --pfile 1000G_${POPNAME}_hg38/1000G.${POPNAME}

    ## 1KG

    ## gnomAD
    echo {1..22} X | tr " " "\n" | head -n 6 | tail -n+1 | parallel -j 6 bash bin/gnomAD_awk_af.sh {}
    cat input/gnomAD/SNP/tmp/extracted_snps_chr*.txt > input/gnomAD/SNP/gnomAD_SNP_AF.txt

    ## Annovar
    /media/bora_A/zhangt/src/bin/annovar/table_annovar.pl output/data_plot/GMTiP.avinput \
    /media/bora_A/zhangt/2023-09-01-sc-aQTL-Project/2023-11-04-xQTL_characteristics/2025-02-27-torus/Annotation/humandb/ \
    -buildver hg38 \
    -out output/data_plot/anno/GMTiP_anno \
    -protocol refGeneWithVer,encRegTfbsClustered,encTfChipPkENCFF003VDB,gwasCatalog,wgRna,rmsk \
    -operation g,r,r,r,r,r \
    -remove -nastring . -csvout -polish
    ##
    for i in ENCFF213AWVH3K27ac ENCFF354JNBH3K4me1 GSE118189ATAC UCSCpolyDB # enhanceratlasPBMC ADAR1binding
    do
    /media/bora_A/zhangt/src/bin/annovar/annotate_variation.pl output/data_plot/GMTiP.avinput /media/bora_A/zhangt/2023-09-01-sc-aQTL-Project/2023-11-04-xQTL_characteristics/2025-02-27-torus/Annotation/humandb/ -bedfile hg38_${i}.bed -buildver hg38 -dbtype bed -regionanno -colsWanted all -out output/data_plot/anno/GMTiP_anno.hg38_${i}
    done
    ## 
    /media/bora_A/zhangt/src/bin/annovar/annotate_variation.pl -out output/data_plot/anno/GMTiP_anno.hg38_gene -build hg38 -neargene 2000 output/data_plot/GMTiP.avinput /media/bora_A/zhangt/2023-09-01-sc-aQTL-Project/2023-11-04-xQTL_characteristics/2025-02-27-torus/Annotation/humandb/
}


AF_SV(){
    ## SV --------------------------
    ## GMTiP frequency ###
    SV_file=/media/london_B/zouxudong/2024-10-21-aGTEx-main/Genotype/SV/GMTiP_SV_LRS.maf05_131INDS.format.sorted.vcf.gz
    ln -s ${SV_file} input/GMTiP/SV/GMTiP_SV.maf05.vcf.gz
    bcftools annotate --remove ^FORMAT/GT -Oz -o input/GMTiP/SV/GMTiP_SV.vcf.gz input/GMTiP/SV/GMTiP_SV.maf05.vcf.gz
    bcftools index -t input/GMTiP/SV/GMTiP_SV.vcf.gz

    plink2 --max-alleles 2 --update-sex input/GMTiP/sex_info.txt --vcf-half-call missing --make-bed --vcf input/GMTiP/SV/GMTiP_SV.vcf.gz --out input/GMTiP/SV/GMTiP_SV
    /media/bora_A/zhangt/src/bin/plink2 --max-alleles 2 --update-sex input/GMTiP/sex_info.txt --vcf-half-call missing  --make-pgen --vcf input/GMTiP/SV/GMTiP_SV.vcf.gz --out input/GMTiP/SV/GMTiP_SV

    plink2 --max-alleles 2 --update-sex input/GMTiP/sex_info.txt --freq --vcf-half-call missing --out input/GMTiP/SV/SV_EAS --vcf input/GMTiP/SV/GMTiP_SV.vcf.gz

    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%REF\t%ALT\n' input/GMTiP/SV/GMTiP_SV.vcf.gz | awk '{OFS="\t";$2=$2-1;if($3=="."){$3=$2+1;print $0}else{print $0}}' > input/GMTiP/SV/GMTiP_LRS_SV.bed

    # 1KG frequency ###
    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.0/delly-unfiltered-hg38/all.delly.hg38.1kGP.ont.bcf
    bcftools view -G GT input/1KG/SV/all.delly.hg38.1kGP.ont.bcf > input/1KG/SV/tmp.vcf
    bcftools annotate --set-id '%CHROM\_%POS\_%SVTYPE\_%SVLEN' -x QUAL,FILTER,INFO,FORMAT/^GT -Oz -o input/1KG/SV/1KG_SV.vcf.gz input/1KG/SV/tmp.vcf
    rm input/1KG/SV/tmp.vcf
    Generate sub-population fileset
    for POPNAME in EUR AFR EAS
    do
        plink2 --vcf input/1KG/SV/1KG_SV.vcf.gz --allow-extra-chr \
                --keep input/1KG/SV/1000G_${POPNAME}_hg38/inds.txt \
                --geno 0.01 --maf 0.05 --hwe 1e-6 \
                --make-bed --out input/1KG/SV/1000G_${POPNAME}_hg38/1000G.${POPNAME}
    done
    cat input/1KG/SV/1000G_*_hg38/1000G*.bim | cut -f 2 | sort -u > input/1KG/SV/selected_variants.txt

    for POPNAME in EUR AFR EAS
    do
        plink2 --vcf input/1KG/SV/1KG_SV.vcf.gz --allow-extra-chr \
                --keep input/1KG/SV/1000G_${POPNAME}_hg38/inds.txt \
                --extract input/1KG/SV/selected_variants.txt \
                --make-bed --out input/1KG/SV/1000G_${POPNAME}_hg38/1000G.${POPNAME}
        plink2 --freq --allow-extra-chr --out input/1KG/SV/1000G_${POPNAME} --bfile 1000G_${POPNAME}_hg38/1000G.${POPNAME}
    done

    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%REF\t%ALT\n' input/1KG/SV/1KG_SV.vcf.gz | awk '{OFS="\t";$2=$2-1;if($3=="."){$3=$2+1;print $0}else{print $0}}' > input/1KG/SV/1KG_SV.bed

    bedtools intersect -a input/GMTiP/SV/GMTiP_LRS_SV.bed -b input/1KG/SV/1KG_SV.bed -wa -wb -f 0.5 -F 0.5 | awk '{split($4,a,"_");split($10,b,"_");if(a[3]==b[3]){print $0}}' > input/GMTiP/SV/GMTiP_LRS_inter_1KG_SV.txt
    cut -f 4 input/GMTiP/SV/GMTiP_LRS_inter_1KG_SV.txt | sort -u | wc # 23703

    ## gnomAD
    bcftools annotate --set-id '%CHROM\_%POS\_%SVTYPE\_%SVLEN' -x QUAL,FILTER -Oz -o input/gnomAD/SV/gnomAD_SV.vcf.gz input/gnomAD/SV/rawdata/gnomad.v4.1.sv.sites.vcf.gz
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%REF\t%ALT\t%INFO/AC_afr\t%INFO/AN_afr\t%INFO/AF_afr\t%INFO/AC_eas\t%INFO/AN_eas\t%INFO/AF_eas\t%INFO/AC_nfe\t%INFO/AN_nfe\t%INFO/AF_nfe\n' input/gnomAD/SV/gnomAD_SV.vcf.gz | awk '{OFS="\t";$2=$2-1;if($6=="<INS>"){$3=$2+1};print $0}' > input/gnomAD/SV/gnomAD_SV.bed

    bedtools intersect -a input/GMTiP/SV/GMTiP_LRS_SV.bed -b input/gnomAD/SV/gnomAD_SV.bed -wo -f 0.5 -F 0.5 | awk '{split($4,a,"_");split($10,b,"_");if(a[3]==b[3]){print $0}}' > input/GMTiP/SV/GMTiP_LRS_inter_gnomAD_SV.txt
    cut -f 4 input/GMTiP/SV/GMTiP_LRS_inter_gnomAD_SV.txt | sort -u | wc # 7742

## AnnotSV  ###########
# AnnotSV
    
    # SVDIR=/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/input/GMTiP/SV
    # cp /media/london_B/zouxudong/2024-10-21-aGTEx-main/Genotype/SV/old/GMTiP_SV_LRS.maf05_131INDS.format.sorted.vcf.gz  ${SVDIR}/GMTiP_SV.maf05.vcf.gz
    # gunzip ${SVDIR}/GMTiP_SV.maf05.vcf.gz
    # cp /media/london_B/zouxudong/2024-10-21-aGTEx-main/Genotype/SV/old/GMTiP_SV_LRS.maf05_131INDS.format.sorted.vcf.gz  ${SVDIR}/GMTiP_SV.maf05.vcf.gz
    # /media/bora_A/zhangt/src/bin/SURVIVOR/Debug/SURVIVOR vcftobed ${SVDIR}/GMTiP_SV.maf05.vcf -99999999 99999999 ${SVDIR}/GMTiP_SV.maf05.vcf.bed

    # cat ${SVDIR}/GMTiP_SV.maf05.vcf.bed | cut -f1,2,5,7,11 | awk -F '\t' '{ $3 = ($3 == "0" ? $2+1 : $3) } 1' OFS='\t' > ${SVDIR}/GMTiP_SV.maf05.vcf.clean.bed

    # /media/bora_A/zhangt/src/bin/AnnotSV/bin/AnnotSV -SVinputFile ${SVDIR}/GMTiP_SV.maf05.vcf.clean.bed -outputDir ${SVDIR}/AnnotSV -SVinputInfo 1 -reciprocal 1 -svtBEDcol 5

    /media/bora_A/zhangt/src/bin/AnnotSV/bin/AnnotSV -SVinputFile ${SVDIR}/GMTiP_SV.maf05.vcf.gz -outputFile ${SVDIR}/AnnotSV

    ${svpipeline_r_lib}/complementAnnotation.R ${SVDIR}/AnnotSV/../SV_AnnotSV.txt.tsv ${runDir}/AnnotSV/samples_merged_ALL.Final.clean.annotated.plus.tsv h38


    zcat ${SV_file} | awk 'BEGIN {OFS="\t"} /^##/ {print} !/^##/ {print $1,$2,$3,$4,$5,$6,$7,$8}' | awk 'NR==FNR{a[$1]=$1;next}{if($0~/^#/ || a[$3]){print $0}}' output/data_plot/SV_EAS_specific.txt - > tmp.sv.vcf

    zcat ${SV_file} | awk 'BEGIN {OFS="\t"} /^##/ {print} !/^##/ {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > input/GMTiP/SV_AnnotSV.vcf
    /media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/bin/AnnotSV/bin/AnnotSV -SVinputFile input/GMTiP/SV_AnnotSV.vcf -outputFile input/GMTiP/SV_AnnotSV.txt

    indNam=AK221
    longshot \
    --bam /media/iceland/share/Datasets/Asian_GTEx/WGS_Long_reads/Mapping/${indNam}-DMSO.align.bam \
    --ref /media/london_B/lixing/2024-08-29-TR-AsianGTEX/2024-12-11-TR-LongReads/BAM/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
    --out_bam ${indNam}-DMSO_chr16_1238154_1247284.haplotype.bam \
    --out ${indNam}-DMSO_chr16_1238154_1247284.vcf \
    -r chr16:1128154-1357284
    samtools index ${indNam}-DMSO_chr16_1238154_1247284.haplotype.bam

    # indNam=AK111
    # samtools view -b input.bam "chr1:1000000-2000000" > region.bam
    # whatshap haplotag -o ${indNam}-DMSO_chr16_1238154_1247284.haplotype.bam \
    # --reference /media/london_B/lixing/2024-08-29-TR-AsianGTEX/2024-12-11-TR-LongReads/BAM/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
    # ${indNam}-DMSO_chr16_1238154_1247284.vcf \
    # /media/iceland/share/Datasets/Asian_GTEx/WGS_Long_reads/Mapping/${indNam}-DMSO.align.bam
}



lixing_SV(){
    ## Lixing final-vcf.phased.vcf.gz
    for POPNAME in EUR AFR EAS
    do
        plink2 --vcf input/1KG/SV/final-vcf.phased.vcf.gz --allow-extra-chr \
                --keep input/1KG/SV/1000G_${POPNAME}_hg38/inds.txt \
                --geno 0.01 --maf 0.05 --hwe 1e-6 \
                --make-bed --out input/1KG/SV/1000G_${POPNAME}_hg38/1000G.${POPNAME}
    done
    cat input/1KG/SV/1000G_*_hg38/1000G*.bim | cut -f 2 | sort -u > input/1KG/SV/selected_variants.txt

    for POPNAME in EUR AFR EAS
    do
        plink2 --vcf input/1KG/SV/final-vcf.phased.vcf.gz --allow-extra-chr \
                --keep input/1KG/SV/1000G_${POPNAME}_hg38/inds.txt \
                --extract input/1KG/SV/selected_variants.txt \
                --make-bed --out input/1KG/SV/1000G_${POPNAME}_hg38/1000G.${POPNAME}
        plink2 --freq --allow-extra-chr --out input/1KG/SV/1000G_${POPNAME} --bfile input/1KG/SV/1000G_${POPNAME}_hg38/1000G.${POPNAME}
    done
}


SNP_SV_LD(){
    CUR_DIR=/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL

    snp_file=/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/input/GMTiP/SNP/GMTiP.maf05.reheader.vcf.gz
    SV_file=/media/london_B/zouxudong/2024-10-21-aGTEx-main/Genotype/SV/GMTiP_SV_LRS.maf05_131INDS.format.sorted.vcf.gz
    # bcftools query -l ${snp_file} > $CUR_DIR/output/GWAS_LD/SNP_samples.txt
    # bcftools query -l ${SV_file} > $CUR_DIR/output/GWAS_LD/SV_samples.txt
    # diff $CUR_DIR/output/GWAS_LD/SNP_samples.txt $CUR_DIR/output/GWAS_LD/SV_samples.txt

    bcftools view -i "ID=@$CUR_DIR/output/GWAS_LD/SNP_SV_ID.txt" ${snp_file} -Oz -o $CUR_DIR/output/GWAS_LD/SNP_variants.vcf.gz  --threads 32
    bcftools view -i "ID=@$CUR_DIR/output/GWAS_LD/SNP_SV_ID.txt" ${SV_file} -Oz -o $CUR_DIR/output/GWAS_LD/SV_variants.vcf.gz

    bcftools index $CUR_DIR/output/GWAS_LD/SNP_variants.vcf.gz
    bcftools index $CUR_DIR/output/GWAS_LD/SV_variants.vcf.gz

    bcftools concat -a $CUR_DIR/output/GWAS_LD/SNP_variants.vcf.gz $CUR_DIR/output/GWAS_LD/SV_variants.vcf.gz -oZ -o $CUR_DIR/output/GWAS_LD/merge_SNP_SV.vcf.gz

    /media/bora_A/zhangt/src/bin/plink2 --update-sex input/GMTiP/sex_info.txt --vcf $CUR_DIR/output/GWAS_LD/merge_SNP_SV.vcf.gz --r2-unphased --ld-window-kb 1000 --ld-window-r2 0.8 --out $CUR_DIR/output/GWAS_LD/merge_SNP_SV_ld_results --split-par 'hg38'
}


function get_GWAS_variant_genotype(){
    host_cur=`hostname`
    if [[ "${host_cur}" == 'mtcook' || "${host_cur}" == 'madeira' ]]; then
        module load bioinformatics/plink_1.90beta
    elif [[ "${host_cur}" == 'alps1' || "${host_cur}" == 'alps2' ]]; then
        module load biosoft/plink_1.90beta
    fi

    snp_file=/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/input/GMTiP/SNP/GMTiP.maf05.reheader.vcf.gz

    cat $CUR_DIR/output/GWAS_LD/get_sig_GWAS.txt |awk 'BEGIN{FS=OFS="\t"}{print $1,$2-1,$2}'  |awk 'NR==1 {print; next} /^chr/ {print; next} {print "chr"$0}'| sort -k1,1 -k2,2n |uniq > $CUR_DIR/output/GWAS_LD/get_sig_GWAS.bed

    # 1. Extract the genotypes of the GWAS significant loci in the GMTiP LRS samples
    bcftools view -R $CUR_DIR/output/GWAS_LD/get_sig_GWAS.bed -c 1 $snp_file -Oz -o  $CUR_DIR/output/GWAS_LD/GMTiP_SNV_INDEL_LRS.get_sig_GWAS.vcf.gz --threads 32

    # 2. Replace rsid with chr:pos to facilitate the subsequent LD analysis in matching the positions
    /media/bora_A/zhangt/miniconda3/bin/python /media/london_B/lixing/2024-08-29-TR-AsianGTEX/2025-04-22-GWAS-TR/bin/replace_rsid_to_chrpos.py $CUR_DIR/output/GWAS_LD/GMTiP_SNV_INDEL_LRS.get_sig_GWAS.vcf.gz $CUR_DIR/output/GWAS_LD/GMTiP_SNV_INDEL_LRS.get_sig_GWAS.changeID.vcf.gz

    # 3. Split chromosomes convert SNV genotypes into dosage format
    mkdir -p $CUR_DIR/output/GWAS_LD/GMTiP_SNV_INDEL_dosage_chr
    for CHROM in {1..22} X Y
    do
    plink --vcf $CUR_DIR/output/GWAS_LD/GMTiP_SNV_INDEL_LRS.get_sig_GWAS.changeID.vcf.gz --keep-allele-order \
        --double-id --recode A --make-bed --chr $CHROM \
        --out $CUR_DIR/output/GWAS_LD/GMTiP_SNV_INDEL_dosage_chr/GMTiP_SNV_INDEL_LRS.$CHROM

    awk 'NR==1 {for(i=1;i<=NF;i++) sub(/_.*/,"",$i)} 1' $CUR_DIR/output/GWAS_LD/GMTiP_SNV_INDEL_dosage_chr/GMTiP_SNV_INDEL_LRS.$CHROM.raw > $CUR_DIR/output/GWAS_LD/GMTiP_SNV_INDEL_dosage_chr/GMTiP_SNV_INDEL_LRS.$CHROM.raw.final
    done
    wait

    rm -rf $CUR_DIR/output/GWAS_LD/GMTiP_SNV_INDEL_LRS.get_sig_GWAS.vcf.{gz,gz.tbi}
    rm -rf $CUR_DIR/output/GWAS_LD/GMTiP_SNV_INDEL_dosage_chr/GMTiP_SNV_INDEL_LRS.*.{bed,bim,fam,log,nosex,raw}
}


function get_TR_genotype(){
    mkdir -p $CUR_DIR/output/GWAS_LD/GMTiP_TR_LRS_dosage_chr

    TR_file=/media/london_B/lixing/2024-08-29-TR-AsianGTEX/2024-12-11-TR-LongReads/output/QTL_Filter/GMTip_LRS_TR_131INDs.Miss85.AF95.AC1.dosage.impute.simpleTR.filter.txt

    awk 'NR==FNR{a[$1]=$1;next}{if(a[$1]||$1=="TRID"){print $0}}' $CUR_DIR/output/GWAS_LD/selected_TR.txt $TR_file | sed '1d' | awk 'BEGIN{FS="_";OFS="\t"}{print $1,$2,$3}' | sort -k1,1 -k2,2n | bedtools window -a $CUR_DIR/output/GWAS_LD/get_sig_GWAS.bed -b stdin -w 1000000 > $CUR_DIR/output/GWAS_LD/get_sig_GWAS.window1M_TRs.txt

    # Split chromosomes convert TR genotypes into dosage format
    awk 'NR==FNR{a[$1]=$1;next}{if(a[$1]||$1=="TRID"){print $0}}' $CUR_DIR/output/GWAS_LD/selected_TR.txt $TR_file |awk -F'\t' 'BEGIN{OFS="\t"} {split($1, a, "_"); $1 = a[1]":"a[2]; print}' |awk 'BEGIN{FS=OFS="\t"} NR==1{$1="site"} {print}'> $CUR_DIR/output/GWAS_LD/GMTip_LRS_TR_131INDs.dosage.txt

    rm $CUR_DIR/output/GWAS_LD/GMTiP_TR_LRS_dosage_chr/*
    chroms=( {1..22} )
    for CHROM in "${chroms[@]}"; do
        awk -v chrom="$CHROM" '
            NR==1 {print} 
            $1 ~ "^chr"chrom":" {print}
        ' "$CUR_DIR/output/GWAS_LD/GMTip_LRS_TR_131INDs.dosage.txt" > "$CUR_DIR/output/GWAS_LD/GMTiP_TR_LRS_dosage_chr/GMTip_LRS_TR_131INDs.dosage.chr$CHROM"
    done
}

function compute_LD(){
#time python get_SNP_STR_LD.v2.py --pair gwas_snps.window250k_pSTRs.txt --chrom chr21 

    rm $CUR_DIR/output/GWAS_LD/TR_SNV_LD/* $CUR_DIR/output/GWAS_LD/slurm/*

    chroms=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 )
    for CHROM in ${chroms[@]}
    do 
    num=`echo "$CHROM" |awk -F"chr" '{print $2}'`
    if [ ! -f $CUR_DIR/output/GWAS_LD/TR_SNV_LD/GMTiP.SNP_TR.LD.$CHROM.txt ]; then
    echo '#!/bin/bash

    ################
    echo "process will start at:"
    date

    Rscript '$CUR_DIR'/bin/get_SNP_STR_LD_split_chrom.R \
    --pair '$CUR_DIR'/output/GWAS_LD/get_sig_GWAS.window1M_TRs.txt \
    --chrom  '$CHROM' \
    --snp '${CUR_DIR}'/output/GWAS_LD/GMTiP_SNV_INDEL_dosage_chr/GMTiP_SNV_INDEL_LRS.'$num'.raw.final \
    --str '${CUR_DIR}'/output/GWAS_LD/GMTiP_TR_LRS_dosage_chr/GMTip_LRS_TR_131INDs.dosage.'$CHROM' \
    --output '$CUR_DIR'/output/GWAS_LD/TR_SNV_LD/ \
    --outLD '$CUR_DIR'/output/GWAS_LD/TR_SNV_LD/GMTiP.SNP_TR.LD.'$CHROM'.txt

    wait
    echo "process will end at:"
    date
    ' > $CUR_DIR/output/GWAS_LD/slurm/compute_ld.$CHROM.sh
    # cd $CUR_DIR/output/GWAS_LD/slurm
    # nohup bash compute_ld.$CHROM.sh > compute_ld.$CHROM.log &
    fi
    done

    ls /media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/output/GWAS_LD/slurm/* | xargs -I {} bash {}

    # neat results
    for file in $CUR_DIR/output/GWAS_LD/TR_SNV_LD/GMTiP.SNP_TR.LD.*.txt
    do
    chr=$(basename $file |awk -F"." '{print $4}')
    cat $file |sed '1d' |awk '$7!="NA" && $7>=0.5' >> $CUR_DIR/output/GWAS_LD/GWAS_TR.High_LD_05.list
    done

    for file in $CUR_DIR/output/TR_SNV_LD/GMTiP.SNP_TR.LD.*.txt
    do
    chr=$(basename $file |awk -F"." '{print $4}')
    cat $file |sed '1d' |awk '$7!="NA"' >> $CUR_DIR/output/GMTiP.SNP_TR.total_LD.list
    done
}



a(){
    /media/bora_A/zhangt/src/bin/plink2 --pfile /media/bora_A/zhangt/single_cell_datasets/seven_genotypes_summary/population/1000G_EAS_hg38/1000G.EAS \
    --extract output/GWAS_LD/tmp_variants.txt \
    --r2-unphased \
    --ld-window-kb 1000 \
    --ld-window-r2 0.6 \
    --out output/GWAS_LD/tmp_variants
    # 585 samples

    /media/bora_A/zhangt/src/bin/plink2 --pfile /media/bora_A/zhangt/single_cell_datasets/seven_genotypes_summary/population/1000G_EAS_hg38/1000G.EAS     --extract output/GWAS_LD/tmp_variant_EAS.txt     --r2-unphased     --ld-window-kb 1000     --ld-window-r2 0.8     --out output/GWAS_LD/tmp_variant_EAS

    /media/bora_A/zhangt/src/bin/plink2 --pfile /media/bora_A/zhangt/single_cell_datasets/seven_genotypes_summary/population/1000G_EUR_hg38/1000G.EUR    --extract output/GWAS_LD/tmp_variant_EUR.txt     --r2-unphased     --ld-window-kb 1000     --ld-window-r2 0.8     --out output/GWAS_LD/tmp_variant_EUR
}
