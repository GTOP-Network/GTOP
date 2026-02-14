#!/bin/bash

dir_nominal=${1}
dir_output=${2}

mkdir ${dir_output}

############################################################################################
#   1: combine_signif_pairs_tjy.py   2: extract_pairs_tjy.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 3. prepare strong SNP-gene pairs for MashR
# strong_pairs_files=(`ls ${dir_output}/*.extracted_pairs.txt.gz | grep ${study} | grep ${state}`)
strong_pairs_files=(`ls ${dir_output}/*.extracted_pairs.txt.gz`)
rm -f ${dir_output}/strong_pairs_files.txt
for l in ${strong_pairs_files[*]}
do
{
    echo ${l} >> ${dir_output}/strong_pairs_files.txt
}
done

# MashR format file (z-score)
python /media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-09-30-mash/bin/MashR/mashr_prepare_input.py ${dir_output}/strong_pairs_files.txt strong_pairs -o ${dir_output} --only_zscore
python /media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-09-30-mash/bin/MashR/mashr_prepare_input.py ${dir_output}/strong_pairs_files.txt strong_beta_se -o ${dir_output} --output_zscore
zcat ${dir_output}/strong_pairs.MashR_input.txt.gz | sed -e 's/_zval//g' | gzip > ${dir_output}/strong_pairs.temp.txt.gz
rm -f ${dir_output}/strong_pairs.MashR_input.txt.gz
mv ${dir_output}/strong_pairs.temp.txt.gz ${dir_output}/strong_pairs.MashR_input.txt.gz
