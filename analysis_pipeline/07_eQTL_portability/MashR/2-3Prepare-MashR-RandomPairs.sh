#!/bin/bash

dir_output=${1}
study=${2}
state=${3}

mkdir ${dir_output}

############################################################################################
#   1: combine_signif_pairs_tjy.py   2: extract_pairs_tjy.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 3. prepare random SNP-gene pairs for MashR
# nominal_pairs_files=(`ls ${dir_output}/*nominal_pairs.extracted_pairs.txt.gz | grep ${study} | grep ${state}`)
nominal_pairs_files=(`ls ${dir_output}/*nominal_pairs.extracted_pairs.txt.gz`)
rm ${dir_output}/nominal_pairs_list_file.txt
for l in ${nominal_pairs_files[*]}
do
{
    echo ${l} >> ${dir_output}/nominal_pairs_list_file.txt
}
done

python /media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-09-30-mash/bin/MashR/mashr_prepare_input.py ${dir_output}/nominal_pairs_list_file.txt nominal_pairs -o ${dir_output} --only_zscore --dropna
