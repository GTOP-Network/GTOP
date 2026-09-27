#!/bin/bash

dir_permutation=${1}
dir_nominal=${2}
dir_output=${3}
tis_i=${4}
# study=${5}
# state=${6}

mkdir ${dir_output}

############################################################################################
#   1: combine_signif_pairs_tjy.py   2: extract_pairs_tjy.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 2. extract all nominal pairs from nominal results for each tissue
# permutation results
# perm_files=(`find ${dir_permutation} -name "*.txt" | grep ${study} | grep ${state}`)
perm_files=(`find ${dir_permutation} -name "*.txt"`)
NAMEs=(${perm_files[@]/*\//})
NAMEs=(${NAMEs[@]/.txt/})

name=${NAMEs[tis_i]}
nominal_files=(`find ${dir_nominal} -name "${name}.txt.gz"`)
rm -f ${dir_output}/${name}.nominal_files.txt
for l in ${nominal_files[*]}
do
{
    echo ${l} >> ${dir_output}/${name}.nominal_files.txt
}
done
# extract_pairs
if [[ ! -f ${dir_output}/${name}.nominal_pairs.extracted_pairs.txt.gz ]];then
python /media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-09-30-mash/bin/MashR/extract_pairs_tjy.py ${dir_output}/${name}.nominal_files.txt ${dir_output}/nominal_pairs.combined_signifpairs.txt.gz ${name}.nominal_pairs -o ${dir_output}
fi
# > output file: *_nominal_pairs.extracted_pairs.txt.gz
rm -f ${dir_output}/${name}.nominal_files.txt
