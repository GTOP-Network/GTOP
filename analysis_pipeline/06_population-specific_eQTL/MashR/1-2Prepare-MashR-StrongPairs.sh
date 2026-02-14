#!/bin/bash

dir_permutation=${1}
dir_nominal=${2}
dir_output=${3}

mkdir -p ${dir_output}

############################################################################################
#   1: combine_signif_pairs_tjy.py   2: extract_pairs_tjy.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 2. extract top pairs from nominal results for each tissue
nFile=`cat ${dir_output}/permutation_files.txt | wc -l`
# echo $nFile # 0..34
# perm_files=(`find ${dir_permutation} -name "*.txt" | grep ${study} | grep ${state}`)
perm_files=(`find ${dir_permutation} -name "*.txt"`)
NAMEs=(${perm_files[@]/*\//})
NAMEs=(${NAMEs[@]/.txt/})
for tis_i in `seq 0 $[nFile-1]`
do
{
    name=${NAMEs[tis_i]}
    # nominal_files=(`ls ${dir_nominal}/${tissue}/${tissue}.cis_qtl_pairs.*.txt.gz`)
    nominal_files=(`find ${dir_nominal} -name "${name}.txt.gz"`)
    rm -f ${dir_output}/${name}.nominal_files.txt
    for l in ${nominal_files[*]}
    do
    {
        echo ${l} >> ${dir_output}/${name}.nominal_files.txt
    }
    done
    # extract_pairs
    python /media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-09-30-mash/bin/MashR/extract_pairs_tjy.py ${dir_output}/${name}.nominal_files.txt ${dir_output}/strong_pairs.combined_signifpairs.txt.gz ${name} -o ${dir_output}
    ### output file: *.extracted_pairs.txt.gz
    # rm -f ${dir_output}/${name}.nominal_files.txt
} &
done
wait

