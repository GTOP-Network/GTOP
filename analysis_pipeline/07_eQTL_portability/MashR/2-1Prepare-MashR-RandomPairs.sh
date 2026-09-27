#!/bin/bash

dir_nominal=${1}
dir_output=${2}
# study=${3}
# state=${4}

mkdir -p ${dir_output}

############################################################################################
#   1: combine_signif_pairs_tjy.py   2: extract_pairs_tjy.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 1. output all nominal SNP-gene pairs information from permutation results to a .txt file
# (colnames: phenotype_id,variant_id,chr,pos)
# file list of permutation results
rm -f -r ${dir_output}/nominal_combined_files.txt
# nominal_combined_files=(`find ${dir_nominal} -name "*.txt.gz" | grep ${study} | grep ${state}`)
nominal_combined_files=(`find ${dir_nominal} -name "*.txt.gz"`)
for l in ${nominal_combined_files[*]}
do
{
    echo ${l} >> ${dir_output}/nominal_combined_files.txt
}
done
/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-09-30-mash/bin/MashR/combine_signif_pairs_tjy.py ${dir_output}/nominal_combined_files.txt nominal_pairs -o ${dir_output}
### output file: nominal_pairs.combined_signifpairs.txt.gz
rm -f ${dir_output}/nominal_combined_files.txt
