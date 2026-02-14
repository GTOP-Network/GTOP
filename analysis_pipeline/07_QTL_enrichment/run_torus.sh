#!/bin/bash

# Bash Script Pipeline for xQTL Data Processing and torus Analysis

# Description:
# This script processes the eQTL-input files, preparation of input files for torus and run torus.

# Example Usage:
# ./script.sh <tensorqtl_nominal_file> <torus_input_file> <anno_file> <torus_res_file>

tensorqtl_nominal_file=$1
torus_input_file=$2
anno_file=$3
torus_res_file=$4

main(){
prepare
torus
}

prepare(){

#torus input file need: phentoype_id variant_id distance p slope slope_se
pigz  -dc ${tensorqtl_nominal_file} |awk -F"\t" -v OFS="\t" '{print $1, $2, $3, $7, $8, $9}' |pigz -c >  ${torus_input_file}

}
torus(){

torus.static -est -d ${torus_input_file} -annot ${anno_file} --fastqtl > ${torus_res_file}


}



main

