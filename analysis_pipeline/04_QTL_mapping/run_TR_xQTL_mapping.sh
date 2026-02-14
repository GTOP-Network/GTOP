#!/bin/bash

TR_dosage=$1 # TR dosages (the averaged TR repeat counts over the two alleles)
phenotype=$2 # Bed file with normalized phenotype expression for samples
covFile=$3 # Formatted covariates file
outDir=$4 # Directory to write output files to
outPrefix=$5 # Prefix of output files


mkdir -p ${outDir}/output/cis_eQTL/independent
mkdir -p ${outDir}/output/cis_eQTL/nominal
mkdir -p ${outDir}/output/cis_eQTL/permutation
mkdir -p ${outDir}/output/trans_eQTL

conda activate tensorqtl_env
#TR eQTL
python3 run_tensorqtl_TR_eQTL.py -t ${tissue} -g ${TR_dosage} -p ${phenotype} -c ${covFile} \
-N ${outDir}/output/cis_eQTL/nominal/ \
-P ${outDir}/output/cis_eQTL/permutation/${tissue}.perm.cis_qtl.txt \
-I ${outDir}/output/cis_eQTL/independent/${tissue}.significant_pairs.txt \
-T ${outDir}/output/trans_eQTL/${tissue}


#TR sQTL(junction usage)
python3 run_tensorqtl_TR_sQTL.py -t ${tissue} -g ${TR_dosage} -p ${phenotype} -c ${covFile} \
-G ${phenotype_groups} \
-N ${outDir}/output/cis_juQTL/nominal/ \
-P ${outDir}/output/cis_juQTL/permutation/${tissue}.perm.cis_qtl.txt \
-I ${outDir}/output/cis_juQTL/independent/${tissue}.significant_pairs.txt \
-T ${outDir}/output/trans_juQTL/${tissue}

#TR sQTL(transcript usage)
python3 run_tensorqtl_TR_sQTL.py -t ${tissue} -g ${TR_dosage} -p ${phenotype} -c ${covFile} \
-G ${phenotype_groups} \
-N ${outDir}/output/cis_tuQTL/nominal/ \
-P ${outDir}/output/cis_tuQTL/permutation/${tissue}.perm.cis_qtl.txt \
-I ${outDir}/output/cis_tuQTL/independent/${tissue}.significant_pairs.txt \
-T ${outDir}/output/trans_tuQTL/${tissue}

