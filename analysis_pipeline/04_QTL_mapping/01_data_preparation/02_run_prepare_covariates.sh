#!/bin/bash

main(){

	s1_run_pca
	s2_prepare_covariates
}

function s1_run_pca(){
	curr_dir=`pwd`
	outdir=$curr_dir/output/pca_res

	if [ ! -d "$outdir" ]
	then
		mkdir -p $outdir
	fi

	for tissue in `cat $curr_dir/input/selected_tissues_for_qtl.txt`
	do
		echo $tissue
		Rscript ./scripts/pca_on_phenotype.R -t $tissue
	done
}

function s2_merge_rnaseqc(){
	curr_dir=`pwd`
	outdir=$curr_dir/output/covariates

	mkdir -p $outdir

	for tissue in `cat $curr_dir/input/selected_tissues_for_qtl.txt` 
	do
		Rscript ./scripts/prepare_covariates.R -t $tissue
	done
}

# >> main
main
