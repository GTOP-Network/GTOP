#!/bin/bash 

dir=`pwd`
outdir=$dir/input/tissue_gene_joint
mkdir -p $outdir
eGenedir=/lustre/home/xdzou/2024-10-21-GTBMap/2026-02-03-GTOP_eQTL_Finemap_susie/input

for tissue in `cat $dir/selected_tissues_11.txt|cut -f1`
do
	echo $tissue
	if [ ! -d "$outdir/$tissue" ]
	then
		mkdir -p $outdir/$tissue
	fi
	cat $eGenedir/tissue_gene_snv/${tissue}_gene_list.txt $eGenedir/tissue_gene_sv/${tissue}_gene_list.txt $eGenedir/tissue_gene_tr/${tissue}_gene_list.txt |sort|uniq > $outdir/${tissue}_gene_list.txt
	split -l 100 $outdir/${tissue}_gene_list.txt -d $outdir/${tissue}/${tissue}.task_
done

