#!/bin/bash 

dir=`pwd`
snv_DIR=/lustre/home/xdzou/2024-10-21-GTBMap/2026-02-02-GTOP_eQTL_mapping/output/QTL_mapping/sQTL/SNV
sv_DIR=/lustre/home/xdzou/2024-10-21-GTBMap/2026-02-02-GTOP_eQTL_mapping/output/QTL_mapping/sQTL/SV
tr_DIR=/lustre/home/xdzou/2024-10-21-GTBMap/2026-02-02-GTOP_eQTL_mapping/output/QTL_mapping/sQTL/TR

outdir_ju=$dir/input/tissue_ju_snv
outdir_tu=$dir/input/tissue_tu_snv
for tissue in `cat $dir/selected_tissues_11.txt|cut -f1`
do
	echo $tissue
	cat $snv_DIR/${tissue}.snv_juqtl.txt |tail -n+2|awk '{print $11"\t"$1}' > $outdir_ju/${tissue}_gene_list.txt
	cat $snv_DIR/${tissue}.snv_tuqtl.txt |tail -n+2|awk '{print $11"\t"$1}' > $outdir_tu/${tissue}_gene_list.txt
done

outdir_ju=$dir/input/tissue_ju_sv
outdir_tu=$dir/input/tissue_tu_sv

for tissue in `cat $dir/selected_tissues_11.txt|cut -f1`
do
	echo $tissue
	cat $sv_DIR/${tissue}.sv_juqtl.txt|tail -n+2|awk '{print $11"\t"$1}' > $outdir_ju/${tissue}_gene_list.txt
	cat $sv_DIR/${tissue}.sv_tuqtl.txt |tail -n+2|awk '{print $11"\t"$1}' > $outdir_tu/${tissue}_gene_list.txt
done

outdir_ju=$dir/input/tissue_ju_tr
outdir_tu=$dir/input/tissue_tu_tr
for tissue in `cat $dir/selected_tissues_11.txt|cut -f1`
do
	echo $tissue
	cat $tr_DIR/${tissue}.tr_juqtl.txt|tail -n+2|awk '{print $11"\t"$1}' > $outdir_ju/${tissue}_gene_list.txt
	cat $tr_DIR/${tissue}.tr_tuqtl.txt|tail -n+2|awk '{print $11"\t"$1}' > $outdir_tu/${tissue}_gene_list.txt
done


outdir_ju=$dir/input/tissue_ju_joint
outdir_tu=$dir/input/tissue_tu_joint

eGenedir=/lustre/home/xdzou/2024-10-21-GTBMap/2026-02-03-GTOP_eQTL_Finemap_susie/input

for tissue in `cat $dir/selected_tissues_11.txt|cut -f1`
do
	echo $tissue
	mkdir -p $outdir_ju/$tissue
	mkdir -p $outdir_tu/$tissue
	cat $eGenedir/tissue_ju_snv/${tissue}_gene_list.txt $eGenedir/tissue_ju_sv/${tissue}_gene_list.txt $eGenedir/tissue_ju_tr/${tissue}_gene_list.txt |sort|uniq > $outdir_ju/${tissue}_gene_list.txt
	cat $eGenedir/tissue_tu_snv/${tissue}_gene_list.txt $eGenedir/tissue_tu_sv/${tissue}_gene_list.txt $eGenedir/tissue_tu_tr/${tissue}_gene_list.txt |sort|uniq > $outdir_tu/${tissue}_gene_list.txt
	split -l 100 $outdir_ju/${tissue}_gene_list.txt -d $outdir_ju/${tissue}/${tissue}.task_
	split -l 100 $outdir_tu/${tissue}_gene_list.txt -d $outdir_tu/${tissue}/${tissue}.task_
done

