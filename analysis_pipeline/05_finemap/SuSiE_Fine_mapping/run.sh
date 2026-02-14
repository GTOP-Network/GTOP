#!/bin/bash
# main function

main(){
	prepare_variants_gene_coordinates
	generate_vcf_by_gene
	get_tissue_gene_list
	run_prepare_phenotype_by_eGene
	add_genotype_to_each_gene
	run_susie_analysis
	run_summarize_susie
}

function prepare_variants_gene_coordinates(){
	currDir=`pwd`
	phenotype_dir=/lustre/home/xdzou/2024-10-21-GTBMap/2026-02-02-GTOP_eQTL_mapping/output/phenotype
	genotype_dir=/lustre/home/xdzou/2024-10-21-GTBMap/2026-02-02-GTOP_eQTL_mapping/output/genotype

# prepare TSS coordinates
	cat $phenotype_dir/*.bed|cut -f1-4|grep -v "phenotype_id"|sort -k1,1 -k2,2n|uniq > $currDir/input/TSS.GTOP_tested_genes.sorted.bed

# prepare small variants coordinates
	cat $genotype_dir/GTOP_LRS.small_variants.maf05.GT.pvar |grep -v "##"|tail -n+2|awk '{print $1"\t"$2-1"\t"$2"\t"$3}'|sort -k1,1 -k2,2n > $currDir/input/GTOP_SNV.coord.sorted.bed
	cat $genotype_dir/GTOP.Phase_I.LRS_160_SV.GT.pvar |grep -v "##"|tail -n+2|awk '{print $1"\t"$2-1"\t"$2"\t"$3}'|sort -k1,1 -k2,2n > $currDir/input/GTOP_SV.coord.sorted.bed


}
# summarize susie
function run_summarize_susie(){
	currDir=`pwd`
	for tissue in `cat $currDir/selected_tissues_11.txt|cut -f1`
	do
		echo $tissue
		echo "#!/bin/bash
#SBATCH --job-name=tr_$tissue
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --error=tr_${tissue}.err
#SBATCH --output=tr_${tissue}.out

TISSUE=$tissue
DIR=$currDir
WKDIR=/flashfs1/scratch.global/xdzou/Fine_map_susie
" > $currDir/submit_summarize_susie.${tissue}.slurm
		echo '
cd $SLURM_SUBMIT_DIR

Rscript $DIR/src/summarize_susie_results_by_tissue.R -d $WKDIR -t $TISSUE -g $DIR/input/tissue_gene_tr/${TISSUE}_gene_list.txt

echo "process end at:"
date
' >> $currDir/submit_summarize_susie.${tissue}.slurm
		sbatch $currDir/submit_summarize_susie.${tissue}.slurm
	done
}


function run_prepare_phenotype_by_eGene(){
	currDir=`pwd`
	for tissue in `cat $currDir/selected_tissues_11.txt|cut -f1`
	do
		echo $tissue
		outDir=/flashfs1/scratch.global/xdzou/Fine_map_susie/output/TR/$tissue
		geneList=$currDir/input/tissue_gene_tr/${tissue}_gene_list.txt
		if [ ! -d "$outDir" ]
		then
			mkdir -p $outDir
		fi
		echo "#!/bin/bash
#SBATCH --job-name=tr_$tissue
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --error=tr_${tissue}.err
#SBATCH --output=tr_${tissue}.out

TISSUE=$tissue
GeneList=$geneList
DIR=$currDir
" > $currDir/submit_prepare_pheno.${tissue}.slurm
		echo '
cd $SLURM_SUBMIT_DIR
Rscript $DIR/src/prepare_pheno_by_gene.residual.tr.R -t $TISSUE -g $GeneList

echo "process end at:"
date
' >> $currDir/submit_prepare_pheno.${tissue}.slurm
		sbatch $currDir/submit_prepare_pheno.${tissue}.slurm
	done
}

function run_prepare_phenotype_by_eGene_2(){
	currDir=`pwd`
	for tissue in `cat $currDir/input/SampleSize_by_tissue.txt|cut -f1`
	do
		echo $tissue
		for gene in `cat $currDir/input/tissue_gene2/${tissue}_gene_list.txt|cut -f2`
		do
			echo $gene
			outdir=$currDir/output/Finemapping_eQTL_sv_snv/${tissue}/${gene}
			if [ ! -d "$outdir" ]
			then
				mkdir -p $outdir
			fi

			Rscript $currDir/src/prepare_pheno_by_gene_2.R -t $tissue -g $gene
		done
	done
}

function generate_vcf_by_gene(){
	currDir=`pwd`
	wkdir=/flashfs1/scratch.global/xdzou/Fine_map_susie
	outdir=$wkdir/input/SNV_by_egenes
	mkdir -p $outdir

	for f in `ls $wkdir/input/task/remained_sv_egenes.75.txt`
	do
		task=`basename $f`
		echo $task
		echo "#!/bin/bash
#SBATCH --job-name=extract_gt_$task
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --error=${task}.err
#SBATCH --output=${task}.out

OUTDIR=$outdir
TASK=$task
DIR=$wkdir

" > $currDir/submit_generate_vcf.${task}.slurm
		echo $'
TSS=$DIR/input/TSS.GTOP_tested_genes.sorted.bed
GTcoord=$DIR/input/GTOP_SV.coord.sorted.bed

#module load bedtools2/2.29.2


#for gene in `cat $DIR/input/task/$TASK`
#do
#	echo $gene
#	mkdir -p $OUTDIR/$gene
#	cat $TSS|grep $gene|awk \'BEGIN{OFS="\\t"}{print $1,$2,$3,$4}\'|bedtools window -w 1000000 -a stdin -b $GTcoord|cut -f8 > $OUTDIR/$gene/sv_list.txt
#done

Rscript $DIR/src/prepare_genotype_by_gene.by_task.R -t $TASK
' >>  $currDir/submit_generate_vcf.${task}.slurm
		sbatch  $currDir/submit_generate_vcf.${task}.slurm
	done
}

function add_genotype_to_each_gene(){
	wkdir=/flashfs1/scratch.global/xdzou/Fine_map_susie
	currDir=`pwd`
	for tissue in `cat $currDir/selected_tissues_11.txt|cut -f1`
	do
		echo $tissue
		outDir=/flashfs1/scratch.global/xdzou/Fine_map_susie/output/SV/$tissue
		geneList=$currDir/input/tissue_gene_sv/${tissue}_gene_list.txt
		for gene in `cat $geneList|cut -f2`
		do
#			if [ -f $wkdir/input/variants_by_tr_egene/$gene/eVar_GT.vcf ]
#			then
#				cp $wkdir/input/variants_by_tr_egene/$gene/eVar_GT.vcf $outDir/$gene/
#			else
#				echo "eVar_GT.vcf not exist for $gene"
#			fi
			if [ ! -f $outDir/$gene/eSV_GT.vcf ]
			then
				cp $wkdir/input/SNV_by_egenes/$gene/eSV_GT.vcf $outDir/$gene/
			fi
		done
	done
}

function run_prepare_genotype_by_eGene(){
	currDir=`pwd`
	for tissue in `cat $currDir/SampleSize_by_tissue.txt|cut -f1`
	do
		echo $tissue
		outDir=/flashfs1/scratch.global/xdzou/Fine_map_susie/output/$tissue
		geneList=/flashfs1/scratch.global/xdzou/Fine_map_susie/input/tissue_gene_snv/${tissue}_gene_list.txt
		if [ ! -d "$outDir" ]
		then
			mkdir -p $outDir
		fi
		echo "#!/bin/bash
#SBATCH --job-name=pre_SNV_$tissue
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --error=${tissue}_SNV.err
#SBATCH --output=${tissue}_SNV.out

TISSUE=$tissue
GeneList=$geneList
DIR=/flashfs1/scratch.global/xdzou/Fine_map_susie
OUTdir=$outDir
" > $currDir/submit_prepare_SNV.${tissue}.slurm
		echo $'

TSS=$DIR/input/TSS.GTOP_tested_genes.sorted.bed
GTcoord=$DIR/input/GTOP_SNV.coord.sorted.bed

module load bedtools2/2.29.2
echo "get snp list..."
for gene in `cat $GeneList|cut -f2`
do
	echo $gene
	mkdir -p $OUTdir/$gene
	cat $TSS|grep $gene|awk \'BEGIN{OFS="\\t"}{print $1,$2,$3,$4}\'|bedtools window -w 1000000 -a stdin -b $GTcoord|cut -f8 > $OUTdir/$gene/snp_list.txt
done

sleep 30
echo "generate genotype matrix for each gene..."
Rscript $DIR/src/prepare_genotype_by_gene.R -t $TISSUE

' >> $currDir/submit_prepare_SNV.${tissue}.slurm
		sbatch $currDir/submit_prepare_SNV.${tissue}.slurm
	done

}
# susie
function run_susie_analysis(){
	currDir=`pwd`
	wkdir=/flashfs1/scratch.global/xdzou/Fine_map_susie
	for tissue in `cat $currDir/selected_tissues_11.txt|cut -f1`
	do
		echo $tissue
		outDir=$wkdir/output/TR/$tissue
#		geneList=$currDir/input/tissue_gene_tr/${tissue}_gene_list.txt
		geneList=$currDir/input/${tissue}.remained_genes.TR.txt
		if [ ! -d "$outDir" ]
		then
			mkdir -p $outDir
		fi
		echo "#!/bin/bash
#SBATCH --job-name=tr_$tissue
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --error=tr_${tissue}.err
#SBATCH --output=tr_${tissue}.out
" > $currDir/submit_susieR_${tissue}.slurm
		echo "

TISSUE=$tissue
DIR=$wkdir
OUTdir=$outDir
GENElist=$geneList
" >> $currDir/submit_susieR_${tissue}.slurm
		echo $'
for gene in `cat $GENElist`
do
	echo $gene
	if [ ! -f "$OUTdir/$gene/eTR_GT.SuSiE.rds" ]
	then
		Rscript $DIR/src/finemapping_tr.R $OUTdir/$gene 10 0.2 0 
	fi
done

echo "+++++++++++++++++++++"
echo "process will end at:"
date
' >> $currDir/submit_susieR_${tissue}.slurm 
		sbatch $currDir/submit_susieR_${tissue}.slurm
	done
}

# get the eGene list in all tissues, permutation FDR<0.05
function get_tissue_gene_list(){
	dir=/lustre/home/xdzou/2024-10-21-GTBMap/2026-02-02-GTOP_eQTL_mapping/output/QTL_mapping
	eGenes=GTOP.tr_egenes.all_tissues.FDR.05.txt
	currdir=`pwd`

	if [ ! -d "$currdir/input/tissue_gene_tr" ]
	then
		mkdir -p $currdir/input/tissue_gene_tr
	fi

	for tissue in `cat $currdir/selected_tissues_11.txt|cut -f1`
	do
		echo $tissue
		cat $dir/$eGenes|tail -n+2|awk -v TISSUE=$tissue '{if($8 == TISSUE) print TISSUE"\t"$1}'|sort|uniq > ${currdir}/input/tissue_gene_tr/${tissue}_gene_list.txt &
		wait
		wc -l $currdir/input/tissue_gene_tr/${tissue}_gene_list.txt
	done
}

main
