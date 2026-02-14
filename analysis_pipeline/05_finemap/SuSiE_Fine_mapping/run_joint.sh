#!/bin/bash
# main function

main(){
	get_tissue_gene_list
	run_prepare_phenotype_by_eGene
	add_genotype_to_each_gene
	run_susie_analysis
	run_summarize_susie
}

function get_var_by_gene(){
	currDir=`pwd`
	for f in `ls $currDir/output/task_04_*`
	do
		task=`basename $f`
		echo $task
		echo "#!/bin/bash
#SBATCH --job-name=$task
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --error=${task}.err
#SBATCH --output=${task}.out


module load bedtools2/2.29.2
TASK=$task
DIR=$currDir
TSS=$currDir/input/TSS.GTOP_tested_genes.sorted.bed
GTcoord=$currDir/input/GTOP_SNV.coord.sorted.bed
GeneList=$currDir/output/$task
" > $currDir/submit_get_SNP_by_gene.${task}.slurm
		echo $'

for gene in `cat $GeneList`
do
	echo $gene
	cat $TSS|grep $gene|awk \'BEGIN{OFS="\\t"}{print $1,$2,$3,$4}\'|bedtools window -w 1000000 -a stdin -b $GTcoord|cut -f8|awk -v GENE=$gene \'{print $0"\\t"GENE}\' > $DIR/output/Finemapping_eQTL_tr/ALL_gene_1Mb_SNV.tr_${TASK}.txt
done

echo "++++++++++++++++"
date
' >> $currDir/submit_get_SNP_by_gene.${task}.slurm
		sbatch $currDir/submit_get_SNP_by_gene.${task}.slurm
	done
}

# summarize susie
function run_summarize_susie(){
	currDir=`pwd`
	for tissue in `cat $currDir/selected_tissues_11.txt|cut -f1`
	do
		echo $tissue
		echo "#!/bin/bash
#SBATCH --job-name=joint_$tissue
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --error=joint_${tissue}.err
#SBATCH --output=joint_${tissue}.out

TISSUE=$tissue
DIR=$currDir
WKDIR=/flashfs1/scratch.global/xdzou/Fine_map_susie
" > $currDir/submit_summarize_susie.${tissue}.slurm
		echo '
cd $SLURM_SUBMIT_DIR

Rscript $DIR/src/summarize_susie_results_by_tissue.R -d $WKDIR -t $TISSUE -g $DIR/input/tissue_gene_joint/${TISSUE}_gene_list.txt

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
		outDir=/flashfs1/scratch.global/xdzou/Fine_map_susie/output/Joint/${tissue}
		geneList=$currDir/input/tissue_gene_joint/${tissue}_gene_list.txt
		if [ ! -d "$outDir" ]
		then
			mkdir -p $outDir
		fi
		echo "#!/bin/bash
#SBATCH --job-name=pre_pheno_$tissue
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --error=pre_pheno_${tissue}.err
#SBATCH --output=pre_pheno_${tissue}.out

TISSUE=$tissue
GeneList=$geneList
DIR=$currDir
" > $currDir/submit_prepare_pheno.${tissue}.slurm
			echo '
cd $SLURM_SUBMIT_DIR
Rscript $DIR/src/prepare_pheno_by_gene.residual.joint.R -t $TISSUE -g $GeneList

echo "process end at:"
date
' >> $currDir/submit_prepare_pheno.${tissue}.slurm
			sbatch $currDir/submit_prepare_pheno.${tissue}.slurm
	done
}

function add_genotype_to_each_gene(){
	wkdir=/flashfs1/scratch.global/xdzou/Fine_map_susie
	currDir=`pwd`
	for tissue in `cat $currDir/selected_tissues_11.txt|cut -f1|tail -n+3`
	do
		echo $tissue
		echo "#!/bin/bash
#SBATCH --job-name=add_gt_${tissue}
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --error=add_gt_${tissue}.err
#SBATCH --output=add_gt_${tissue}.out

TISSUE=$tissue
GeneList=$currDir/input/tissue_gene_joint/${tissue}_gene_list.txt
DIR=$currDir
" > $currDir/submit_add_genotype.${tissue}.slurm
		echo $'

for gene in `cat $GeneList|cut -f2`
do
	echo $gene
	Rscript $DIR/src/generate_joint_GT_by_gene.R -g $gene -t $TISSUE
done
' >> $currDir/submit_add_genotype.${tissue}.slurm
		sbatch $currDir/submit_add_genotype.${tissue}.slurm
	done
}

function run_prepare_genotype_by_eGene(){
	currDir=`pwd`
	tissue_list=( Spleen Whole_Blood )
#	tissue_list=( Adrenal_Gland Gallbladder Liver Muscle Pancreas_Body Pancreas_Head Pancreas_Tail Spleen Whole_Blood )
	for tissue in ${tissue_list[@]}
	do
		echo $tissue
		for f in `ls $currDir/input/tissue_gene_joint/$tissue/${tissue}.task_*`
		do
			task=`basename $f`
#			tissue=${task%.*}
#			echo $tissue
			outDir=$currDir/output/Finemapping_eQTL_joint/$tissue
			geneList=$currDir/input/tissue_gene_joint/${tissue}/$task
			if [ ! -d "$outDir" ]
			then
				mkdir -p $outDir
			fi
			echo "#!/bin/bash
#SBATCH --job-name=pre_gt_$task
#SBATCH --partition=fat-1
#SBATCH --nodes=1
#SBATCH --error=${task}_gt.err
#SBATCH --output=${task}_gt.out

TISSUE=$tissue
GeneList=$geneList
DIR=$currDir
OUTdir=$outDir
" > $currDir/submit_prepare_genotype.${task}.slurm
			echo $'

TSS=$DIR/input/TSS.GTOP_tested_genes.sorted.bed
GTcoord=$DIR/input/GTOP_all_var_joint.coord_sorted.bed

module load bedtools2/2.29.2
echo "get snp list..."
for gene in `cat $GeneList|cut -f2`
do
	echo $gene
	mkdir -p $OUTdir/$gene
	cat $TSS|grep $gene|awk \'BEGIN{OFS="\\t"}{print $1,$2,$3,$4}\'|bedtools window -w 1000000 -a stdin -b $GTcoord|cut -f8 > $OUTdir/$gene/snp_list.txt
done

echo "generate genotype matrix for each gene..."
Rscript $DIR/src/prepare_genotype_by_gene.tr.R -g $GeneList -t $TISSUE

' >> $currDir/submit_prepare_genotype.${task}.slurm
			sbatch $currDir/submit_prepare_genotype.${task}.slurm
		done
	done
}
# susie
function run_susie_analysis(){
	currDir=`pwd`
	wkdir=/flashfs1/scratch.global/xdzou/Fine_map_susie
	for tissue in `cat $currDir/selected_tissues_11.txt|cut -f1`
	do
		echo $tissue
		outDir=$wkdir/output/Joint/$tissue
		geneList=$currDir/input/tissue_gene_joint/${tissue}_gene_list.txt
		if [ ! -d "$outDir" ]
		then
			mkdir -p $outDir
		fi
		echo "#!/bin/bash
#SBATCH --job-name=susie_$tissue
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --error=susie_${tissue}.err
#SBATCH --output=susie_${tissue}.out
" > $currDir/submit_susieR_${tissue}.slurm
		echo "
TISSUE=$tissue
DIR=$wkdir
OUTdir=$outDir
GENElist=$geneList
" >> $currDir/submit_susieR_${tissue}.slurm
		echo $'
for gene in `cat $GENElist|cut -f2`
do
	echo $gene
	Rscript $DIR/src/finemapping.R $OUTdir/$gene 10 0.2 0 
done

echo "+++++++++++++++++++++"
echo "process will end at:"
date
' >> $currDir/submit_susieR_${tissue}.slurm 
		sbatch $currDir/submit_susieR_${tissue}.slurm
	done
}

# get the atssGene list in all tissues, here we use conditional pass QTLs
function get_tissue_gene_list(){
	currdir=`pwd`
	if [ ! -d "$currdir/input/tissue_gene_joint" ]
	then
		mkdir -p $currdir/input/tissue_gene_joint
	fi
	bash $currdir/src/split_tissue_egenes.sh
}

function split_files(){
	curr_dir=`pwd`
	dir=/lustre/home/xdzou/2024-10-21-GTBMap/2025-06-15-joint-split-gt-QTL/output/QTL_mapping/all_joint/cis_QTL/text_format

	if [ ! -d "$curr_dir/input/split_gene" ]
	then
		mkdir -p $curr_dir/input/split_gene/
	fi

	for f in `ls $dir/*.cis_eQTL.all_pairs.add_tstat.txt`
	do
		tissue=`basename $f .cis_eQTL.all_pairs.add_tstat.txt`
		echo $tissue
		if [ ! -d "$curr_dir/input/split_gene/$tissue" ]
		then
			mkdir -p $curr_dir/input/split_gene/$tissue
		fi
		echo "#!/bin/bash" > ${curr_dir}/submit_split_files.${tissue}.slurm
		echo "
#SBATCH --job-name=splitFiles_$tissue
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=splitFile_${tissue}.err
#SBATCH --output=splitFile_${tissue}.out

baseDir=$curr_dir
TISSUE=$tissue
">> ${curr_dir}/submit_split_files.${tissue}.slurm
		echo $' 
/lustre/home/xdzou/src/anaconda2/bin/python $baseDir/src/Split_file.py $TISSUE

echo "+++++++++++++++++++++++++++++"
echo "process will end at:"
date
' >> ${curr_dir}/submit_split_files.${tissue}.slurm &
	wait
	sbatch ${curr_dir}/submit_split_files.${tissue}.slurm
	done
}


# split qtl by gene
function run_split_qtl_by_gene(){
	dir=/lustre/home/xdzou/2022-08-05-altTSS_QTL-Project
	current_dir=`pwd`
	for line in `cat ${dir}/input/gtex_tissues_and_sampleSize.txt`
	do
		tissue=${line%:*}
		if [ ! -d "$dir/input/atssQTL_byGene/$tissue" ]
		then
			mkdir -p $dir/input/atssQTL_byGene/$tissue
		fi

		echo "#!/bin/bash" > ${current_dir}/submit_splitQTL_${tissue}.slurm
		echo "
#SBATCH --job-name=$tissue
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=%j.err 
#SBATCH --output=%j.out " >> ${current_dir}/submit_splitQTL_${tissue}.slurm

		echo $'
echo "process will start at:"
date
echo "+++++++++++++++++++++++++++++"

cd $SLURM_SUBMIT_DIR' >> ${current_dir}/submit_splitQTL_${tissue}.slurm

		echo "
baseDir=$dir
TISSUE=$tissue" >> ${current_dir}/submit_splitQTL_${tissue}.slurm

		echo $'
python $baseDir/src/split_qtl_by_gene.py --qtl_file $baseDir/input/atssQTL/Cis_atssQTL_all.${TISSUE}.formatted.txt --outprefix $baseDir/input/atssQTL_byGene/${TISSUE}/

echo "++++++++++++++++++++++++++++++"
echo "process will sleep 10s"
sleep 10
echo "process end at:"
date
' >> ${current_dir}/submit_splitQTL_${tissue}.slurm &
#		wait
#		sbatch ${current_dir}/submit_splitQTL_${tissue}.slurm
	done
}






function generate_vcf_by_gene(){
	currDir=`pwd`
	wkdir=/flashfs1/scratch.global/xdzou/Fine_map_susie
	outdir=$wkdir/input/variants_by_joint_egene

	for f in `ls $wkdir/input/task/task_*`
	do
		task=`basename $f`
		echo $task
		echo "#!/bin/bash
#SBATCH --job-name=extract_gt_$task
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --error=${task}.err
#SBATCH --output=${task}.out

OUTDIR=$outdir
TASK=$task
DIR=$wkdir

" > $currDir/submit_generate_vcf.${task}.slurm
		echo $'
TSS=$DIR/input/TSS.GTOP_tested_genes.sorted.bed
GTcoord=$DIR/input/GTOP_all_var_joint.coord_sorted.bed

#module load bedtools2/2.29.2


#for gene in `cat $DIR/input/task/$TASK`
#do
#	echo $gene
#	mkdir -p $OUTDIR/$gene
#	cat $TSS|grep $gene|awk \'BEGIN{OFS="\\t"}{print $1,$2,$3,$4}\'|bedtools window -w 1000000 -a stdin -b $GTcoord|cut -f8 > $OUTDIR/$gene/snp_list.txt
#done

Rscript $DIR/src/prepare_genotype_by_gene.by_task.R -t $TASK
' >>  $currDir/submit_generate_vcf.${task}.slurm
		sbatch  $currDir/submit_generate_vcf.${task}.slurm
	done

}


main
