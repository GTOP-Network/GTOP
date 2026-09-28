#!/bin/bash
# main function

main(){
	get_tissue_gene_list
	run_prepare_phenotype_by_eGene
	add_genotype_to_each_gene
	run_susie_analysis
	run_summarize_susie
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
		outDir=/path/to/Fine_map_susie/output/Joint/${tissue}
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
	wkdir=/path/to/Fine_map_susie
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


# susie
function run_susie_analysis(){
	currDir=`pwd`
	wkdir=/path/to/Fine_map_susie
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



main
