#!/bin/bash


main(){

sort_pregroup
run_tensorQTL

}


sort_pregroup(){
	basedir=`pwd`
	slurm=$basedir/slurm/$inputtype/sort_pregroup
	mkdir -p $slurm


	for covfile in `ls $basedir/output/covariates/*.covariates.txt`
	do
		tissue=`basename $covfile .covariates.txt`

	less $basedir/input/ori_phe/${tissue}/leafcutter2.filter.leafcutter.bed.gz|awk 'NR > 1 {split($4, a, ":"); print $4 "\t" a[6]}' > $basedir/input/ori_phe/${tissue}/group.list


		# tissue="Adipose"
		echo $tissue
		echo "#!/bin/bash" > ${slurm}/${tissue}.slurm
		echo "
#SBATCH --job-name=sorted_${tissue}
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error=sorted_${tissue}.err
#SBATCH --output=sorted_${tissue}.out

DIR=$basedir
OutDir=$outdir
Tissue=$tissue
		" >> ${slurm}/${tissue}.slurm
		echo $'

echo "process start at:"
date

Rscript ./01_data_preparation/scripts/sorted.group.r -p leafcutter2.filter.leafcutter.bed.gz -g group.list -w $DIR/input/ori_phe/$Tissue


echo "process end at:"
date
' >> ${slurm}/${tissue}.slurm
cd $slurm
		sbatch ${slurm}/${tissue}.slurm
	done
}


function run_tensorQTL(){
	basedir=`pwd`
	outdir=$basedir/output/nominal
	per=$basedir/output/permutation
	slurm=$basedir/slurm/tensorqtl_cis_${type}


		mkdir -p $outdir $slurm $per

	for covfile in `ls $basedir/output/covariates/*.covariates.txt`
	do
		tissue=`basename $covfile .covariates.txt`

		echo $tissue
		echo "#!/bin/bash" > ${slurm}/submit_tensorqtl_${tissue}.slurm
		echo "
#SBATCH --job-name=ts_${tissue}_${type}
#SBATCH --partition=gpu-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --error=${tissue}.err
#SBATCH --output=${tissue}.out

DIR=$basedir
OutDir=$outdir
Tissue=$tissue
		" >> ${slurm}/submit_tensorqtl_${tissue}.slurm
		echo $'

echo "process start at:"
date




# run nominal
python3 -m tensorqtl $DIR/output/genotype/GTOP.GT $DIR/input/ori_phe/${Tissue}/phenotype.sorted.bed $OutDir/${Tissue} \
 	--covariates $DIR/output/covariates/${Tissue}.covariates.txt \
 	--mode cis_nominal


# run permutated
python3 -m tensorqtl $DIR/output/genotype/GTOP.GT $DIR/input/ori_phe/${Tissue}/phenotype.sorted.bed '$per'/${Tissue} \
  	--covariates $DIR/output/covariates/${Tissue}.covariates.txt \
	--phenotype_groups $DIR/input/ori_phe/${Tissue}/phenotype_groups.sorted.txt \
  	--mode cis --chunk_size chr





echo "process end at:"
date
' >> ${slurm}/submit_tensorqtl_${tissue}.slurm
cd $slurm
		sbatch submit_tensorqtl_${tissue}.slurm
	done
}



main
