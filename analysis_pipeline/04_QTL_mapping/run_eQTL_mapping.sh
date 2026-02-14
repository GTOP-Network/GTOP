#!/bin/bash

main(){
	run_tensorQTL_nominal
	run_tensorQTL_permutation
}


function run_tensorQTL_nominal(){
	basedir=`pwd`
	outdir=$basedir/output/nominal

	if [ ! -d "$outdir" ]
	then
		mkdir -p $outdir
	fi

	for covfile in `ls $basedir/output/covariates/*.covariates.txt`
	do
		tissue=`basename $covfile .covariates.txt`
		echo $tissue
		echo "#!/bin/bash" > ${basedir}/submit_tensorqtl_${tissue}.nominal.slurm
		echo "
#SBATCH --job-name=tensor_${tissue}_nominal
#SBATCH --partition=gpu-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --error=tensor_${tissue}.nominal.err
#SBATCH --output=tensor_${tissue}.nominal.out


DIR=$basedir
OutDir=$outdir
Tissue=$tissue
		" >> ${basedir}/submit_tensorqtl_${tissue}.nominal.slurm
		echo $'

echo "process start at:"
date

cd $SLURM_SUBMIT_DIR

python3 -m tensorqtl $DIR/output/genotype/GTOP.GT $DIR/output/phenotype/${Tissue}.phenotype.bed $OutDir/${Tissue}.nominal \
	--covariates $DIR/output/covariates/${Tissue}.covariates.txt \
	--mode cis_nominal


echo "process end at:"
date
' >> ${basedir}/submit_tensorqtl_${tissue}.nominal.slurm
		sbatch ${basedir}/submit_tensorqtl_${tissue}.nominal.slurm
	done
}


function run_tensorQTL_perm(){
	basedir=`pwd`
	outdir=$basedir/output/QTL_mapping/permutation

	if [ ! -d "$outdir" ]
	then
		mkdir -p $outdir
	fi

	for covfile in `ls $basedir/output/covariates/*.covariates.txt`
	do
		tissue=`basename $covfile .covariates.txt`
		echo $tissue
		echo "#!/bin/bash" > ${basedir}/submit_tensorqtl_${tissue}.perm.slurm
		echo "
#SBATCH --job-name=tensor_${tissue}_perm
#SBATCH --partition=gpu-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --error=tensor_${tissue}.perm.err
#SBATCH --output=tensor_${tissue}.perm.out


DIR=$basedir
OutDir=$outdir
Tissue=$tissue
		" >> ${basedir}/submit_tensorqtl_${tissue}.perm.slurm
		echo $'

echo "process start at:"
date

cd $SLURM_SUBMIT_DIR

python3 -m tensorqtl $DIR/output/genotype/GTOP.GT $DIR/output/phenotype/${Tissue}.phenotype.bed $OutDir/${Tissue}.perm\
	--covariates $DIR/output/covariates/${Tissue}.covariates.txt \
	--mode cis


echo "process end at:"
date
' >> ${basedir}/submit_tensorqtl_${tissue}.perm.slurm
		sbatch ${basedir}/submit_tensorqtl_${tissue}.perm.slurm
	done
}

#main
main
