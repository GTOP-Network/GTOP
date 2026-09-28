#!/bin/bash

main(){

	s1_run_rnaseqc
	s2_merge_rnaseqc
	s3__generate_expMat_gct_by_tissue
	s4_prepare_phenotype
	s5_generate_phenotype_bed
}
rawbam_dir=$curr_dir/input/bam
curr_dir=$(pwd)

function s1_run_rnaseqc(){
	curr_dir=`pwd`
	outdir=$curr_dir/output/RNASeQC_out

	if [ ! -d "$outdir" ]
	then
		mkdir -p $outdir
	fi
# generate gene models for RNASeQC2 quantification
	python ./gtex-pipeline-master/gene_model/collapse_annotation.py --collapse_only ./input/gtop_gencode47.gtf ./output/gtop_gencode47.genes.gtf

	for f in `ls $curr_dir/task/task_rnaseqc_*`
	do
		task=`basename $f`
		echo $task
		echo "#!/bin/bash" > ${curr_dir}/submit_rnaseqc_${task}.slurm
		echo "
#SBATCH --job-name=$task
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=rnaseqc_${task}.err
#SBATCH --output=rnaseqc_${task}.out

DIR=$curr_dir
BamDir=$rawbam_dir
OutDir=$outdir
TASK=$task

" >> ${curr_dir}/submit_rnaseqc_${task}.slurm
		echo $'
echo "+++++++++++++++++++++"
echo "process start at:"
date

cd $SLURM_SUBMIT_DIR

for SAMP in `cat $DIR/task/$TASK`
do
	echo $SAMP
	rnaseqc --sample=$SAMP --gene-length=100 --coverage --stranded="RF" $DIR/input/gtop_gencode47.genes.gtf $BamDir/${SAMP}.Aligned.sortedByCoord.out.bam $OutDir

done

echo "++++++++++++++++++++"
echo "process end at:"
date
' >> ${curr_dir}/submit_rnaseqc_${task}.slurm
		sbatch ${curr_dir}/submit_rnaseqc_${task}.slurm
	done
}

function s2_merge_rnaseqc(){
	RNASeQCdir=$base_dir/output/RNASeQC_out
	curr_dir=`pwd`
	Rscript $curr_dir/script/merge_rnaseqc.R -p $RNASeQCdir
}


function s3_generate_expMat_by_tissue(){
	curr_dir=`pwd`
	mkdir -p $curr_dir/output/reads_gct
	mkdir -p $curr_dir/output/tpm_gct

	Rscript ./script/prepare_exp_matrix_by_tissue.R
}

function s4_prepare_phenotype(){
	curr_dir=`pwd`
	outdir=$curr_dir/output/phenotype

	mkdir -p $outdir

	for tissue in `cat $curr_dir/input/selected_tissues_for_qtl.txt`
	do
		echo $tissue
		echo "#!/bin/bash" > $curr_dir/submit_prepare_phenotype_${tissue}.slurm
		echo "
#SBATCH --job-name=$tissue
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=pre_pheno.${tissue}.err
#SBATCH --output=pre_pheno.${tissue}.out

DIR=$curr_dir
TISSUE=$tissue
OutDir=$outdir
" >> $curr_dir/submit_prepare_phenotype_${tissue}.slurm
		echo $'
python ./gtex-pipeline-master/qtl/src/eqtl_prepare_expression.py $DIR/output/tpm_gct/${TISSUE}.tpm.gct $DIR/output/reads_gct/${TISSUE}.reads_count.gct --tpm_threshold 0.1 --count_threshold 6 --sample_frac_threshold 0.2 --normalization_method tmm --output $OutDir/${TISSUE}.phenotype.txt

' >> $curr_dir/submit_prepare_phenotype_${tissue}.slurm
		sbatch $curr_dir/submit_prepare_phenotype_${tissue}.slurm
	done
}

function s5_generate_phenotype_bed(){
	curr_dir=`pwd`
	outdir=$curr_dir/output/phenotype

	mkdir -p $outdir

	for f in `ls $outdir/*.txt`
	do
		tissue=`basename $f .phenotype.txt`
		echo $tissue
		Rscript $curr_dir/scripts/prepare_phenotype.bed.R -t $tissue
		cat $outdir/${tissue}.bed|sort -k1,1 -k2,2n > $curr_dir/output/phenotype/${tissue}.phenotype.bed
	done
}

# >> main
main
