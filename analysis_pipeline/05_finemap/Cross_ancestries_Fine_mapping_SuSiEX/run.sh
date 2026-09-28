#!/bin/bash

main(){
	generate_gene_list
	extract_asso_by_gene
	extract_asso_by_gene_gtex
	run_format_gtop_asso
	run_susieX
	summarize_susieX
}

# global variables
WKDIR=`pwd`


function generate_gene_list(){
	Tissue=$1
	annotation=$WKDIR/input/GTOP.gene_info.txt

	Rscript $WKDIR/src/generate_egene_list.R -a $annotation -t $Tissue
}

function run_susieX(){
	Refdir=$WKDIR/input/reference
	tissue=Adipose
	susiedir=$WKDIR/output/${tissue}/SuSiEx
	outdir=$WKDIR/output/$tissue/SuSiEx/results
	mkdir -p $outdir
	for f in `ls $WKDIR/task/task_*`
	do
		task=`basename $f`
		echo $task
		echo "#!/bin/bash
#SBATCH --job-name=susieX_${tissue}_$task
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error=susieX_${tissue}.${task}.err
#SBATCH --output=susieX_${tissue}.${task}.out


REFDIR=$Refdir
TISSUE=$tissue
OUTDIR=$outdir
TASK=$task
basedir=$WKDIR
SuSiExDIR=$susiedir
" > $WKDIR/submit_susieX.${tissue}.${task}.slurm
		echo $'

echo "start at:"
date
#cd $SLURM_SUBMIT_DIR

cd $SuSiExDIR/examples

while read line
do
	GENE=`echo $line|awk \'{print $1}\'`
	CHROM=`echo $line|awk \'{print $2}\'`
	CHROM=${CHROM#chr}
	Nsig=`cat GTEX_${GENE}.sumstat.txt|wc -l`
	TSS=`echo $line|awk \'{print $3}\'`
	START=`expr $TSS - 1000000`
	END_p=`expr $TSS + 1000000`
	if [ $START -lt 0 ]
	then
		START=0
	fi

	if [ $Nsig -gt 5 ]
	then
		../bin/SuSiEx --sst_file=GTOP_${GENE}.sumstat.txt,GTEX_${GENE}.sumstat.txt \
			--n_gwas=122,258 \
			--ref_file=$REFDIR/GTOP,$REFDIR/GTEx \
			--ld_file=${GENE}_GTOP,${GENE}_GTEx \
			--out_dir=$OUTDIR \
			--out_name=susiex_${GENE} \
			--chr=$CHROM \
			--bp=${START},${END_p} \
			--chr_col=1,1 \
			--snp_col=5,5 \
			--bp_col=2,2 \
			--a1_col=4,4 \
			--a2_col=3,3 \
			--eff_col=7,7 \
			--se_col=8,8 \
			--pval_col=6,6 \
			--plink=../utilities/plink \
			--threads=8
	fi


done < $basedir/task_blood/$TASK

echo "end at:"
date
' >> $WKDIR/submit_susieX.${tissue}.${task}.slurm
		sbatch $WKDIR/submit_susieX.${tissue}.${task}.slurm
	done
}

function run_format_gtop_asso(){
	Tissue=Adipose
	Reference=$WKDIR/input/reference/GTOP.SNP_list.txt

	echo "#!/bin/bash
#SBATCH --job-name=format_gtop_asso.${Tissue}
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --error=${Tissue}.format.err
#SBATCH --output=${Tissue}.format.out


TISSUE=$Tissue
REF=$Reference
baseDIR=$WKDIR
" > $WKDIR/submit_format_${Tissue}.slurm
	echo '
echo "start at:"
date

cd $SLURM_SUBMIT_DIR
Rscript $baseDIR/src/add_snp_meta.R -r $REF -t $TISSUE

echo "end at:"
date
' >> $WKDIR/submit_format_${Tissue}.slurm
	sbatch $WKDIR/submit_format_${Tissue}.slurm
}

function extract_asso_by_gene(){
	eqtlDIR=$WKDIR/input/eQTL_GTOP
	tissue=Adipose
	outdir=$WKDIR/output/${tissue}/SuSiEx/examples
	mkdir -p $WKDIR/output/${tissue}

	if [ ! -d "$outdir" ]
	then
		cp -r $WKDIR/output/SuSiEx $WKDIR/output/$tissue
	fi

	for f in `ls $WKDIR/task/task_*`
	do
		task=`basename $f`
		echo $task
		echo "#!/bin/bash
#SBATCH --job-name=${tissue}_$task
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --error=${tissue}_${task}.err
#SBATCH --output=${tissue}_${task}.out


baseDir=$WKDIR
QTLdir=$eqtlDIR
TISSUE=$tissue
OUTdir=$outdir
TASK=$task
" > $WKDIR/submit_${tissue}_${task}.slurm
		echo $'

echo "task start:"
date

cd $SLURM_SUBMIT_DIR

for gene in `cat $baseDir/task/$TASK|tail -n+2|cut -f1`
do
	echo $gene
	genePRE=${gene%.*}
	echo -e "SNP\\tpvalue\\tbeta\\tse" > $OUTdir/GTOP_${genePRE}.txt
	zcat $QTLdir/cis_QTL/text_format/${TISSUE}.cis_eQTL.all_pairs.txt.gz|grep "$gene" |awk \'{print $2"\\t"$7"\\t"$8"\\t"$9}\' >> $OUTdir/GTOP_${genePRE}.txt
done
' >> $WKDIR/submit_${tissue}_${task}.slurm
		sbatch $WKDIR/submit_${tissue}_${task}.slurm
	done
}


function extract_asso_by_gene_gtex(){
	eqtlDIR=$WKDIR/input/eQTL_GTEx
	tissue=Adipose
	outdir=$WKDIR/output/${tissue}/SuSiEx/examples
	if [ ! -d "$outdir" ]
	then
		cp -r $WKDIR/output/SuSiEx $WKDIR/output/$tissue
	fi

	for f in `ls $WKDIR/task/task_*`
	do
		task=`basename $f`
		echo $task
		echo "#!/bin/bash
#SBATCH --job-name=gtex_${tissue}_$task
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --error=gtex_${tissue}_${task}.err
#SBATCH --output=gtex_${tissue}_${task}.out


baseDir=$WKDIR
QTLdir=$eqtlDIR
TISSUE=$tissue
OUTdir=$outdir
TASK=$task
" > $WKDIR/submit_${tissue}_${task}.gtex.slurm
		echo $'

echo "task start:"
date

cd $SLURM_SUBMIT_DIR

for gene in `cat $baseDir/task/$TASK|tail -n+2|cut -f1`
do
	echo $gene
	genePRE=${gene%.*}
	echo -e "chrom\\tbp\\tA2\\tA1\\tSNP\\tpvalue\\tbeta\\tse" > $OUTdir/GTEX_${genePRE}.sumstat.txt
	zcat $QTLdir/${TISSUE}.allpairs_addrsid.txt.gz|grep "$genePRE"|cut -f1-5,12,13,14|sed 's/^chr//'|grep -v "nan" >> $OUTdir/GTEX_${genePRE}.sumstat.txt
done


echo "task end:"
date
' >> $WKDIR/submit_${tissue}_${task}.gtex.slurm
		sbatch $WKDIR/submit_${tissue}_${task}.gtex.slurm
	done
}
function summarize_susieX(){
	WKDIR=`pwd`
	tissue_list=( Adipose Adrenal_Gland Liver Muscle Pancreas_Body Spleen Skin Whole_Blood )
	for Tissue in ${tissue_list[@]}
	do
		echo $Tissue
		bash $WKDIR/src/summarize_susiex_res.sh $WKDIR/output/$Tissue/SuSiEx/results $Tissue
	done

}


main
