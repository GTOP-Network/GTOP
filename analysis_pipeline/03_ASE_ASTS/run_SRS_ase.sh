#!/bin/bash

# Bash Script Pipeline for SRS ASE analysis

# Description:
# This script processes short-read sequencing data for ASE analysis



# Example Usage:
# ./script.sh <pwd> <WASP_PASS_BAM_dir> <PHASE_VCF_dir> <gtfref> <genomeref> <genebed>


pwd=$1
WASP_PASS_BAM_dir=$2
PHASE_VCF_dir=$3
gtfref=$4
genomeref=$5
genebed=$6


SNP_ASE_dir=$pwd/SNP_ASE 
mkdir -p $SNP_ASE_dir


main(){

 ase_snp_level #Uses GATK to calculate reference and alternate allele counts at heterozygous sites for each sample
 ase_snp_cal_lamp #Calculates the global foreign allele frequency (lamp value) per individual, based on GTEx's ase_calculate_lamp.py
 ase_snp_sum #Integrates ASE data for each donor, performs statistical tests, and applies quality filters, inspired by GTEx's

}

ase_snp_level(){
slurm=$pwd/slurm/ase_snp_level

mkdir -p $slurm
cd $SNP_ASE_dir

cd $pwd
for i in `ls ${WASP_PASS_BAM_dir}/*.Aligned.sortedByCoord.WASPpass.out.bam`;do
sample=`basename $i .Aligned.sortedByCoord.WASPpass.out.bam`
vcffilename=`echo $sample| cut -d'-' -f1-2`

chrlist=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

for chr in ${chrlist[@]};do


echo '#!/bin/bash
#SBATCH --job-name='${sample}'.chr'$chr'_gatk
#SBATCH --partition=cu-debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error='${sample}'.chr'$chr'.err
#SBATCH --output='${sample}'.chr'$chr'.out
##########################################

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
vcffilename=`echo '$sample'| cut -d'-' -f1-2`
echo $vcffilename
module load gatk/4.1.0.0-gcc-9.2.0


tabix -p vcf -f '$PHASE_VCF_dir'/${vcffilename}_het_chr'${chr}'.vcf.gz

gatk --java-options "-Xmx24g -Xms12g -XX:ParallelGCThreads=8" ASEReadCounter -R '$genomeref' -L chr'${chr}' --min-base-quality 10 -min-depth 1  --min-mapping-quality 255 -I '${WASP_PASS_BAM_dir}'/'$sample'.Aligned.sortedByCoord.WASPpass.out.bam -V '$PHASE_VCF_dir'/${vcffilename}_het_chr'${chr}'.vcf.gz -O '$SNP_ASE_dir'/'$sample'.chr'$chr'.csv



module unload gatk/4.1.0.0-gcc-9.2.0


echo "processs will sleep 30s"
sleep 5
echo "process end at : "
date
#module unload biosoft/torus
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
wait' > $slurm/$sample.chr$chr.slurm

wait
cd ${slurm}
sbatch $sample.chr$chr.slurm
echo $sample $chr" is finished"
done
done

}


ase_snp_cal_lamp(){
	slurm=$pwd/slurm/cal_lamp
	lampDir=$pwd/SNP_ASE_lamp/

	mkdir -p $slurm $lampDir
	echo "#!/bin/bash
#SBATCH --job-name=cal_ase_lamp
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error=cal_ase_lamp.err
#SBATCH --output=cal_ase_lamp.out

echo "+++++++++++++++++++++"
echo "process start at:"
date


# prepare inputlist
ls $SNP_ASE_dir/*.csv > $pwd/SNP_ASE.txt

# cal lamp
python $pwd/bin/ase_calculate_lamp.py $pwd/SNP_ASE.txt my_studyGTOP -o $lampDir

echo "++++++++++++++++++++"
echo "process end at:"
date
" > $slurm/cal_ase_lamp.slurm

cd $slurm
sbatch cal_ase_lamp.slurm

}


ase_snp_sum(){

	slurm=$pwd/slurm/ase_snp_sum

	outpath=$pwd/SNP_ASE_sum/

	mkdir -p $slurm $outpath

	# #combine ase res for per sample
	ls $SNP_ASE_dir|cut -d. -f1|sort|uniq| \
	while read line;
	do
	files=(`realpath $SNP_ASE_dir/*.csv|grep chr|grep $line|sort -k1,1V`)
	cat ${files[@]}|grep -v contig|sed '1icontig\tposition\tvariantID\trefAllele\taltAllele\trefCount\taltCount\ttotalCount\tlowMAPQDepth\tlowBaseQDepth\trawDepth\totherBases\timproperPairs' > $SNP_ASE_dir/$line.all.tsv
	done
	## pre list output input/SNP_ASE.files.frosum.list
	Rscript bin/pre_SNP_ASE_listforsum.r


	
	ls $SNP_ASE_dir|cut -d. -f1|cut -d- -f1-2|sort|uniq|\
	while read ind;do


	echo "#!/bin/bash
#SBATCH --job-name=ase_sum_$ind
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error=$ind.err
#SBATCH --output=$ind.out

echo "+++++++++++++++++++++"
echo "process start at:"
date

python $pwd/bin/ase_aggregate_by_individual.py \
  $pwd/input/SNP_ASE.files.list \
  $PHASE_VCF_dir/${ind}_het.vcf.gz \
  $genebed \
  $pwd/input/wgEncodeCrgMapabilityAlign100mer.hg38.bigWig.bw \
  $pwd/SNP_ASE_lamp/my_study.lamp_values.txt  \
  $ind \
  -o $outpath


echo "++++++++++++++++++++"
echo "process end at:"
date
" > $slurm/$ind.slurm

cd $slurm
sbatch $ind.slurm
echo $ind " is finished!!!!"
done


}







main
