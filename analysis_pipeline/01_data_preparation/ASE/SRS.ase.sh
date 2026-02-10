main(){
# SNP level ASE
 ase_snp_level
 ase_snp_cal_lamp
 ase_snp_sum
# Haplotype level ASE
 ase_haplotype_level
 ase_haplotype_level_gene
 ase_haplotype_exp_matrix
 exp_matrix_splittissue
 ase_haplotype_cis

}

pwd=/flashfs1/scratch.global/ynwang/ynwang/2024-11-15-GTBMap/2025-02-15-ASE/2025-12-21-SRS-LRSvcf
WASP_PASS_BAM_dir=/flashfs1/scratch.global/ynwang/ynwang/2024-11-15-GTBMap/2025-03-25-RNA-seq-mapping/160LRS_match
PHASE_VCF_dir=/flashfs1/scratch.global/ynwang/ynwang/2024-11-15-GTBMap/2025-02-15-ASE/2025-10-07-lorals/input/processed_vcf
SNP_ASE_dir=$pwd/SNP_ASE
HAP_ASE_allvcf_dir=$pwd/HAP_ASE_allvcf
HAP_ASE_gene_dir=$pwd/HAP_ASE_gene
gtfref=/lustre/home/xdzou/2024-10-21-GTBMap/2025-02-11-RNA-mapping/gencode.v47.annotation.gtf
genomeref=/lustre/home/ynwang/2024-11-15-GTBMap/2025-05-04-editing-quant/GRCh38.primary_assembly.genome.fa
genebed=/flashfs1/scratch.global/ynwang/ynwang/2024-11-15-GTBMap/2025-02-15-ASE/2025-12-21-phaser-LRSvcf/GTOP_novel-GENCODE_v47.gene.forphaser.bed

mkdir -p $HAP_ASE_gene_dir
mkdir -p $WASP_PASS_BAM_dir
mkdir -p $HAP_ASE_allvcf_dir $SNP_ASE_dir

ase_snp_sum(){

	slurm=$pwd/slurm/ase_snp_sum

	outpath=$pwd/SNP_ASE_sum/

	mkdir -p $slurm $outpath

	# #combine ase res for per sample
	# ls $SNP_ASE_dir|cut -d. -f1|sort|uniq| \
	# while read line;
	# do
	# files=(`realpath $SNP_ASE_dir/*.csv|grep chr|grep $line|sort -k1,1V`)
	# cat ${files[@]}|grep -v contig|sed '1icontig\tposition\tvariantID\trefAllele\taltAllele\trefCount\taltCount\ttotalCount\tlowMAPQDepth\tlowBaseQDepth\trawDepth\totherBases\timproperPairs' > $SNP_ASE_dir/$line.all.tsv
	# done
	## pre list output input/SNP_ASE.files.frosum.list
	# Rscript bin/pre_SNP_ASE_listforsum.r


	
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

python $pwd/bin/ase_aggregate_by_individual.me.py \
  $pwd/input/SNP_ASE.files.frosum.list \
  $PHASE_VCF_dir/${ind}_het.vcf.gz \
  /flashfs1/scratch.global/ynwang/ynwang/2024-11-15-GTBMap/2025-02-15-ASE/2025-10-07-lorals/input/GTOP_novel-GENCODE_v47_all_exons.bed \
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

export FONTCONFIG_PATH=/etc/fonts
export FONTCONFIG_FILE=/etc/fonts/fonts.conf
# prepare inputlist
ls $SNP_ASE_dir/*.csv > $pwd/SNP_ASE.txt

# cal lamp
python $pwd/bin/ase_calculate_lamp.py $pwd/SNP_ASE.txt my_study -o $lampDir

echo "++++++++++++++++++++"
echo "process end at:"
date
" > $slurm/cal_ase_lamp.slurm

cd $slurm
sbatch cal_ase_lamp.slurm

}


exp_matrix_splittissue(){

bed_dir="$pwd/phaser_gw_bed_final"
inputfile=(GTOP.bed.gz GTOP.gw_phased.bed.gz)


for t in "${inputfile[@]}"; do


less GMTiP_tissue_code_and_colors.csv|tail -n+2|tr ',' '\t' | while read tis abv tiscode color ; do

    if [[ "$t" == *"gw_phased"* ]]; then
    outfile=$bed_dir/$tis.gw_phased.bed.gz
	else
	outfile=$bed_dir/$tis.bed.gz
    fi

echo $t $tis $tiscode

less $bed_dir/GTOP.bed.gz | awk '
NR==1 {
    for(i=1;i<=NF;i++) {
        header[i]=$i
    }
    printf "%s\t%s\t%s\t%s", $1, $2, $3, $4
    for(i=5;i<=NF;i++) {
        if($i ~ /'$tiscode'/) {
            col_idx[i]=1
            printf "\t%s", $i
        }
    }
    print ""
}
NR>1 {
    printf "%s\t%s\t%s\t%s", $1, $2, $3, $4
    for(i=5;i<=NF;i++) {
        if(i in col_idx) {
            printf "\t%s", $i
        }
    }
    print ""
}' | bgzip -c >  $outfile
tabix -p bed $outfile

done

done


}


ase_snp_level(){
slurm=$pwd/slurm/ase_snp_level
#index 
gatk CreateSequenceDictionary   -R $genomeref   -O /lustre/home/ynwang/2024-11-15-GTBMap/2025-05-04-editing-quant/GRCh38.primary_assembly.genome.dict

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


ase_haplotype_level_gene(){


#new and old in this folder
HAP_ASE_gene_dir=$pwd/HAP_ASE_gene_allvcf_combinebam
mkdir -p $HAP_ASE_gene_dir

HAP_ASE_dir=$pwd/HAP_ASE_allvcf_combinebam

cd $HAP_ASE_dir



slurm=$pwd/slurm/ase_haplotype_level_gene_allvcf_combinebam
mkdir -p $slurm



ls ${WASP_PASS_BAM_dir}/*.Aligned.sortedByCoord.WASPpass.out.bam| \
awk -F/ '{print $NF}' | \
cut -d- -f1-2|sort|uniq| \
while read vcffilename ; do



if [ ! -f ${HAP_ASE_dir}/$vcffilename.vcf.gz  ];then
	echo "$vcffilename step1 no finished,skipping..."
	continue
fi

if [ -e $HAP_ASE_gene_dir/$vcffilename.gene_ae.txt  ];then
	echo "$vcffilename exist,skipping..."
	continue
fi

echo "$vcffilename is starting!!!"
	# sample=AGTEX-AJ141-0087-SN-2ALA
	echo '#!/bin/bash
#SBATCH --job-name='${vcffilename}'_phsgene
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --error='${vcffilename}'.err
#SBATCH --output='${vcffilename}'.out
##########################################
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"

if [ ! -f '$HAP_ASE_gene_dir'\'$sample'.gene_ae.txt ] ;then
	# gene level
	python3 /lustre/home/ynwang/software/phaser/phaser_gene_ae/phaser_gene_ae.py \
		--haplotypic_counts '${HAP_ASE_dir}'/'$vcffilename'.haplotypic_counts.txt \
		--features '$genebed' \
		--o '$HAP_ASE_gene_dir'/'$vcffilename'.gene_ae.txt #\
		# --min_haplo_maf 0.5
fi

# module unload python
echo "processs will sleep 30s"
sleep 5
echo "process end at : "
date
#module unload biosoft/torus
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
wait' > $slurm/$sample.slurm

wait
cd ${slurm}
sbatch $slurm/$sample.slurm
echo $sample" is finished"
  process=`squeue | grep ynwang |wc -l`
		while (( process >= 300))
		do
			echo "Current number of jobs is larger than 120"
			echo "Wait another 1 minutes"
			sleep 1m
			process=`squeue | grep ynwang|wc -l`
		done
done


}

ase_haplotype_exp_matrix(){


slurm=$pwd/slurm/ase_haplotype_matrix/
mkdir -p $slurm
output=$pwd/phaser_gw_bed_final/
mkdir -p $output

# if [ -f $output/$tissue.gw_phased.bed.gz ] ;then
# echo "$tissue exist, skipping..."
# continue
# fi

HAP_ASE_gene_dir=$pwd/HAP_ASE_gene_allvcf_combinebam

cd $pwd/HAP_ASE_gene_allvcf_combinebam

	echo '#!/bin/bash
#SBATCH --job-name=ase_haplotype_matrix
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error=ase_haplotype_matrix.err
#SBATCH --output=ase_haplotype_matrix.out
##########################################
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"


	python3 /lustre/home/ynwang/software/phaser/phaser_pop/phaser_expr_matrix.py \
		--gene_ae_dir '${HAP_ASE_gene_dir}' \
		--features '$genebed' \
		--t 8 \
		--o '$output'/GTOP

		module load pigz

		pigz -dc '$output'/GTOP.gw_phased.bed.gz | awk '\''NR == 1 {print $0} NR > 1 {for (i = 5; i <= NF; i++) {if ($i != "0|0") {print $0
		break}}}'\'' | bgzip -c > '$output'/GTOP.gw_phased.filter.bed.gz

tabix -p bed '$output'/GTOP.gw_phased.filter.bed.gz

		module unload pigz


# module unload python
echo "processs will sleep 30s"
sleep 5
echo "process end at : "
date
#module unload biosoft/torus
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
wait' > $slurm/$tissue.slurm

wait
cd ${slurm}
sbatch $tissue.slurm
# echo $sample" is finished"
  process=`squeue | grep ynwang |wc -l`
		while (( process >= 120))
		do
			echo "Current number of jobs is larger than 120"
			echo "Wait another 1 minutes"
			sleep 1m
			process=`squeue | grep ynwang|wc -l`
		done




}

ase_haplotype_level(){

HAP_ASE_allvcf_dir=$pwd/HAP_ASE_allvcf_combinebam

slurm=$pwd/slurm/ase_haplotype_level_allvcf_combinebam
mkdir -p $slurm $HAP_ASE_allvcf_dir


cd $HAP_ASE_allvcf_dir

ls ${WASP_PASS_BAM_dir}/*.Aligned.sortedByCoord.WASPpass.out.bam| \
awk -F/ '{print $NF}' | \
cut -d- -f1-2|sort|uniq| \
while read vcffilename ; do

# 
files=(`ls ${WASP_PASS_BAM_dir}/*.Aligned.sortedByCoord.WASPpass.out.bam | grep $vcffilename`)

# 
bamfile=$(printf ",%s" "${files[@]}")
bamfile=${bamfile:1} 

echo "$vcffilename is starting!!!"
	echo '#!/bin/bash
#SBATCH --job-name=combie_'${vcffilename}'_phaser
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error='${vcffilename}'.err
#SBATCH --output='${vcffilename}'.out
##########################################
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"


python3 /lustre/home/ynwang/software/phaser/phaser/phaser.py \
--vcf /flashfs1/scratch.global/ynwang/ynwang/2024-11-15-GTBMap/2025-02-15-ASE/2025-12-21-phaser-LRSvcf/GMTiP_SNV_INDEL_LRS.71inds.forphaser.vcf.gz \
--bam '$bamfile' \
--paired_end 1 --mapq 255 \
--baseq 10 --sample '$vcffilename' \
--blacklist /lustre/home/ynwang/software/phaser/phaser_test/hg38_hla.chr.bed --haplo_count_blacklist /lustre/home/ynwang/software/phaser/phaser_test/hg38_haplo_count_blacklist.chr.bed --threads 8 --o '${HAP_ASE_allvcf_dir}'/'$vcffilename' --id_separator _ --gw_phase_method 1 --include_indels 0 --gw_phase_vcf 1 --gw_af_field AF

echo "processs will sleep 30s"
sleep 5
echo "process end at : "
date
#module unload biosoft/torus
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
wait' > $slurm/$vcffilename.slurm

wait
cd ${slurm}
sbatch $vcffilename.slurm
echo $vcffilename" is finished"
done


}



ase_haplotype_cis(){
slurm=$pwd/slurm/ase_haplotype_cis_LRSasepair_tis
mkdir -p $slurm


bed_dir="$pwd/phaser_gw_bed_final"
vcf=/flashfs1/scratch.global/ynwang/ynwang/2024-11-15-GTBMap/2025-02-15-ASE/2025-12-21-phaser-LRSvcf/GMTiP_SNV_INDEL_LRS.71inds.forphaser.vcf.gz
map_dir="$pwd/sample_map/"
phaser_cis_var_out_dir=$pwd/phaser_cis_var_LRSasepair_tis
mkdir -p $phaser_cis_var_out_dir $map_dir


less GMTiP_tissue_code_and_colors.csv|tail -n+2|tr ',' '\t' | while read tis abv tiscode color ; do


# make sample map
if [ ! -f  $map_dir/$tis.sample_map.txt ];then
less ${bed_dir}/$tis.gw_phased.bed.gz|head -1|tr '\t' '\n' |tail -n +5 > $map_dir/$tis.bed.sample
less $map_dir/$tis.bed.sample|awk -v OFS="\t" -F "-" '{print $1"-"$2, $0}'|sed '1i vcf_sample\tbed_sample' > $map_dir/$tis.sample_map.txt
fi
#
chrlist=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

for chr in ${chrlist[@]};do
	echo '#!/bin/bash
#SBATCH --job-name='$tis'_'$chr'_cisourvcf_nogw
#SBATCH --partition=cu-1,cu-debug,cu-short,cpuPartition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --error='$tis'_'$chr'.err
#SBATCH --output='$tis'_'$chr'.out
##########################################
echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"


python3 /lustre/home/ynwang/software/phaser/phaser_pop/phaser_cis_var.qym.py --bed '${bed_dir}'/'$tis'.bed.gz \
					--vcf '$vcf' \
					--pairs /lustre/home/ynwang/2024-11-15-GTBMap/2025-02-15-ASE/allLRSsig.genesnppairs.for.phasercis.txt \
					--map  '${map_dir}'/'$tis'.sample_map.txt \
					--o '$phaser_cis_var_out_dir'/'$tis'.'${chr}'.results.txt --ignore_v 0 --t 16 \
					--chr chr'${chr}'

# module unload python
echo "processs will sleep 30s"
sleep 5
echo "process end at : "
date
#module unload biosoft/torus
rm -rf /tmp/nodefile.$$
rm -rf /tmp/nodes.$$/
wait' > $slurm/$tis.$chr.slurm
cd ${slurm}
sbatch $tis.$chr.slurm
echo $tis.$chr " is finished"

done
done
}
main
