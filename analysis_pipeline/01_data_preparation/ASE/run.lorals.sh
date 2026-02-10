#!/bin/bash

# Bash Script Pipeline for lorals

# Description:
# This script processes the eQTL-input files, preparation of input files for torus and run torus.

# Example Usage:
# ./script.sh <workpath> <ori_vcf> <refgenome> <samplelist> <inputbamdir> <reftrans> <refbed> <genetsv>


workpath=$1
ori_vcf=$2
refgenome=$3
samplelist=$4  #
inputbamdir=$5
reftrans=$6
refbed=$7 #genebed
genetsv=$8


provcfdir=$workpath/output/processed_vcf
mapgenodir=$workpath/output/map_genome
maptransdir=$workpath/output/map_trans


mkdir -p $provcfdir $mapgenodir $maptransdir $wkpath/output/ase $wkpath/output/asts

main(){

process_vcf
hap_map
hap_map_trans
ase_cal_chr
asts_cal_quant
process_asts

}
wkpath=`pwd`
# lorals_dir=/media/bora_A/wangyn/software/lorals
lorals_dir=/lustre/home/ynwang/software/lorals

process_vcf(){

  
./scripts/process_vcf.sh -v $ori_vcf -f $refgenome -s $samplelist -o $provcfdir

}

hap_map(){


for i in `ls $inputbamdir/*.bam`; do

ind=`basename $i .bam| cut -d- -f1-2`
sample=`basename $i .bam| cut -d- -f1-3`

mkdir -p $mapgenodir/$sample
bam2fasta -o $mapgenodir/$sample/$sample $i
cd $mapgenodir/$sample

./scripts/hap_aligner_pacbio.sh -f $mapgenodir/$sample/$sample.fasta.gz -G $provcfdir/$ind -o $mapgenodir/$sample


fi


}



hap_map_trans(){


for i in `ls $inputbamdir/*.bam`; do

ind=`basename $i .bam| cut -d- -f1-2`
sample=`basename $i .bam| cut -d- -f1-3`


mkdir -p $maptransdir/$sample


 cd $maptransdir/$sample


minimap2 -t 16 -ax map-hifi --MD $reftrans $mapgenodir/$sample/$sample.fasta.gz > $maptransdir/$sample/$sample.sam
samtools view -F 2304 -hSb $maptransdir/$sample/$sample.sam | \
samtools sort -@ 16 -o $maptransdir/$sample/$sample.sorted.bam
samtools index -b $maptransdir/$sample/$sample.sorted.bam


done

}


ase_cal_chr(){


mkdir -p $wkpath/output/ase

for sample in `ls $wkpath/output/map_genome/`; do
for chr in chr{1..22};do 

ind=`echo $sample| cut -d- -f1-2`


calc_ase -b $wkpath/output/map_genome/$sample/${sample}_reads_aln_sorted.merged.bam \
-f $provcfdir/${ind}_het.vcf.gz \
-o $wkpath/output/ase/$sample.ase.tsv



annotate_ase -i $wkpath/output/ase/$sample.${chr}.ase.tsv \
-b $refbed \
-f $provcfdir/${ind}_het.vcf.gz \
-o $wkpath/output/ase/$sample.ase.annotated.tsv


awk -F'\t' '!(\$15==1 || \$16==1 || \$17==1 || \$18==1 || \$19==1)' $wkpath/output/ase/$sample.ase.annotated.tsv > $wkpath/output/ase/$sample.ase.annotated.clean.tsv


}


asts_cal_quant(){


for sample in `ls $wkpath/output/map_genome/`; do


ind=`echo $sample| cut -d- -f1-2`

cd $wkpath/output/asts


calc_asts -m quant \
-b $mapgenodir/$sample/${sample}_reads_aln_sorted.merged.bam \
-i $wkpath/output/ase/$sample.ase.annotated.clean.tsv \
-o $wkpath/output/asts/$sample \
-x $maptransdir/$sample/$sample.sorted.bam


done




}

process_asts(){

for sample in `ls $wkpath/output/map_genome/`; do

ind=`echo $sample| cut -d- -f1-2`


fileall=`realpath $wkpath/output/asts/*.clean.tsv|grep $sample|grep clean`


mkdir -p  $wkpath/output/asts/${sample}_g10t10
process_asts -i ${fileall[@]} -g $wkpath/input/GTOP_novel-GENCODE_v47.genes.tsv \
-o $wkpath/output/${folder}/${sample}_g10t10 --multiple-snps


mkdir -p $wkpath/output/asts/${sample}_g20t20
process_asts -i ${fileall[@]} -g $wkpath/input/GTOP_novel-GENCODE_v47.genes.tsv \
-o $wkpath/output/asts/${sample}_g20t20 --multiple-snps --min-reads-gene 20 --min-reads-transcript 20

done
}

main
