#!/bin/bash

# Bash Script Pipeline for S-LDSC Analysis

# Description:
# This script processes input files and run S-LDSC.

# Example Usage:
# ./script.sh <workpath> <sigVariantdir> <finemapdir> <tissueName> <genodir> <gwassumdir> <weightsdir> <baselinedir>

workpath=$1             # Directory: work path
sigVariantdir=$2        # Directory: contains significant variant/gene lists
finemapdir=$3           # Directory: fine-mapping results (e.g., SuSiE, FINEMAP outputs)
tissueName=$4           # String: tissue name (e.g., Liver, Brain, etc.)
genodir=$5            # Directory: genotype data (PLINK format files) download from https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays
gwassumdir=$6           # Directory: GWAS summary statistics
weightsdir=$7           # Directory: LD score regression weights files download from https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays
baselinedir=$8          # Directory: baseline annotation files download from https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays



snplistdir=$workpath/input/generate_snp
ldscanndir=$workpath/input/ldsc_annot
annofiledir=$workpath/input/paste_all
ldscresdir=$workpath/output/res
tmpscript=$workpath/script_tmp

mkdir -p $snplistdir $ldscanndir $annofiledir $ldscresdir $tmpscript

main(){
generate_snp_f
ldsc_annot_f
make_annot_f
enrichment_f
}



generate_snp_f(){



#print sigVariant list
less $sigVariantdir/${tissueName}.significant_pairs.txt | awk '{print $2}' | sort | uniq > $snplistdir/$tissueName.allsig.snp
#print cs variant with pip > 0.01
less $finemapdir/susieR_res.snv_egenes.${tissueName}.txt |awk 'NR > 1 && $3 > 0.01 && $4 !~ /NA/ {print $2}'|sort|uniq > $snplistdir/${tissueName}.finemap.snp

#print cs variant with max pip value
less  $finemapdir/susieR_res.snv_egenes.${tissueName}.txt| awk 'NR > 1 && $3 > 0.01 && $4 !~ /NA/ {print $2, $3}'|sort -k1,1 -k2,2nr | awk '!seen[$1]++ {print $1, $2}' > $snplistdir/${tissueName}.finemappip.snp 



}


ldsc_annot_f(){

filetype=( allsig finemappip finemap )
for file_t in ${filetype[@]};do

		tissueName=`echo "$file" |awk -F"/" '{print $NF;exit}'|awk -F"." '{print $1;exit}'`

for chr in $(seq 1 22); do
if [ $file_t != "finemappip" ];then

echo '
library(data.table)

QTL_snp=fread(paste0("'${sigVariantdir}'","'${tissueName}'.'$file_t'.snp"),head=F,stringsAsFactors=F,data.table=F)

snp_top=unique(QTL_snp[,1])

bim=fread(paste0("'${genodir}'","'${chr}'",".bim"),head=F,stringsAsFactors=F,data.table=F)
bim$snp_top=0; 
index3=match(snp_top,bim$V2,nomatch=0)
bim$snp_top[index3]=1;

anno=bim[,c(1,4,2,3,7)]
colnames(anno)=c("CHR","BP","SNP","CM","'${tissueName}'")

write.table(anno[,5],paste0("'${ldscanndir}'/'${tissueName}'.","'${chr}'",".annot"),row=F,col=T,quo=F,sep="\t")

'> ${tmpscript}/${tissueName}.${chr}.r
else
echo '
library(data.table)

QTL_snp <- fread(paste0("'${sigVariantdir}'","'${tissueName}'.'$file_t'.snp"),
                head=F, stringsAsFactors=F, data.table=F)
bim <- fread(paste0("'${genodir}'","'${chr}'",".bim"),
             head=F, stringsAsFactors=F, data.table=F)


match_values <- QTL_snp$V2[match(bim$V2, QTL_snp$V1)]
bim$V7[!is.na(match_values)] <- match_values[!is.na(match_values)]
bim$V7[is.na(match_values)] <- 0
anno <- bim[, c(1,4,2,3,7)]
colnames(anno) <- c("CHR", "BP", "SNP", "CM", "'${tissueName}'")
write.table(anno[, 5], 
            paste0("'${ldscanndir}'/'${tissueName}'.","'${chr}'",".annot"),
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

'> ${tmpscript}/${tissueName}.${chr}.r
fi
wait


R CMD BATCH --no-save ${tmpscript}/${tissueName}.${chr}.r

done
done

}

make_annot_f(){

filetype=(finemappip finemap allsig)


for file_t in ${filetype[@]};do
for i in $(seq 1 22); do

gzip -c $ldscanndir/$tissueName.${i}.annot > ${annofiledir}/$tissueName.${i}.annot.gz


ldsc.py --l2 --bfile ${genodir}/1000G.EAS.QC.${i} --ld-wind-cm 1 \
--annot ${annofiledir}/$tissueName.${i}.annot.gz \
--thin-annot \
--out ${annofiledir}/$tissueName.${i} \
--print-snps weights.EAS.hm3_noMHC.list


done
done

}



enrichment_f(){

filetype=(allsig finemappip finemap)
for file_t in ${filetype[@]};do
for file in ${gwassumdir}/*.sumstats.gz;do
		trait=`echo "$file" |awk -F"/" '{print $NF;exit}'|awk -F".sumstats.gz" '{print $1;exit}'`
		cd ${ldscresdir}
python ldsc.py --h2 ${gwassumdir}/${trait}.sumstats.gz \
--ref-ld-chr ${baselinedir}/baseline.,$annofiledir/$tissueName. \
--w-ld-chr ${weightsdir}/weights.EAS.hm3_noMHC. \
--overlap-annot --print-coefficients \
--frqfile-chr ${genodir}/1000G.EAS.QC. \
--out $ldscresdir/${trait} --print-delete-vals

done
done

}








main
