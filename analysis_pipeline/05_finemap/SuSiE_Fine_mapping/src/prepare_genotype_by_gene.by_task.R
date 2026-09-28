library(optparse)

option_list <- list(
	make_option(c("-t","--task"),type="character",default="NA",action="store",help="specify a snp list")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)

setwd("/path/to/dir")
genotype_bed <- "./input/GTOP_LRS_SRS.small_variants.autosome_X.maf_05.gt_dosage.txt"

genelist <- paste0("./input/task/",opt$task)

df_geno <- fread(genotype_bed,header=T,sep="\t")
setkey(df_geno,ID)

df_genes <- fread(genelist,header=F,sep="\t")
names(df_genes) <- c("Gene")

for(i in 1:dim(df_genes)[1]){
	gene <- df_genes$Gene[i]
	snp_file <- paste0("./input/SNV_by_egenes/",gene,"/snv_list.txt")
	df_snp <- fread(snp_file,header=F,sep="\t")
	snp_list <- df_snp$V1
	rm(df_snp)
#	df_gt <- df_geno %>% filter(variant_id %in% snp_list)
	df_gt <- df_geno[J(snp_list),nomatch = 0]
	outdir <- paste0("./input/SNV_by_egenes/",gene)

	if(file.exists(outdir)){
		cat(outdir," already exist!\n")
		fwrite(df_gt,file=paste0(outdir,"/eSNV_GT.vcf"),quote=F,sep="\t",row.names=F,na=".")
	}else{
		dir.create(outdir)
		fwrite(df_gt,file=paste0(outdir,"/eSNV_GT.vcf"),quote=F,sep="\t",row.names=F,na=".")
	}
}



