library(optparse)

option_list <- list(
	make_option(c("-g","--gene"),type="character",default="NA",action="store",help="specify a snp list"),
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a snp list")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)

setwd("/flashfs1/scratch.global/xdzou/Fine_map_susie")
GTdir <- paste0("./input/SNV_by_egenes/",opt$gene)

gt_files <- list.files(path=GTdir,pattern=".vcf")

df_gt <- data.frame()

df_snv <- fread(paste0(GTdir,"/",gt_files[1]),header=T,sep="\t")
header <- names(df_snv)
df_snv %<>% distinct(variant_id,.keep_all=T)

for(i in 2:length(gt_files)){
	df_var <- fread(paste0(GTdir,"/",gt_files[i]),header=T,sep="\t")
	df_var %<>% distinct(variant_id,.keep_all=T)
	df_var %<>% select(all_of(header))
	df_snv <- rbind(df_snv,df_var)
}

outfile <- paste0("./output/Joint/",opt$tissue,"/",opt$gene,"/eVar_GT.vcf")

fwrite(df_snv,outfile,quote=F,sep="\t",row.names=F)
