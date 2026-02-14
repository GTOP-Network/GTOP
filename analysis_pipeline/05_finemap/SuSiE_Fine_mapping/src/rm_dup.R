library(optparse)

option_list <- list(
	make_option(c("-d","--dir"),type="character",default="NA",action="store",help="specify a tissue")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)

setwd(opt$dir)

df_gt <- fread("eTR_GT.vcf",header=T,sep="\t")

df_gt %<>% distinct(variant_id,.keep_all=T)

write.table(df_gt,file="eTR_GT.vcf",quote=F,sep="\t",row.names=F)
