library(optparse)

option_list <- list(
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)

setwd("/path/to/dir")

# functions 
rm_version <- function(x){
	return(strsplit(x,split=".",fixed=T)[[1]][1])
}


df_gtex <- fread(paste0("./input/GTEx_eQTL/eGene/",opt$tissue,".egene.txt"),header=F)

df_gtop <- fread(paste0("./input/",opt$tissue,".egene_list.txt"),header=F)

df_gtex$V1 <- sapply(df_gtex$V1,rm_version)

df_gtop %<>% filter(V1 %in% df_gtex$V1)

write.table(df_gtop,file=paste0("./input/",opt$tissue,".shared_egene.txt"),quote=F,sep="\t",row.names=F,col.names=F)
