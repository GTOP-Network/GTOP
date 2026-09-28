library(optparse)

option_list <- list(
	make_option(c("-m","--map"),type="character",default="NA",action="store",help="specify a tissue"),
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue"),
	make_option(c("-g","--genes"),type="character",default="NA",action="store",help="specify a tissue")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)

setwd("/path/to/2025-12-09-cross-pop")

# functions 
# load reference
df_ref <- fread(opt$map,header=T,sep="\t")
names(df_ref) <- c("SNP","SNP_156")

# load gene list
genelist_file <- paste0("./input/",opt$tissue,"/",opt$genes)
df_genes <- fread(genelist_file,header=F,sep="\t")
names(df_genes) <- c("Gene","chrom","pos")

selected_genes <- unique(df_genes$Gene)
rm(df_genes)

for(i in 1:length(selected_genes)){
	cat(i,selected_genes[i]," \n")
	infile <- paste0("./output/",opt$tissue,"/SuSiEx/examples/GTEX_",selected_genes[i],".txt")
	if(file.exists(infile)){
		df_gene <- fread(infile,header=T,sep="\t")
		df_ref.sub <- df_ref %>% filter(SNP %in% df_gene$SNP)
		df_gene <- merge(df_gene,df_ref.sub,by="SNP",all.x=T)
		df_gene$SNP_156[is.na(df_gene$SNP_156)] <- df_gene$SNP[is.na(df_gene$SNP_156)]
		df_gene %<>% select(chrom,bp,A2,A1,SNP_156,pvalue,beta,se)
		names(df_gene) <- c("chrom","bp","A2","A1","SNP","pvalue","beta","se")
		output_file <- paste0("./output/",opt$tissue,"/SuSiEx/examples/GTEX_",selected_genes[i],".sumstat.txt")
		fwrite(df_gene,file=output_file,quote=F,sep="\t",row.names=F)
	}else{
		cat(infile," not exists!!\n")
	}
}

