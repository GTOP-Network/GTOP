# format alt TSS quant results

setwd("/lustre/home/xdzou/2024-10-21-GTBMap/2025-07-03-Finemap_susie")
library(data.table)
library(dplyr)
library(magrittr)


dir <- "/lustre/home/xdzou/2024-10-21-GTBMap/2025-02-13-expression-quant"
# -- functions
extract_tss <- function(x){
	chrom <- x[1]
	pos1 <- x[2]
	pos2 <- x[3]
	strand <- x[4]
	if(strand=="+"){
		tss = pos1
	}else{
		tss = pos2
	}
	return(tss)
}

# -- load tested genes
dat_test_gene <- fread("./input/TSS_location.GMTiP.sorted.uniq.bed",header=F,sep="\t")
names(dat_test_gene) <- c("chrom","pos1","pos2","Gene")

# -- load gene annotation
df_gene <- fread(paste0(dir,"/input/gencode.v47.gene_info.txt"),header=F,sep="\t")
names(df_gene) <- c("chrom","start","end","strand","Gene","Name","GeneType")
df_gene$Name <- NULL
df_gene$GeneType <- NULL

df_gene$TSS <- apply(df_gene[,c(1,2,3,4)],1,extract_tss)
df_gene %<>% select(Gene,chrom,TSS)
df_gene$pos0 <- as.numeric(df_gene$TSS) - 1
df_gene %<>% select(chrom,pos0,TSS,Gene) %>% filter(Gene %in% dat_test_gene$Gene)


write.table(df_gene,file="./input/TSS.GTOP_tested_genes.bed",quote=F,sep="\t",row.names=F)


