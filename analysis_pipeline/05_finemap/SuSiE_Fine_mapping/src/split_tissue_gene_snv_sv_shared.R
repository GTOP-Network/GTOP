library(optparse)

option_list <- list(
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)

setwd("/lustre/home/xdzou/2024-10-21-GTBMap/2025-07-03-Finemap_susie")
egenelist <- "/lustre/home/xdzou/2024-10-21-GTBMap/2025-02-13-expression-quant/output/QTL_mapping_131/covariates/PC_0_no_Batch/ALL_varTypes.eGenes.FDR_0.1.separate.txt"

df_egenes_all <- fread(egenelist,header=T,sep="\t")

df_tissue.snv <- df_egenes_all %>% filter(varType!="SV",Tissue==opt$tissue)
df_tissue.sv <- df_egenes_all %>% filter(varType=="SV",Tissue==opt$tissue)


shared_egenes <- intersect(df_tissue.snv$Gene,df_tissue.sv$Gene)

df_out <- data.frame(Tissue=opt$tissue,Gene=shared_egenes)


outfile <- paste0("./input/tissue_gene2/",opt$tissue,"_gene_list.txt")

write.table(df_out,file=outfile,quote=F,sep="\t",row.names=F,col.names=F)



