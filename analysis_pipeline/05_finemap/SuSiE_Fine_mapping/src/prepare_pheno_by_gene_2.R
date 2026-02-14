library(optparse)

option_list <- list(
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue"),
	make_option(c("-g","--gene"),type="character",default="NA",action="store",help="specify a gene id")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)

setwd("/lustre/home/xdzou/2024-10-21-GTBMap/2025-07-03-Finemap_susie")
pheno_bed <- paste0("/lustre/home/xdzou/2024-10-21-GTBMap/2025-02-13-expression-quant/output/phenotype_new/",opt$tissue,".phenotype.bed")

df_pheno <- fread(pheno_bed,header=T,sep="\t")
df_pheno <- df_pheno[,c(-1,-2,-3)]

df_pheno <- as.data.frame(df_pheno,check.names=F,stringsAsFactors=F)

rownames(df_pheno) <- df_pheno$phenotype_id
df_pheno$phenotype_id <- NULL
df_mat <- as.matrix(df_pheno)
inds <- names(df_pheno)

#df_genes <- fread(genelist,header=F,sep="\t")
#names(df_genes) <- c("Tissue","Gene")

gene <- opt$gene
exp_gene <- as.numeric(df_mat[gene,])
outdir <- paste0("./output/Finemapping_eQTL_sv_snv/",opt$tissue,"/",gene)
df_exp_gene <- data.frame(PID=inds,IID=inds,Exp=exp_gene)
fwrite(df_exp_gene,file=paste0(outdir,"/expr.phen"),quote=F,sep="\t",row.names=F,col.names=F)



