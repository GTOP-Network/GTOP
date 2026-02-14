library(optparse)

option_list <- list(
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue"),
	make_option(c("-g","--genes"),type="character",default="NA",action="store",help="specify a gene list")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)

setwd("/lustre/home/xdzou/2024-10-21-GTBMap/2026-02-03-GTOP_eQTL_Finemap_susie")
#pheno_dir <- "/lustre/home/xdzou/2024-10-21-GTBMap/2026-02-02-GTOP_eQTL_mapping/output/Exp_residual/"
pheno_dir <- "/lustre/home/xdzou/2024-10-21-GTBMap/2026-02-02-GTOP_eQTL_mapping/output/sQTL_residual/ju/" # juQTL
#pheno_dir <- "/lustre/home/xdzou/2024-10-21-GTBMap/2026-02-02-GTOP_eQTL_mapping/output/sQTL_residual/tu/" # tuQTL
pheno_bed <- paste0(pheno_dir,opt$tissue,".residual_ztrans.txt")

df_pheno <- fread(pheno_bed)

df_pheno <- as.data.frame(df_pheno,check.names=F,stringsAsFactors=F)

rownames(df_pheno) <- df_pheno$V1
df_pheno$V1 <- NULL
df_mat <- as.matrix(df_pheno)
inds <- names(df_pheno)

df_genes <- fread(opt$genes,header=F,sep="\t")
names(df_genes) <- c("Tissue","Gene")

for(i in 1:dim(df_genes)[1]){
	gene <- df_genes$Gene[i]
	exp_gene <- as.numeric(df_mat[gene,])
	outdir <- paste0("/flashfs1/scratch.global/xdzou/Fine_map_susie_ju/output/Joint/",opt$tissue,"/",gene)
	df_exp_gene <- data.frame(PID=inds,IID=inds,Exp=exp_gene)

	if(file.exists(outdir)){
		cat(outdir," already exist!\n")
		fwrite(df_exp_gene,file=paste0(outdir,"/expr.phen"),quote=F,sep="\t",row.names=F,col.names=F)
	}else{
		dir.create(outdir)
		fwrite(df_exp_gene,file=paste0(outdir,"/expr.phen"),quote=F,sep="\t",row.names=F,col.names=F)
	}
}



