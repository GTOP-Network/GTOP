library(optparse)

option_list <- list(
	make_option(c("-g","--genes"),type="character",default="NA",action="store",help="specify a gene list"),
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue"),
	make_option(c("-v","--varType"),type="character",default="SNV",action="store",help="specify variant type: SNV SV")

)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)

setwd("/lustre/home/xdzou/2024-10-21-GTBMap/2025-07-03-Finemap_susie")
dir <- paste0("./output/Finemapping_eQTL_",opt$varType,"/",opt$tissue,"/")
df_genes <- fread(opt$genes,header=F,sep="\t")
genelist <- unique(df_genes$V2)

for(i in 1:length(genelist)){
	infile <- paste0(dir,genelist[i],"/",snp_list.txt)

}

rownames(df_pheno) <- df_pheno$phenotype_id
df_pheno$phenotype_id <- NULL
df_mat <- as.matrix(df_pheno)
inds <- names(df_pheno)

df_genes <- fread(genelist,header=F,sep="\t")
names(df_genes) <- c("Tissue","Gene")

for(i in 1:dim(df_genes)[1]){
	gene <- df_genes$Gene[i]
	exp_gene <- as.numeric(df_mat[gene,])
	outdir <- paste0("./output/Finemapping_eQTL/",opt$tissue,"/",gene)
	df_exp_gene <- data.frame(PID=inds,IID=inds,Exp=exp_gene)

	if(file.exists(outdir)){
		cat(outdir," already exist!\n")
		fwrite(df_exp_gene,file=paste0(outdir,"/expr.phen"),quote=F,sep="\t",row.names=F,col.names=F)
	}else{
		dir.create(outdir)
		fwrite(df_exp_gene,file=paste0(outdir,"/expr.phen"),quote=F,sep="\t",row.names=F,col.names=F)
	}
}



