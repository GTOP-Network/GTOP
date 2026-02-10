library(optparse)

option_list <- list(
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)

setwd("/lustre/home/xdzou/2024-10-21-GTBMap/2025-05-07-GCTA-eQTL")

dir <- "/lustre/home/xdzou/2024-10-21-GTBMap/2025-02-13-expression-quant/output/"

input_file <- paste0(dir,"phenotype_new/",opt$tissue,".phenotype.bed")
#output_file <- paste0("./cis_QTL_SE/",opt$tissue,"cis_QTL_all.add_se.txt")


# load raw cis associations
dat <- as.data.frame(fread(input_file,header=T,sep="\t"),check.names=F,stringsAsFactors=F)
dat <- dat[,c(-1,-2,-3)]
rownames(dat) <- dat$phenotype_id

expr <- t(dat[,-1])


# load covariates
cov_file <- paste0(dir,"covariates_new/",opt$tissue,".covariates.txt")
df_cov <- as.data.frame(fread(cov_file,header=T,sep="\t"),check.names=F,stringsAsFactors=F)
rownames(df_cov) <- df_cov$ID
df_cov$ID <- NULL


inds_kept <- intersect(colnames(df_cov),rownames(expr))


df_cov <- df_cov[,inds_kept]
expr <- expr[inds_kept,]

resids <- matrix(, ncol = ncol(expr), nrow = nrow(expr))
rownames(resids) = rownames(expr)
colnames(resids) = colnames(expr)

#genes <- colnames(expr)
for(i in 1:ncol(expr)){
	data = as.data.frame(cbind(expr[, i], t(df_cov)))
	colnames(data) = c('Exp', rownames(df_cov))
	model <- lm(Exp ~ ., data=data)
	resids[,i] = model$residuals
}


resids <- t(scale(resids))
#rownames(resids) <- genes
output_file <- paste0("./output/Exp_residual/",opt$tissue,".residual_ztrans.txt")
write.table(resids,file=output_file,quote=F,sep="\t",row.names=T)
