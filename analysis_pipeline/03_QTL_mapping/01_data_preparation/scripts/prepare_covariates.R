library(optparse)

option_list <- list(
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)
library(PCAForQTL)


prcompRes <- readRDS(paste0("./output/pca_res/",opt$tissue,".prcompRes.RDS"))

PCs <- prcompRes$x
dim(PCs)
PCs[1:5,1:5]


# choose PCs by runElbow
K_elbow <- runElbow(prcompResult=prcompRes)



# load known covariates
df_known <- as.data.frame(fread("./input/Known_covariates.txt",header=T,sep="\t"),
						  check.names=F,stringsAsFactors=F)

#df_known <- df_known %>% filter(ID %in% known_covariates)
df_known <- as.data.frame(df_known,check.names=F)
rownames(df_known) <- df_known$ID
df_known$ID <- NULL
knownCovariates <- t(df_known)

# intersected inds
inds_comm <- intersect(rownames(PCs),rownames(knownCovariates))

colnames(PCs) <- paste0("InferredCov_",1:dim(PCs)[2])
PCs_top <- PCs[inds_comm,1:K_elbow]
knownCovariates <- knownCovariates[inds_comm,]


knownCov_filtered <- filterKnownCovariates(knownCovariates,PCs_top,unadjustedR2_cutoff=0.9)
PCs_top <- scale(PCs_top)
final_covs <- t(cbind(knownCov_filtered,PCs_top))

df_c <- cbind(data.frame(ID=rownames(final_covs)),as.data.frame(final_covs,check.names=F))

write.table(df_c,file=paste0("./output/covariates/",opt$tissue,".covariates.txt"),quote=F,sep="\t",row.names=F)



