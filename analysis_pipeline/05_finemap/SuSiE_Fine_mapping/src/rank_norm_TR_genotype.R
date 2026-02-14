setwd("/lustre/home/xdzou/2024-10-21-GTBMap/2025-07-03-Finemap_susie")

library(data.table)
library(dplyr)
library(magrittr)

library(scales)
# function
rank_norm_scale <- function(x){
	r <- rank(x,ties.method="average")
	return(round(rescale(r,to=c(0,2)),digits=4))
}

# load data
dat <- fread("./input/GMTip_LRS_TR_147INDs.Miss85.AF95.AC1.dosage.simpleTR.filter.final.txt",header=T,sep="\t")

dat <- as.data.frame(dat,check.names=F)
dat$MOTIFS <- NULL

dat_out <- cbind(data.frame(variant_id=dat$TRID),lapply(dat[,-1],rank_norm_scale))

cat(dim(dat_out),"\n")

write.table(dat_out,file="./input/GTOP_LRS_131inds.TR.genotype_norm_scaled.V2.txt",quote=F,sep="\t",row.names=F)
