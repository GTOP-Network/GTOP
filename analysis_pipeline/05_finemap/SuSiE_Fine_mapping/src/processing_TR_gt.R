setwd("/lustre/home/xdzou/2024-10-21-GTBMap/2025-07-03-Finemap_susie/input/")
library(data.table)
library(dplyr)



dat <- fread("GMTip_LRS_TR_131INDs.Miss85.AF95.AC1.dosage.impute.simpleTR.filter.txt",header=T,sep="\t")

names(dat) <- c("variant_id",names(dat)[-1])

saveRDS(dat,file="GTOP_TR_131inds.GT_dosage.RDS")
