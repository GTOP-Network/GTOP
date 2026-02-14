setwd("/lustre/home/xdzou/2024-10-21-GTBMap/2025-07-03-Finemap_susie")
rm(list=ls())

library(data.table)
library(dplyr)
library(magrittr)

# load SV list

df_snv <- fread("./input/GTOP_SNV.coord.sorted.bed",header=F,sep="\t")
df_sv <- fread("./input/GTOP_SV.coord.sorted.bed",header=F,sep="\t")
snv_list <- unique(df_snv$V4)
sv_list <- unique(df_sv$V4)
rm(df_snv,df_sv)

# load fine map summarized file
df_fm <- fread("./output/GTOP_eQTL.finemap_Joint.all_tissues.txt",header=T,sep="\t")

df_fm$varType <- "TR"
df_fm$varType[df_fm$variant_id %in% snv_list] <- "SNV"
df_fm$varType[df_fm$variant_id %in% sv_list] <- "SV"

fwrite(df_fm,file="./output/GTOP_eQTL.finemap_Joint.all_tissues.anno_vartype.txt",quote=F,sep="\t",row.names=F)

