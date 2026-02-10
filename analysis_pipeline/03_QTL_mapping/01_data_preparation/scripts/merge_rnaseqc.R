library(optparse)
option_list <- list(
	make_option(c("-p","--path"),type="character",default="NA",action="store",help="specify a tissue")
	)

setwd(opt$path)

library(data.table)
library(dplyr)
library(magrittr)


# 01. merge reads count of each sample into a multi-smaple count matrix
infile_list <- list.files(path="./output/RNASeQC_out/", pattern=".gene_reads.gct")

df_m <- fread(paste0("./output/RNASeQC_out/",infile_list[1]),skip=2,header=T,sep="\t")
df_m <- df_m[,c(1,3)]
colnames(df_m)[1] <- "Gene"
head(df_m)

for(i in 2:length(infile_list)){
	cat(i,": ", infile_list[i],"\n")
	dat <- fread(paste0("./output/RNASeQC_out/",infile_list[i]),skip=2,header=T,sep="\t")
	cat(dim(dat),"\n")
	dat <- dat[,c(1,3)]
	colnames(dat)[1] <- "Gene"

	df_m <- merge(df_m,dat,by="Gene")
	cat(dim(df_m),"\n")
}

dim(df_m)
df_m[1:5,1:5]
write.table(df_m,file="./output/All_samples.readscount.txt",quote=F,sep="\t",row.names=F)


# 02. merge TPM of each sample into a multi-smaple TPM matrix
infile_list <- list.files(path="./output/RNASeQC_out/", pattern=".gene_tpm.gct")
df_m <- fread(paste0("./output/RNASeQC_out/",infile_list[1]),skip=2,header=T,sep="\t")
df_m <- df_m[,c(1,3)]
colnames(df_m)[1] <- "Gene"

for(i in 2:length(infile_list)){
	cat(i,": ", infile_list[i],"\n")
	dat <- fread(paste0("./output/RNASeQC_out/",infile_list[i]),skip=2,header=T,sep="\t")
	cat(dim(dat),"\n")
	dat <- dat[,c(1,3)]
	colnames(dat)[1] <- "Gene"

	df_m <- merge(df_m,dat,by="Gene")
	cat(dim(df_m),"\n")
}

write.table(df_m,file="./output/All_samples.TPM.txt",quote=F,sep="\t",row.names=F)
