library(data.table)
library(dplyr)
library(magrittr)

meta_file <- "./input/tissue_code_and_colors.csv"
raw_tpm_file <- "./output/All_samples.TPM.txt"
raw_reads_file <- "./output/All_samples.readscount.txt"

# functions

# annotate tissue by code map in meta. meta is a data.frame with two columns: Tissue + Tissue_Code
annotate_tissue <- function(x,meta){
	tissueCode <- strsplit(x,split="-",fixed=T)[[1]][3]
	return(as.character(meta$Tissue[meta$Tissue_Code==tissueCode]))
}

extract_inds <- function(x){
	w <- strsplit(x,split="-",fixed=T)[[1]]
	return(paste(w[1:2],collapse="-"))
}

# load meta:

df_tissue_meta <- read.csv(meta_file,header=T,stringsAsFactors=F,colClasses=rep("character",4))
names(df_tissue_meta) <- c("Tissue","Abbreviation","Tissue_Code","Tissue_Color")
# load TPM and reads count matrix

df_tpm <- as.data.frame(fread(raw_tpm_file,header=T,sep="\t"),check.names=F,stringsAsFactors=F)
df_reads <- as.data.frame(fread(raw_reads_file,header=T,sep="\t"),check.names=F,stringsAsFactors=F)

df_tpm[1:5,1:5]
df_reads[1:5,1:5]
# extract all samples from above matrix
df_all_samples <- data.frame(SampleName=colnames(df_tpm)[-1])
head(df_all_samples)
dim(df_all_samples)

df_all_samples$Tissue <- as.character(sapply(df_all_samples$SampleName,annotate_tissue,meta=df_tissue_meta[,c(1,3)]))
head(df_all_samples)
dim(df_all_samples)

class(df_all_samples$Tissue)
unique_tissue_list <- unique(df_all_samples$Tissue)
head(unique_tissue_list)

unique_tissue_list

for(i in 1:length(unique_tissue_list)){
	tissue <- unique_tissue_list[i]
	cat(tissue,"\n")
	sample_by_tissue <- df_all_samples$SampleName[df_all_samples$Tissue == tissue]
	kept_cols <- c("Gene",sample_by_tissue)
	df_tpm_tissue <- df_tpm[,kept_cols]
	df_reads_tissue <- df_reads[,kept_cols]
	inds <- sapply(sample_by_tissue,extract_inds)
	nRow <- dim(df_tpm_tissue)[1]
	nCol <- length(inds)
	names(df_tpm_tissue) <- c("Gene",inds)
	names(df_reads_tissue) <- c("Gene",inds)
	line1 <- "#1.2"
	line2 <- paste(nRow,nCol,sep="\t")
	outFile_1 <- paste0("./output/reads_gct/",tissue,".reads_count.gct")
	outFile_2 <- paste0("./output/tpm_gct/",tissue,".tpm.gct")
	write.table(line1,file=outFile_1,row.names=F,sep="\t",col.names=F,quote=F)
	write.table(line2,file=outFile_1,append=T,row.names=F,sep="\t",col.names=F,quote=F)
	write.table(df_reads_tissue,file=outFile_1,append=T,row.names=F,sep="\t",col.names=T,quote=F)

	write.table(line1,file=outFile_2,row.names=F,sep="\t",col.names=F,quote=F)
	write.table(line2,file=outFile_2,append=T,row.names=F,sep="\t",col.names=F,quote=F)
	write.table(df_tpm_tissue,file=outFile_2,append=T,row.names=F,sep="\t",col.names=T,quote=F)
}

