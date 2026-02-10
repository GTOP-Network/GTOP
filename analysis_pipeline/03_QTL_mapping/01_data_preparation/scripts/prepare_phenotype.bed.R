# format alt TSS quant results
library(optparse)

option_list <- list(
		    make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)

# -- functions
extract_tss <- function(x){
	chrom <- x[1]
	pos1 <- x[2]
	pos2 <- x[3]
	strand <- x[4]
	if(strand=="+"){
		tss = pos1
	}else{
		tss = pos2
	}
	return(tss)
}


# -- load gene annotation
#df_gene <- fread("./input/gencode.v47.gene_info.txt",header=F,sep="\t")
df_gene <- fread("./input/GTOP.gene_info.txt",header=F,sep="\t")
names(df_gene) <- c("chrom","start","end","strand","Gene","Name","GeneType")
df_gene$Name <- NULL
df_gene$GeneType <- NULL

df_gene$TSS <- apply(df_gene[,c(1,2,3,4)],1,extract_tss)
df_gene %<>% select(Gene,chrom,TSS)
df_gene$pos0 <- as.numeric(df_gene$TSS) - 1

# -- load genotype individuals
#df_gt <- fread("/lustre/home/xdzou/2024-10-21-GTBMap/2025-02-13-expression-quant/output/genotype/GMTiP.LRS131_SV.GT.psam",skip=1,header=F,sep="\t")
#names(df_gt) <- c("INDS","SEX")
df_gt <- fread("./input/new_header.txt",header=F,sep="\t")
names(df_gt) <- c("INDS")

# -- load file
infile <- paste0("./output/phenotype/",opt$tissue,".phenotype.txt")
dat <- as.data.frame(fread(infile,header=T,sep="\t"),stringsAsFactors=F,check.names=F)
inds <- intersect(colnames(dat)[-1],df_gt$INDS)
cat("Data shape (before): ",dim(dat),"\n")
dat <- merge(dat,df_gene,by="Gene")
cat("Data shape (after): ",dim(dat),"\n")

dat %<>% select(all_of(c("chrom","pos0","TSS","Gene",inds)))
kept_chrs <- paste0("chr",1:22)
kept_chrs <- c(kept_chrs,"chrX")
dat.f <- dat %>% filter(chrom %in% kept_chrs)

names(dat.f) <- c("#chr","start","end","phenotype_id",inds)

# -- output
write.table(dat.f,file=paste0("./output/phenotype/",opt$tissue,".bed"),quote=F,sep="\t",row.names=F)


