library(optparse)

option_list <- list(
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)


input_file <- ""
output_file <- ""
indir <- "/lustre/home/xdzou/2024-10-21-GTBMap/2025-06-15-joint-split-gt-QTL/output/QTL_mapping/all_joint/cis_QTL/text_format/"
input_file <- paste0(indir,opt$tissue,".cis_eQTL.all_pairs.txt.gz")
output_file <- paste0(indir,opt$tissue,".cis_eQTL.all_pairs.add_tstat.txt")

# load sample count in each tissue, and atGene in each tissue
currDir <- "/lustre/home/xdzou/2024-10-21-GTBMap/2025-06-15-joint-split-gt-QTL/output/QTL_mapping/all_joint/2025-06-18-Finemap_susie/"
sample_n <- read.table(paste0(currDir,"input/SampleSize_by_tissue.txt"),header=F,sep="\t",stringsAsFactors=F)
names(sample_n) <- c("Tissue","Size")

#atGene_df <- read.table(paste0(currDir,"input/tissue_gene/",opt$tissue,"_gene_list.txt"),header=F,sep="\t",stringsAsFactors=F)
#atGenes <- as.character(atGene_df$V2)
#rm(atGene_df)

# load raw cis associations
dat <- as.data.frame(fread(input_file,header=T,sep="\t"),stringsAsFactors=F)
cat("dat dimension:",dim(dat),"\n")
dat$pval_nominal <- as.numeric(dat$pval_nominal)


# filter chromosome X Y MT
N <- sample_n[sample_n$Tissue==opt$tissue,2]

# two groups by slope sign
dat.plus <- dat %>% filter(slope>0)
dat.minor <- dat %>% filter(slope<=0)
rm(dat)
dat.plus$N <- N
dat.minor$N <- N
dat.plus$zscore <- -qnorm(dat.plus$pval_nominal/2)
dat.minor$zscore <- qnorm(dat.minor$pval_nominal/2)

dat.plus$tstat <- qt(dat.plus$pval_nominal,dat.plus$N,lower=FALSE)
dat.minor$tstat <- -qt(dat.minor$pval_nominal,dat.minor$N,lower=FALSE)

dat.plus$N <- NULL;dat.plus$zscore <- NULL
dat.minor$N <- NULL;dat.minor$zscore <- NULL
dat.out <- rbind(dat.plus,dat.minor)
rm(dat.plus,dat.minor)
cat("Fields in final mat:",dim(dat.out),"\n")

fwrite(dat.out,file=output_file,quote=F,sep="\t",row.names=F)
