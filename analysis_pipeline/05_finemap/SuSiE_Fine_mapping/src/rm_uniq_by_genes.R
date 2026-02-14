library(optparse)

option_list <- list(
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)

dir1 <- "/flashfs1/scratch.global/xdzou/Fine_map_susie/output/TR"
setwd(dir1)
gene_out_dir <- "/lustre/home/xdzou/2024-10-21-GTBMap/2026-02-03-GTOP_eQTL_Finemap_susie/input/"
file_suc <- paste0(opt$tissue,"_success.txt")
file_all <- paste0(opt$tissue,".all.txt")
df_suc <- fread(file_suc,header=F,sep="\t")
df_all<- fread(file_all,header=F,sep="\t")

genes <- setdiff(df_all$V1,df_suc$V1)

write.table(data.frame(Gene=genes),file=paste0(gene_out_dir,opt$tissue,".remained_genes.TR.txt"),quote=F,row.names=F,col.names=F)

for(i in 1:length(genes)){
	cat(genes[i],'\n')
	file_gt <- paste0(opt$tissue,"/",genes[i],"/eTR_GT.vcf")
	df_gt <- fread(file_gt,header=T,sep="\t")
	dim(df_gt)
	df_gt %<>% distinct(variant_id,.keep_all=T)
	dim(df_gt)

	write.table(df_gt,file=paste0(opt$tissue,"/",genes[i],"/eTR_GT.nr.vcf"),quote=F,sep="\t",row.names=F)
}
