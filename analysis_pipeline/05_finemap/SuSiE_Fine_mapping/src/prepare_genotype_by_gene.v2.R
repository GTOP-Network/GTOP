library(optparse)

option_list <- list(
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue"),
	make_option(c("-g","--genes"),type="character",default="NA",action="store",help="specify a gene list")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)

setwd("/lustre/home/xdzou/2024-10-21-GTBMap/2025-07-03-Finemap_susie")
genotype_bed <- "./input/GMTiP_LRS.joint_var.dosage.RDS"
#genelist <- paste0("./input/tissue_gene2/",opt$tissue,"_gene_list.txt")
genelist <- opt$genes

df_geno <- readRDS(genotype_bed)

df_genes <- fread(genelist,header=F,sep="\t")
names(df_genes) <- c("Tissue","Gene")

for(i in 1:dim(df_genes)[1]){
	gene <- df_genes$Gene[i]
	snp_file <- paste0("./output/Finemapping_eQTL_sv_snv/",opt$tissue,"/",gene,"/snp_list.txt")
	df_snp <- fread(snp_file,header=F,sep="\t")
	snp_list <- df_snp$V1
	rm(df_snp)
	df_gt <- df_geno %>% filter(variant_id %in% snp_list)
	outdir <- paste0("./output/Finemapping_eQTL_sv_snv/",opt$tissue,"/",gene)

	if(file.exists(outdir)){
		cat(outdir," already exist!\n")
		fwrite(df_gt,file=paste0(outdir,"/eVar_GT.vcf"),quote=F,sep="\t",row.names=F,na=".")
	}else{
		dir.create(outdir)
		fwrite(df_gt,file=paste0(outdir,"/eVar_GT.vcf"),quote=F,sep="\t",row.names=F,na=".")
	}
}



