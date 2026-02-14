library(optparse)

option_list <- list(
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)


# global settings
fields <- c("Gene","gene_chr","gene_start","gene_end","strand","N_snp","dist","SNP","SNP_chr","SNP_pos1","SNP_pos2","noml_pval","r_squared","slope","best_hit")

extract_snpid <- function(x){
	return(strsplit(x,split="_",fixed=T)[[1]][1])
}
basedir <- "/lustre/home/xdzou/2022-08-05-altTSS_QTL-Project/"
qtl_file <- paste0(basedir,"output/QTL_mapping/cis_QTLs/",opt$tissue,".cis_QTL_all.txt")
susie_file <- paste0(basedir,"2022-10-05-heritability/output/susieR_res.all_gene.",opt$tissue,".txt")

output_file <- paste0(basedir,"2022-10-05-heritability/output/susieR_res.add_dist.",opt$tissue,".txt")

cat("QTL file: ",qtl_file,"\n")
cat("SuSiE file: ",susie_file,"\n")
# load qtl file
df_qtl <- fread(qtl_file,header=F,sep=" ")
names(df_qtl) <- fields

df_qtl %<>% select(Gene,SNP,dist)
df_qtl$pair <- paste(df_qtl$Gene,df_qtl$SNP,sep=":")
df_qtl %<>% select(pair,dist)
# load susieR results

df_susie <- fread(susie_file,header=T,sep="\t")
df_susie %<>% select(locus_id,variant_id,pip) %>% filter(!is.na(pip))
df_susie$SNP <- sapply(df_susie$variant_id,extract_snpid)
df_susie$pair <- paste(df_susie$locus_id,df_susie$SNP,sep=":")
df_susie %<>% select(pair,pip)

df_susie <- merge(df_susie,df_qtl,by="pair",all.x=T)



cat("output file: ",output_file,"\n")
write.table(df_susie,file=output_file,quote=F,sep="\t",row.names=F)
