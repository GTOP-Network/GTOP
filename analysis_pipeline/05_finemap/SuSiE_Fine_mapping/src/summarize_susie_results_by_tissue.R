library(optparse)

option_list <- list(
                    make_option(c("-d","--directory_finemap"),type="character",default="",action="store",help="the base location of running fine-mapping, default is ./FineMapping"),
                    make_option(c("-t","--tissue"),type="character",default="",action="store",help="the base location of running fine-mapping, default is ./FineMapping"),
                    make_option(c("-g","--gene_list"),type="character",default="",action="store",help="the file contains gene list used for running fine-mapping, default is picked_asso_list.loc_1000000.txt"))

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

basedir <- opt$directory_finemap # the basic path of susieR analysis
aGene_file <- opt$gene_list
#task <- basename(aGene_file)
setwd(basedir)
cat('Options:\n','basedir:',basedir,'\naGene_file:',aGene_file,'\n')
gene_list <- read.table(opt$gene_list,header=F,sep="\t")
gene_list <- as.character(gene_list[,2])
cs_size <- c()
N_cs <- c()
susie_df <- data.frame(locus_id=c(),variant_id=c(),pip=c(),cs=c(),cs_size=c(),cs_purity=c())
for(idx in 1:length(gene_list)){
        file_name <- paste0("./output/Joint/",opt$tissue,"/",gene_list[idx],"/eVar_GT.SuSiE.txt")
        cat(file_name,"\n")
        if (file.exists(file_name)){
                df <- read.table(file_name,header=T,sep=" ")
				df <- df[!is.na(df$cs_size),]

                if (dim(df)[1]>0){
					cs_size[idx] <- unique(df$cs_size)
					N_cs[idx] <- unique(df$cs)[1]
                    df$locus_id <- gene_list[idx]
                    susie_df <- rbind(susie_df,df)
                }else{
					cs_size[idx] <- NA
					N_cs[idx] <- NA
				}
		}else{
            cs_size[idx] <- NA
			N_cs[idx] <- NA
	}
}

summary_susie <- data.frame(Gene=gene_list,N_CS=N_cs,CS_Size=cs_size)

susie_df$Tissue <- opt$tissue
write.table(susie_df,file=paste0("./SuSiE_Summary/susieR_res.joint_egenes.",opt$tissue,".txt"),quote=F,row.names=F,sep="\t")
