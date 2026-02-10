library(optparse)

option_list <- list(
	make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a tissue")
)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(data.table)
library(dplyr)
library(magrittr)
library(PCAForQTL)


input_file <- paste0("./output/phenotype/",opt$tissue,".phenotype.bed")


# load raw cis associations
dat <- as.data.frame(fread(input_file,header=T,sep="\t"),check.names=F,stringsAsFactors=F)
dat <- dat[,c(-1,-2,-3)]
#names(dat) <- c("Gene",colnames(dat)[-1])

# remove genes with variance of zero
dat$Var <- apply(dat[,-1],1,var)
dat <- dat[order(-dat$Var),]
dat <- dat[1:8000,]
dat$Var <- NULL
rownames(dat) <- dat$Gene

expr <- t(dat[,-1])

# perform PCA
prcompRes <- prcomp(expr,center=TRUE,scale.=TRUE)
PCs <- prcompRes$x
dim(PCs)

# save PCA results
saveRDS(prcompRes,file=paste0("./output/pca_res/",opt$tissue,".prcompRes.RDS"))


# choose PCs by runElbow
K_pc_elbow <- runElbow(prcompResult=prcompRes)

# choose PCs by runBE
RNGkind("L'Ecuyer-CMRG")
set.seed(100)
resultRunBE<-PCAForQTL::runBE(expr,B=20,alpha=0.05)
K_pc_BE <- resultRunBE$numOfPCsChosen


makeScreePlot(prcompRes,labels=c("Elbow","BE"),values=c(K_pc_elbow,K_pc_BE),titleText=as.character(opt$tissue))
ggplot2::ggsave(file=paste0("./output/pca_res/",opt$tissue,".pdf"),width=16,height=11,unit="cm")



