#==============================================#
# phenotype-pca-tuQTL #
# Supp-Figure-18#
#==============================================#
library(data.table)
library(dplyr)
library(magrittr)
library(PCAForQTL)
library(cowplot)

setwd("/path/to/GTOP_code/supp/supp_fig18")
# Supp.Fig.18 transcript usage PCA ----------------------------------------

reslist <- readRDS("./input/Figure S18.rds")
figlist <- vector("list", length(reslist))
for (i in 1:length(reslist)) {
  K_pc_elbow <- runElbow(prcompResult=reslist[[i]])
  figlist[[i]] <- makeScreePlot(reslist[[i]] ,labels=c("Elbow"),values=c(K_pc_elbow),titleText=names(reslist)[i])
}

cowplot::plot_grid(plotlist = figlist, ncol = 3)





