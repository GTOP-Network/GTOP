#==============================================#
# phenotype-pca-sQTL #
# Supp-Figure-23#
#==============================================#
library(data.table)
library(dplyr)
library(magrittr)
library(PCAForQTL)
library(cowplot)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig23.FigureS17.phenotype-pca-sQTL")
# Supp.Fig.23 splicing PCA ------------------------------------------------


reslist <- readRDS("./input/Figure.S23.rds")

figlist <- list()

for (i in 1:length(reslist)) {
  K_pc_elbow <- runElbow(prcompResult=reslist[[i]])
  figlist[[i]] <- makeScreePlot(reslist[[i]] ,labels=c("Elbow"),values=c(K_pc_elbow),titleText=names(reslist)[i])
}

cowplot::plot_grid(plotlist = figlist, ncol = 3)
