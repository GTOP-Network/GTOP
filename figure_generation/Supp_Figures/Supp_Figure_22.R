#==============================================#
# expression PCA #
# Supp-Figure-22#
#==============================================#

library(data.table)
library(dplyr)
library(magrittr)
library(PCAForQTL)
library(cowplot)

setwd("/path/to/GTOP_code/supp/supp_fig22")

# Supp.Fig.22 expression PCA ----------------------------------------------


reslist <- readRDS("./input/Figure.S22.rds.gz")

figlist <- list()

for (i in 1:length(reslist)) {
  K_pc_elbow <- runElbow(prcompResult=reslist[[i]])
  figlist[[i]] <- makeScreePlot(reslist[[i]] ,labels=c("Elbow"),values=c(K_pc_elbow),titleText=names(reslist)[i])
}

cowplot::plot_grid(plotlist = figlist, ncol = 3)


