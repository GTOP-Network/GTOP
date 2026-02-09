#==============================================#
# ASE-eQTL-correlation #
# Supp-Figure-22#
#==============================================#

library(ggplot2)
library(ggpubr)
library(data.table)
library(tidyverse)
library(cowplot)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig22")


# Supp.Fig.22 ASE-eQTL-correlation ----------------------------------------


df_plot <- fread("./input/Figure S22.txt", sep = "\t")


figlist <- list()
for (tis in unique(df_plot$tissue)) {
  df_tmp <- df_plot %>% dplyr::filter(tissue==tis)
  cor_test <- cor.test(df_tmp$slope, df_tmp$logafc)
  rho <- cor_test$estimate[1]
  p_value <- cor_test$p.value[1]

  
  
  figlist[[tis]] <- ggplot(df_tmp)+
    geom_point(aes(x=logafc,y=slope),
               fill="grey", size=2.5, shape=21)+
    geom_text(
      aes(x=0, y=max(slope)-0.5),
      label = paste("rho =", round(rho, 2), ", p =", format(p_value, scientific = TRUE, digits = 3)),
      hjust = 1, vjust = 1, size = 3, color = "black")+
    geom_smooth(aes(x = logafc, y = slope),
                method = "lm", se = FALSE,
                color = "black", size = 0.5, linetype = "dashed") +
    xlab("ASE aFC") + 
    ylab("Tensorqtl slope")+
    theme_pubr()+
    ggtitle(sprintf("%s", tis))
  
}

plot_grid(plotlist = figlist, cols = 3)
