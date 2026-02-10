#==============================================#
# GTOP-GTEX-eQTL-correlation-SNV #
# Supp-Figure-23#
#==============================================#
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(data.table)

setwd("/path/to/GTOP_code/supp/supp_fig23")
# Supp.Fig.24a correlation of eQTL effects between GTOP and GTEx -----------


plotdf <- fread("./input/Figure S23a.txt")
figlist <- list()
for (tis in unique(plotdf$tissue_gtop)) {
  plotdf_t <- plotdf %>% dplyr::filter(tissue_gtop==tis)
  
  
  figlist[[tis]] <- ggplot(plotdf_t, aes(x=slope_gtop, y=slope_gtex)) + #zscore
    geom_point(color="#cb9b5e", size=1)+xlab("Slope of eQTL from GTOP")+ylab("Slope of eQTL from GTEx")+
    geom_vline(xintercept = 0, "black")+
    geom_hline(yintercept = 0, "black")+
    stat_cor(method = "pearson", label.x =-3, label.y = 3)+
    scale_x_continuous(limits = c(-4, 4))+
    scale_y_continuous(limits = c(-4, 4))+
    ggtitle(sprintf("%s", tis))+
    theme_bw()+
    theme(panel.grid = element_blank(), aspect.ratio = 1)
}


cowplot::plot_grid(plotlist = figlist, ncol = 3)


# Supp.Fig.24b rb of eQTL between GTOP and GTEx -----------

df <- readRDS("input/Figure S23b.rds")

ggplot(df, aes(x = tis1, y = tis2, fill = r_b)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), 
            size = 3, 
            color = "black") +
  scale_fill_gradient2(
    high = "#2171B5",
    mid = "#9ECAE1",
    low = "white",
    midpoint = median(df$r_b, na.rm = TRUE),
    name = "Rb"
  )+
  labs(
    x = "Discovery",
    y = "Replication",
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  coord_fixed(ratio = 1)

