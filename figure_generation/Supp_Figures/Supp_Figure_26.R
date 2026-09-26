#==============================================#
# MAJIQTL repliaction #
# Supp-Figure-26#
#==============================================#
library(data.table)
library(ggplot2)
library(stringi)
library(stringr)
library(dplyr)
library(ggsci)
library(tidyverse)
library(ggpubr)
library(magrittr)
library(scales)
library(ggrastr)
library(ggupset)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig26. R3.2.MAJIQTL-replication")

# Supp.Fig.26a: sGene overlapping --------------------------------------------------------
df <- fread("./input/Supp_Fig26a.txt")
df$category <- factor(
  df$category,
  levels = c( "Tensor-QTL Unique", "majiQTL Unique","Overlap")
)
ggplot(df, aes(x = tissue, y = proportion, fill = category)) +
    geom_bar(stat = "identity", width = 0.75) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_y_continuous(labels = function(x) abs(x),name = "Proportion (%)") +
    scale_fill_manual(
      values = c("Overlap" = "#C48A6B","Tensor-QTL Unique" = "#7A8B8B","majiQTL Unique" = "#8B9A6E")) +
    facet_wrap(.~ gtype, ncol = 1)+
    labs(x = "Tissue",fill = "Category") +
    theme_pubr() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14)
    )
# Supp.Fig.26b: scatter of slope --------------------------------------------------------
df <- fread("./input/Supp_Fig26b.txt")
df$gtype<-factor(df$gtype,levels = c("SNP","TR","SV"))
ggplot(df,aes(x = slope_tensor, y = slope_majiqtl, color = class)) +
    geom_point(size = 0.8, alpha = 0.3) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed") +
    stat_cor(aes(color = class),method = "spearman",label.x.npc = "left",label.y.npc = c(0.98, 0.90),size = 4,r.digits = 3,show.legend = FALSE) +
    scale_color_manual(values=c("overlap"="#C48A6B", "tensorqtlunique"="#7A8B8B"))+
    facet_wrap(.~ gtype, ncol = 1)+
    theme_pubr() +
    labs(x = "TensorQTL slope",y = "MAJIQTL slope",color = "sGene class")
