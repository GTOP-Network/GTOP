#==============================================#
# Proportion of expression variance captured by the PCs  #
# Supp-Figure-25#
#==============================================#

library(data.table)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(ggplot2)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig25.Figure.S19-PCs-explain")

# Supp.Fig.25a juQTL ------------------------------------------------------
tissue_factor_summary<-readRDS("input/Fig_S25a.input.RDS")

ggplot(
  tissue_factor_summary,aes(x = Factor,y = Tissue,fill = Variance_fraction)) +
  geom_tile() + geom_text(aes(label = sprintf("%.1f%%",Variance_fraction * 100)),size = 3) +
  scale_fill_gradientn(
    colours = c("#feebe2","#fbb4b9","#f768a1","#c51b8a","#7a0177")) +
  theme_bw() +
  labs(x = NULL,y = NULL,fill = "Variance\nfraction") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
# Supp.Fig.25b eQTL--------------------------------------------------------
adjR2_mat<-readRDS("input/supp25b.adjR2_mat.RDS")
row_ha<-readRDS("input/supp25b.row_ha.RDS")


col_fun <- colorRamp2(
  c(0, 0.25, 0.5, 0.75, 1),
  c("#feebe2","#fbb4b9","#f768a1","#c51b8a","#7a0177")
)
p1<-Heatmap(
  adjR2_mat,
  name = "Adj.R²",
  col = col_fun,
  
  na_col = "grey80",
  left_annotation = row_ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 10)
);p1


# Supp.Fig.25c Cumulative PVE-------------------------------------------------


df <- fread("./input/Figure.S25.txt")
figlist <- list()
type <- unique(df$qtltype)[c(2, 1, 3)]

for (qtltype_t in type) {
  plotdf <- df %>% dplyr::filter(qtltype==qtltype_t)
  
  plotdf$tissue <- factor(plotdf$tissue, levels = plotdf$tissue)
  
  
  scale_factor <- max(plotdf$num_PCs) / max(plotdf$Cumulative_PVE)
  
  figlist[[qtltype_t]] <- ggplot(plotdf, aes(x = tissue)) +
    geom_bar(aes(y = Cumulative_PVE, fill = "pro_explain"), 
             stat = "identity") +
    geom_point(aes(y = num_PCs / scale_factor, color = "PCnumber"),
               shape = 19) +
    scale_y_continuous(
      name = "Proportion of expression variance captured by the PCs",
      limits = c(0, 0.8),
      sec.axis = sec_axis(~. * scale_factor, 
                          name = "Number of PCs in xQTL mapping")
    ) +
    scale_fill_manual(values = c("pro_explain" = "#c9c9c9")) +
    scale_color_manual(values = c("PCnumber" = "#9d3929")) +
    theme_pubr()+
    theme(axis.text.x = element_blank(),
          legend.position = "none")+
    ggtitle(qtltype_t)
}

cowplot::plot_grid(plotlist = figlist[type], ncol = 1)


