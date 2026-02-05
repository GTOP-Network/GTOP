#==============================================#
# GWAS-fdQTL #
# Supp-Figure-30#
#==============================================#


library(ggplot2)
library(tidyverse)
library(ggpubr)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig30")


# Supp.Fig.30a fd-QTL- GWAS association -----------------------------------------
fd_EAS_GWAS_res3 <- fread("./input/supp_fig30a_data.txt")

pheatmap::pheatmap(fd_EAS_GWAS_res3[, 3:ncol(fd_EAS_GWAS_res3)], border_color = NA, 
                   treeheight_row = 15, treeheight_col = 15,
                   color = colorRampPalette(c("white", "red"))(100),
                   labels_row=paste0(fd_EAS_GWAS_res3$variant_id, "_", 
                                     fd_EAS_GWAS_res3$gene_symbol) )


# Supp.Fig.30b fd-QTL- trait number --------------------------------------------
fd_EAS_GWAS_res2 <- fread("./input/supp_fig30b_data.txt")

fd_EAS_GWAS_res2$count <- factor(as.character(fd_EAS_GWAS_res2$count),
                                 levels = as.character(unique(fd_EAS_GWAS_res2$count)))

ggplot(fd_EAS_GWAS_res2, aes(qtl_count, count)) +
  geom_col() +
  theme_classic() +
  labs(x="Number of fd-QTLs", y="Number of disease/traits")


# Example Disease associated-fd-QTL --------------------------------------------
gene_qtl_subdf <- fread("./input/supp_fig30c_data.txt")

gene_qtl_subdf$cs <- factor(gene_qtl_subdf$cs, levels = unique(gene_qtl_subdf$cs))
ggplot(gene_qtl_subdf[-log10(gene_qtl_subdf$p_value) > 1, ], aes(x=maf_GTEx, y=maf_GTOP) ) +
  ggpointdensity::geom_pointdensity( alpha=1, aes(size=-log10(p_value), fill=cs), adjust=1, shape=21 )+
  scale_color_viridis_b()+
  ggdensity::geom_hdr_lines( linetype="dashed", linewidth=0.5)+
  labs(x="maf_GTEx", y="maf_GTOP")+
  theme_classic()+
  theme(
    axis.line = element_line(color="black"),
    axis.ticks = element_line(color="black"),
    axis.text.x = element_text(color="black"),
    axis.title = element_text(size=rel(1.7)),
    legend.title = element_blank(),
    legend.text = element_text(size=rel(1.2)),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm")
  ) +
  ggh4x::coord_axes_inside(labels_inside = F) +
  scale_fill_manual(values = c("-"="white", "L1"="#f0c479", "L2"="#e19b64"))
