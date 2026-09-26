#==============================================#
# GTOP-GTEX-eQTL-correlation-SNV #
# Supp-Figure-31#
#==============================================#
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(data.table)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig31.Figure.S23-GTOP-GTEX-eQTL-correlation-SNV")
# Supp.Fig.31a (GTOP lead to GTEx)correlation of eQTL effects between GTOP and GTEx -----------


plotdf <- fread("./input/Figure S31a.txt.gz")
figlist <- list()
for (tispair in unique(plotdf$tissuepairs)) {
  plotdf_t <- plotdf %>% dplyr::filter(tissuepairs==tispair)

    figlist[[tispair]] <- ggplot(plotdf_t ) +
    geom_point(aes(x = slope_gtop, y = slope_gtex, color = match_type), size = 1, alpha = 0.75) +
    geom_vline(xintercept = 0, color = "black", linetype="dashed") +
    geom_hline(yintercept = 0, color = "black", linetype="dashed") +
    scale_color_manual(values = c("direct" = "#cb9b5e","ld_proxy" = "#6B8E9F"),drop = FALSE) +
    stat_cor(aes(x = slope_gtop, y = slope_gtex), method = "spearman",label.x = .1,label.y = -1.5) +
    scale_x_continuous(limits = c(-2, 2)) +
    scale_y_continuous(limits = c(-2, 2)) +
    xlab("Slope of eQTL from GTOP") +
    ylab("Slope of eQTL from GTEx / LD-proxy aligned GTEx") +
    ggtitle(tispair) +
    theme_pubr() +
    theme(
      panel.grid = element_blank(),
      aspect.ratio = 1,
      legend.title = element_blank()
    )
}


cowplot::plot_grid(plotlist = figlist, ncol = 3)


# Supp.Fig.31a (GTEx lead to GTOP)correlation of eQTL effects between GTOP and GTEx -----------


plotdf <- fread("./input/Figure S31b.txt.gz")
figlist <- list()
for (tispair in unique(plotdf$tissuepairs)) {
  plotdf_t <- plotdf %>% dplyr::filter(tissuepairs==tispair)

    figlist[[tispair]] <- ggplot(plotdf_t) +
    geom_point(aes(x = slope_gtex, y = slope_gtop, color = match_type), size = 1, alpha = 0.75) +
    geom_vline(xintercept = 0, color = "black", linetype="dashed") +
    geom_hline(yintercept = 0, color = "black", linetype="dashed") +
    scale_color_manual(values = c("direct" = "#cb9b5e","ld_proxy" = "#6B8E9F"),drop = FALSE) +
    stat_cor(aes(x = slope_gtex, y = slope_gtop), method = "spearman",label.x = .1,label.y = -1.5) +
    scale_x_continuous(limits = c(-2, 2)) +
    scale_y_continuous(limits = c(-2, 2)) +
    xlab("Slope of eQTL from GTEx") +
    ylab("Slope of eQTL from GTOP / LD-proxy aligned GTOP") +
    ggtitle(tispair) +
    theme_pubr() +
    theme( panel.grid = element_blank(), aspect.ratio = 1, legend.title = element_blank())

}


cowplot::plot_grid(plotlist = figlist, ncol = 3)

# Supp.Fig.31c rb of eQTL between GTOP and GTEx -----------

df <- fread("./input/Figure S31c.txt")
df$class<-factor(df$class,levels = c("GTOP2GTEx","GTEx2GTOP"))
p <- ggplot(
  df,
  aes(x = xlabel, y = r, color = class)) +
  geom_errorbar( aes(ymin = ci_low, ymax = ci_high), width = 0.15, position = position_dodge(width = 0.45), linewidth = 0.5) +
  geom_point( position = position_dodge(width = 0.45), size = 2.6) +
  scale_color_manual(
    values = c( "GTOP2GTEx" = "#A66A5B", "GTEx2GTOP" = "#4E6E8E")) +
  coord_cartesian(ylim = c(0.78, 0.92)) +
  theme_pubr(base_size = 12) +
  labs( x = NULL, y = expression(R[b]), color = NULL) +
  theme( axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top");p

