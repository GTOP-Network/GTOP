#==============================================#
# Extended Fig.5 #
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
library(pheatmap)
setwd("/media/london_A/mengxin/GTOP_code/extend/extend_145/input")


# Extended Fig.6a:  Visualization of the effect size of the same QTL pair in different tissues ------------------------------------------------------------------

cor <- fread("ExtendFig5a.txt",data.table = FALSE)
rownames(cor) <- cor[, 1]
cor <- cor[, -1]
rownames(cor) <- gsub(".mashr.txt", "", rownames(cor))
colnames(cor) <- gsub(".mashr.txt", "", colnames(cor))
pheatmap(cor,
         color = colorRampPalette(c("#F7FBFF", "#08306B"))(10),
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = NA,
         cellwidth = 15,
         cellheight = 8,
         fontsize_row = 7,
         fontsize_col = 9)


# Extended Fig.6b: tissue specific QTL number  --------------------------------------------------------

plot_df <- fread("ExtendFig5b.txt") %>% mutate(Specific_Tissue = reorder(Specific_Tissue, -total))
plot_df$QTL <- factor(plot_df$QTL, levels = c( "SV","TR", "SNV"))

ggplot(plot_df, aes(x = Specific_Tissue, y = n, fill = QTL)) +
  geom_col(width = 0.7) +
  labs(x = "Tissue",y = "Number of tissue-specific QTL",fill = "QTL") +
  theme_classic() +
  scale_fill_manual(values = c("SV"="#c88565","TR"="#931e2a","SNV"="#ebd1bf"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title = element_text(color = "black"))

# Extended Fig.6c: tissue-sharing with TSS  --------------------------------------------------------

dat <- fread("ExtendFig5c.txt") %>% dplyr::filter(QTL %in% "SNV")
ggplot(dat, aes(x = TSS_dis, group = tissue_group, color = tissue_group)) +
  geom_density(adjust = 1.5, linewidth = 1) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(x = "TSS distance (Mb)",y = "Density",
       color = "# tissues\neQTLs active in") +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.7),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12),
        legend.position = c(0.9,0.7),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) +
  scale_x_continuous(labels = function(x) x/1e6,limits = c(-1e6,1e6))+
  facet_wrap(~QTL,scales = "free_y")


# tissue sharing QTL with TSS  --------------------------------------------

df_binned <- fread("ExtendFig5d.txt")
ggplot(df_binned,aes(x = tissue_median, y = TSS_median, color = QTL)) +
  geom_smooth(se = T, linewidth = 1.1, method = "loess", span = 0.8) +
  scale_color_manual(values = c("SV"="#227a7e","TR"="#931e2a","SNV"="#7c8bad","Control"="grey"))+
  scale_x_continuous(breaks = seq(floor(min(df_binned$tissue_median, na.rm = TRUE)),
                                  ceiling(max(df_binned$tissue_median, na.rm = TRUE)), by = 1) )+
  labs(x = "#Tissues eGenes expressed in", y = "Distance of QTL to TSS (Mb)") +
  theme_classic() +
  theme( axis.text = element_text(color = "black", size = 10),legend.position = "top")

