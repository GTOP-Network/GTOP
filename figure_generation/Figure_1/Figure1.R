#==============================================#
# Figure-1#
#==============================================#
setwd("/media/london_A/mengxin/GTOP_code/fig-1/input")

library(data.table)
library(ggplot2)
library(stringi)
library(stringr)
library(dplyr)
library(ggsci)
library(tidyverse)
library(gg.gap)
library(ggbreak)
library(ggpubr)


# Fig.1b:  geno PCA ------------------------------------------------------------------

df_pca <- read.table("Fig1b.txt",header=TRUE,check.names=FALSE)
title <- colnames(df_pca)
outpre <- as.character(args$outpre)
pop_tar <- 'GTOP'
df_gtex <- df_pca[df_pca$superpop == 'GTEX',]
df_agtex <- df_pca[df_pca$superpop == 'GTOP',]
df_1kg <- df_pca[!df_pca$superpop %in% c('GTEX','GTOP'),]
df_pca <- rbind(df_1kg,df_gtex,df_agtex)
df_pca$dataset <- df_pca$superpop
df_pca[! df_pca$superpop %in% c('GTEX','GTOP'),]$dataset <- '1KGP'
colors_pop <- c('gray',"#e7c798","#bd5a45")
pops <- c('1KGP','GTEX','GTOP')
p <- ggplot(df_pca,aes(x=-PC1,y=-PC2)) +
  geom_point(aes(color=dataset),size = 3) +
  xlab('PC1 (7.6%)') + ylab('PC2 (4.8%)') +
  scale_color_manual(breaks=pops,labels=c('1KGP','GTEx','GTOP'),values=colors_pop,
                     guide=guide_legend(override.aes=list(size=3))) +
  theme_classic() + 
  theme(axis.line=element_line(color='black'),
        axis.text=element_text(face='bold'),
        legend.position=c(0.3,0.4),
        legend.text=element_text(face='bold',size=rel(1.1)),
        legend.key=element_blank())


p


# Fig.1c Multidimensional Scaling -----------------------------------------



suppressPackageStartupMessages({
  library(data.table)    
  library(DESeq2)        
  library(matrixStats)   
  library(ggplot2)       
  library(dendextend)    
  library(ggdendro)      
  library(patchwork)     
  library(dplyr)         
})
set.seed(1202)
options(stringsAsFactors = FALSE)

# Theme settings
THEME_SIZE <- 12
TITLE_SIZE <- 14
POINT_SIZE <- 1.8
LINE_WIDTH <- 0.6
ALPHA_NORMAL <- 0.7
ALPHA_OUTLIER <- 0.9

unified_theme <- theme_bw(base_size = THEME_SIZE) +
  theme(
    plot.title = element_text(hjust = 0.5, size = TITLE_SIZE, face = "bold"),
    axis.title = element_text(size = THEME_SIZE, face = "bold"),
    axis.text = element_text(size = THEME_SIZE - 1),
    legend.title = element_text(size = THEME_SIZE, face = "bold"),
    legend.text = element_text(size = THEME_SIZE - 1),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )



# Load data ---------------------------------------------------------------
load("sample_annot_full_1586.RData")

final_annot <- sample_annot_full[, c("sample_id", "Subject", "Tissue","Batch", "Tissue_Color_Code")]
colnames(final_annot) <- c("sample", "individual", "tissue", "batch", "color")
rownames(final_annot) <- final_annot$sample

meta <- final_annot[, c("sample", "individual", "tissue","batch", "color")]
rownames(meta) <- meta$sample

# Color setup -------------------------------------------------------------

all_tissue <- sort(unique(meta$tissue))
tissue_colors <- unique(meta[, c("tissue", "color")])
tissue_colors <- setNames(paste0("#", tissue_colors$color), tissue_colors$tissue)
load("vst_matrix.RData")
#mds_full_corr <- cmdscale(dist(t(vst_matrix)), k=2, eig=TRUE)
load("mds_full_corr.RData")

eigvals_corr <- mds_full_corr$eig
prop1_corr <- round(100 * eigvals_corr[1] / sum(eigvals_corr[eigvals_corr > 0]), 1)
prop2_corr <- round(100 * eigvals_corr[2] / sum(eigvals_corr[eigvals_corr > 0]), 1)

mds_df_corrected <- data.frame(
  MDS1 = mds_full_corr$points[,1],
  MDS2 = mds_full_corr$points[,2],
  sample = colnames(vst_matrix)
)
mds_df_corrected <- merge(mds_df_corrected, meta, by="sample", all.x=TRUE)

# MDS plot  ---------------------------------------------------------------

gg_mds_corr <- ggplot(mds_df_corrected, aes(x=MDS1, y=MDS2, color=tissue)) +
  geom_point(size=POINT_SIZE, alpha=ALPHA_NORMAL) +
  scale_color_manual(values=tissue_colors) +
  unified_theme +
  theme(legend.position="none") +
  labs(
    x=paste0("Coordinate 1 (", prop1_corr, "%)"),
    y=paste0("Coordinate 2 (", prop2_corr, "%)"),
    title="Multidimensional Scaling Analysis",
    color="Tissue"
  )
print(gg_mds_corr)


# Fig.1d: tissue cluster ---------------------------------------------------

expr_by_tissue <- sapply(all_tissue, function(grp) {
  idx <- meta$tissue == grp
  if(sum(idx) > 0) {
    rowMedians(vst_matrix[, idx, drop=FALSE], na.rm = TRUE)
  } else {
    rep(NA, nrow(vst_matrix))
  }
})
expr_by_tissue <- expr_by_tissue[, !apply(is.na(expr_by_tissue), 2, any), drop=FALSE]

if(ncol(expr_by_tissue) > 1) {
  dist_mat_tissue <- as.dist(1 - cor(expr_by_tissue,method="spearman"))
  hc <- hclust(dist_mat_tissue, method="ward.D2")
  dend <- as.dendrogram(hc)
  
  labels_cex(dend) <- 0.7
  labels_colors(dend) <- tissue_colors[labels(dend)]
  labels_dend <- labels(dend)
  label_col <- tissue_colors[labels_dend]
  
  # pdf(paste0(output_prefix,"_test.pdf"), width=10, height=8)
  # par(mar=c(7,4,2,2))
  plot(dend, horiz=FALSE, main="Tissue Hierarchical Clustering",
       xlab="", ylab="Distance (1 - Pearson correlation)", axes=TRUE,
       cex.main=1.2, cex.lab=1.1, cex.axis=0.9, lwd=1.5)
  
  xticks <- seq_along(labels_dend)
  ypos_start <- par("usr")[3] + 0.02 * diff(par("usr")[3:4])
  ypos_end <- ypos_start + 0.03 * diff(par("usr")[3:4])
  
  for(i in seq_along(xticks)){
    segments(xticks[i], ypos_start, xticks[i], ypos_end,
             lwd=6, col=label_col[i], xpd=TRUE)
  }
  
  #dev.off()
}

# === combine plot ===
dend <- as.dendrogram(hc)
dend_data <- dendro_data(dend)
labels_df <- dend_data$labels
labels_dend <- labels(dend)
label_col <- tissue_colors[labels_dend]
labels_df$color <- label_col

segments_df <- dend_data$segments
N <- nrow(labels_df)

seg_length <- 0.005
text_offset <- 0.01

tree_plot <- ggplot() +
  geom_segment(data = segments_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               size = 1.2, color = "grey30", alpha = 0.8
  ) +
  geom_segment(data = labels_df,
               aes(x = x, y = y , xend = x, yend = y ),
               color = labels_df$color,
               size = 4, lineend = "round"
  ) +
  geom_text(data = labels_df,
            aes(x = x, y = y , label = label),
            color = labels_df$color,
            size = 3.8,
            fontface = "bold",
            hjust = 0.5, vjust = 1
  ) +
  labs(title = "Tissue Hierarchical Clustering",
       x = "Tissue",
       y = "Cluster Distance") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(size = 13, face = "bold", color = "black", hjust = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  coord_flip()

gg_mds_corr_adjusted <- gg_mds_corr +
  theme(
    plot.title = element_text(size = 13, face = "bold", color = "black", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 9)
  )

library(cowplot)
combined_plot <- plot_grid(
  gg_mds_corr_adjusted, tree_plot,
  ncol = 2, rel_widths = c(1.2, 1.2)
)


# Fig.1e: variant number and length ----------------------------------------

variant_counts <- data.frame(SNV = c(7247578, 12434256),SV = c(36735, 63708),TR = c(38615,1056782),
                             row.names = c("common", "rare")) %>% 
  tibble::rownames_to_column(var = "variant_type") %>% 
  pivot_longer(cols = -variant_type, names_to = "category",values_to = "count") %>%
  mutate(category=factor(category,levels=c("SNV","TR","SV")))
variant_counts$variant_type <- factor(variant_counts$variant_type,levels = c("rare","common"))

variant_counts <- variant_counts %>%
  mutate(range_group = case_when(count < 70000               ~ "Low",
    count >= 900000 & count < 5000000 ~ "Mid",
    count >= 5000000              ~ "High"))


p1 <- ggplot(variant_counts, aes(y = category, x = count, fill = variant_type)) +
  geom_col(position = "dodge") +
  facet_wrap(~factor(range_group, levels = c("Low", "Mid", "High")),scales = "free_x", nrow = 1) +
  theme_classic() +
  scale_fill_grey() +
  labs(x = "Variant Number (1e4)") +
  theme(axis.text = element_text(color = "black", size = 10),
    axis.ticks = element_line(color = "black"),legend.position = "none")

variant_length <- data.frame(SNV = c(5774020+6375135, 10428404+12348887),
                             SV = c(29293416, 96216228),TR = c(3214629, 17818302),
                             row.names = c("common", "rare")) %>% 
  tibble::rownames_to_column(var = "variant_type") %>% 
  pivot_longer(cols = -variant_type,  names_to = "category", values_to = "count") %>%
  mutate(category=factor(category,levels=c("SNV","TR","SV")))
variant_length$variant_type <- factor(variant_length$variant_type,levels = c("rare","common"))
p2 <- ggplot(variant_length, aes(y = category, x = count, fill = variant_type)) +
  geom_col(position="dodge") +
  scale_x_continuous(labels = function(x) ifelse(x >= 1000000, paste0(x/1000000), x)) +
  theme_classic()+
  scale_fill_grey()+
  labs(x="Varinat length(Mb)")+
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"),
        legend.position = "none")


# Fig.1f: LRS/SRS per Genome variant number --------------------------------
source("geom_boxplot2.R")
plot <- fread("Fig1f.txt") %>% dplyr::filter(Major_Type %in% c("STR","VNTR","SV","Total"))
plot$Major_Type <- factor(plot$Major_Type,levels = c("Total","SV","STR","VNTR"))
ggplot(plot, aes(x = tec, y = Total_count, fill = tec)) +
  geom_boxplot2(width = 0.7, width.errorbar = 0.5) +
  scale_fill_manual(values = c("LRS"="#913628","SRS"="#227e85")) +
  facet_wrap(~ Major_Type, nrow = 1, scales = "free_y") +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = function(x) { max_val <- max(x, na.rm = TRUE)
    magnitude <- 10^floor(log10(max_val))
    max_break <- ceiling(max_val / magnitude) * magnitude
    step <- max_break / 4
    seq(0, max_break, by = step)},labels = scales::comma) +
  labs(x = "Technology", y = "Total Count per Sample") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold"))



# Fig.1g: compare small variants number between LRS and SRS ---------------


df_plot <- readRDS("Fig.1g.RDS")
df_plot$ShareStatus <- factor(df_plot$ShareStatus,levels = c("Specific","Share"))
p1 <- ggplot(df_plot,aes(x=Tech,y=Count_variants,fill=ShareStatus)) + geom_bar(stat = "identity",width=.8,position = position_stack()) + theme_pubr() + 
  scale_fill_manual(breaks = c("Share","Specific"),values = c("#8da0cb","#9d3929"));p1


# Fig.1h: compare total SV count in LRS and SRS ---------------------------


df_sv.lrs_vs_srs <- readRDS("Fig.1h.SV.LRS_vs_SRS.RDS")
p2 <- ggplot(df_sv.lrs_vs_srs,aes(x=Var1,y=Freq,fill=Var2)) + geom_bar(stat = "identity",width=.8) + theme_pubr() + 
  scale_fill_manual(breaks = c(0,1),values = c("#983628","#8da0cb"));p1
df_sv.srs_vs_lrs <- readRDS("Fig.1h.SV.SRS_vs_LRS.RDS")
p3 <- ggplot(df_sv.srs_vs_lrs,aes(x=Var1,y=Freq,fill=Var2)) + geom_bar(stat = "identity",width=.8) + theme_pubr() + 
  scale_fill_manual(breaks = c(0,1),values = c("#247E85","#8da0cb"));p3


cowplot::plot_grid(p2,p3,align="h")


# Fig.1i: compare total TR count in LRS and SRS ------------------------------------------------

col_overlap  <- "#8b9dc5"
col_lrs_spec <- "#933628"
col_srs_spec <- "#247e87"
dat <- fread("Fig1i.txt")%>%
  mutate( group = factor(group, levels = c("Specific","Overlap")),
          TR_type = factor(TR_type, levels = c("2","3","4","5","6","VNTR")))

fill_map <- c("Overlap.LRS"  = col_overlap,"Overlap.SRS"  = col_overlap,
              "Specific.LRS" = col_lrs_spec,"Specific.SRS" = col_srs_spec)
ggplot() +
  geom_col(data = dat %>% filter(tec == "LRS"),
           aes(y = count, x = TR_type, fill = interaction(group, tec)),
           position = "stack" ) +
  geom_col(data = dat %>% filter(tec == "SRS"),
           aes(y = -count, x = TR_type, fill = interaction(group, tec)),
           position = "stack") +
  scale_fill_manual(values = fill_map, name = "Group") +
  scale_y_continuous(labels = function(x) { lx <- abs(x); ifelse(lx >= 1000, paste0(lx/1000, "k"), lx) },
                     limits = c(-500000, 500000)
                     ) +
  geom_hline(yintercept = 0, color = "black") +
  theme_classic(base_size = 13) +
  theme(axis.text.y = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top")

ggplot() +
  geom_col(data = dat %>% filter(tec == "LRS",TR_type=="VNTR"),
           aes(y = count, x = TR_type, fill = interaction(group, tec)),
           position = "stack" ) +
  geom_col(data = dat %>% filter(tec == "SRS",TR_type=="VNTR"),
           aes(y = -count, x = TR_type, fill = interaction(group, tec)),
           position = "stack") +
  scale_fill_manual(values = fill_map, name = "Group") +
  scale_y_continuous(labels = function(x) { lx <- abs(x); ifelse(lx >= 100, paste0(lx/1000, "k"), lx) },
                     limits = c(-50000, 50000)) +
  geom_hline(yintercept = 0, color = "black") +
  theme_classic(base_size = 13) +
  theme(axis.text.y = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top")






