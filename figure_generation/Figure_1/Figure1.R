#==============================================#
# Figure-1#
#==============================================#
setwd("/path/to/GTOP_code/fig-1/input")

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
set.seed(2026)
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


load("sample_annot_full_1613.RData")

meta <- sample_annot_full[, c("sample_id", "Subject", "Tissue","Batch", "Tissue_Color_Code")]
colnames(meta) <- c("sample", "individual", "tissue", "batch", "color")
rownames(meta) <- meta$sample

all_tissue <- sort(unique(meta$tissue))
tissue_colors <- unique(meta[, c("tissue", "color")])
tissue_colors <- setNames(paste0("#", tissue_colors$color), tissue_colors$tissue)
load("expr_by_tissue.RData")
load("mds_full_corr.RData")

eigvals_corr <- mds_full_corr$eig
prop1_corr <- round(100 * eigvals_corr[1] / sum(eigvals_corr[eigvals_corr > 0]), 1)
prop2_corr <- round(100 * eigvals_corr[2] / sum(eigvals_corr[eigvals_corr > 0]), 1)

mds_df_corrected <- data.frame(
  MDS1 = mds_full_corr$points[,1],
  MDS2 = mds_full_corr$points[,2],
  sample = rownames(mds_full_corr$points)
)
mds_df_corrected <- merge(mds_df_corrected, meta, by="sample", all.x=TRUE)

gg_mds_corr <- ggplot(mds_df_corrected, aes(x=MDS1, y=MDS2, color=tissue)) +
  geom_point(size=3, alpha=ALPHA_NORMAL) +
  scale_color_manual(values=tissue_colors) +
  unified_theme +
  theme(legend.position="none") +
  labs(
    x=paste0("Coordinate 1 (", prop1_corr, "%)"),
    y=paste0("Coordinate 2 (", prop2_corr, "%)"),
    #title="Multidimensional Scaling Analysis",
    color="Tissue"
  )+
  theme_classic() +   
  theme(
    legend.position = "right"   
  )
print(gg_mds_corr)


# Fig.1d: tissue cluster ---------------------------------------------------


if(ncol(expr_by_tissue) > 1) {
  dist_mat_tissue <- as.dist(1 - cor(expr_by_tissue, method = "spearman"))
  hc <- hclust(dist_mat_tissue, method = "average")
  dend <- as.dendrogram(hc)
  
  labels_cex(dend) <- 0.7
  labels_colors(dend) <- tissue_colors[labels(dend)]
  labels_dend <- labels(dend)
  label_col <- tissue_colors[labels_dend]
  
  plot_dend <- function(){
    par(mar = c(4, 10, 2, 2)) 
    
    max_h <- attr(dend, "height")
    dend_noLabels <- dend
    labels(dend_noLabels) <- rep("", length(labels_dend))
    
    plot(dend_noLabels, horiz = TRUE, main = "",
         xlab = "Cluster distance", ylab = "", axes = TRUE,
         xlim = c(0, max_h),
         cex.main = 1.2, cex.lab = 1.1, cex.axis = 0.9, lwd = 1.5)
    
    mtext("b", side = 3, line = 0.5, at = par("usr")[1], adj = 1.5, cex = 1.4, font = 2)
    
    yticks <- seq_along(labels_dend)
    dot_x <- par("usr")[1] + 0.01 * diff(par("usr")[1:2])
    
    for(i in seq_along(yticks)){
      points(dot_x, yticks[i], pch = 19, col = label_col[i], cex = 1.1, xpd = TRUE)
      text(dot_x, yticks[i], labels_dend[i], col = label_col[i], cex = 0.7,
           pos = 2, offset = 0.3, xpd = TRUE)
    }
  }
}
plot_dend()
# Fig.1e: variant number and length ----------------------------------------

variant_counts <- data.frame(SNV = c(7748122, 14044743),SV = c(33219,53753),TR = c(38615,1056776),
                             row.names = c("common", "rare")) %>% 
  tibble::rownames_to_column(var = "variant_type") %>% 
  pivot_longer(
    cols = -variant_type,  
    names_to = "category", 
    values_to = "count") %>%
  mutate(category=factor(category,levels=c("TR","SNV","SV")))
variant_counts$variant_type <- factor(variant_counts$variant_type,levels = c("rare","common"))

p1 <- ggplot(variant_counts, aes(y = category, x = count/1000, fill = variant_type)) +
  geom_col(position = "dodge") + 
  scale_x_log10() +
  theme_classic() +scale_fill_manual(values = c("#c9c9cb","#939eb2"))+
  labs(x = "Variant Number (1e3)") +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.ticks = element_line(color = "black"),
    legend.position = "none")

variant_length <- data.frame(SNV = c(6031603+7320975, 10958697+16907410),
                             SV = c(25353740, 88516072),TR = c(3214629, 17818221),
                             row.names = c("common", "rare")) %>% 
  tibble::rownames_to_column(var = "variant_type") %>% 
  pivot_longer(
    cols = -variant_type,  
    names_to = "category", 
    values_to = "count") %>%
  mutate(category=factor(category,levels=c("TR","SNV","SV")))
variant_length$variant_type <- factor(variant_length$variant_type,levels = c("rare","common"))

p2 <- ggplot(variant_length, aes(y = category, x = count, fill = variant_type)) +
  geom_col(position="dodge") +
  #geom_text(aes(label = ifelse(count >= 1e6, 
  #                            paste0(round(count/1e6, 1), "Mb"),count))) +
  scale_x_continuous(labels = function(x) ifelse(x >= 1000000, paste0(x/1000000), x)) +
  theme_classic()+
  scale_fill_manual(values = c("#c9c9cb","#939eb2"))+
  labs(x="Varinat length(Mb)")+
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"),
        legend.position = "none")

cowplot::plot_grid(p2,p1)

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



# Fig.1g: compare SV between GTOP LRS and external and internal datasets ---------------

plot <- fread("Fig.1g.txt")
plot$svtype<-factor(plot$svtype,levels = c("INS","DEL"))
ggplot( plot,aes(x=svtype,y = value,fill = category)) +
  geom_col( width = 0.7) +
  facet_wrap(~AF_group, nrow = 1) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c( "Reported" = "#bfbebe", "Novel" = "#66869b" )) +
  labs(x = NULL, y = "Proportion",fill = NULL ) +
  theme_classic() +
  theme(axis.line = element_line(color = "black"),
        legend.position = "bottom",
        axis.text = element_text( color = "black",size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12,face = "bold"))


# Fig.1h: compare total TR count in LRS and SRS ------------------------------------------------

col_overlap  <- "#8b9dc5"
col_lrs_spec <- "#933628"
col_srs_spec <- "#247e87"
dat <- fread("Fig1h.txt")%>%
  mutate( group = factor(group, levels = c("Specific","Overlap")),
          TR_type = factor(TR_type, levels = c("2","3","4","5","6","VNTR")))

fill_map <- c("Overlap.LRS"  = col_overlap,"Overlap.SRS"  = col_overlap,
              "Specific.LRS" = col_lrs_spec,"Specific.SRS" = col_srs_spec)
ggplot() +
  geom_col(data = dat %>% dplyr::filter(tec == "LRS"),
           aes(y = count, x = TR_type, fill = interaction(group, tec)),
           position = "stack" ) +
  geom_col(data = dat %>% dplyr::filter(tec == "SRS"),
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
  geom_col(data = dat %>% dplyr::filter(tec == "LRS",TR_type=="VNTR"),
           aes(y = count, x = TR_type, fill = interaction(group, tec)),
           position = "stack" ) +
  geom_col(data = dat %>% dplyr::filter(tec == "SRS",TR_type=="VNTR"),
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






