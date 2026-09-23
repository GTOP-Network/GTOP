#==================================#
# long & short RNA #
# Figure-2 #
#==================================#
library(data.table)
library(ggplot2)
setwd("/media/london_A/mengxin/GTOP_code/fig-2/input/")
# Fig.2a LRS_RNA_MDS + hc------------------------------------------------------

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

mds_df_corrected<-fread("LR_Sample_MDS_coor.txt")
mds_coor<-fread("LR_Sample_MDS_var.txt")

prop1_corr<-mds_coor[1,1]
prop2_corr<-mds_coor[1,2]

gg_mds_corr <- ggplot(mds_df_corrected, aes(x = MDS1, y = MDS2, color = Tissue)) +
  geom_point(size = 2, alpha = 1) +
  scale_color_manual(values = tissue_colors) +
  unified_theme +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme_classic()+
  labs(
    x = paste0("Coordinate 1 (", prop1_corr, ")"),
    y = paste0("Coordinate 2 (", prop2_corr, ")"),
    color = "Tissue"
  )
print(gg_mds_corr)

library(dendextend)
load("LRS_expr_by_tissue.RData")

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

# Fig.2b LR_transcript_novel_stat_Transcript --------------------------------------------
library(data.table)
library(ggplot2)

theme_gtop <- function(base_size = 10) {
    theme_classic(base_size = base_size) +
      theme(
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_blank()
      )
  }
gtop_ref_colors <- c(
    "Annotated" = "#cf928f",
    "Novel" = "#ad3b2b",
    "GENCODE v47" = "#9aaac2",
    "GTOP + GENCODE" = "#6f86ad"
  )
data <- read.delim("Fig_2b.txt", check.names = FALSE)
data$component <- ifelse(data$component == "Merged reference", "GTOP + GENCODE", data$component)
data$group <- factor(data$group, levels = c("GTOP", "GENCODE", "GTOP + GENCODE"))
data$component <- factor(data$component, levels = c("Annotated", "Novel", "GENCODE v47", "GTOP + GENCODE"))
ggplot(data, aes(group, count / 1000, fill = component)) +
    geom_col(width = 0.65, color = NA, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = gtop_ref_colors, breaks = c("Annotated", "Novel"), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Number of transcripts (x10^3)") +
    theme_gtop()

# Fig.2c LR_SQANTI3_annot_coding_status -----------------------------------
library(patchwork)
  
theme_gtop <- function(base_size = 10) {
    theme_classic(base_size = base_size) +
      theme(
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA)
      )
  }
  
data <- read.delim("Fig_2c.txt", check.names = FALSE)
long <- reshape(
    data,
    direction = "long",
    varying = c("Protein-coding", "Noncoding", "NMD-sensitive"),
    v.names = "count",
    timevar = "coding_status",
    times = c("Protein-coding", "Noncoding", "NMD-sensitive")
  )
long$structural_category <- factor(long$structural_category, levels = data$structural_category)
long$coding_status <- factor(long$coding_status, levels = c("Protein-coding", "Noncoding", "NMD-sensitive"))
  
total_counts <- rowSums(data[, setdiff(names(data), "structural_category"), drop = FALSE])
y_cut <- 3.8
y_upper <- max(total_counts, na.rm = TRUE) * 1.15
  
coding_colors <- c("Protein-coding" = "#c0392b", "Noncoding" = "#f39c12", "NMD-sensitive" = "#16a085")
  

breaks_top    <- c(20, 40, 60,80)     
breaks_bottom <- c(0, 2)           

  
p_top <- ggplot(long, aes(structural_category, count, fill = coding_status)) +
    geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
    coord_cartesian(ylim = c(y_cut, y_upper), clip = "on") +
    scale_y_continuous(
      breaks = breaks_top,
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_fill_manual(values = coding_colors, drop = FALSE) +
    labs(x = NULL, y = NULL, fill = "Coding status") +
    theme_gtop() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = margin(5.5, 5.5, 2, 22)
    )
  
p_bottom <- ggplot(long, aes(structural_category, count, fill = coding_status)) +
    geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
    coord_cartesian(ylim = c(0, y_cut), clip = "on") +
    scale_y_continuous(
      breaks = breaks_bottom,
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_fill_manual(values = coding_colors, drop = FALSE, guide = "none") +
    labs(x = NULL, y = "Number of transcripts (x10^3)") +
    theme_gtop() +
    theme(plot.margin = margin(2, 5.5, 5.5, 22))
  
(p_top / p_bottom) +
  plot_layout(heights = c(2, 1), guides = "collect") &
  theme(legend.position = "right")
  
# Fig2.d LR_suppa_splicing_events_stat ------------------------------------

data<-fread("Fig_2d.txt")
theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
data$event <- factor(data$event, levels = unique(data$event))
data$component <- factor(data$component, levels = c("GENCODE v47", "Novel"))
ggplot(data, aes(event, count / 1000, fill = component)) +
  geom_col(width = 0.65, color = NA, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("GENCODE v47" = "#9aaac2", "Novel" = "#ad3b2b")) +
  labs(x = "Splicing events", y = "Number of splicing events (x10^3)") +
  theme_gtop()
# Fig2.e tissue breadth distribution for annotated and novel transcripts using the min-10-samples table.-----------------------------
library(tidyr)
data<-fread("Fig_2e.txt")

data$x_num <- suppressWarnings(as.numeric(data$x))
data$category <- factor(data$category, levels = c("Annotated", "Novel"))
ggplot(data, aes(x_num, proportion, color = category)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  scale_color_manual(values = c("Annotated" = "#cf928f", "Novel" = "#ad3b2b")) +
  labs(x = "Number of tissues", y = "Proportion of transcripts") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


# Fig2.f tissue-level peptide support for novel transcripts. ------------------------------------------
library(dplyr)
pretty_label <- function(x) {
  x <- gsub("_", " ", as.character(x))
  x <- tolower(x)
  ifelse(nchar(x) > 0, paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))), x)
}
first_existing_col <- function(data, cols, default = NA_character_) {
  for (col in cols) {
    if (col %in% names(data)) return(data[[col]])
  }
  rep(default, nrow(data))
}
data <- read.delim("Fig_2f.txt", check.names = FALSE)
data <- data[data$validate_status %in% c("Unique peptide", "Shared peptide"), ]
data$tissue_label <- first_existing_col(data, c("tissue_label"), NA_character_)
data$tissue_label[is.na(data$tissue_label)] <- pretty_label(data$tissue[is.na(data$tissue_label)])
data$tissue_color <- first_existing_col(data, c("tissue_color"), "#7f8080")
data$tissue_label <- factor(data$tissue_label, levels = sort(unique(as.character(data$tissue_label))))
data$validate_status <- factor(data$validate_status, levels = c("Unique peptide", "Shared peptide"))
tissue_points <- data[!duplicated(data$tissue_label), c("tissue_label", "tissue_color")]
ggplot(data, aes(tissue_label, ratio, fill = validate_status)) +
  geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
  geom_point(
    data = tissue_points,
    aes(x = tissue_label, y = -0.03, color = tissue_color),
    inherit.aes = FALSE,
    size = 2
  ) +
  scale_fill_manual(values = c("Unique peptide" = "#5b8fc7", "Shared peptide" = "#a3c1e1")) +
  scale_color_identity() +
  scale_y_continuous(breaks = seq(0,1,by=0.2))+
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(x = "Tissues", y = "Proportion of isoform with peptides") +
  theme_gtop() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Fig2.g WGCNA ------------------------------------------------------------

library(tidyverse)  
library(limma)  
library(RColorBrewer) 
library(vroom)
library(stringr)

load("heart_go_use.RData")
load("merged_expr.RData")
load("sample_annot_full_1613.RData")

heart<-sample_annot_full[str_detect(sample_annot_full$Tissue,"Heart"),]

tissue <- heart$Tissue
names(tissue) <- heart$sample_id


merged_expr_norm <- t(scale(t(merged_expr), center = TRUE, scale = TRUE))
merged_expr_norm <- as.data.frame(merged_expr_norm)

my_colors <- c(colorRampPalette(c("#8089a9", "white"))(50),
               colorRampPalette(c("white", "#bb4633"))(50))
my_breaks <- c(seq(-1.5, 0, length.out = 51), seq(0.01, 1.5, length.out = 50))

library(pheatmap)
p<-pheatmap(merged_expr_norm,
            color = my_colors,
            breaks = my_breaks,
            cluster_rows = T,
            cluster_cols = T,
            border_color = NA,
            show_rownames = T,
            show_colnames = T
)

p

row_order <- p$tree_row$order
row_names_sorted <- rownames(merged_expr_norm)[row_order]

go_result_bar$log10p<- -log10(go_result_bar$p.adjust)
library(ggplot2)
go_result_bar<-go_result_bar[match(row_names_sorted,go_result_bar$Module),]
go_result_bar$Description <- factor(go_result_bar$Description, levels = rev(go_result_bar$Description))

ggplot(go_result_bar, aes(x = Description, y = log10p)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#56B4E9") +  
  coord_flip() +  
  labs(x = "GO Term", y = "-log10(p.adjust)", title = "GO Enrichment") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),          
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )


library(data.table)

load("Fig2g.p_novel.RData")
p_novel <- p_novel %>%
  group_by(Var1) %>%
  mutate(Proportion = Freq / sum(Freq)) %>%
  ungroup()
my_colors <- c("annotated" = "white", "novel" = "#CCB1CC")
p_novel$Var1<-factor(p_novel$Var1,levels = rev(row_names_sorted))

ggplot(p_novel, aes(x = Var1, y = Proportion, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = my_colors) +  
  labs(
    x = "Var1",
    y = "Proportion",
    fill = "Type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank()
  )+
  coord_flip() 


# Fig2.h ASE/ASTS gene number  ---------------------------------------------------------
library(ggpubr)
library(dplyr)
df_plot <- fread("Fig 2h.txt")
ggplot(df_plot, aes(x=number, y=reorder(class, number)))+
  geom_bar(stat = "identity", fill="#5c86af")+
  #geom_text(aes(x=number+2, label=number))+
  theme_pubr()+
  ylab("")

# Fig.2i the number of significant ase/asts gene for per tissue -----------


df_count.ase_tissue <- fread("Fig 2i.txt")

ggplot(df_count.ase_tissue,aes(x=reorder(Tissue, -Freq),y=Freq)) + geom_bar(stat = "identity",fill="#a7b2c5") + theme_pubr() + 
  xlab("")+
  ylab("# ASE(gene)")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))

df_plot <- fread("Fig 2i.asts.txt")
df_plotsum <- df_plot %>% 
  group_by(Tissue) %>% 
  dplyr::summarise(total=sum(value)) %>% 
  arrange(desc(total))


df_plot$Tissue <- factor(df_plot$Tissue,levels = df_plotsum$Tissue)

ggplot(df_plot,aes(x=Tissue,y=value,fill=variable)) + geom_bar(stat = "identity") + theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) + 
  scale_fill_manual(breaks = c("A","B","C"),values = c("#fee5d9","#fc9272","#de2d26"),
                    labels = c("no novel transcript", "with novel transcript", "with novel transcript aFC>1.5")) + 
  ylab("# ASTS(gene)")




