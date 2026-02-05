#==================================#
# long & short RNA #
# Extend Data Figure 2 #
#==================================#
library(data.table)

setwd("/media/london_A/mengxin/GTOP_code/extend/extend_2/input")

# Extend.Data Fig.2a  ------------------------------------------------------
library(data.table)
library(ggplot2)

dat<-fread("ext2a.LR_transcript_novel_stat_Gene.txt")
dat_long <- melt(
  dat,
  id.vars = "index",
  variable.name = "Type",
  value.name = "Count"
)
dat_long$Type <- factor(
  dat_long$Type,
  levels = c( "Novel","Annotated", "GENCODE v47")
)
dat_long$index <- factor(
  dat_long$index,
  levels = c("GTOP", "GENCODE", "GTOP + GENCODE")
)
type_colors <- c(
  "GENCODE v47"        = "#7d8bad",
  "Annotated"  = "#c3968f",
  "Novel"          = "#9d3929"
)
dat_long$Count<-dat_long$Count/1000
ggplot(dat_long, aes(x = index, y = Count, fill = Type)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = type_colors) +
  labs(
    x = NULL,
    y = "Number of genes (x10³)",
    fill = "Category"
  ) +
  theme_classic(base_size = 13)

# Extend Data Fig.2b  ------------------------------------------------------
library(data.table)
library(UpSetR)

dat <- fread("ext2b.LR_transcript_compare_database_upset.txt", sep = "\t")
setDT(dat)

cols <- c("CHESS3", "GTEx long-read", "GTOP")
dat[, (cols) := lapply(.SD, as.logical), .SDcols = cols]

count_col <- setdiff(names(dat), cols)
dat[, (count_col) := as.integer(get(count_col))]

dat_expand <- dat[rep(seq_len(.N), get(count_col))]
dat_expand[, (cols) := lapply(.SD, as.integer), .SDcols = cols]
dat_expand[, (count_col) := NULL]
upset(
  dat_expand,
  sets = cols,
  order.by = "freq",
  decreasing = TRUE,
  sets.bar.color = "grey60",
  main.bar.color = "grey60",
  point.size = 3.5,
  line.size = 1.2,
  show.numbers = FALSE,
  text.scale = c(1.6, 1.6, 1.2, 1.2, 1.4, 1.2)
)

# Extend Data Fig.2c  ------------------------------------------------------
df <- fread("ext2c.LR_transcript_compare_database_bar.txt", sep = "\t")
df_long <- melt(df, id.vars = "index", variable.name = "Type", value.name = "Count")
df_long$Type <- factor(df_long$Type, levels = rev(c("Annotated", "Novel")))
type_colors <- c(
  "Annotated" = "#c3968f",
  "Novel" = "#9d3929"
)
ggplot(df_long, aes(x = index, y = Count, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = type_colors) +
  theme_classic() +
  labs(
    x = "",
    y = "Number of novel isoforms",
    fill = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

# Extend.Data.Fig.2d  ------------------------------------------------------

dat <- fread("ext2d.LR_transcript_tissue_contribute.txt", sep = "\t")
dat$tissue <- factor(dat$tissue, levels = dat[order(-number_transcript)]$tissue)
tissue_colors <- setNames(dat$color, dat$tissue_matched)
tissue_colors <- tissue_colors[!duplicated(names(tissue_colors))]

ggplot(dat, aes(x = tissue, y = number_transcript, fill = tissue_matched)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = tissue_colors) +
  labs(
    x = "",
    y = "Number of novel isoforms",
    fill = ""
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  guides(fill = guide_legend(ncol = 1))


# Extend.Data.Fig.2e  ------------------------------------------------------

dat <- fread("ext2e.LR_transcript_tissue_group_contribute.txt", sep = "\t")
dat$tissue_group <- factor(dat$tissue_group, levels = dat$tissue_group)

tissue_colors <- setNames(dat$color, dat$tissue_group)

ggplot(dat, aes(x = tissue_group, y = number_transcript, fill = tissue_group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = tissue_colors) +
  labs(
    x = "",
    y = "Number of novel isoforms",
    fill = ""
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
  ) +
  guides(fill = guide_legend(ncol = 1))


# Extend.Data.Fig.2fg -----------------------------------------------------

load("tpm_matrix.RData")
load("LRS_RNA_metadata.RData")

# novel&annotated -------------------------------------------------------

novel<-tpm_matrix[grepl("GTOP",rownames(tpm_matrix)),]
annotated<-tpm_matrix[grepl("ENS",rownames(tpm_matrix)),]

# remove batch (使用ComBat) -------------------------------------------------
library(sva)
metadata<-metadata[metadata$sample %in% colnames(tpm_matrix),]
identical(metadata$sample,colnames(novel))

batch <- metadata$batch  
mod <- model.matrix(~Tissue, data = metadata)

#ComBat
novel <- ComBat(
  dat = novel,
  batch = batch,
  mod = mod,  
  par.prior = TRUE,
  prior.plots = FALSE
)

annotated <- ComBat(
  dat = annotated,
  batch = batch,
  mod = mod,  
  par.prior = TRUE,
  prior.plots = FALSE
)


library(pheatmap)

cor_mat_novel <- cor(novel, method = "spearman")  
cor_mat_annotated <- cor(annotated, method = "spearman")  

identical(colnames(cor_mat_novel),colnames(cor_mat_annotated))


# # Convert correlation to distance --------------------------------------------------------------
dist_novel <- as.dist(1 - cor_mat_novel)
dist_annotated <- as.dist(1 - cor_mat_annotated)

anno_row<-metadata[,c(5,7)]
rownames(anno_row)<-metadata$sample
identical(rownames(anno_row),colnames(cor_mat_novel))
head(anno_row)

tissue_color<-unique(anno_row[,1:2])

tissue_colors <-list(Tissue = setNames(paste0("#",tissue_color$Tissue_Color_Code),tissue_color$Tissue))
ann_colors <- list(
  Tissue = tissue_colors
)
anno_row<-anno_row[,1,drop=F]

pheatmap(
  cor_mat_novel,
  clustering_distance_rows = dist_novel,
  clustering_distance_cols = dist_novel,
  annotation_colors = tissue_colors,
  clustering_method = "average",
  annotation_row = anno_row,
  show_rownames = FALSE,
  show_colnames = FALSE,
  display_numbers = FALSE,
  color = colorRampPalette(c("#305089", "white", "#be4e2f"))(100),
  main = "novel"
)

pheatmap(
  cor_mat_annotated,
  annotation_colors = tissue_colors,
  clustering_distance_rows = dist_annotated,
  clustering_distance_cols = dist_annotated,
  clustering_method = "average",
  annotation_row = anno_row,
  show_rownames = FALSE,
  show_colnames = FALSE,
  display_numbers = FALSE,
  color = colorRampPalette(c("#305089", "white", "#be4e2f"))(100),
  main = "annotated"
)

# Extend Data Fig.2h  ------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(gridExtra)
library(grid)

df <- fread("ext2h.SR_transcript_tissue_express_count.txt", sep = "\t")
df <- df %>% arrange(transcript_enhanced)
tissue_order <- df %>% 
  group_by(tissue) %>% 
  summarise(sum_val = sum(transcript_enhanced, na.rm = TRUE)) %>% 
  arrange(desc(sum_val)) %>% 
  pull(tissue)
df$tissue <- factor(df$tissue, levels = tissue_order)

df_transcript <- df %>% 
  select(tissue, transcript_gencode, transcript_enhanced, color) %>% 
  pivot_longer(cols = -c(tissue, color), 
               names_to = "type", 
               values_to = "value") %>% 
  mutate(group = gsub("transcript_", "", type))

df_gene <- df %>% 
  select(tissue, gene_gencode, gene_enhanced, color) %>% 
  pivot_longer(cols = -c(tissue, color), 
               names_to = "type", 
               values_to = "value") %>% 
  mutate(group = gsub("gene_", "", type))

y_offset_transcript <- -0.05 * max(df_transcript$value, na.rm = TRUE)
y_offset_gene <- -0.05 * max(df_gene$value, na.rm = TRUE)

p1 <- ggplot(df_transcript, aes(x = tissue, y = value, group = group)) +
  geom_area(aes(fill = group), alpha = 0.3, position = "identity") +
  geom_line(aes(color = group), linewidth = 1) +
  geom_point(data = df, aes(x = tissue, y = y_offset_transcript, color = I(color)), 
             size = 3, shape = 19, inherit.aes = FALSE) +
  scale_fill_manual(values = c("gencode" = "#1f77b4", "enhanced" = "#d62728")) +
  scale_color_manual(values = c("gencode" = "#1f77b4", "enhanced" = "#d62728")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 0, 10),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(y = "Number of transcripts")

p2 <- ggplot(df_gene, aes(x = tissue, y = value, group = group)) +
  geom_area(aes(fill = group), alpha = 0.3, position = "identity") +
  geom_line(aes(color = group), linewidth = 1) +
  geom_point(data = df, aes(x = tissue, y = y_offset_gene, color = I(color)), 
             size = 3, shape = 19, inherit.aes = FALSE) +
  scale_fill_manual(values = c("gencode" = "#1f77b4", "enhanced" = "#d62728")) +
  scale_color_manual(values = c("gencode" = "#1f77b4", "enhanced" = "#d62728")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 10, 10, 10)
  ) +
  labs(x = "Tissue", y = "Number of genes")

title <- textGrob("Gene & Transcript Values by Tissue", gp = gpar(fontsize = 14, fontface = "bold"))
final_plot <- grid.arrange(title, p1, p2, 
                           ncol = 1, 
                           heights = c(0.1, 0.45, 0.45))

print(final_plot)


# Extend Data Fig.2i  ------------------------------------------------------
library(data.table)
library(pheatmap)
df <- fread("ext2i.SR_novel_isoform_example_expression_heatmap.txt", sep = "\t")
expr_mat <- as.matrix(df[, -1, with = FALSE])
rownames(expr_mat) <- df$Transcript
expr_mat <- log2(expr_mat + 1)
my_colors <- colorRampPalette(c("#F2F2F2", "#B2182B"))(256)

pheatmap(
  expr_mat,
  color = my_colors,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 8,
  border_color = 'white',
  angle_col = 90,
  treeheight_row = 10,
  treeheight_col = 10
)




