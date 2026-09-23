#==================================#
# long & short RNA #
# Extend Data Figure 2 #
#==================================#
library(data.table)
library(ggplot2)
setwd("/media/london_A/mengxin/GTOP_code/extend/extend_2/input")

# Extended.Fig.2a  ------------------------------------------------------
gtop_ref_colors <- c(
  "Annotated" = "#cf928f",
  "Novel" = "#ad3b2b",
  "GENCODE v47" = "#9aaac2",
  "GTOP + GENCODE" = "#6f86ad"
)
data <- read.delim("Ext_2a.txt", check.names = FALSE)
data$component <- ifelse(data$component == "Merged reference", "GTOP + GENCODE", data$component)
data$group <- factor(data$group, levels = c("GTOP", "GENCODE", "GTOP + GENCODE"))
data$component <- factor(data$component, levels = c("Annotated", "Novel", "GENCODE v47", "GTOP + GENCODE"))
ggplot(data, aes(group, count / 1000, fill = component)) +
  geom_col(width = 0.65, color = NA, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = gtop_ref_colors, breaks = c("Annotated", "Novel"), drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 80, by = 10),expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Number of genes (x10^3)") +
  theme_gtop()
# Extended.Fig.2b  ------------------------------------------------------

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
theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
data <- read.delim("Ext_2b.txt", check.names = FALSE)

data$tissue_group <- factor(data$tissue_group, levels = c("Novel tissue", "GTEx-matched tissue"))
data$tissue_label <- factor(data$tissue_label, levels = data$tissue_label[order(data$all_count, decreasing = TRUE)])
group_colors <- c("Novel tissue" = "#47689e", "GTEx-matched tissue" = "#b5b5b6")
group_labels <- c("Novel tissue" = "Uncharacterized tissue", "GTEx-matched tissue" = "GTEx-matched tissue")
ggplot(data, aes(tissue_label, all_count, fill = tissue_group)) +
  geom_col(width = 0.8, color = NA) +
  scale_fill_manual(values = group_colors, labels = group_labels, drop = FALSE) +
  scale_y_continuous(limits = c(0, 35000),breaks = seq(0,35000,by=5000),expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Number of novel transcripts", fill = NULL) +
  theme_gtop() +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1))
# Extended.Fig.2c  ------------------------------------------------------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
data <- read.delim( "Ext_2c.txt", check.names = FALSE)
data$component <- factor(data$component, levels = c("Exclusive transcripts", "Shared transcripts"))
ggplot(data, aes(group, count, fill = component)) +
  geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("Exclusive transcripts" = "#3f6f9f", "Shared transcripts" = "#a7bfd8")) +
  labs(x = NULL, y = "Number of novel transcripts") +
  theme_gtop()

# Extended.Fig.2d  ------------------------------------------------------
library(data.table)
library(UpSetR)
library(patchwork)

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
as_bool <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}
data <- read.delim("Ext_2d.txt", check.names = FALSE)
data$gtop <- as_bool(data$gtop)
data$chess <- as_bool(data$chess)
data$combination <- ifelse(data$gtop & data$chess, "GTOP|CHESS3", ifelse(data$gtop, "GTOP", "CHESS3"))
data <- data[order(data$count, decreasing = TRUE), ]
data$combination <- factor(data$combination, levels = data$combination)
matrix <- data.frame(
  combination = rep(data$combination, each = 2),
  method = rep(c("GTOP", "CHESS3"), times = nrow(data)),
  present = as.vector(t(as.matrix(data[, c("gtop", "chess")]))),
  stringsAsFactors = FALSE
)
matrix$method <- factor(matrix$method, levels = c("GTOP", "CHESS3"))
totals <- data.frame(
  method = c("GTOP", "CHESS3"),
  count = c(sum(data$count[data$gtop]), sum(data$count[data$chess]))
)
totals$method <- factor(totals$method, levels = c("GTOP", "CHESS3"))
matrix$method <- factor(matrix$method, levels = c("GTOP", "CHESS3"))
p_top <- ggplot(data, aes(combination, count)) +
  geom_col(fill = "grey60", width = 0.72, color = "black", linewidth = 0.2) +
  labs(x = NULL, y = "Intersection size") +
  theme_gtop() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p_matrix <- ggplot(matrix, aes(combination, method)) +
  geom_point(color = "#DDDDDD", size = 2.2) +
  geom_line(data = subset(matrix, present), aes(group = combination), color = "black", linewidth = 0.4) +
  geom_point(data = subset(matrix, present), color = "black", size = 2.2) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line = element_blank(), axis.ticks.y = element_blank())
p_left <- ggplot(totals, aes(method, count)) +
  geom_col(fill = "grey60", width = 0.65, color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_y_reverse() +
  labs(x = NULL, y = "Set size") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
((plot_spacer() / p_left + plot_layout(heights = c(3, 1))) |
    (p_top / p_matrix + plot_layout(heights = c(3, 1)))) +
  plot_layout(widths = c(0.75, 4.2))
# Extended.Fig.2e  ------------------------------------------------------
theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
data <- read.delim("Ext_2e.txt", check.names = FALSE)
data <- data[, c("label", "Annotated", "Novel")]
long <- reshape(data, direction = "long", varying = c("Annotated", "Novel"), v.names = "count", timevar = "annotation", times = c("Annotated", "Novel"))
long$label <- factor(long$label, levels = unique(data$label))
long$annotation <- factor(long$annotation, levels = c("Annotated", "Novel"))
ggplot(long, aes(label, count, fill = annotation)) +
  geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("Annotated" = "#cf928f", "Novel" = "#ad3b2b")) +
  labs(x = NULL, y = "Number of transcripts") +
  theme_gtop()

# Extended.Fig.2f-g -----------------------------------------------------


load("LRS_RNA_metadata.RData")
load("cor_mat_novel.RData")
load("cor_mat_annotated.RData")

dist_novel <- as.dist(1 - cor_mat_novel)
dist_annotated <- as.dist(1 - cor_mat_annotated)


anno_row<-metadata[,c(5,7)]
rownames(anno_row)<-metadata$sample
anno_row<-anno_row[colnames(cor_mat_novel),]
identical(rownames(anno_row),colnames(cor_mat_novel))
head(anno_row)

tissue_color<-unique(anno_row[,1:2])

tissue_colors <-list(Tissue = setNames(paste0("#",tissue_color$Tissue_Color_Code),tissue_color$Tissue))
ann_colors <- list(
  Tissue = tissue_colors
)
anno_row<-anno_row[,1,drop=F]

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
# Extended.Fig.2h  ------------------------------------------------------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
first_existing_col <- function(data, cols, default = NA_character_) {
  for (col in cols) {
    if (col %in% names(data)) return(data[[col]])
  }
  rep(default, nrow(data))
}
data <- read.delim("Ext_2h.txt", check.names = FALSE)
data <- data[order(data$transcript_enhanced, decreasing = TRUE), ]
data$rank <- seq_len(nrow(data))
data[, c("gene_gencode", "gene_enhanced", "transcript_gencode", "transcript_enhanced")] <-
  data[, c("gene_gencode", "gene_enhanced", "transcript_gencode", "transcript_enhanced")] / 1000
scale_factor <- max(data$transcript_enhanced, na.rm = TRUE) / max(data$gene_enhanced, na.rm = TRUE)
data$gene_gencode_mirror <- -data$gene_gencode * scale_factor
data$gene_enhanced_mirror <- -data$gene_enhanced * scale_factor
data$tissue_color <- first_existing_col(data, c("tissue_color"), "#7f8080")
color_gencode <- "#1f77b4"
color_enhanced <- "#d62728"

gene_max_true <- ceiling(max(data$gene_enhanced, na.rm = TRUE) / 10) * 10   
gene_breaks_true <- seq(0, gene_max_true, by = 10)          
gene_breaks_scaled <- -gene_breaks_true * scale_factor      

trans_max <- ceiling(max(data$transcript_enhanced, na.rm = TRUE) / 50) * 50
trans_breaks <- seq(0, trans_max, by = 50)                   

y_breaks <- c(rev(gene_breaks_scaled), trans_breaks[-1])     
y_labels <- c(rev(as.character(gene_breaks_true)), as.character(trans_breaks[-1]))

ggplot(data, aes(rank)) +
  geom_ribbon(aes(ymin = 0, ymax = transcript_gencode), fill = "#d0e1f2", alpha = 0.55) +
  geom_ribbon(aes(ymin = transcript_gencode, ymax = transcript_enhanced), fill = "#f4b6b6", alpha = 0.45) +
  geom_line(aes(y = transcript_gencode, color = "GENCODE v47"), linewidth = 0.75) +
  geom_line(aes(y = transcript_enhanced, color = "GTOP + GENCODE v47"), linewidth = 0.75) +
  geom_ribbon(aes(ymin = gene_gencode_mirror, ymax = 0), fill = "#d0e1f2", alpha = 0.55) +
  geom_ribbon(aes(ymin = gene_enhanced_mirror, ymax = gene_gencode_mirror), fill = "#f4b6b6", alpha = 0.45) +
  geom_line(aes(y = gene_gencode_mirror, color = "GENCODE v47"), linewidth = 0.75) +
  geom_line(aes(y = gene_enhanced_mirror, color = "GTOP + GENCODE v47"), linewidth = 0.75) +
  geom_point(aes(y = min(gene_enhanced_mirror, na.rm = TRUE) * 1.05, fill = tissue_color), shape = 21, color = "black", size = 1.8, stroke = 0.25) +
  annotate("text", x = nrow(data) * 0.4, y = max(data$transcript_enhanced, na.rm = TRUE) * 0.96, label = "Transcript", size = 3.5) +
  annotate("text", x = nrow(data) * 0.4, y = min(data$gene_enhanced_mirror, na.rm = TRUE) * 0.75, label = "Gene", size = 3.5) +
  scale_color_manual(
    values = c("GENCODE v47" = color_gencode, "GTOP + GENCODE v47" = color_enhanced),
    na.value = "#7f8080",
    breaks = c("GENCODE v47", "GTOP + GENCODE v47")
  ) +
  scale_fill_identity() +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  labs(x = NULL, y = "Number of expressed features (x10^3)", color = NULL) +
  theme_gtop() +
  theme(axis.text.x = element_blank(), legend.position = c(0.98, 0.98), legend.justification = c(1, 1))

# Extended.Fig.2i  ------------------------------------------------------
theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
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
isoform_structure_colors <- c("Novel" = "#C44E52", "Annotated" = "#4C72B0")
expr <- read.delim("Ext_2i_expr.txt", check.names = FALSE)
expr$log2_tpm <- as.numeric(expr$log2_tpm)
transcript_order <- unique(expr$transcript_id)
tissue_labels <- first_existing_col(expr, c("tissue_label"), NA_character_)
tissue_labels[is.na(tissue_labels)] <- pretty_label(expr$tissue[is.na(tissue_labels)])
tissue_map <- unique(data.frame(tissue = expr$tissue, tissue_label = tissue_labels, stringsAsFactors = FALSE))
expr$tissue_label <- factor(tissue_labels, levels = tissue_map$tissue_label)
expr$transcript_id <- factor(expr$transcript_id, levels = rev(transcript_order))
p1 <- ggplot(expr, aes(tissue_label, transcript_id, fill = log2_tpm)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "#f2f2f2", high = "#b2182b", na.value = "white") +
  labs(x = NULL, y = NULL, fill = "log2(TPM+1)") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

structure <- read.delim("Ext_2i_structure.txt", check.names = FALSE)
if (!("chrom" %in% names(structure))) structure$chrom <- "chrom"
if (!("strand" %in% names(structure))) structure$strand <- "."
transcript_ids <- c("GTOPT000378535", "ENST00000228841.15", "ENST00000548438.1")
transcript_ids <- transcript_ids[transcript_ids %in% unique(structure$transcript_id)]
if (length(transcript_ids) == 0) transcript_ids <- unique(structure$transcript_id)
plot_tx <- rev(transcript_ids)
exons <- structure[structure$feature == "exon" & structure$transcript_id %in% transcript_ids, ]
exons <- exons[order(exons$chrom, exons$strand, exons$start, exons$end), ]
merged <- data.frame()
for (key in unique(paste(exons$chrom, exons$strand, sep = "\r"))) {
  sub <- exons[paste(exons$chrom, exons$strand, sep = "\r") == key, ]
  if (nrow(sub) == 0) next
  cur_start <- sub$start[1]
  cur_end <- sub$end[1]
  rows <- list()
  for (i in seq_len(nrow(sub))[-1]) {
    if (sub$start[i] <= cur_end) {
      cur_end <- max(cur_end, sub$end[i])
    } else {
      rows[[length(rows) + 1]] <- data.frame(start = cur_start, end = cur_end)
      cur_start <- sub$start[i]
      cur_end <- sub$end[i]
    }
  }
  rows[[length(rows) + 1]] <- data.frame(start = cur_start, end = cur_end)
  merged <- rbind(merged, do.call(rbind, rows))
}
merged <- merged[order(merged$start, merged$end), ]
intron_spacing <- 100
merged$meta_start <- c(0, head(cumsum((merged$end - merged$start) + intron_spacing), -1))
map_to_meta <- function(x) {
  hit <- which(merged$start <= x & x <= merged$end)
  if (length(hit) > 0) return(merged$meta_start[hit[1]] + (x - merged$start[hit[1]]))
  before <- which(x < merged$start)
  if (length(before) > 0) return(merged$meta_start[before[1]] - intron_spacing)
  max(merged$meta_start + (merged$end - merged$start), na.rm = TRUE)
}
structure <- structure[structure$transcript_id %in% transcript_ids, ]
structure$x_start <- vapply(structure$start, map_to_meta, numeric(1))
structure$x_end <- vapply(structure$end, map_to_meta, numeric(1))
introns <- do.call(rbind, lapply(plot_tx, function(tx) {
  tx_exons <- structure[structure$transcript_id == tx & structure$feature == "exon", ]
  tx_exons <- tx_exons[order(tx_exons$x_start), ]
  if (nrow(tx_exons) < 2) return(NULL)
  data.frame(
    transcript_id = tx,
    x_start = head(tx_exons$x_end, -1),
    x_end = tail(tx_exons$x_start, -1),
    stringsAsFactors = FALSE
  )
}))
structure$transcript_id <- factor(structure$transcript_id, levels = plot_tx)
if (!is.null(introns)) introns$transcript_id <- factor(introns$transcript_id, levels = plot_tx)
p2 <- ggplot() +
  {if (!is.null(introns)) geom_segment(data = introns, aes(x = x_start, xend = x_end, y = transcript_id, yend = transcript_id), color = "black", linewidth = 0.3)} +
  geom_segment(data = subset(structure, feature == "exon"), aes(x = x_start, xend = x_end, y = transcript_id, yend = transcript_id, color = annotation), linewidth = 4, lineend = "butt") +
  geom_segment(data = subset(structure, feature == "CDS"), aes(x = x_start, xend = x_end, y = transcript_id, yend = transcript_id, color = annotation), linewidth = 6, lineend = "butt") +
  scale_color_manual(values = isoform_structure_colors, drop = FALSE) +
  labs(x = "Transcript relative position (intron spacing = 100)", y = NULL) +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
list(expression = p1, structure = p2)
# Extended.Fig.2j  ------------------------------------------------------
# generated using Adobe Illustrator.



