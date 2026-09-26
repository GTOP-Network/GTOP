#===============================================#
# LRS read QC #
# Supp Figure 4   #
#===============================================#
library(ggplot2)
library(ggpubr)


setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig4/input")


# Supp.Figure.4a all-read versus aligned-read counts per sample --------------


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
data <- read.delim("Fig_S4a.txt", check.names = FALSE)
data$tissue_color <- first_existing_col(data, c("tissue_color"), "#7f8080")
ggplot(data, aes(total_reads_million, mapped_reads_million)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.6, color = "lightgray") +
  geom_point(aes(fill = tissue_color), shape = 21, size = 2.6, alpha = 1, color = "black", stroke = 0.3) +
  scale_fill_identity() +
  labs(x = "Number of all reads (million)", y = "Number of aligned reads (million)") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Supp.Figure.4b read Phred quality-score density from sampled per-read values-------------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_blank())
}

first_existing_col <- function(data, cols, default = NA_character_) {
  for (col in cols) if (col %in% names(data)) return(data[[col]])
  rep(default, nrow(data))
}

data <- read.delim("Fig_S4b.txt", check.names = FALSE)
data$tissue_color <- first_existing_col(data, "tissue_color", "#7f8080")
data <- data[order(data$line_order, data$sample, data$qscore), ]
data$sample <- factor(data$sample, levels = unique(data$sample))
y_breaks <- seq(0, ceiling(max(data$density, na.rm = TRUE) / 0.2) * 0.2, by = 0.2)

ggplot(data, aes(qscore, density, group = sample)) +
  geom_line(aes(color = tissue_color), alpha = 0.4, linewidth = 0.35) +
  scale_color_identity() +
  scale_y_continuous(breaks = y_breaks) +
  labs(x = "Phred quality score", y = "Density") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Supp.Figure.4c sample averaged aligned-read length across tissues--------------------------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_blank())
}

data <- read.delim("Fig_S4c.txt", check.names = FALSE)
data$tissue_label <- factor(data$tissue_label, levels = sort(unique(data$tissue_label)))
fill_values <- setNames(data$tissue_color[match(levels(data$tissue_label), data$tissue_label)],
                        levels(data$tissue_label))

ggplot(data, aes(tissue_label, mapped_mean_length, fill = tissue_label)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.35, show.legend = FALSE) +
  geom_jitter(width = 0.22, size = 0.8, alpha = 0.55, color = "black") +
  scale_fill_manual(values = fill_values) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = NULL, y = "Averaged aligned read length") +
  theme_gtop() +
  theme(legend.position = "none")
