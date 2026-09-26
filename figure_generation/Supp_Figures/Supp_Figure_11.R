#===============================================#
# LRS RNA transcript-assemble #
# Supp-Figure-11                         #
#===============================================#
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig11.R1.5-LRS-RNA-transcript-assemble/input")


# Supp.Fig.11b structural category composition by merge method----------------

structural_palette <- c(
  "FSM" = "#4E79A7",
  "ISM" = "#59A14F",
  "NIC" = "#F28E2B",
  "NNC" = "#E15759",
  "Fusion" = "#B07AA1",
  "Full-splice match" = "#4E79A7",
  "Incomplete-splice match" = "#59A14F",
  "Novel in catalog" = "#F28E2B",
  "Novel not in catalog" = "#E15759",
  "Genic" = "#76B7B2",
  "Antisense" = "#EDC948",
  "Intergenic" = "#9C755F"
)
scale_count_k <- function(values) {
  max_value <- max(values, na.rm = TRUE)
  if (max_value > 100000) {
    list(scale = 1000, suffix = " (x10^3)")
  } else {
    list(scale = 1, suffix = "")
  }
}
data <- read.delim("Fig_S11b.txt", check.names = FALSE)
data$method <- factor(data$method, levels = unique(data$method[order(data$method_order)]))
data$category <- factor(data$category, levels = unique(data$category[order(data$category_order)]))
count_scale <- scale_count_k(data$total)
ggplot(data, aes(method, count / count_scale$scale, fill = category)) +
  geom_col(width = 0.9, color = "white", linewidth = 0.3, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = structural_palette, drop = FALSE) +
  labs(x = "Method", y = paste0("Number of transcripts", count_scale$suffix), fill = "Structural category") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(1.05, 0.65),
    legend.justification = c(0, 0.5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    plot.margin = margin(5.5, 65, 30, 5.5)
  )

# Supp.Fig.11c UpSet plot of transcript intersections supported by at least three methods-------

library(patchwork)
structural_palette <- c(
  "FSM" = "#4E79A7",
  "ISM" = "#59A14F",
  "NIC" = "#F28E2B",
  "NNC" = "#E15759",
  "Fusion" = "#B07AA1",
  "Full-splice match" = "#4E79A7",
  "Incomplete-splice match" = "#59A14F",
  "Novel in catalog" = "#F28E2B",
  "Novel not in catalog" = "#E15759",
  "Genic" = "#76B7B2",
  "Antisense" = "#EDC948",
  "Intergenic" = "#9C755F"
)
scale_count_k <- function(values) {
  max_value <- max(values, na.rm = TRUE)
  if (max_value > 100000) {
    list(scale = 1000, suffix = " (x10^3)")
  } else {
    list(scale = 1, suffix = "")
  }
}
intersections <- read.delim("Fig_S11c_intersections.txt", check.names = FALSE)
matrix <- read.delim("Fig_S11c_matrix.txt", check.names = FALSE)
totals <- read.delim("Fig_S11c_totals.txt", check.names = FALSE)

intersection_levels <- unique(intersections$intersection[order(intersections$intersection_order)])
method_levels <- unique(matrix$method[order(matrix$method_order)])
category_levels <- unique(intersections$category[order(intersections$category_order)])

intersections$intersection <- factor(intersections$intersection, levels = intersection_levels)
intersections$category <- factor(intersections$category, levels = category_levels)
matrix$intersection <- factor(matrix$intersection, levels = intersection_levels)
matrix$method <- factor(matrix$method, levels = rev(method_levels))
matrix$present <- tolower(as.character(matrix$present)) %in% c("true", "1")
totals$method <- factor(totals$method, levels = rev(method_levels))
totals$category <- factor(totals$category, levels = category_levels)

inter_sum <- aggregate(count ~ intersection, intersections, sum)
total_sum <- aggregate(count ~ method, totals, sum)
inter_scale <- scale_count_k(inter_sum$count)
total_scale <- scale_count_k(total_sum$count)

p_intersections <- ggplot(intersections, aes(intersection, count / inter_scale$scale, fill = category)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.25, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = structural_palette, drop = FALSE) +
  labs(x = NULL, y = paste0("Intersection size", inter_scale$suffix), fill = "Structural category") +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(5.5, 5.5, 0, 5.5)
  )

p_matrix <- ggplot(matrix, aes(intersection, method)) +
  geom_point(color = "#DDDDDD", size = 1.7) +
  geom_line(
    data = subset(matrix, present),
    aes(group = intersection),
    color = "black",
    linewidth = 0.35
  ) +
  geom_point(data = subset(matrix, present), color = "black", size = 1.7) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0, 5.5, 5.5, 5.5)
  )

p_totals <- ggplot(totals, aes(method, count / total_scale$scale, fill = category)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.25, position = position_stack(reverse = TRUE)) +
  coord_flip() +
  scale_y_reverse() +
  scale_fill_manual(values = structural_palette, drop = FALSE, guide = "none") +
  labs(x = NULL, y = paste0("Set size", total_scale$suffix)) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 0, 5.5, 5.5)
  )

((plot_spacer() / p_totals + plot_layout(heights = c(3, 1.2))) |
    (p_intersections / p_matrix + plot_layout(heights = c(3, 1.2)))) +
  plot_layout(widths = c(0.65, 4.8))
