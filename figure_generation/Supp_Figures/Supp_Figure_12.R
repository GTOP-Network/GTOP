#===============================================#
# Transcript discovery, annotation, and filtering #
# Supp-Figure-12                         #
#===============================================#
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig12. Figure S9/input")


# Supp.Fig.12a alternative-splicing count QC distribution------------------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
python_status_colors <- c("passed" = "#4682B4", "filtered" = "#FFA07A")
data <- read.delim("Fig_S12a.txt", check.names = FALSE)
data$status <- factor(data$status, levels = c("passed", "filtered"))
ggplot(data, aes(x, count, fill = status)) +
  geom_col(width = 1, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = python_status_colors, drop = FALSE) +
  labs(x = "Number of adenines in 20bp downstream TTS", y = "Number of transcripts") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Supp.Fig.12b RT-switching filter summary------------------------------------
theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
python_status_colors <- c("passed" = "#4682B4", "filtered" = "#FFA07A")
data <- read.delim("Fig_S12b.txt", check.names = FALSE)
data$plot_status <- factor(data$plot_status, levels = c("passed", "filtered"))
ggplot(data, aes(x_label, count, fill = plot_status)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = python_status_colors, drop = FALSE) +
  labs(x = "RT-switching status", y = "Number of transcripts") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Supp.Fig.12c noncanonical splice-junction filter summary------------------------
theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
python_status_colors <- c("passed" = "#4682B4", "filtered" = "#FFA07A")
data <- read.delim("Fig_S12c.txt", check.names = FALSE)
data$plot_status <- factor(data$plot_status, levels = c("passed", "filtered"))
ggplot(data, aes(x_label, count, fill = plot_status)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = python_status_colors, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 500000, by = 100000)
  ) +
  labs(x = "Non-canonical splice junctions", y = "Number of transcripts") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Supp.Fig.12d FLNC read-count support distribution------------------------------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
python_status_colors <- c("passed" = "#4682B4", "filtered" = "#FFA07A")
data <- read.delim("Fig_S12d.txt", check.names = FALSE)
data$status <- factor(data$status, levels = c("passed", "filtered"))
ggplot(data, aes(x_center, y, fill = status)) +
  geom_col(aes(width = x_width), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = python_status_colors, drop = FALSE) +
  scale_y_continuous( breaks = seq(0, 50000, by = 10000))+
  labs(x = "Supported FLNC read count (log10(x + 1))", y = "Number of transcripts") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_vline(xintercept = log10(11), color = "red", linetype = "dashed", linewidth = 0.5)

# Supp.Fig.12e support-sample distribution for novel transcripts-----------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
data <- read.delim("Fig_S12e.txt", check.names = FALSE)
data$status <- "all"
ggplot(data, aes(x_center, y, fill = status)) +
  geom_col(aes(width = x_width), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("all" = "#7f8080"), drop = FALSE) +
  scale_y_continuous(breaks = seq(2.0, 5.5, by = 0.5))+
  coord_cartesian(ylim = c(2.0, 5.5),expand = c(0, 0) ) +
  labs(x = "Supported long-read RNA-seq samples for novel transcripts (FLNC >= 5)", y = "Number of transcripts (log10)") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Supp.Fig.12f -novel first-exon length distribution---------------------------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}

python_status_colors <- c(
  "passed" = "#4682B4",
  "filtered" = "#FFA07A"
)

data <- read.delim("Fig_S12f.txt", check.names = FALSE)

data$status <- factor(
  data$status,
  levels = c("passed", "filtered")
)

x_values <- c(data$x_left, data$x_right)
x_values <- x_values[is.finite(x_values)]

x_min <- min(x_values, na.rm = TRUE)
x_max <- max(x_values, na.rm = TRUE)
x_pad <- (x_max - x_min) * 0.05

x_lower <- if (x_min > 0) {
  x_min - min(x_pad, x_min * 0.05)
} else {
  x_min
}

x_upper <- x_max + x_pad
y_upper <- max(data$y, na.rm = TRUE) * 1.05

ggplot(data, aes(fill = status)) +
  geom_rect(
    aes(
      xmin = x_left,
      xmax = x_right,
      ymin = 0,
      ymax = y
    ),
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = python_status_colors,
    breaks = c("passed", "filtered"),
    labels = c("Passed filters", "Filtered out"),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(x_lower, x_upper),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, y_upper),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "First exon length (bp, log10)",
    y = "Number of transcripts (log10)"
  ) +
  theme_gtop() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = c(1, 1),
    legend.justification = c(1, 1)
  ) +
  geom_vline(
    xintercept = log10(30),
    color = "red",
    linetype = "dashed",
    linewidth = 0.5
  )

# Supp.Fig.12g -mismatch per 30 nt distribution for first exons longer than 30 nt---------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}

python_status_colors <- c(
  "passed" = "#4682B4",
  "filtered" = "#FFA07A"
)

data <- read.delim("Fig_S12g.txt", check.names = FALSE)

data$status <- factor(
  data$status,
  levels = c("passed", "filtered")
)

x_values <- c(data$x_left, data$x_right)
x_values <- x_values[is.finite(x_values)]

x_min <- min(x_values, na.rm = TRUE)
x_max <- max(x_values, na.rm = TRUE)
x_pad <- (x_max - x_min) * 0.05

x_lower <- if (x_min > 0) {
  x_min - min(x_pad, x_min * 0.05)
} else {
  x_min
}

x_upper <- x_max + x_pad
y_upper <- max(data$y, na.rm = TRUE) * 1.05

ggplot(data, aes(fill = status)) +
  geom_rect(
    aes(
      xmin = x_left,
      xmax = x_right,
      ymin = 0,
      ymax = y
    ),
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = python_status_colors,
    breaks = c("passed", "filtered"),
    labels = c("Passed filters", "Filtered out"),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(x_lower, x_upper),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, y_upper),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Mismatches per 30 nt (log10(x + 1))",
    y = "Number of transcripts (log10)"
  ) +
  theme_gtop() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = c(1, 1),
    legend.justification = c(1, 1)
  ) +
  geom_vline(
    xintercept = log10(1 + 1),
    color = "red",
    linetype = "dashed",
    linewidth = 0.5
  )

# Supp.Fig.9h unique short-read support for novel junctions----------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
data <- read.delim("Fig_S12h.txt", check.names = FALSE)
data$status <- "all"
ggplot(data, aes(x_center, y, fill = status)) +
  geom_col(aes(width = x_width), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("all" = "#7f8080"), drop = FALSE) +
  labs(x = "Total unique reads", y = "Count") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
# Supp.Fig.9i number of samples with >=5 reads for novel junctions----------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
data <- read.delim("Fig_S12i.txt", check.names = FALSE)
data$status <- "all"
ggplot(data, aes(x_center, y, fill = status)) +
  geom_col(aes(width = x_width), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("all" = "#7f8080"), drop = FALSE) +
  scale_y_continuous(breaks = seq(3, 5.5, by = 0.5))+
  coord_cartesian(ylim = c(2.5, 5.5)) +
  labs(x = "Samples with >=5 reads", y = "log10(count)") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


# Supp.Fig.9j -------------------------------------------------------------

library(patchwork)
theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
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
data <- read.delim("Fig_S12j.txt", check.names = FALSE)
tier_order <- c("Catalogue expansion", "High confidence", "With unique peptide")
data$tier <- factor(data$tier, levels = tier_order)
totals <- aggregate(count ~ tier, data, sum)
category_order <- names(sort(tapply(data$count, data$category, sum), decreasing = TRUE))
data$category <- factor(data$category, levels = category_order)
tier3_total <- totals$count[match("With unique peptide", as.character(totals$tier))]
left_xlim_max <- max(tier3_total * 1.4, 10, na.rm = TRUE)
right_xlim_min <- left_xlim_max * 1.1
right_xlim_max <- max(totals$count, na.rm = TRUE) * 1.05

base <- ggplot(data, aes(tier, count, fill = category)) +
  geom_col(width = 0.65, position = position_stack(reverse = TRUE)) +
  scale_x_discrete(limits = rev(tier_order)) +
  scale_fill_manual(values = structural_palette, drop = FALSE) +
  labs(x = NULL, y = "Number of transcripts", fill = "Structural category") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

left_breaks  <- c(0, 2500)
right_breaks <- c(100000, 200000)

p_left <- base +
  coord_flip(ylim = c(0, left_xlim_max), clip = "on") +
  scale_y_continuous(
    breaks = left_breaks,
    labels = function(x) format(x, big.mark = ",", trim = TRUE),
    expand = c(0, 0)
  ) +
  theme(legend.position = "none")

p_right <- base +
  coord_flip(ylim = c(right_xlim_min, right_xlim_max), clip = "on") +
  scale_y_continuous(
    breaks = right_breaks,
    labels = function(x) format(x, big.mark = ",", trim = TRUE),
    expand = c(0, 0)
  ) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

(p_left | p_right) + plot_layout(widths = c(1.2, 2.2))
