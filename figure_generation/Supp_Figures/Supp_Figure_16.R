#===============================================#
# LRS transcript  #
#Supp-Figure-16 #
#===============================================#
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig16. R2.6.LRS_saturation/input")

# Fig S16a: paired annotated/novel transcript saturation increment at cutoff 10 --------


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
data <- read.delim("Fig_S16a.txt", check.names = FALSE)
data$transcript_label <- factor(data$transcript_label, levels = c("Annotated", "Novel"))
data$group_label <- pretty_label(data$group_label)
color_col <- first_existing_col(data, c("group_color"), NA_character_)
color_values <- setNames(color_col[match(unique(data$group_label), data$group_label)], unique(data$group_label))
color_values[is.na(color_values)] <- "#7f8080"
ggplot(data, aes(sample_size, incremental_detection_ratio, color = group_label, group = group_label)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.2) +
  geom_hline(yintercept = 0, color = "red", linewidth = 0.3,linetype = "dashed" ) +
  facet_wrap(~transcript_label, nrow = 1) +
  scale_color_manual(values = color_values, drop = FALSE) +
  scale_y_continuous(breaks = seq(0,0.6,by=0.1))+
  coord_cartesian(ylim = c(0, 0.6))+
  labs(x = "Sample size", y = "Added transcript proportion") +
  theme_gtop() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.background = element_blank()
  )
#Fig S16b: tissue-specific transcript proportion by sample size, Tau > 0.85.---------------------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
data <- read.delim("Fig_S16b.txt", check.names = FALSE)
data$transcript_label <- factor(data$transcript_label, levels = c("Annotated", "Novel"))
ggplot(data, aes(sample_size, mean_specificity_ratio, color = transcript_label, group = transcript_label)) +
  geom_errorbar(
    aes(
      ymin = mean_specificity_ratio - sd_specificity_ratio,
      ymax = mean_specificity_ratio + sd_specificity_ratio
    ),
    width = 0.18,
    linewidth = 0.35
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Annotated" = "#4C72B0", "Novel" = "#C44E52")) +
  coord_cartesian(ylim = c(0.1, 0.7), xlim = c(0.3, 10.7)) +
  scale_x_continuous(breaks = sort(unique(data$sample_size))) +
  scale_y_continuous(breaks = seq(0.1,0.7,by=0.1))+
  coord_cartesian(ylim = c(0.1, 0.7))+
  labs(x = "Sample size per tissue", y = "Proportion of transcripts (Tau > 0.85)") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

