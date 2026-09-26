#===============================================#
# PacBio versus Illumina  #
#Supp-Figure-15 #
#===============================================#
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig15. Figure S11/input")


# Supp.Fig.15a PacBio versus Illumina sample Spearman correlations---------------------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
data <- read.delim("Fig_S15a.txt", check.names = FALSE)
data$level <- factor(data$level, levels = c("gene", "transcript"), labels = c("Gene", "Transcript"))
ggplot(data, aes(level, correlation, fill = level)) +
  geom_violin(width = 0.65, color = "black", alpha = 0.8) +
  geom_jitter(width = 0.08, size = 0.9, alpha = 0.4, color = "black") +
  stat_summary(fun.min = min, fun.max = max, fun = median,
               geom = "errorbar", width = 0.25, linewidth = 0.5, color = "#4470a4") +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.2)) +
  scale_fill_manual(values = c("Gene" = "#B0C4DE", "Transcript" = "#B0C4DE")) +
  coord_cartesian(ylim = c(0, 0.95)) +
  labs(x = NULL, y = "\u03c1 between PacBio and Illumina") +
  theme_gtop() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

# Supp.Fig.15b  example gene-level LR/SR correlation-------------------------------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
data <- read.delim("Fig_S15b.txt", check.names = FALSE)
rho <- suppressWarnings(cor(data$`short-read`, data$`long-read`, method = "spearman", use = "complete.obs"))
ggplot(data, aes(`short-read`, `long-read`)) +
  geom_point(size = 0.5, color="grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_y_continuous(breaks = seq(0, 17.5, by = 2.5)) +
  coord_cartesian(ylim = c(0, 17.5))+
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, label = sprintf("Spearman rho = %.3f", rho)) +
  labs(x = "log2(TPM + 1) of Illumina", y = "log2(TPM + 1) of PacBio", title = "Gene") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Supp.Fig.15c gene expression correlation between Illumina and Pacbio --------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
data <- read.delim("Fig_S15c.txt", check.names = FALSE)
rho <- suppressWarnings(cor(data$`short-read`, data$`long-read`, method = "spearman", use = "complete.obs"))
ggplot(data, aes(`short-read`, `long-read`)) +
  geom_point(size = 0.5, color="grey50") +
  scale_y_continuous(breaks = seq(0, 15, by = 2.5)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, label = sprintf("Spearman rho = %.3f", rho)) +
  labs(x = "log2(TPM + 1) of Illumina", y = "log2(TPM + 1) of PacBio", title = "Transcript") +
  theme_gtop() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
