#===============================================#
# MS peptide  #
#Supp-Figure-11 #
#===============================================#
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)

setwd("/path/to/GTOP_code/supp/supp_fig11")


# Supp.Fig.11a rho between Salmon and RSEM ----------------------------------

plot_df <- fread("./input/supp11a.SR_transcript_expr_two_method_corr.txt")


ggplot(plot_df)+
  geom_histogram(aes(x=salmon_rsem), color="white", fill  = "#aec1d9")+
  geom_vline(xintercept = median(plot_df$salmon_rsem),
             linetype="dashed", color="red")+
  xlab("Spearman rho between Salmon and RSEM")+
  ylab("Number of samples")+
  theme_pubr()



# Supp.Fig.11b  -------------------------------------------------------------

plot_df <- fread("./input/supp11b.LR_SR_corr.txt")

ggplot(plot_df, aes(x = level, y = correlation)) +
  geom_violin(alpha = 1, 
              width = 0.5,
              fill = "#aec1d9",
              trim = TRUE) +
  geom_boxplot(width = 0.15, 
               outlier.shape = NA,
               fill = "#aec1d9",
               alpha = 1) +
  geom_jitter(width = 0.1,
              size = 1.5,
              alpha = 0.5) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "",
       y = "Spearman rho") +
  theme_pubr()



# Supp.Fig.11c gene expression correlation between Illumina and Pacbio --------

plot_df <- fread("./input/supp11c.LR_SR_expr_corr_example.Gene.txt")

plot_df_filtered <- plot_df %>%
  filter(`short-read` > 0 & `long-read` > 0)

ggplot(plot_df_filtered) +
  geom_point(
    aes(x = log2(`short-read` + 1),
        y = log2(`long-read` + 1)),
    color = "grey60",
    size = 1
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  xlab("log2(TPM + 1) - Illumina") +
  ylab("log2(TPM + 1) - PacBio") +
  theme_pubr()


# Supp.Fig.11d transcript expression correlation between Illumina  --------


plot_df <- fread("./input/supp11d.LR_SR_expr_corr_example.Transcript.txt")

plot_df_filtered <- plot_df %>%
  filter(`short-read` > 0 & `long-read` > 0)

ggplot(plot_df_filtered) +
  geom_point(
    aes(x = log2(`short-read` + 1),
        y = log2(`long-read` + 1)),
    color = "grey60",
    size = 1
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  xlab("log2(TPM + 1) - Illumina") +
  ylab("log2(TPM + 1) - PacBio") +
  theme_pubr()

