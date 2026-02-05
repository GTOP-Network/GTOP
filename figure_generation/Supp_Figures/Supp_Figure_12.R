#===============================================#
# MS peptide  #
#Supp-Figure-11 #
#===============================================#
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig12")


# Supp.Fig.12a rho between Salmon and RSEM ----------------------------------

plot_df <- fread("./input/supp12a.SR_transcript_expr_two_method_corr.txt")


ggplot(plot_df)+
  geom_histogram(aes(x=salmon_rsem), color="white")+
  geom_vline(xintercept = median(plot_df$salmon_rsem),
             linetype="dashed", color="red")+
  xlab("Spearman rho between Salmon and RSEM")+
  ylab("Number of samples")+
  theme_pubr()



# Supp.Fig.12b  -------------------------------------------------------------

plot_df <- fread("./input/supp12b.LR_SR_corr.txt")

ggplot(plot_df, aes(x = level, y = correlation)) +
  geom_violin(alpha = 0.3, 
              width = 0.7,
              trim = TRUE) +
  geom_boxplot(width = 0.15, 
               outlier.shape = NA,
               fill = "white",
               alpha = 0.7) +
  geom_jitter(width = 0.1,
              size = 1.5,
              alpha = 0.5) +
  labs(x = "",
       y = "Spearman rho") +
  theme_pubr()



# Supp.Fig.12c gene expression correlation between Illumina and Pacbio --------

plot_df <- fread("./input/supp12c.LR_SR_expr_corr_example.Gene.txt")

ggplot(plot_df)+
  geom_point(aes(x=log2(`short-read`+1),
                 y=log2(`long-read`+1)))+
  xlab("log2(TPM + 1) of Illumina")+
  ylab("log2(TPM + 1) of PacBio")+
  theme_pubr()


# Supp.Fig.12d transcript expression correlation between Illumina  --------


plot_df <- fread("./input/supp12d.LR_SR_expr_corr_example.Transcript.txt")

ggplot(plot_df)+
  geom_point(aes(x=log2(`short-read`+1),
                 y=log2(`long-read`+1)))+
  xlab("log2(TPM + 1) of Illumina")+
  ylab("log2(TPM + 1) of PacBio")+
  theme_pubr()
