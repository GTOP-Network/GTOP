#==============================================#
# SV-eQTL-example IRGM #
# Supp-Figure-24#
#==============================================#
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(data.table)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig24")
# Supp.Fig.24b SV-eQTL-example IRGM -----------

dat <- fread("./input/Figure S24.txt",sep = "\t")
df_long <- dat %>%pivot_longer(cols = c(sv_geno),names_to = "geno_type",values_to = "geno") %>%
  mutate( geno = factor(geno, levels = c(0,1,2)),geno_type = factor(geno_type, levels = c("sv_geno", "snp_geno"))) %>% filter(!is.na(geno)) 

ggplot() +
  #geom_violin(data = df_long,aes(x = geno, y = pheno, group = geno), fill = NA,color = "black", width = 0.8) +
  geom_boxplot(data = df_long,aes(x = geno, y = pheno, group = geno),width = 0.5,fill = NA,color = "black") +
  geom_jitter(data = df_long,aes(x = geno, y = pheno, color = geno_type),width = 0.15,size = 2 ) +
  labs(x = "chr5_150823598_DEL_0M137_20103",y = "Normalized IRGM expression" ) +
  scale_color_manual(values = c("#dd9128")) +
  theme_classic(base_size = 14)+theme(legend.position = "none")




