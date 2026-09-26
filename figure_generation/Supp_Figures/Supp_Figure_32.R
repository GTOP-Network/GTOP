#===============================================#
# SV eQTL example  #
#Supp-Figure-32 #
#===============================================#
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(data.table)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig32.Figure.S24-SV-eQTL-example/input")


# supp_fig32.b ------------------------------------------------------------

data<-readRDS("supp32b.RDS")
df_long <- data %>%pivot_longer(cols = c(sv_geno),names_to = "geno_type",values_to = "geno") %>%
  mutate( geno = factor(geno, levels = c(0,1,2)),geno_type = factor(geno_type, levels = c("sv_geno", "snp_geno"))) %>% filter(!is.na(geno)) 

p<-ggplot() +
  #geom_violin(data = df_long,aes(x = geno, y = pheno, group = geno), fill = NA,color = "black", width = 0.8) +
  geom_boxplot(data = df_long,aes(x = geno, y = IRGM, group = geno),width = 0.5,fill = NA,color = "black",outlier.shape = NA) +
  geom_jitter(data = df_long,aes(x = geno, y = IRGM, color = geno_type),width = 0.15,size = 2 ) +
  labs(x = "chr5_150823598_DEL_0M137_20103",y = "Normalized IRGM expression" ) +
  scale_color_manual(values = c("grey50")) +
  theme_classic(base_size = 14)+theme(legend.position = "none");p

df_long$geno_num <- as.numeric(as.character(df_long$geno))

fit <- lm(IRGM ~ geno_num, data = df_long)
summary(fit)
