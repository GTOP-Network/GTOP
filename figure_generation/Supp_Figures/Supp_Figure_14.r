#==============================================#
# LRS vs SRS phased SNP comparation #
# Supp-Figure-14#
#==============================================#
library(data.table)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyverse)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig14")

# Supp.Fig.14a ------------------------------------------------------------
dat.m <- fread("./input/Figure S4.txt") %>% 
  mutate(chromosome=factor(chromosome, levels=seq(1:22)))

# phased counts
p1 <- ggplot(dat.m,aes(x=chromosome,y=phased,color=group)) + geom_boxplot(width=.8,outlier.shape = NA) + theme_pubr() + 
  scale_color_manual(breaks = c("SRS","LRS"),values = c("#7888a4","#a63b2a"));p1

# blocks
p2 <- ggplot(dat.m,aes(x=chromosome,y=blocks,color=group)) + geom_boxplot(width=.8,outlier.shape = NA) + theme_pubr() + 
  scale_color_manual(breaks = c("SRS","LRS"),values = c("#7888a4","#a63b2a"));p2

# variant_per_block_avg
p3 <- ggplot(dat.m,aes(x=chromosome,y=variant_per_block_avg,color=group)) + geom_boxplot(width=.8,outlier.shape = NA) + theme_pubr() + 
  scale_color_manual(breaks = c("SRS","LRS"),values = c("#7888a4","#a63b2a"));p3

# bp_per_block_avg
p4 <- ggplot(dat.m,aes(x=chromosome,y=bp_per_block_avg,color=group)) + geom_boxplot(width=.8,outlier.shape = NA) + theme_pubr() +
  scale_color_manual(breaks = c("SRS","LRS"),values = c("#7888a4","#a63b2a"));p4

# heterozygous_variants
p5 <- ggplot(dat.m,aes(x=chromosome,y=heterozygous_variants,color=group)) + geom_boxplot(width=.8,outlier.shape = NA) + theme_pubr() +
  scale_color_manual(breaks = c("SRS","LRS"),values = c("#7888a4","#a63b2a"));p5


# phased fraction
p6 <- ggplot(dat.m,aes(x=chromosome,y=phased_fraction,color=group)) + geom_boxplot(width=.8,outlier.shape = NA) + theme_pubr() +
  scale_color_manual(breaks = c("SRS","LRS"),values = c("#7888a4","#a63b2a")) ;p6


cowplot::plot_grid(p1, p6, p2,p3,p4,p5,ncol=1,align = "v")


# Supp.Fig.14a plot phased variants by individuals ------------------------


df_plot <- dat.m %>% dplyr::select(Indis, chromosome, variants, phased, unphased, group) %>% 
  dplyr::group_by(Indis, group) %>% 
  dplyr::summarise(total_var=sum(variants),total_phased=sum(phased),total_unphased=sum(unphased)) %>% ungroup()
  # rbind(df_plot.l,df_plot.s)

p8 <- ggplot(df_plot,aes(x=group,y=total_phased,color=group)) + geom_boxplot() + geom_jitter(width = .2,aes(color=group)) + theme_pubr() + 
  scale_color_manual(breaks = c("SRS","LRS"),values = c("#7888a4","#a63b2a"));p8



df_plot <- dat.m %>% dplyr::select(Indis, chromosome, variants, blocks, group) %>% 
  dplyr::group_by(Indis, group) %>% 
  dplyr::summarise(total_block=sum(blocks)) %>% ungroup()
  # rbind(df_plot.l,df_plot.s)
p9 <- ggplot(df_plot,aes(x=group,y=total_block,color=group)) + geom_boxplot() + geom_jitter(width = .2,aes(color=group)) + theme_pubr() + 
  scale_color_manual(breaks = c("SRS","LRS"),values = c("#7888a4","#a63b2a")) + scale_y_log10();p9



df_plot <- dat.m %>% dplyr::select(Indis, chromosome, variants, heterozygous_variants, group) %>% 
  dplyr::group_by(Indis, group) %>% 
  dplyr::summarise(total_het=sum(heterozygous_variants)) %>% ungroup()
  # rbind(df_plot.l,df_plot.s)
p10 <- ggplot(df_plot,aes(x=group,y=total_het,color=group)) + geom_boxplot() + geom_jitter(width = .2,aes(color=group)) + theme_pubr() + 
  scale_color_manual(breaks = c("SRS","LRS"),values = c("#7888a4","#a63b2a"));p10


cowplot::plot_grid(p8,p9,p10,ncol=1,align="v")

