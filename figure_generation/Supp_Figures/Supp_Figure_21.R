#==============================================#
#  isolaser resault #
# Supp-Figure-21#
#==============================================#
library(data.table)
library(dplyr)
library(magrittr)
library(tidyverse)
library(magrittr)
library(dplyr)
library(ggpubr)
library(patchwork)
setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig21.Figure.S15")


# Supp.Fig.21a  Number of unique cis-directed events and the numbe --------

df_count_by_sample <- fread("./input/Figure S21a.txt")
df_count_long <- melt(df_count_by_sample,id.vars = c("sample","donorID","Tissue","Tissue_Color_Code"))
df_count_long$sample <- factor(df_count_long$sample,levels = df_count_by_sample$sample)

p1 <- ggplot(df_count_long,aes(x=sample,y=value,fill=variable)) + 
  geom_bar(stat = "identity",width=.8,position = position_dodge()) + 
  geom_hline(yintercept = c(median(df_count_by_sample$eventCount),median(df_count_by_sample$geneCount)),linetype="dashed") +
  xlab("Number of events/genes")+
  theme_pubr() + 
  scale_fill_manual(
    values = c(eventCount = "#b9aece", geneCount = "#1d76b3")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) + ylab("Count");p1

df_count_by_sample$sample <- factor(df_count_by_sample$sample,levels = df_count_by_sample$sample)
p2 <- ggplot(df_count_by_sample,aes(x=sample,y=1,fill=Tissue)) + geom_tile() +
  theme_void() + scale_fill_manual(breaks = df_count_long$Tissue,values = paste0("#",df_count_long$Tissue_Color_Code)) + 
  theme(legend.position = "none");p2
library(patchwork)

p1 / p2  + plot_layout(heights = c(5, 0.5))


###Fig. S21b  Number of variants, exonic parts, and genes detected by mergingsamples of the same tissue type using isoLASER-joint


df_count_m <- fread("./input/Figure S21b.txt")

df_count_m$Tissue <- factor(df_count_m$Tissue,levels=unique(df_count_m$Tissue))

p3 <- ggplot(df_count_m,aes(x=Tissue,y=value,fill=variable)) + 
  geom_bar(stat = "identity",width=.8,position = position_dodge()) + 
  geom_hline(yintercept = c(median(df_count_m$Count_Variants),median(df_count_m$Count_exon),median(df_count_m$Count_gene)),
             linetype="dashed")  + theme_pubr() + 
  scale_fill_manual(
    values = c(Count_Variants = "#deebf7", Count_exon = "#9ecae1",Count_gene="#3182bd")) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) + ylab("Count");p3


color <- df_count_m %>% distinct(Tissue, Tissue_Color_Code)


p4 <- ggplot(color,aes(x=Tissue,y=1,fill=Tissue)) + geom_tile() +
  theme_void() + 
  scale_fill_manual(breaks = color$Tissue,
                    values = setNames(paste0("#",color$Tissue_Color_Code),
                                      color$Tissue)) + 
  theme(legend.position = "none");p4


p3 / p4  + plot_layout(heights = c(5, 0.5))

