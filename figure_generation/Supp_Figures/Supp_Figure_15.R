#==============================================#
#  ASE ASTS#
# Supp-Figure-15#
#==============================================#

library(dplyr)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrastr)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig15/input")


#Supp.Fig.15a the number of significant ase/asts -------------------------------

df_plot.m <- fread("Figure S15a.txt")
ggplot(df_plot.m,aes(x=group,y=Count,fill=Significance)) + geom_bar(stat="identity",width=.8) + theme_pubr() + 
  scale_fill_manual(breaks = c("no","yes"),values = c("#D3DCE5","#C0B5D7"))


# Supp.Fig.15b comparation bewteen SRS ASE and LRS ASE
df_plot <- fread("Figure S15b.txt")
ggplot(df_plot)+
  # geom_point(size=3, aes(x=logafc, y=as.numeric(log2_aFC)), alpha=0.6, color="grey80") +
  geom_point_rast(size=3, aes(x=logafc, y=as.numeric(log2_aFC)), alpha=0.6, color="#5d86af") +
  xlab("LRS log2aFC") + ylab("SRS log2aFC") +
  theme_pubr()

# Supp.Fig.15c concordance of ASE direction between LRS and SRS across thresholds

df_plot <- fread("Figure S15c.txt")

ggplot(df_plot)+
  geom_point(aes(x=thresh, y=ratio, fill=dir, color=dir))+
  geom_line(aes(x=thresh, y=ratio, fill=dir, color=dir))+
  scale_fill_manual(values=c("diff"="grey", "same"="#d8bbb3"))+
  scale_color_manual(values=c("diff"="grey", "same"="#d8bbb3"))+
  scale_y_continuous(limits = c(0, 1))+
  theme_pubr()+
  xlab("Threshold of ASE")+
  ylab("Proportion of significanct ASE")


# Supp.Fig.15d the distribution of ASE/ASTS events for donors
df_plot <- fread("Figure S15d.txt")
ggplot(df_plot,aes(x=count,y=Freq,color=group)) + 
  scale_color_manual(
  values = c(
    "ASE" = "#bcb3d1",
    "ASTS" = "#a4c48a"
  ))+
    geom_point(size=3) + theme_pubr()

# Supp.Fig.15e log2aFC in ASE events for donor 1 and donor 2

df_plot <- fread("Figure S15e.txt")
ggplot(df_plot,aes(x=logafc.x,y=logafc.y)) + geom_point(color="steelblue") + 
  xlab("log aFC in donor-1")+
  ylab("log aFC in donor-2")+
  theme_pubr()


# S15 f-g the distribution of ASE/ASTS genes

df_plot <- fread("Figure S15f-g.txt")

p1 <- ggplot(df_plot,aes(x=ase_tissues))+
  geom_bar(aes(y = after_stat(count / 1000)), stat = "count")+
  xlab("Number of tissues")+
  ylab("Number of ASE genes (x10³)")+
  theme_pubr()

p2 <- ggplot(df_plot,aes(x=asts_tissues))+
  geom_bar(aes(y = after_stat(count / 1000)), stat = "count")+
  xlab("Number of tissues")+
  ylab("Number of ASTS genes (x10³)")+
  theme_pubr()

cowplot::plot_grid(p1, p2, ncol = 2)
