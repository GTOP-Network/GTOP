#==============================================#
# SV_Finemap #
# Supp-Figure-39#
#==============================================#
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(data.table)
library(dplyr)
library(magrittr)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig39.Figure.S27-SV-finemap-example-SIGIRR/input")


# Supp.Fig.39a --SV-TR enrichment----------------------------------------------------

df_m<-readRDS("SVTR_enrich_causal_CS.compare_to_SNV.RDS")

df_m$PIP_num <- df_m$PIP
df_m$PIP <- factor(df_m$PIP)
df_tissue_color <- read.csv("GMTiP_tissue_code_and_colors.csv",header = T)
df_tissue_color %<>% filter(Tissue %in% unique(df_m$Tissue))
df_m$PIP_num <- as.numeric(as.character(df_m$PIP_num))

p1 <- ggplot(df_m,aes(x=PIP,y=OR)) + geom_boxplot(outlier.shape = NA,aes(fill = PIP_num),alpha=.8) +
  geom_jitter(aes(color=Tissue),width=.2,size=1) + theme_pubr() +
  scale_color_manual(breaks = df_tissue_color$Tissue,values = paste0("#",df_tissue_color$Tissue_Color_Code)) +
  scale_fill_gradient(low = "#fff7ec",high="#b30000")
print(p1)

res <- lm(OR ~ PIP_num, data=df_m)
summary(res)


# Supp.Fig.39b --SV-TR enrichment filter AF -------------------------------

df_m <- readRDS("SVTR_enrich_causal_CS.compare_to_SNV.filter_by_AF_10.RDS")
df_m$PIP_num <- df_m$PIP
df_m$PIP <- factor(df_m$PIP)
df_tissue_color <- read.csv("GMTiP_tissue_code_and_colors.csv",header = T)
df_tissue_color %<>% filter(Tissue %in% unique(df_m$Tissue))

p1 <- ggplot(df_m,aes(x=PIP,y=OR)) + geom_boxplot(outlier.shape = NA,aes(fill = PIP_num),alpha=.8) +
  geom_jitter(aes(color=Tissue),width=.2,size=1) + theme_pubr() +
  scale_color_manual(breaks = df_tissue_color$Tissue,values = paste0("#",df_tissue_color$Tissue_Color_Code)) +
  scale_fill_gradient(low = "#fff7ec",high="#b30000")
print(p1)

res <- lm(OR ~ PIP_num, data=df_m)
summary(res)


# Supp.Fig.39d-e example NFATC1 ---------------------------------------------
load("supp39d.RData")
df_long <- data %>%
  pivot_longer(
    cols = c(chr18_79400832_INS_5CM493_172, rs2596606, rs657693),
    names_to = "variant_id",
    values_to = "genotype"
  ) %>%
  mutate(
    VarType = case_when(
      grepl("chr18_79400832_INS_5CM493_172", variant_id) ~ "SV",  
      TRUE ~ "SNV"
    ),
    genotype = factor(genotype, levels = c(0, 1, 2))
  )

df_long$variant_id <- factor(
  df_long$variant_id,
  levels = c("rs2596606", "rs657693", "chr18_79400832_INS_5CM493_172")
)
df_long<-df_long[!is.na(df_long$genotype),]
p<-ggplot(df_long, aes(x = genotype, y = NFATC1)) +
  geom_boxplot(
    aes(group = interaction(genotype, variant_id), fill = VarType),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_point(
    aes(group = interaction(genotype, variant_id)),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 1.3, color = "black", alpha = 0.7
  ) +
  scale_fill_manual(values = c("SNV" = "#4A7C7C", "SV" = "#B8C4E0")) +
  labs(
    title = expression(italic("NFATC1")),
    x = "Genotype (0/1/2)",
    y = "Normalized expression",
    fill = ""
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, hjust = 0),
    legend.position = "top",
    legend.justification = "right",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  );p

dat<-readRDS("supp39e.RDS")
p<-ggplot(dat,aes(x=pos1,y=pip,color=cs)) + geom_point(size=1.5) + theme_pubr() + 
  scale_color_manual(breaks = c("L1","not_cs_var"),values=c("red","grey"));p

