#==============================================#
# QTL #
# Figure-3#
#==============================================#
library(data.table)
library(ggplot2)
library(stringi)
library(stringr)
library(dplyr)
library(ggsci)
library(tidyverse)
library(gg.gap)
library(ggbreak)
library(ggpubr)
library(magrittr)
library(ggpointdensity)
library(ggdensity)
library(viridis)
library(ggh4x)

setwd("/path/to/GTOP_code/fig-3/input")
# Fig.3a: summary of eQTL/sQTL ------------------
rm_version <- function(x){
  return(strsplit(x,split = ".",fixed = T)[[1]][1])
}

extract_sgenes <- function(x){
  return(strsplit(x,split=":",fixed=T)[[1]][6])
}

extract_transcript <- function(x){
  return(strsplit(x,split="_",fixed = T)[[1]][1])
}
#load eqtl data
## input files
df_eQTL <- readRDS("Fig3a.df_eQTL.RDS")
df_sQTL <- readRDS("Fig3a.df_sQTL.RDS")

df_juQTL <- df_sQTL %>% filter(QTLtype=="juQTL")
df_tuQTL <- df_sQTL %>% filter(QTLtype=="tuQTL")
df_tuQTL$transcript <- sapply(df_tuQTL$phenotype_id,extract_transcript)
df_tuQTL$novel <- "NO"
df_tuQTL$novel[grepl("PB",df_tuQTL$transcript)] <- "YES" # 22.37%

##count QTLs
df_plot.eqtl <- as.data.frame(table(df_eQTL$Tissue,df_eQTL$VarType))

df_snv_eqtl <- df_plot.eqtl %>% filter(Var2=="SNV") %>% arrange(desc(Freq))
tissue_orders_1 <- df_snv_eqtl$Var1
rm(df_snv_eqtl)
df_plot.eqtl$Var1 <- factor(df_plot.eqtl$Var1,levels = tissue_orders_1)
p1 <- ggplot(df_plot.eqtl,aes(x=Var1,y=Freq)) + geom_bar(stat="identity",width=.8,fill="#aac79c") + theme_pubr() + 
  facet_wrap(~Var2,scales = "free_y",ncol=1);p1


df_plot.juqtl <- as.data.frame(table(df_juQTL$Tissue,df_juQTL$VarType))
df_snv_juqtl <- df_plot.juqtl %>% filter(Var2=="SNV") %>% arrange(desc(Freq))
tissue_orders_2 <- df_snv_juqtl$Var1
rm(df_snv_juqtl)

df_plot.juqtl$Var1 <- factor(df_plot.juqtl$Var1,levels = tissue_orders_2)
p2 <- ggplot(df_plot.juqtl,aes(x=Var1,y=Freq)) + geom_bar(stat="identity",width=.8,fill="#7d8bad") + theme_pubr() + 
  facet_wrap(~Var2,scales = "free_y",ncol=1);p2

df_plot.tuqtl <- as.data.frame(table(df_tuQTL$Tissue,df_tuQTL$VarType,df_tuQTL$novel))
df_snv_tuqtl <- df_tuQTL %>%  filter(VarType=="SNV")
df_count_snv_tu <- as.data.frame(table(df_snv_tuqtl$Tissue))
df_count_snv_tu <- df_count_snv_tu[order(-df_count_snv_tu$Freq),]
tissue_orders_3 <- df_count_snv_tu$Var1
rm(df_count_snv_tu)

df_plot.tuqtl$Var1 <- factor(df_plot.tuqtl$Var1,levels = tissue_orders_3)
df_plot.tuqtl$Var3 <- factor(df_plot.tuqtl$Var3,levels = c("YES","NO"))
p3 <- ggplot(df_plot.tuqtl,aes(x=Var1,y=Freq,fill = Var3)) + geom_bar(stat="identity",width=.8,position = position_stack()) + theme_pubr() + 
  scale_fill_manual(breaks = c("NO","YES"),values = c("#CF928F","#9d3929")) +
  facet_wrap(~Var2,scales = "free_y",ncol=1);p3

cowplot::plot_grid(p1,p2,p3,ncol=3,align = "vh")


# Fig.3b: distance of finemapped eQTL to TSS ------------------
# load data
df_plot_snv <- readRDS("Fig3b.snv_to_tss.dist.RDS")
df_plot_sv <- readRDS("Fig3b.sv_to_tss.dist.RDS")
df_plot_tr <- readRDS("Fig3b.tr_to_tss.dist.RDS")

p1 <- ggplot(df_plot_snv,aes(x=distance,color=group)) + geom_density(size=1.5) + theme_pubr() + 
  scale_color_manual(values = c("#fee0d2","#fc9272","#de2d26","grey"));p1

p2 <- ggplot(df_plot_sv,aes(x=distance,color=group)) + geom_density(size=1.5) + theme_pubr() + 
  scale_color_manual(values = c("#fee0d2","#fc9272","#de2d26","grey"));p2


p3 <- ggplot(df_plot_tr,aes(x=distance,color=group)) + geom_density(size=1.5) + theme_pubr() + 
  scale_color_manual(values = c("#fee0d2","#fc9272","#de2d26","grey"));p3


cowplot::plot_grid(p1,p2,p3,ncol = 1,align = "v")



# Fig.3c: effect compare between GTEx and GTOP ----------------------------


df_plot <- fread("Fig3c.txt")

ggplot(df_plot, aes(x=slope_gtop, y=slope_gtex) ) +
  scale_fill_continuous(type = "viridis") +
  ggpointdensity::geom_pointdensity( alpha=.8, size=2, adjust=1, shape=21 )+
  scale_color_viridis_b()+
  ggdensity::geom_hdr_lines( linetype="dashed", linewidth=0.5)+
  labs(x="GTOP effect size", y="GTEx effect size")+
  theme_classic()+
  stat_cor(aes(x=slope_gtop, y=slope_gtex),method = "pearson", label.x =-3, label.y = 3)+
  theme(
    axis.line = element_line(color="black", linewidth=1),
    axis.ticks = element_line(color="black", linewidth=1),
    axis.text.x = element_text(size=rel(2), color="black"),
    axis.text.y=element_text(size=rel(2), color="black"),
    axis.title = element_text(size=rel(1.7)),
    legend.title = element_blank(),
    legend.text = element_text(size=rel(1.2)),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm")
  ) +
  ggh4x::coord_axes_inside(labels_inside = F)


# Fig.3d: venn plot- eGenes sharing by different variant types ----------------------

library(ggVennDiagram)
egene_list <- readRDS("Fig3d.eGenes_by_VarType.RDS")
set.seed(101)
ggVennDiagram(egene_list) + scale_fill_gradient(low="grey90",high = "red")


# Fig.3e:  Pathogenic TR QTL ------------------------------------------------------------------

count <- fread("Fig3e.pathogenic_TR-xQTL.txt") %>%dplyr::select(QTL, Tissue, TR_GeneName) %>%
  mutate(QTL=ifelse(QTL=="eQTL","eQTL","sQTL")) %>%
  distinct() %>%group_by(Tissue, TR_GeneName) %>%
  summarise(Class = case_when(all(QTL == "eQTL") ~ "eQTL",all(QTL == "sQTL") ~ "sQTL",
                              any(QTL == "eQTL") & any(QTL == "sQTL") ~ "eQTL & sQTL",
                              TRUE ~ "Other"),.groups = "drop")
dosage <- fread("Fig3e.pathogenic_TR_dosage.txt")
dosage <- dosage %>% mutate(TR_GeneName=factor(TR_GeneName,levels=c(dosage %>% group_by(TR_GeneName) %>% 
                                                   summarise(cnv=median(CNV,na.rm = TRUE))%>% 
                                                   arrange(cnv) %>% pull(TR_GeneName))))
tissue_order <- count %>%
  group_by(Tissue) %>%
  summarise(total = n(), .groups = "drop") %>%
  arrange(total) %>%
  pull(Tissue)

count$Class <- factor(count$Class, levels = c("eQTL", "sQTL", "eQTL & sQTL"))
count$TR_GeneName <- factor(count$TR_GeneName, levels = levels(dosage$TR_GeneName))
count$Tissue  <- factor(count$Tissue,  levels = rev(tissue_order))

ggplot(count, aes(x = TR_GeneName, y = Tissue, fill = Class)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("eQTL" = "#aac79c","sQTL" = "#808bab","eQTL & sQTL" ="#c7b3cf")) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid = element_blank(),legend.position = "bottom")+
  coord_flip()

ggplot(dosage, aes(x = TR_GeneName, y = CNV)) +
  geom_boxplot(aes(fill=GeneRegion)) +
  scale_fill_npg() +
  theme_bw() +
  scale_y_log10()+
  theme(axis.text = element_text(colour = "black")) +
  coord_flip()


# Fig.3F: pathogenic TR-example --------------------------------------------------------------


dat <- fread("./Fig3f.txt") %>% mutate(CNV_f = factor(as.character(CNV), levels = sort(unique(CNV))))

ggplot(data = dat, aes(x = CNV_f, y = pheno)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1, na.rm = TRUE) +
  geom_smooth(aes(x = as.numeric(factor(CNV_f)), y = pheno, group = Tissue),
              method = "lm", formula = y ~ x,se = T, size = 1, linetype = "solid") +
  labs( x= "MUC1 60bp VNTR length (number of repeats)",
        color = "Tissue") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        strip.text = element_text(size = 10, face = "bold")) +
  facet_grid( ~ Tissue, scales = "free_x", space = "free_x")


# Fig.3g:  MASH -------------------------------------------------------------

plot_data_all <- fread("Fig3g.MASH.txt")
ggplot(plot_data_all, aes(x = factor(Number, levels = c("1","2","3","4","5","6","7","8","9","10","11")),
                                y = Density_sum, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width=0.9), color="black") +
  scale_x_discrete() +
  labs(x="Number of Tissues", y="Proportion of eQTLs", fill="Variant") +
  scale_y_continuous(limits = c(0, max(plot_data_all$Density_sum)+0.1)) +
  scale_fill_manual(values = c("sv_eQTL"="#c88565","tr_eQTL"="#931e2a","snv_eQTL"="#ebd1bf"))+
  theme_classic(base_size = 12) +
  theme(legend.position = "top")





