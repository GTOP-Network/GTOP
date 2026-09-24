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
df_plot.eqtl <- readRDS("df_plot.eqtl.RDS")
df_plot.juqtl <- readRDS("df_plot.juqtl.RDS")
df_plot.tuqtl <- readRDS("df_plot.tuqtl.RDS")


p1 <- ggplot(df_plot.eqtl,aes(x=Var1,y=Freq)) + geom_bar(stat="identity",width=.8,fill="#aac79c") + theme_pubr() + 
  facet_wrap(~Var2,scales = "free_y",ncol=1);p1

p2 <- ggplot(df_plot.juqtl,aes(x=Var1,y=Freq)) + geom_bar(stat="identity",width=.8,fill="#7d8bad") + theme_pubr() + 
  facet_wrap(~Var2,scales = "free_y",ncol=1);p2

p3 <- ggplot(df_plot.tuqtl,aes(x=Var1,y=Freq,fill = Var3)) + geom_bar(stat="identity",width=.8,position = position_stack()) + theme_pubr() +
  scale_fill_manual(breaks = c("NO","YES"),values = c("#CF928F","#9D3929")) +
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


df_plot <- fread("Fig 3c.txt")

ggplot(df_plot, aes(x=slope_gtop, y=slope_gtex) ) +
  scale_fill_continuous(type = "viridis") +
  ggpointdensity::geom_pointdensity( alpha=.8, size=2, adjust=1, shape=21 )+
  scale_color_viridis_b()+
  ggdensity::geom_hdr_lines( linetype="dashed", linewidth=0.5)+
  labs(x="GTOP effect size", y="GTEx effect size")+
  theme_classic()+
  stat_cor(aes(x=slope_gtop, y=slope_gtex),method = "spearman", label.x =-3, label.y = 3)+
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



# Fig.3e:  SV effect size correlation with SV length ------------------------------------------------------------------


sv_qtl <- fread("Fig3e.sv_length_effect.txt")
sv_qtl$group <- factor(sv_qtl$group,levels = c("50-100bp","100-200bp","200-500bp","500-10kb",">10kb"))
groups2 <- levels(sv_qtl$group)
comparisons2 <- list( c("50-100bp", "100-200bp"), c("100-200bp", "200-500bp"),c("200-500bp", "500-10kb"),c("500-10kb", ">10kb"))
ggplot(sv_qtl, aes(x = group, y = abs(slope), fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, color = "black") +
  scale_fill_manual(values = c("#f1eef6","#bdc9e1","#74a9cf","#2b8cbe","#045a8d")) +
  stat_summary(aes(group = 1),fun = median,geom = "line",linewidth = 1,color = "black") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(colour = "black"),legend.position = "none") +
  labs(y = "|Effect Size|") +
  stat_compare_means( comparisons = comparisons2,method = "wilcox.test", label = "p",step.increase = 0.08 )


# Fig.3f:  Pathogenic TR QTL ------------------------------------------------------------------

count <- fread("Fig3f.pathogenic_TR-xQTL.txt") %>%dplyr::select(QTL, Tissue, TR_GeneName) %>%
  mutate(QTL=ifelse(QTL=="eQTL","eQTL","sQTL")) %>%
  distinct() %>%group_by(Tissue, TR_GeneName) %>%
  summarise(Class = case_when(all(QTL == "eQTL") ~ "eQTL",all(QTL == "sQTL") ~ "sQTL",
                              any(QTL == "eQTL") & any(QTL == "sQTL") ~ "eQTL & sQTL",
                              TRUE ~ "Other"),.groups = "drop")
dosage <- fread("Fig3f.pathogenic_TR_dosage.txt")
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


# Fig.3g: pathogenic TR-example --------------------------------------------------------------

dat2 <- fread("./Fig3g.junction2.txt") %>% mutate(CNV_f = factor(as.character(CNV), levels = sort(unique(CNV))))
ggplot(data = dat2, aes(x = CNV_f, y = pheno,colour = Tissue)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, na.rm = TRUE) +
  theme_classic()+theme(legend.position = "none")+
  scale_color_manual(values = c("Pancreas_Head"="#AB5F2C")) +
  geom_smooth(aes(x = as.numeric(factor(CNV_f)), y = pheno, group = Tissue),
              method = "lm", formula = y ~ x,se = T, size = 1, linetype = "solid")+
  labs(y = paste0( "Normalized Junction Usage", " | ", unique(dat2$pheno)))

dat1 <- fread("./Fig3g.junction1.txt") %>% mutate(CNV_f = factor(as.character(CNV), levels = sort(unique(CNV))))
ggplot(data = dat1, aes(x = CNV_f, y = pheno,colour = Tissue)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, na.rm = TRUE) +
  theme_classic()+theme(legend.position = "none")+
  scale_color_manual(values = c("Pancreas_Head"="#AB5F2C")) +
  geom_smooth(aes(x = as.numeric(factor(CNV_f)), y = pheno, group = Tissue),
              method = "lm", formula = y ~ x,se = T, size = 1, linetype = "solid")+
  labs(y = paste0( "Normalized Junction Usage", " | ", unique(dat1$pheno)))


# Fig.3g:  MASH -------------------------------------------------------------

plot_data_all <- fread("Fig3h.MASH.txt")
plot_data_all$Group <- factor(plot_data_all$Group,
                              levels = c("sv_eQTL", "tr_eQTL", "snv_eQTL"))
ggplot(plot_data_all, aes(x = factor(Number, levels = c("1","2","3","4","5","6","7","8","9","10","11")),
                                y = Density_sum, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width=0.9), color="black") +
  scale_x_discrete() +
  labs(x="Number of Tissues", y="Proportion of eQTLs", fill="Variant") +
  scale_y_continuous(limits = c(0, max(plot_data_all$Density_sum)+0.1)) +
  scale_fill_manual(values = c("sv_eQTL"="#c88565","tr_eQTL"="#931e2a","snv_eQTL"="#ebd1bf"))+
  theme_classic(base_size = 12) +
  theme(legend.position = "top")





