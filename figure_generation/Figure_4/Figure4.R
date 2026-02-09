#==================================#
# Population-specific QTLs #
# Figure-4 #
#==================================#

setwd("/media/london_A/mengxin/GTOP_code/fig-4")

library(data.table)
library(tidyverse)
library(reticulate)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(patchwork)


## Fig.4a, Geographic frequencies of eQTL -------------------------------------

if (!py_module_available("matplotlib")) {
  py_install("matplotlib")
}

if (!py_module_available("geovar")) {
  py_install("git+https://github.com/aabiddanda/geovar")
}

py_run_string('
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import pkg_resources
from geovar import *

plt.rcParams[\'pdf.fonttype\'] = 42

geovar_test = GeoVar()
geovar_test.add_freq_mat("./input/Fig4a.txt")
geovar_test.geovar_binning()

geovar_plot = GeoVarPlot()
geovar_plot.add_data_geovar(geovar_test)
geovar_plot.filter_data()
geovar_plot.add_cmap()

fig, ax = plt.subplots(1,1,figsize=(3,6))
geovar_plot.plot_geovar(ax)
ax.set_xticklabels(geovar_plot.poplist)
plt.savefig("./Fig4a.pdf", dpi=300, bbox_inches=\'tight\')
')


## Fig.4b, Allele frequency difference -----------------------------------------

freq_plot_data <- fread("./input/Fig4b.txt")

ggplot(freq_plot_data, aes(x = af_GTEx, y = af_GTOP)) +
  geom_point(aes(color = proxy_delta_group), alpha = 0.6) +
  scale_color_manual(
    name = "Frequency Difference",
    values =  c("#f8e4da", "#f1cf9c", "#e6986d", "#cf5e47", "#a43331", "#5c95d6") # , "#475c92"
  ) +
  labs(
    x = "Allele frequency (1KGP European)",
    y = "Allele frequency (GTOP)"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  ) +
  geom_hline(yintercept=c(0.05, 0.95), linetype='dashed', color='grey', size=0.5) +
  geom_abline(linetype='dashed', color='darkgrey', size=0.5)

length(freq_plot_data$variant_proxy[freq_plot_data$proxy_delta_group%in%
                                      c("|ΔAF| > 0.2", "|ΔAF| > 0.3", 
                                        "|ΔAF| > 0.4", "|ΔAF| > 0.5")])/
  length(freq_plot_data$variant_proxy)


## Fig.4c, Distribution of fd-QTL variants -------------------------------------

freq_for_type <- fread("./input/Fig4c.txt")

ggplot(data = freq_for_type, aes(x=proxy_delta_group, y=value, fill=variable)) +
  geom_col(position = 'dodge') +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 30, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_fill_manual(values = c( "#7d8cad","#0f3d7b")) +
  labs(x="", y="Proportion of fine-mapped variants")


## Fig.4d, Distribution of fd-QTL genes ----------------------------------------

ff_gene_summary <- fread("./input/Fig4d.txt")
color_vec <- readRDS("./input/tissue_color.RDS")

ggplot(ff_gene_summary, aes(x = Gene_count, y = new_two_group)) +
  geom_boxplot(outlier.colour = NA, aes(fill = new_two_group), width=0.8) +
  ggbeeswarm::geom_quasirandom(aes(group=tissue, col=tissue), size=2, width = 0.4) +
  scale_color_manual(values = color_vec) +
  scale_fill_manual(values = c("#c0c3d1", "#959db7")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x="", y="Gene number")


## Fig.4e, fd-eQTL example -----------------------------------------------------
load("./input/Fig4e.RData")

SNP_name <- GWAS_data$rsid[which.min(GWAS_data$p)]
SV_name <- GWAS_data$rsid[which.min(GWAS_data$p)]
loci_start <- min(GWAS_data$pos) + 200000
loci_end <- max(GWAS_data$pos) - 400000
GWAS_name <- "Urticaria"
chr_name <- GWAS_data$chrom[1]

GWAS_locus <- locus(data = as.data.frame(GWAS_data),
                    xrange = c(loci_start, loci_end),
                    seqname = chr_name, index_snp = SNP_name,
                    ens_db = "EnsDb.Hsapiens.v86")
GWAS_locus <- link_LD(GWAS_locus, token = "b6336e5da5d3", pop = "EAS")
GWAS_plot <- gg_scatter(GWAS_locus, pcutoff = FALSE, yzero=T,  labels = SNP_name, legend_pos = "right",
                        LD_scheme = c("#e5e5e5", "#e5e5e5", "#3e70b4", "#3f7d1d",
                                      "orange", "red", "red")) + 
  annotate("text",x=loci_end/10^6, y=max(-log10(GWAS_data$p)),label=GWAS_name, hjust="right") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  ylim(c(0, max(-log10(GWAS_data$p)+1))) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", size = 1) 

SNV_eQTL_locus <- locus(data = as.data.frame(SNV_eQTL_data),
                        xrange = c(loci_start, loci_end),
                        seqname = chr_name, index_snp = SNP_name,
                        ens_db = "EnsDb.Hsapiens.v86")
SNV_eQTL_locus <- link_LD(SNV_eQTL_locus, token = "b6336e5da5d3", pop = "EAS")
SNV_eQTL_plot <- gg_scatter(SNV_eQTL_locus, pcutoff = FALSE, yzero=T,  labels = SNP_name, legend_pos = "right",
                            LD_scheme = c("#e5e5e5", "#e5e5e5", "#3e70b4", "#3f7d1d",
                                          "orange", "red", "red")) + 
  annotate("text",x=loci_end/10^6, y=max(-log10(SNV_eQTL_data$p)),label="eQTL_Whole Blood", hjust="right") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

gene_plot <- gg_genetracks(SNV_eQTL_locus, highlight = "STIM1",
                           filter_gene_name = c("STIM1"),
                           filter_gene_biotype = c("protein_coding"))

wrap_plots(list(GWAS_plot, SNV_eQTL_plot, gene_plot), ncol = 1, heights = c(3,3,1))


## Fig.4f he-QTL ---------------------------------------------------------------
count_df <- fread("./input/Fig4f.txt")
color_vec <- readRDS("./input/tissue_color.RDS")

ggplot(count_df, aes(count, reorder(tissue, count), color=tissue)) +
  ylab('') + xlab('Number of he-QTLs') +
  geom_segment(aes(x=0, xend=count, y=tissue, yend=tissue)) +
  geom_point(shape=16, size=4) +
  scale_color_manual(values = color_vec)+
  theme_classic(base_size = 12) +
  theme(text=element_text(size=12),
        legend.position = "none",
        legend.title = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))


## Fig.4g Z-score of he-QTL between GTOP and GTEx ------------------------------
count_df <- fread("./input/Fig4g.txt")

ggplot(data = count_df, aes(x = type1, y = abs(value), fill=type1)) + 
  geom_violin() +
  geom_boxplot(width=0.1, outliers = FALSE) +
  theme_classic() +
  scale_fill_manual(values = c("GTEx"="#e7c798", "GTOP"="#bd5a45")) +
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x="", y="|Z-score|")


# Fig.4h cs size compare----------------------------------------
library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)

dat.gtex <- readRDS("./input/dat.gtex.RDS")
dat.gtop <- readRDS("./input/dat.gtop.RDS")
dat.cran <- readRDS("./input/dat.cran.RDS")
# examine CS number per gene
df.gtex <- dat.gtex %>% select(locus_id,cs,Tissue) %>% distinct(.keep_all = T)
df.gtop <- dat.gtop %>% select(locus_id,cs,Tissue) %>% distinct(.keep_all = T)
df.cran <- dat.cran %>% select(Gene,CS_ID,Tissue) %>% distinct(.keep_all = T)


df_count.cran <- df.cran %>% group_by(Gene,Tissue) %>% summarise(cs_count=n())
df_count.gtex <- df.gtex %>% group_by(locus_id,Tissue) %>% summarise(cs_count=n())
df_count.gtop <- df.gtop %>% group_by(locus_id,Tissue) %>% summarise(cs_count=n())

df_count.cran$group <- "GTEx+GTOP"
df_count.gtex$group <- "GTEx"
df_count.gtop$group <- "GTOP"

names(df_count.gtex) <- c("Gene","Tissue","cs_count","group")
names(df_count.gtop) <- c("Gene","Tissue","cs_count","group")
df_count <- rbind(df_count.cran,df_count.gtex,df_count.gtop)

# extract genes with single cs
df_count.cran$tissue_gene <- paste(df_count.cran$Gene,df_count.cran$Tissue,sep=":")
df_count.gtex$tissue_gene <- paste(df_count.gtex$Gene,df_count.gtex$Tissue,sep = ":")
df_count.gtop$tissue_gene <- paste(df_count.gtop$Gene,df_count.gtop$Tissue,sep = ":")
x.cran <- as.character(df_count.cran$tissue_gene[df_count.cran$cs_count==1])
x.gtex <- as.character(df_count.gtex$tissue_gene[df_count.gtex$cs_count==1])
x.gtop <- as.character(df_count.gtop$tissue_gene[df_count.gtop$cs_count==1])

overlap_tissuegenes <- intersect(intersect(x.gtex,x.gtop),x.cran)

df.cran <- dat.cran %>% filter(tissue_gene %in% overlap_tissuegenes) %>% select(Gene,CS_ID,MAX_PIP,CS_LENGTH,Tissue,tissue_gene)
df.gtex <- dat.gtex %>% filter(tissue_gene %in% overlap_tissuegenes) %>% select(locus_id,cs,pip,cs_size,Tissue,tissue_gene) %>% 
  group_by(locus_id,cs,Tissue) %>% mutate(max_pip = max(pip)) %>% ungroup() %>% filter(pip==max_pip) %>% 
  select(locus_id,cs,max_pip,cs_size,Tissue,tissue_gene) %>% distinct(.keep_all = T)
df.gtop <- dat.gtop %>% filter(tissue_gene %in% overlap_tissuegenes) %>% select(locus_id,cs,pip,cs_size,Tissue,tissue_gene) %>% 
  group_by(locus_id,cs,Tissue) %>% mutate(max_pip = max(pip)) %>% ungroup() %>% filter(pip==max_pip) %>% 
  select(locus_id,cs,max_pip,cs_size,Tissue,tissue_gene) %>% distinct(.keep_all = T)

df.cran$group <- "GTEx+GTOP"
df.gtex$group <- "GTEx"
df.gtop$group <- "GTOP"
names(df.cran) <- c("locus_id","cs","max_pip","cs_size","Tissue","tissue_gene","group")
df_plot1 <- rbind(df.cran,df.gtex,df.gtop)
df_plot1$group <- factor(df_plot1$group,levels = c("GTOP","GTEx","GTEx+GTOP"))
# all tissue combined
p1 <- ggplot(df_plot1,aes(x=group,y=log2(cs_size))) +  geom_violin(aes(fill=group)) + 
  geom_boxplot(width=.2,fill="white") + theme_pubr() + 
  scale_fill_manual(breaks =c("GTOP","GTEx","GTEx+GTOP"), values = c("#B65844","#E2C396","#7784A3") );p1 # n=3242 cs


x <- df.cran %>% mutate(maxPIP_cran=max_pip,CS_size_cran=cs_size) %>% select(tissue_gene,maxPIP_cran,CS_size_cran)
y <- df.gtex %>% mutate(maxPIP_gtex=max_pip,CS_size_gtex=cs_size) %>% select(tissue_gene,maxPIP_gtex,CS_size_gtex) 
z <- df.gtop %>% mutate(maxPIP_gtop=max_pip,CS_size_gtop=cs_size) %>% select(tissue_gene,maxPIP_gtop,CS_size_gtop) 
xy <- merge(x,y,by="tissue_gene")
xyz <- merge(xy,z,by="tissue_gene")

test1 <- wilcox.test(xyz$CS_size_cran,xyz$CS_size_gtex,alternative = "less",paired = T) # p<2.2e-16
test2 <- wilcox.test(xyz$CS_size_cran,xyz$CS_size_gtop,alternative = "less",paired = T) # p<2.2e-16



# Fig.4i plot prop of PIP>0.8 in single population fine-mapping and cross-ancestry fine-mapping------

overlap_tissuegenes <- intersect(dat.cran$tissue_gene,intersect(dat.gtex$tissue_gene,dat.gtop$tissue_gene))
df.cran <- dat.cran %>% filter(tissue_gene %in% overlap_tissuegenes)
df.gtex <- dat.gtex %>% filter(tissue_gene %in% overlap_tissuegenes) %>% group_by(tissue_gene,cs) %>% summarize(max_pip=max(pip))
df.gtop <- dat.gtop %>% filter(tissue_gene %in% overlap_tissuegenes) %>% group_by(tissue_gene,cs) %>% summarize(max_pip=max(pip))


df.cran$PIP_bin <- cut(df.cran$MAX_PIP,breaks = c(0,0.8,1.0),labels = c("0-0.8","0.8-1.0"))
df.gtex$PIP_bin <- cut(df.gtex$max_pip,breaks = c(0,0.8,1.0),labels = c("0-0.8","0.8-1.0"))
df.gtop$PIP_bin <- cut(df.gtop$max_pip,breaks = c(0,0.8,1.0),labels = c("0-0.8","0.8-1.0"))



extract_tissue <- function(x){
  return(strsplit(x,split=":",fixed=T)[[1]][2])
}

df.gtex$Tissue <- sapply(df.gtex$tissue_gene,extract_tissue)
df.gtop$Tissue <- sapply(df.gtop$tissue_gene,extract_tissue)

df.gtex_count <-  as.data.frame.matrix(as.matrix(table(df.gtex$Tissue,df.gtex$PIP_bin)))
df.gtop_count <-  as.data.frame.matrix(as.matrix(table(df.gtop$Tissue,df.gtop$PIP_bin)))
df.cran_count <-  as.data.frame.matrix(as.matrix(table(df.cran$Tissue,df.cran$PIP_bin)))

df.gtex_count$group <- "GTEx"
df.gtop_count$group <- "GTOP"
df.cran_count$group <- "GTEx+GTOP"
names(df.gtex_count) <- c("lowPIP","highPIP","group")
names(df.gtop_count) <- c("lowPIP","highPIP","group")
names(df.cran_count) <- c("lowPIP","highPIP","group")

df.gtex_count %<>% mutate(PRP1=lowPIP/(lowPIP+highPIP),PRP2=highPIP/(lowPIP+highPIP))
df.gtop_count %<>% mutate(PRP1=lowPIP/(lowPIP+highPIP),PRP2=highPIP/(lowPIP+highPIP))
df.cran_count %<>% mutate(PRP1=lowPIP/(lowPIP+highPIP),PRP2=highPIP/(lowPIP+highPIP))

df_plot <- rbind(df.gtex_count[,c(-1,-2)],df.gtop_count[,c(-1,-2)],df.cran_count[,c(-1,-2)])

df_plot.l <- df_plot %>% group_by(group) %>% summarise(meanPRP1=mean(PRP1),meanPRP2=mean(PRP2),sdPRP2=sd(PRP2))

library(reshape2)
df_plot.l2 <- melt(df_plot.l,id.vars = c("group","sdPRP2"))
df_plot.l2$group <- factor(df_plot.l2$group,levels = c("GTOP","GTEx","GTEx+GTOP"))
df_plot.l2$variable <- factor(df_plot.l2$variable,levels = c("meanPRP2","meanPRP1"))
dat_errorbar <- df_plot.l2 %>% filter(variable=="meanPRP1")

p4 <- ggplot(df_plot.l2) + geom_bar(aes(x=group,y=value,fill=variable),stat="identity") + 
  geom_errorbar(data=dat_errorbar,aes(x=group,ymin = value-sdPRP2,ymax = value+sdPRP2),width=.4)+ theme_pubr() + 
  scale_fill_manual(breaks = c("meanPRP2","meanPRP1"),values = c("#913627","#7d8caf"));p4


# Fig.4j ABO --------------------------------------------------------------


df_fm.gtex <- readRDS("./input/ABO.gtex_finemap.RDS")
df_fm.gtop <- readRDS("./input/ABO.gtop_finemap.RDS")
df_fm.susiex <- readRDS("./input/ABO.susiex_finemap.RDS")


susiex.cs <- c("rs2129834802","rs8176719")
df_fm.susiex$PIP <- apply(df_fm.susiex[,c(3,4)],1,max)
df_fm.gtex$in_cs <- 0
df_fm.gtex$in_cs[!is.na(df_fm.gtex$cs)] <- 1

df_fm.gtop$in_cs <- 0
df_fm.gtop$in_cs[!is.na(df_fm.gtop$cs)] <- 1

df_fm.gtex %<>% select(rsid156,chrom,pos,pip,in_cs)
df_fm.gtop %<>% select(variant_id,chrom,pos1,pip,in_cs)

df_fm.susiex$in_cs <- 0
df_fm.susiex$in_cs[df_fm.susiex$SNP %in% susiex.cs] <- 1



df_fm.gtex$in_cs <- factor(df_fm.gtex$in_cs)
df_fm.gtop$in_cs <- factor(df_fm.gtop$in_cs)
df_fm.susiex$in_cs <- factor(df_fm.susiex$in_cs)
p1 <- ggplot(df_fm.gtex,aes(x=pos,y=pip,color=in_cs)) + geom_point(size=2) + theme_pubr() + 
  scale_color_manual(breaks = c(0,1),values = c("grey","red"))+ ylim(0,1);p1
p2 <- ggplot(df_fm.gtop,aes(x=pos1,y=pip,color=in_cs)) + geom_point(size=2) + theme_pubr() + 
  scale_color_manual(breaks = c(0,1),values = c("grey","red")) + ylim(0,1);p2
p3 <- ggplot(df_fm.susiex,aes(x=BP,y=PIP,color=in_cs)) + geom_point(size=2) + theme_pubr() + 
  scale_color_manual(breaks = c(0,1),values = c("grey","red")) + ylim(0,1);p3

#pdf(file="SuSiEX.example.ABO.liver.pdf",height = 6.8,width = 4.5)
cowplot::plot_grid(p1,p2,p3,ncol=1,align = "v")
#dev.off()