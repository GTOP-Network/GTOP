#==================================#
# Population-specific QTLs #
# Figure-4 #
#==================================#

setwd("/media/london_A/mengxin/GTOP_code/fig-4")

library(ggrastr)
library(data.table)
library(tidyverse)
library(patchwork)
library(reticulate)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)


# Fig.4a cs size compare----------------------------------------
library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)

dat.gtex <- readRDS("./input/dat.gtex.RDS")
dat.gtop <- readRDS("./input/dat.gtop.RDS")
dat.cran <- readRDS("./input/dat.cran.RDS")
# examine CS number per gene
df.gtex <- dat.gtex %>% dplyr::select(locus_id,cs,Tissue) %>% distinct(.keep_all = T)
df.gtop <- dat.gtop %>% dplyr::select(locus_id,cs,Tissue) %>% distinct(.keep_all = T)
df.cran <- dat.cran %>% dplyr::select(Gene,CS_ID,Tissue) %>% distinct(.keep_all = T)


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

df.cran <- dat.cran %>% dplyr::filter(tissue_gene %in% overlap_tissuegenes) %>% dplyr::select(Gene,CS_ID,MAX_PIP,CS_LENGTH,Tissue,tissue_gene)
df.gtex <- dat.gtex %>% dplyr::filter(tissue_gene %in% overlap_tissuegenes) %>% dplyr::select(locus_id,cs,pip,cs_size,Tissue,tissue_gene) %>% 
  group_by(locus_id,cs,Tissue) %>% mutate(max_pip = max(pip)) %>% ungroup() %>% dplyr::filter(pip==max_pip) %>% 
  dplyr::select(locus_id,cs,max_pip,cs_size,Tissue,tissue_gene) %>% distinct(.keep_all = T)
df.gtop <- dat.gtop %>% dplyr::filter(tissue_gene %in% overlap_tissuegenes) %>% dplyr::select(locus_id,cs,pip,cs_size,Tissue,tissue_gene) %>% 
  group_by(locus_id,cs,Tissue) %>% mutate(max_pip = max(pip)) %>% ungroup() %>% dplyr::filter(pip==max_pip) %>% 
  dplyr::select(locus_id,cs,max_pip,cs_size,Tissue,tissue_gene) %>% distinct(.keep_all = T)

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


x <- df.cran %>% mutate(maxPIP_cran=max_pip,CS_size_cran=cs_size) %>% dplyr::select(tissue_gene,maxPIP_cran,CS_size_cran)
y <- df.gtex %>% mutate(maxPIP_gtex=max_pip,CS_size_gtex=cs_size) %>% dplyr::select(tissue_gene,maxPIP_gtex,CS_size_gtex) 
z <- df.gtop %>% mutate(maxPIP_gtop=max_pip,CS_size_gtop=cs_size) %>% dplyr::select(tissue_gene,maxPIP_gtop,CS_size_gtop) 
xy <- merge(x,y,by="tissue_gene")
xyz <- merge(xy,z,by="tissue_gene")

test1 <- wilcox.test(xyz$CS_size_cran,xyz$CS_size_gtex,alternative = "less",paired = T) # p<2.2e-16
test2 <- wilcox.test(xyz$CS_size_cran,xyz$CS_size_gtop,alternative = "less",paired = T) # p<2.2e-16


# Fig.4b plot prop of PIP>0.8 in single population fine-mapping and cross-ancestry fine-mapping------

overlap_tissuegenes <- intersect(dat.cran$tissue_gene,intersect(dat.gtex$tissue_gene,dat.gtop$tissue_gene))
df.cran <- dat.cran %>% dplyr::filter(tissue_gene %in% overlap_tissuegenes)
df.gtex <- dat.gtex %>% dplyr::filter(tissue_gene %in% overlap_tissuegenes) %>% group_by(tissue_gene,cs) %>% summarize(max_pip=max(pip))
df.gtop <- dat.gtop %>% dplyr::filter(tissue_gene %in% overlap_tissuegenes) %>% group_by(tissue_gene,cs) %>% summarize(max_pip=max(pip))


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
df_plot.l2 <- reshape2::melt(df_plot.l,id.vars = c("group","sdPRP2"))
df_plot.l2$group <- factor(df_plot.l2$group,levels = c("GTOP","GTEx","GTEx+GTOP"))
df_plot.l2$variable <- factor(df_plot.l2$variable,levels = c("meanPRP2","meanPRP1"))
dat_errorbar <- df_plot.l2 %>% dplyr::filter(variable=="meanPRP1")

p4 <- ggplot(df_plot.l2) + geom_bar(aes(x=group,y=value,fill=variable),stat="identity") + 
  geom_errorbar(data=dat_errorbar,aes(x=group,ymin = value-sdPRP2,ymax = value+sdPRP2),width=.4)+ theme_pubr() + 
  scale_fill_manual(breaks = c("meanPRP2","meanPRP1"),values = c("#913627","#7d8caf"));p4

x <- df.gtex_count$PRP2
y <- df.gtop_count$PRP2
z <- df.cran_count$PRP2
t.test(x,z,paired = T,alternative = "less") # GTOP+GTEx VS GTEX: p=3.942e-09
t.test(y,z,paired = T,alternative = "less") # GTOP+GTEx VS GTOP: p=7.768e-07


## Fig.4c, PIP of SAMD9L -------------------------------------

pipdf <- fread("./input/Fig4c.txt")

pipdf <- pipdf %>%
  mutate(
    cs_label = as.logical(cs_label),
    study = factor(study, levels = rev(c("GTEx", "GTOP", "Cross")))
  )

ggplot(pipdf, aes(x = BP, y = PIP, color = cs_label)) +
  geom_point(size = 1.5, alpha = 0.75) +
  facet_wrap(~ study, ncol = 1) +
  scale_color_manual(values = c(`TRUE` = "red3", `FALSE` = "grey70"), labels = c(`TRUE` = "CS", `FALSE` = "non-CS")) +
  theme_pubr() +
  labs(x = "Position", y = "PIP", color = NULL) +
  theme(strip.text = element_text(face = "bold"), legend.position = "top")


## Fig.4d, fd-QTL genes ----------------------------------------

library(reticulate)

py_require("matplotlib==3.6.3")
py_require("git+https://github.com/aabiddanda/geovar")

py_config()

py_run_string(
  '
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

if not hasattr(np, "row_stack"):
    np.row_stack = np.vstack

if not hasattr(cm, "get_cmap"):
    def get_cmap(name=None, lut=None):
        cmap = plt.colormaps[name if name is not None else "viridis"]
        if lut is not None:
            cmap = cmap.resampled(lut)
        return cmap
    cm.get_cmap = get_cmap

from geovar import *

plt.rcParams["pdf.fonttype"] = 42

geovar_test = GeoVar()
geovar_test.add_freq_mat("input/Fig4d.txt")
geovar_test.geovar_binning()

geovar_plot = GeoVarPlot()
geovar_plot.add_data_geovar(geovar_test)
geovar_plot.filter_data()
geovar_plot.add_cmap()

fig, ax = plt.subplots(1, 1, figsize=(3, 6))
geovar_plot.plot_geovar(ax)

ax.set_xticklabels(geovar_plot.poplist)

plt.savefig(
    "input/freq_data.pdf",
    dpi=300,
    bbox_inches="tight"
)

plt.close()
'
)

## Fig.4e, fd-eQTL gene number -------------------------------------------------
tissue_gene_summary <- fread("input/Fig4e.txt")
color_vec <- readRDS("./input/tissue_color.RDS")
tissue_gene_summary$variable <- factor(tissue_gene_summary$variable, 
                                       levels = c("AFR", "EUR", "SAS", "AMR"))

ggplot(
  tissue_gene_summary,
  aes(x = variable, y = count)
) +
  geom_boxplot(outlier.color = "NA", fill = "#e5e5e4") +
  geom_point(
    aes(group = tissue, color = tissue),
    position = position_dodge(width = 0.4)
  ) +
  scale_color_manual(values = color_vec) +
  theme_classic() +
  labs(x = "", y = "Gene number") +
  scale_fill_manual(values = c("#6874b4", "#4bb9b9"))


## Fig.4f, fd-eQTL example -----------------------------------------------------
library(EnsDb.Hsapiens.v86)

load("./input/Fig4e.RData")

SNP_name <- "rs74971849"
loci_start <- SNV_eQTL_data$pos[which.min(SNV_eQTL_data$p)] - 250000
loci_end <- SNV_eQTL_data$pos[which.min(SNV_eQTL_data$p)] + 250000
GWAS_name <- "Total bilirubin"
t_xqtl_type <- "eQTL"
gene_ensg <- "ENSG00000205754.13"
t_tissue <- "Liver"
gene_symbol <- "SLCO1B7"
chr_name <- SNV_eQTL_data$chrom[1]

## GWAS locuszoom
GWAS_locus <- locuszoomr::locus(
  data = as.data.frame(GWAS_data),
  xrange = c(loci_start, loci_end),
  seqname = chr_name,
  index_snp = SNP_name,
  ens_db = "EnsDb.Hsapiens.v86"
)

GWAS_locus$data$ld <- LD_info[GWAS_locus$data$rsid, "R2"]
GWAS_locus$data$col <- "transparent"

GWAS_plot <- gg_scatter(
  GWAS_locus,
  pcutoff = FALSE,
  yzero = T,
  size = 2,
  ylim = c(0, -log10(2.130846e-26)),
  labels = SNP_name,
  color = "black",
  legend_pos = "right",
  LD_scheme = c(
    "#e5e5e5",
    "#e5e5e5",
    "#3e70b4",
    "#3f7d1d",
    "orange",
    "red",
    "red"
  )
) +
  annotate(
    "text",
    x = loci_end / 10^6,
    y = max(-log10(GWAS_data$p)),
    label = GWAS_name,
    hjust = "right"
  ) +
  guides(
    fill = guide_legend(
      reverse = TRUE,
      override.aes = list(
        colour = "transparent",
        stroke = 0
      )
    )
  ) +
  labs(x = "") +
  theme(axis.text.x = element_blank())

## eQTL locuszoom
SNV_eQTL_locus <- locus(
  data = as.data.frame(SNV_eQTL_data),
  xrange = c(loci_start, loci_end),
  seqname = chr_name,
  index_snp = SNP_name,
  ens_db = "EnsDb.Hsapiens.v86"
)

SNV_eQTL_locus$data$ld <- LD_info[SNV_eQTL_locus$data$rsid, "R2"]
SNV_eQTL_locus$data$col <- "transparent"

SNV_eQTL_plot <- gg_scatter(
  SNV_eQTL_locus,
  pcutoff = FALSE,
  size = 2,
  yzero = T,
  labels = SNP_name,
  color = "black",
  legend_pos = "right",
  LD_scheme = c(
    "#e5e5e5",
    "#e5e5e5",
    "#3e70b4",
    "#3f7d1d",
    "orange",
    "red",
    "red"
  )
) +
  annotate(
    "text",
    x = loci_end / 10^6,
    y = max(-log10(SNV_eQTL_data$p)),
    label = sprintf("%s | %s", t_xqtl_type, t_tissue),
    hjust = "right"
  ) +
  guides(
    fill = guide_legend(
      reverse = TRUE,
      override.aes = list(
        colour = "transparent",
        stroke = 0
      )
    )
  )

## gene structure
gene_plot <- gg_genetracks(
  SNV_eQTL_locus,
  highlight = gene_symbol,
  highlight_col = "#3055a3"
  # filter_gene_name = c(gene_symbol)
) +
  labs(x = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )

## joint-fine-mapping
fm_locus <- locus(
  data = as.data.frame(fm_data),
  yvar = "PIP",
  xrange = c(loci_start, loci_end),
  seqname = chr_name, #index_snp = SNP_name,
  ens_db = "EnsDb.Hsapiens.v86"
)
fm_locus$data$ld <- LD_info[fm_locus$data$rsid, "R2"]
fm_locus$data$col <- "transparent"

fm_plot <- gg_scatter(
  fm_locus,
  pcutoff = FALSE,
  yzero = T,
  labels = SNP_name,
  color = "black",
  legend_pos = "right",
  LD_scheme = c(
    "#e5e5e5",
    "#e5e5e5",
    "#3e70b4",
    "#3f7d1d",
    "orange",
    "red",
    "red"
  )
) +
  annotate(
    "text",
    x = loci_end / 10^6,
    y = max(fm_data$PIP),
    label = "fine-mapping",
    hjust = "right"
  ) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())


wrap_plots(
  list(
    GWAS_plot,
    SNV_eQTL_plot,
    fm_plot,
    gene_plot
  ),
  ncol = 1,
  heights = c(3, 3, 3, 2)
)


## Fig.4g eQTL portability -----------------------------------------------------
library(ggpubr)
summary_df1 <- fread("./input/Fig4g.txt")

summary_df1$type2 <- factor(
  summary_df1$type2,
  levels = c(
    "FDR",
    "ratio_nominal",
    "ratio_mash",
    "ratio_uncertainty-aware",
    "gene"
  )
)

p3 <- ggplot(
  summary_df1[
    summary_df1$type3 == "Lead eVariant",
  ],
  aes(type2, ratio, fill = type1)
) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.7),
    width = 0.6
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.7),
    width = 0.2,
    linewidth = 0.6
  ) +
  lims(y = c(0, 1)) +
  theme_pubr() +
  scale_fill_manual(
    values = c("#e9bd76", "#3e71a6")
  ) +
  labs(x = "", y = "Number of eVariants")


p4 <- ggplot(
  summary_df1[
    summary_df1$type3 != "Lead eVariant",
  ],
  aes(type3, ratio, fill = type1)
) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.7),
    width = 0.6
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.7),
    width = 0.2,
    linewidth = 0.6
  ) +
  lims(y = c(0, 1)) +
  theme_pubr() +
  scale_fill_manual(
    values = c("#e9bd76", "#3e71a6")
  ) +
  labs(x = "", y = "Number of eGenes")

p3 +
  p4 +
  plot_layout(
    guides = "collect",
    widths = c(2, 1)
  )

