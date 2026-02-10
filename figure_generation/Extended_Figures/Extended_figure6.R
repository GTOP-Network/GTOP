#==================================#
# Population-specific QTLs #
# Extended Figure-6 #
#==================================#

setwd("/path/to/GTOP_code/extend/extend_6")

library(data.table)
library(tidyverse)
library(reticulate)

## Extended.Fig.6a, Fine-mapped eQTLs and sQTLs --------------------------

fm_proxy_count <- fread("./input/Extended_fig6a_data.txt")

ggplot(fm_proxy_count, aes(x=cs_num, y=phenotype_count, fill=xQTL_type)) +
    geom_col(position = "dodge") +
    scale_x_continuous(breaks = 1:7, labels = c(1:6,">=7")) +
    theme_classic() +
    labs(x="Number of independent QTLs", y="Gene number") +
    scale_fill_manual(values = c("eQTL"="#9ac294", "sQTL"="#7b86a7"))


## Extended.Fig.6b, Geographic frequencies of eQTL -----------------------

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
geovar_test.add_freq_mat("./input/Extended_fig6b_data.txt")
geovar_test.geovar_binning()

geovar_plot = GeoVarPlot()
geovar_plot.add_data_geovar(geovar_test)
geovar_plot.filter_data()
geovar_plot.add_cmap()

fig, ax = plt.subplots(1,1,figsize=(3,6))
geovar_plot.plot_geovar(ax)
ax.set_xticklabels(geovar_plot.poplist)
plt.savefig("./ExtendFig6b.pdf", dpi=300, bbox_inches=\'tight\')
')


## Extended.Fig.6c, Geographic frequencies of sQTL -----------------------

py_run_string('
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import pkg_resources
from geovar import *

plt.rcParams[\'pdf.fonttype\'] = 42

geovar_test = GeoVar()
geovar_test.add_freq_mat("./input/Extended_fig6c_data.txt")
geovar_test.geovar_binning()

geovar_plot = GeoVarPlot()
geovar_plot.add_data_geovar(geovar_test)
geovar_plot.filter_data()
geovar_plot.add_cmap()

fig, ax = plt.subplots(1,1,figsize=(3,6))
geovar_plot.plot_geovar(ax)
ax.set_xticklabels(geovar_plot.poplist)
plt.savefig("./ExtendedFig6c.pdf", dpi=300, bbox_inches=\'tight\')
')


## Extended.Fig.6d, Compared with GTEx -----------------------------------
GTE_comparison <- fread("./input/Extended_fig6d.txt")

tissue_oder<-c(unique(GTE_comparison$ntissue))
GTE_comparison$ntissue<-factor(GTE_comparison$ntissue,levels=rev(tissue_oder))

ggplot(GTE_comparison, aes(x = count, y=ntissue, fill = Type)) +
    geom_col(position = "fill") +
    # facet_wrap(tissue~., nrow = 1) +
    labs(x = "Proportion of eGenes", y = "") +
    theme_classic() +
    theme(
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12)
    ) +
    scale_fill_manual(values = c("#9d3929", "#fc9272", "#ab889a", "#3578ac","1")) 

## Extende.Fig.6e, Number of tissues sharing fd-QTLs --------------------
EASspecific_genelist <- readRDS("./input/Extended_fig6e_data.RDS")

ggvenn::ggvenn(EASspecific_genelist, fill_color=c("#7d8cad", "#0f3d7b"), show_percentage = F)

ff_tissue_summary <- fread("./input/Extended_fig6e_data.txt")

ff_tissue_summary$sum <- factor(as.character(ff_tissue_summary$sum), 
                                levels = rev(as.character(unique(ff_tissue_summary$sum))))
ff_tissue_summary$new_Type_GTEx <- factor(ff_tissue_summary$new_Type_GTEx, 
                                          levels = c("U", "R"))

ggplot(ff_tissue_summary, aes(x=gene_count, y=sum, fill=new_Type_GTEx)) +
    geom_col(position = "dodge") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black")) +
    labs(x="Number of tissues", y="Gene number", fill="") +
    scale_fill_manual(values = c("#7b89a8", "#0d3b75")) +
    scale_x_continuous(position = "top") 


## Extended.Fig.6f, beta value of he-QTLs --------------------------------
sig_res_df <-  fread("./input/Extended_fig6f_data.txt")

ggplot(data = sig_res_df, aes(x = GTOP_beta, y = GTEx_beta)) + 
    geom_point(color="#0c0d0c", size=0.3) +
    ggpointdensity::geom_pointdensity(size=2, adjust=0.5, shape=19 )+
    scale_color_viridis_c(option = "D", direction = 1) +
    ggdensity::geom_hdr_lines(linetype="solid", linewidth=0.5, 
                              probs = c(0.9, 0.7, 0.5), color="#ff0000") + 
    geom_abline(intercept = 0, slope = 1, color="#282c33", 
                size=0.5, linetype="dashed") +
    geom_abline(intercept = 0, slope = -1, color="#282c33", 
                size=0.5, linetype="dashed")+
    theme_classic(base_size = 12) +
    geom_hline(yintercept = 0, color = "black") +
    labs(x="Effect size (GTOP)", y="Effect size (GTEx)") +
    ggh4x::coord_axes_inside(labels_inside = T)


## Extended.Fig.6g, Example of he-QTLs --------------------------------
GTOP_GTEx_data <- readRDS("./input/Extended_fig6g_data.RDS")
GTOP_GTEx_data$source <- factor(GTOP_GTEx_data$source, levels = c("GTOP", "GTEx"))

ggplot(GTOP_GTEx_data, aes(genotype_label, scale_counts, fill = source)) +
    geom_boxplot(outlier.colour = NA, width = 0.8) +
    ggbeeswarm::geom_quasirandom(
        aes(
            x = genotype_label,
            col = source,
            group = interaction(genotype_label, source)
        ),
        size = 2,
        dodge.width = 0.8,
        width = 0.2
    ) +
    theme_classic(base_size = 12) +
    scale_fill_manual(values = c("white", "white")) +
    scale_color_manual(values = c("#bd5a45", "#e7c798")) +
    theme(
        legend.position = 'top',
        plot.title = element_text(size = 10),
        axis.text = element_text(colour = "black")
    ) 
