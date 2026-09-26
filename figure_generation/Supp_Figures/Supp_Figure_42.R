#==============================================#
# GWAS-fdQTL #
# Supp-Figure-42#
#==============================================#


library(ggplot2)
library(tidyverse)
library(ggpubr)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig42.Figure S28-BRAP")

# --------- Helper functions
theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      ext = element_text(colour = "black"),
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(color = "black", size = base_size),
      axis.line = element_line(
        colour = "black", # Line color
        linewidth = 0.5, # Line thickness
        lineend = "round" # Rounded ends
      ),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.10, "cm"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0),
      legend.title = element_blank(),
      legend.text = element_text(color = "black", size = base_size),
      plot.margin = margin(5, 8, 5, 5)
    )
}


# Supp.Fig.42a Example Disease associated-fd-QTL -------------------------------

## Figure a
load("./input/supp_fig42a_data.RData")

SNP_name <- "rs4646776"
loci_start <- SNV_eQTL_data$pos[which.min(SNV_eQTL_data$p)] - 1000000
loci_end <- SNV_eQTL_data$pos[which.min(SNV_eQTL_data$p)] + 1000000
GWAS_name <- "Gout"
t_xqtl_type <- "eQTL"
t_tissue <- "Adipose"
gene_symbol <- "BRAP"
chr_name <- SNV_eQTL_data$chrom[1]

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
  color = "black",
  labels = SNP_name,
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
  theme_pub() +
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
  ) +
  theme_pub() +
  labs(x = "") +
  theme(axis.text.x = element_blank())

## gene structure
gene_plot <- gg_genetracks(
  SNV_eQTL_locus,
  highlight = gene_symbol,
  highlight_col = "#3055a3",
  filter_gene_name = c(gene_symbol),
  filter_gene_biotype = c("protein_coding")
)

## plot
wrap_plots(list(GWAS_plot, SNV_eQTL_plot, gene_plot), ncol = 1, heights = c(3,3,1))



# Supp.Fig.42b ------------------------------------------------------------
tmp_df <- fread("./input/supp_fig42b_data.txt")

ggplot(tmp_df[tmp_df$p<0.01,], aes(af_other*100, af_GTOP*100)) +
  ggpointdensity::geom_pointdensity( alpha=1, aes(size=-log10(p)), adjust=1, shape=21 )+
  scale_color_viridis_b()+
  ggdensity::geom_hdr_lines( linetype="dashed", linewidth=0.5)+
  labs(x="Alternative allele frequency (Other populations)", 
       y="Alternative allele frequency (EAS)")+
  facet_grid(.~type) +
  theme_classic()+
  theme(
    axis.line = element_line(color="black"),
    axis.ticks = element_line(color="black"),
    axis.text = element_text(color="black"),
    legend.position = "right",
  ) + ylim(c(0,100)) +
  ggh4x::coord_axes_inside(labels_inside = F)
