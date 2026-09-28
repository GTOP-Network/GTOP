#==============================================#
# GWAS-fdQTL #
# Supp-Figure-42#
#==============================================#

library(ggplot2)
library(tidyverse)
library(ggpubr)

setwd("/path/to/GTOP_code/supp/supp_fig42/")

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
plot_df <- fread("./input/supp_fig42b_data.txt")

SNP_name <- "rs4646776"
lead_df <- plot_df %>%
  dplyr::filter(variant_id == SNP_name)

plot_df$ld_bin <- factor(plot_df$ld_bin, levels = c(
  "No LD information",
  "< 0.2",
  "0.2–0.4",
  "0.4–0.6",
  "0.6–0.8",
  "≥ 0.8",
  "Lead variant"
))

ggplot(
  plot_df,
  aes(
    x = af_other * 100,
    y = af_eas * 100
  )
) +
  ## 1. Plot non-lead variants first
  geom_point(
    data = subset(plot_df, variant_id != SNP_name),
    aes(
      size = -log10(p),
      fill = ld_bin
    ),
    shape = 21,
    colour = "grey30",
    stroke = 0.2,
    alpha = 1
  ) +
  ## 2. Plot lead variant
  geom_point(
    data = lead_df,
    aes(
      x = af_other * 100,
      y = af_eas * 100,
      size = -log10(p),
      fill = ld_bin
    ),
    shape = 21,
    colour = "black",
    stroke = 0.2,
    inherit.aes = FALSE
  ) +
  ## Label lead variant
  ggrepel::geom_text_repel(
    data = lead_df,
    aes(
      x = af_other * 100,
      y = af_eas * 100,
      label = variant_id
    ),
    inherit.aes = FALSE,
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0
  ) +
  scale_fill_manual(
    name = expression(LD ~ (r^2)),
    values = c(
      "Lead variant" = "purple",
      "≥ 0.8" = "red",
      "0.6–0.8" = "orange",
      "0.4–0.6" = "#3f7d1d",
      "0.2–0.4" = "#3e70b4",
      "< 0.2" = "#e5e5e5"
    ),
    breaks = c(
      "Lead variant",
      "≥ 0.8",
      "0.6–0.8",
      "0.4–0.6",
      "0.2–0.4",
      "< 0.2"
    ),
    na.value = "#e5e5e5",
    drop = FALSE
  ) +
  scale_size_continuous(
    name = expression(-log[10](italic(P)))
  ) +
  labs(
    x = "Aallele frequency (Other populations)",
    y = "Allele frequency (EAS)"
  ) +
  facet_grid(. ~ type) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "right"
  ) +
  coord_cartesian(
    xlim = c(0, 100),
    ylim = c(0, 100)
  ) +
  ggh4x::coord_axes_inside(labels_inside = FALSE)
