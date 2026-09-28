#!/usr/bin/env Rscript

# ==============================================================================
# GTOP vs. gnomAD / 1KGP Allele Frequency Comparison
# ==============================================================================

#%% ------------------------ 0. prepare files (packages, input files, output files)
PROJECT_DIR <- "path/to/dir"
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(ggrastr)
})

# ------------ input
IN_GTOP_FREQ_FILE <- "/path/to/data/genotype/gtop/gtop_snv.afreq"
IN_GNOMAD_FREQ_FILE <- "/path/to/data/genotype/gnomad/gtop_intersect_gnomad_info.txt"
IN_CHRPOS_RSID_FILE <- "/path/to/data/genotype/gtop/gtop_chrpos_rsid.txt"

# ------------ output
OUT_DATA_FILE <- "output/data/integrated_freq.txt"
OUT_AF_PLOT_FILE <- "output/figures/AF_correlation_GTOP_gnomAD.pdf"
lapply(
    unique(dirname(c(OUT_DATA_FILE, OUT_AF_PLOT_FILE))),
    function(x) {
        if (!dir.exists(x)) {
            dir.create(recursive = TRUE)
        }
    }
)

# ------------ function
theme_pub <- function(base_size = 12) {
    theme_classic(base_size = base_size) +
        theme(
            axis.text = element_text(color = "black", size = base_size),
            axis.title = element_text(color = "black", size = base_size),
            axis.line = element_line(
                colour = "black",
                linewidth = 0.5,
                lineend = "round"
            ),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            axis.ticks.length = unit(0.15, "cm"),
            plot.title = element_text(face = "bold", hjust = 0),
            plot.subtitle = element_text(hjust = 0),
            legend.text = element_text(color = "black", size = base_size)
        )
}


#%% ------------------------ 1. Load GTOP data (GTOP variant, allele frequency)
gtop_freq <- fread(IN_GTOP_FREQ_FILE)
chrpos2rsid <- fread(
    IN_CHRPOS_RSID_FILE,
    col.names = c("chr_pos_ref_alt", "variant_id")
)

## Rename columns and attach chr:pos → rsID mapping
gtop_freq <- merge(
    gtop_freq[, .(variant_id = ID, af_gtop = ALT_FREQS)],
    chrpos2rsid,
    by = "variant_id"
)[, .(chr_pos_ref_alt, variant_id, af_gtop)]

## Keep common variants only (MAF 5%–95%)
gtop_freq <- gtop_freq[af_gtop >= 0.05 & af_gtop <= 0.95]


#%% ------------------------ 2. Load and Filter gnomAD (minimum allele number)
gnomad_freq <- fread(IN_GNOMAD_FREQ_FILE)

gnomad_freq$AC_eur <- gnomad_freq$AC_nfe + gnomad_freq$AC_fin
gnomad_freq$AN_eur <- gnomad_freq$AN_nfe + gnomad_freq$AN_fin
gnomad_freq$AF_eur <- gnomad_freq$AC_eur / gnomad_freq$AN_eur

MIN_AN <- 2000L

gnomad_freq <- gnomad_freq[
    AN_afr > MIN_AN &
        AN_amr > MIN_AN &
        AN_eas > MIN_AN &
        AN_eur > MIN_AN &
        AN_sas > MIN_AN
]


#%% ------------------------ 3. Integrate GTOP and gnomAD (five ancestry groups)
integrated_freq <- merge(
    gnomad_freq,
    gtop_freq[, .(ID = variant_id, chr_pos_ref_alt, af_gtop)],
    by = "ID"
)

## Cast AF columns to numeric (raw gnomAD export may store these as character)
integrated_freq[, `:=`(
    af_afr = as.numeric(AF_afr),
    af_amr = as.numeric(AF_amr),
    af_eas = as.numeric(AF_eas),
    af_eur = as.numeric(AF_eur),
    af_sas = as.numeric(AF_sas)
)]

integrated_freq$af_neas <- (integrated_freq$AC_afr +
    integrated_freq$AC_amr +
    integrated_freq$AC_eur +
    integrated_freq$AC_sas) /
    (integrated_freq$AN_afr +
        integrated_freq$AN_amr +
        integrated_freq$AN_eur +
        integrated_freq$AN_sas)

integrated_freq$type_afr <- case_when(
    integrated_freq$af_afr == 0 | integrated_freq$af_afr == 1 ~ "U",
    integrated_freq$af_afr < 0.01 | integrated_freq$af_afr > 0.99 ~ "R",
    TRUE ~ "C"
)
integrated_freq$type_amr <- case_when(
    integrated_freq$af_amr == 0 | integrated_freq$af_amr == 1 ~ "U",
    integrated_freq$af_amr < 0.01 | integrated_freq$af_amr > 0.99 ~ "R",
    TRUE ~ "C"
)
integrated_freq$type_eur <- case_when(
    integrated_freq$af_eur == 0 | integrated_freq$af_eur == 1 ~ "U",
    integrated_freq$af_eur < 0.01 | integrated_freq$af_eur > 0.99 ~ "R",
    TRUE ~ "C"
)
integrated_freq$type_sas <- case_when(
    integrated_freq$af_sas == 0 | integrated_freq$af_sas == 1 ~ "U",
    integrated_freq$af_sas < 0.01 | integrated_freq$af_sas > 0.99 ~ "R",
    TRUE ~ "C"
)
integrated_freq$type_neas <- case_when(
    integrated_freq$af_neas == 0 | integrated_freq$af_neas == 1 ~ "U",
    integrated_freq$af_neas < 0.01 | integrated_freq$af_neas > 0.99 ~ "R",
    TRUE ~ "C"
)

integrated_freq <- integrated_freq[, .(
    variant_id = ID,
    chr_pos_ref_alt,
    af_gtop,
    af_afr,
    af_amr,
    af_eas,
    af_eur,
    af_sas,
    af_neas,
    type_afr,
    type_amr,
    type_eur,
    type_sas,
    type_neas
)]

nrow(gtop_freq) # 7740074 # 7370055
nrow(integrated_freq) # 7032964 (90.9%) # n = 6,715,823 (91.12%)
table(integrated_freq$af_eas >= 0.05 & integrated_freq$af_eas <= 0.95) # 6701993

integrated_freq <- integrated_freq[af_eas >= 0.05 & af_eas <= 0.95]


#%% ------------------------ 4. Pearson correlation (GTOP vs. gnomAD EAS)
r_pearson <- cor(
    integrated_freq$af_gtop,
    integrated_freq$af_eas,
    method = "pearson"
)
cat(sprintf("Pearson r (GTOP vs. gnomAD EAS): %.4f\n", r_pearson))

ggplot(integrated_freq, aes(x = af_eas, y = af_gtop)) +
    rasterise(geom_hex(bins = 250), dpi = 600) +
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.6) +
    scale_fill_viridis_c(trans = "log10") +
    coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
        x = "AF gnomAD EAS",
        y = "AF GTOP",
        fill = "Variant count",
        subtitle = sprintf("Pearson r = %.4f", r_pearson)
    ) +
    theme_pub()

ggsave(OUT_AF_PLOT_FILE, width = 6, height = 5)


#%% ------------------------ 5. Save integrated frequency table
fwrite(integrated_freq, OUT_DATA_FILE, sep = "\t")
