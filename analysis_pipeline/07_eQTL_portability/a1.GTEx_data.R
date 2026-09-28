#!/usr/bin/env Rscript

# ==============================================================================
# GTEx 8 tissues （QTL and gene information）
# ==============================================================================

#%% -------------------- 0. prepare files (packages, input, output path)
PROJECT_DIR <- "/path/to/dir"
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
})

# --------- input
GTOP_DIR <- "/path/to/gtop/snv_eqtl/07_summary"
GTEX_V8_DIR <- "/path/to/src/data/GTEx/v8/eQTL/GTEx_Analysis_v8_eQTL"

# --------- output
FIG_DIR <- file.path(PROJECT_DIR, "output/figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)


#%% -------------------- 1. Helper functions
strip_gene_version <- function(x) {
    sub("\\..*$", "", as.character(x))
}

clean_genes <- function(x) {
    unique(na.omit(strip_gene_version(x)))
}

safe_pct <- function(n, d) {
    if (length(d) == 0 || is.na(d) || d == 0) {
        return(NA_real_)
    }
    100 * n / d
}

read_gtex_v8 <- function(gtex_tissue) {
    file <- file.path(
        GTEX_V8_DIR,
        sprintf("%s.v8.egenes.txt.gz", gtex_tissue)
    )

    dt <- fread(file)
    dt[, gene_id_clean := strip_gene_version(gene_id)]
    dt
}

list_for_upset <- function(data_list) {
    all_genes <- unique(unlist(data_list))

    upset_input <- map_dfc(data_list, ~ all_genes %in% .x) %>%
        mutate(gene_id = all_genes, .before = 1)

    upset_input
}


#%% -------------------- 2. Load metadata
tissue_raw <- fread("input/gtop_gtex_tissues.txt")

tissue_link <- data.table(
    GTOP_tissue = tissue_raw[["GTOP_tissue"]],
    GTEx_tissue = tissue_raw[["GTEx_tissue"]]
)

gtop_samplesize <- fread(
    "/path/to/gtop/snv_eqtl/07_summary/sample_size.txt",
)

tissue_link <- tissue_link[
    order(factor(
        GTOP_tissue,
        levels = intersect(gtop_samplesize$tissue, GTOP_tissue)
    ))
]
fwrite(tissue_link, "input/gtop_gtex_tissues.txt", sep = "\t")


#%% -------------------- 3. GTOP GTEx overlap
gtex_samplesize <- fread(
    "path/to/gtexv8_eqtl/07_summary/sample_size.txt"
)

lapply(2:nrow(tissue_link), function(row_index) {
    gtop_tissue_name <- tissue_link$GTOP_tissue[row_index]
    gtex_tissue_name <- tissue_link$GTEx_tissue[row_index]

    gtop_eqtl <- fread(sprintf(
        "/path/to/gtop/snv_eqtl/06_nominal_slim/split/%s.xgene_allpairs.txt.gz",
        gtop_tissue_name
    ))
    gtop_eqtl$stable_ensg <- gsub("\\..+", "", gtop_eqtl$ensg)

    gtex_allpairs <- fread(sprintf(
        "/path/to/gtexv8_eqtl/05_nominal/%s.txt.gz",
        gtex_tissue_name
    ))
    gtex_allpairs$stable_ensg <- gsub("\\..+", "", gtex_allpairs$gene_id)
    gtex_allpairs$chr_pos_ref_alt <- gsub("_b38", "", gtex_allpairs$variant_id)

    gtex_gtop_overlap <- merge(
        gtex_allpairs[, .(
            stable_ensg,
            chr_pos_ref_alt,
            gtex_beta = slope,
            gtex_se = slope_se,
            gtex_maf = maf,
            gtex_p_value = pval_nominal,
            gtex_n = gtex_samplesize$N[
                gtex_samplesize$tissue == gtex_tissue_name
            ]
        )],
        gtop_eqtl[, .(
            stable_ensg,
            chr_pos_ref_alt,
            gtop_beta = beta,
            gtop_se = se,
            gtop_maf = 0.5 - abs(0.5 - af),
            gtop_p_value = p_value,
            gtop_n = N
        )],
        by = c("stable_ensg", "chr_pos_ref_alt")
    )

    fwrite(
        gtex_gtop_overlap,
        sprintf(
            "/path/to/gtexv8_eqtl/08_gtop_overlap/%s.txt.gz",
            gtop_tissue_name
        ),
        sep = "\t"
    )
})
