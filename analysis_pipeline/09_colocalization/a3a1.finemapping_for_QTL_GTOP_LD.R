#!/usr/bin/env Rscript

# Load required libraries
library(ieugwasr)
library(data.table)
library(coloc)
library(susieR)

argvs <- commandArgs(trailingOnly = TRUE)

# Check number of arguments
if (length(argvs) != 3) {
    stop(
        "Error: Exactly 3 command-line arguments are required: GWASNAME CHRNAME SENTINELPOS SENTINELSNP"
    )
}

QTLTYPE <- tolower(argvs[1])
TISSUENAME <- argvs[2]
GENENAME <- argvs[3]

# ------------------------------------------------------------
# Read xQTL summary statistics
# ------------------------------------------------------------
message(QTLTYPE, "\t", TISSUENAME, "\t", GENENAME)

QTL_data <- fread(sprintf(
    "/path/to/gtop/%s/06_nominal_slim/split/%s/%s.txt.gz",
    QTLTYPE,
    TISSUENAME,
    GENENAME
))
QTL_data$varbeta <- QTL_data$se * QTL_data$se

QTL_data$variant_id <- gsub("_b38", "", QTL_data$variant_id)
CHRNAME <- strsplit(QTL_data$chr_pos_ref_alt[1], "_")[[1]][1]
setDF(QTL_data)
QTL_data$variant_id[QTL_data$variant_id == ""] <- QTL_data$chr_pos_ref_alt[
    QTL_data$variant_id == ""
]
rownames(QTL_data) <- QTL_data$variant_id

# ------------------------------------------------------------
# Compute LD matrix using 1000G EAS reference panel
# ------------------------------------------------------------
# bfile_path <- sprintf("/path/to/src/data/1000G/five_ancestry_groups/EAS/splitbychr/1000G.EAS.maf01.%s", CHRNAME)
# plink_bin  <- "/path/to/src/bin/plink"

bfile_path <- sprintf(
    "/path/to/src/GTOP_hg38/splitbychr/gtop_snv.maf05.%s",
    CHRNAME
)
plink_bin <- "/lustre/home/tzhang/src/bin/plink"

# Verify that PLINK binary exists
if (!file.exists(plink_bin)) {
    stop("Error: PLINK binary not found at specified path.")
}

# Verify that PLINK bfile prefix exists (check .bed/.bim/.fam)
if (!all(file.exists(paste0(bfile_path, c(".bed", ".bim", ".fam"))))) {
    stop(sprintf(
        "Error: PLINK reference files missing for chromosome: %s",
        bfile_path
    ))
}

# Compute LD matrix
QTL_EAS_ld <- tryCatch(
    {
        ld_matrix(
            variants = QTL_data$variant_id,
            bfile = bfile_path,
            plink_bin = plink_bin,
            with_alleles = FALSE
        )
    },
    error = function(e) {
        stop("Error during LD matrix computation: ", conditionMessage(e))
    }
)

# Ensure LD matrix matches available SNPs
common_snps <- intersect(rownames(QTL_EAS_ld), rownames(QTL_data))
if (length(common_snps) == 0) {
    stop("Error: No overlapping SNPs between GWAS data and LD reference panel.")
}

# Subset GWAS data and LD matrix to common SNPs only
QTL_data_sub <- QTL_data[common_snps, , drop = FALSE]
QTL_EAS_ld_sub <- QTL_EAS_ld[common_snps, common_snps, drop = FALSE]
QTL_EAS_ld_sub[is.na(QTL_EAS_ld_sub)] <- 0

# ------------------------------------------------------------
# Prepare dataset for coloc fine-mapping (SuSiE)
# ------------------------------------------------------------
d2 <- list(
    beta = QTL_data_sub$beta,
    varbeta = QTL_data_sub$varbeta,
    pval = QTL_data_sub$p_value,
    snp = common_snps,
    LD = as.matrix(QTL_EAS_ld_sub),
    N = QTL_data_sub$N[1],
    type = "quant",
    sdY = 1
)

# ------------------------------------------------------------
# Run SuSiE fine-mapping with comprehensive error handling
# ------------------------------------------------------------
fm_prefix <- sprintf(
    "input/finemapping_QTL_GTOP/%s/%s/%s",
    QTLTYPE,
    TISSUENAME,
    GENENAME
)
kriging_prefix <- sprintf(
    "input/finemapping_QTL_GTOP/%s/%s/kriging_rss/%s",
    QTLTYPE,
    TISSUENAME,
    GENENAME
)
if (!dir.exists(dirname(kriging_prefix))) {
    dir.create(dirname(kriging_prefix), recursive = TRUE, showWarnings = FALSE)
}

s2 <- tryCatch(
    {
        runsusie(d2, maxit = 1000, repeat_until_convergence = FALSE)
    },
    error = function(e) {
        message("SuSiE failed to converge: ", conditionMessage(e))
        NULL
    }
)

condz_in <- kriging_rss(d2$beta / sqrt(d2$varbeta), d2$LD, n = d2$N)
condz_in_df <- condz_in$conditional_dist
rownames(condz_in_df) <- d2$snp
fwrite(
    condz_in_df,
    file = sprintf("%s.txt.gz", kriging_prefix),
    sep = "\t",
    row.names = TRUE
)

d2$LD <- NULL

# Save result based on outcome
if (is.null(s2)) {
    # Save input data when SuSiE fails entirely
    save(d2, file = paste0(fm_prefix, ".noConverged.RData"))
} else {
    susie_summary <- tryCatch(summary(s2), error = function(e) NULL)
    if (is.null(susie_summary) || is.null(susie_summary$cs)) {
        # Save input data when no credible sets are returned
        save(d2, file = paste0(fm_prefix, ".noCS.RData"))
    } else {
        # Save successful SuSiE result
        save(s2, d2, file = paste0(fm_prefix, ".finemapping.RData"))
        load(paste0(fm_prefix, ".finemapping.RData"))
    }
}

message("Processing completed for: ", TISSUENAME, " at ", GENENAME)
