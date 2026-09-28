#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(data.table)
library(ieugwasr)
library(susieR)
library(coloc)

argvs <- commandArgs(trailingOnly = TRUE)

# Check number of arguments
if (length(argvs) != 5) {
    stop(
        "Error: Exactly 5 command-line arguments are required: GWASNAME CHRNAME SENTINELSNP LOCISTART LOCIEND"
    )
}

GWASNAME <- as.character(argvs[1])
CHRNAME <- as.character(argvs[2])
SENTINELSNP <- as.character(argvs[3])
LOCISTART <- as.numeric(argvs[4])
LOCIEND <- as.numeric(argvs[5])

# 1. Load GWAS trait type from metadata file
gwas_qtype_file <- "input/GWAS_qtype.txt"
if (!file.exists(gwas_qtype_file)) {
    stop("Error: GWAS_qtype.txt file not found at expected path.")
}
gwas_qtype_df <- fread(gwas_qtype_file, stringsAsFactors = FALSE)
if (!GWASNAME %in% gwas_qtype_df$GWAS_name) {
    stop(sprintf(
        "Error: GWAS name '%s' not found in GWAS_qtype.txt.",
        GWASNAME
    ))
}
gwas_qtype <- gwas_qtype_df[GWAS_name == GWASNAME, Type]


# 2. Read GWAS summary statistics for the specified chromosome
gwas_input_path <- sprintf(
    "/path/to/All_EAS/susie_coloc_split/%s/%s_coloc.%s.txt.gz",
    GWASNAME,
    GWASNAME,
    CHRNAME
)
if (!file.exists(gwas_input_path)) {
    stop(sprintf("Error: GWAS input file not found: %s", gwas_input_path))
}

gwas_data <- fread(gwas_input_path, stringsAsFactors = FALSE)
gwas_loci_data <- gwas_data[position >= LOCISTART & position <= LOCIEND]
gwas_loci_data <- gwas_loci_data[!is.na(se)]

required_cols <- c("position", "se", "beta", "maf", "P", "rsid", "N") # Ensure required columns exist
missing_cols <- setdiff(required_cols, names(gwas_loci_data))
if (length(missing_cols) > 0) {
    stop(sprintf(
        "Error: Missing required columns in GWAS data: %s",
        paste(missing_cols, collapse = ", ")
    ))
}

gwas_loci_data[, varbeta := se^2] # Compute variance of beta (needed for coloc)
gwas_loci_data[,
    rsid_clean := ifelse(
        grepl("^rs\\d+", rsid),
        rsid,
        paste0(
            "chr",
            gsub(":", "_", gsub("^chr", "", rsid, ignore.case = TRUE))
        )
    )
]

setDF(gwas_loci_data)
gwas_loci_data <- gwas_loci_data[!duplicated(gwas_loci_data$rsid_clean), ]
rownames(gwas_loci_data) <- gwas_loci_data$rsid_clean


# 3. Compute LD matrix using 1000G EAS reference panel
# bfile_path <- sprintf("/media/bora_A/zhangt/src/data/1000G/five_ancestry_groups/EAS/splitbychr/1000G.EAS.maf01.%s", CHRNAME)
# plink_bin  <- "/media/bora_A/zhangt/src/bin/plink"
bfile_path <- sprintf(
    "/path/to/splitbychr/1000G.EAS.maf01.%s",
    CHRNAME
)
plink_bin <- "/lustre/home/tzhang/src/bin/plink"

if (!file.exists(plink_bin)) {
    stop("Error: PLINK binary not found at specified path.")
}

if (!all(file.exists(paste0(bfile_path, c(".bed", ".bim", ".fam"))))) {
    stop(sprintf(
        "Error: PLINK reference files missing for chromosome: %s",
        bfile_path
    ))
}

# Compute LD matrix
GWAS_EAS_ld <- tryCatch(
    {
        ld_matrix(
            variants = gwas_loci_data$rsid_clean,
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
common_snps <- intersect(rownames(GWAS_EAS_ld), rownames(gwas_loci_data))
if (length(common_snps) == 0) {
    stop("Error: No overlapping SNPs between GWAS data and LD reference panel.")
}

# Subset GWAS data and LD matrix to common SNPs only
gwas_loci_data_sub <- gwas_loci_data[common_snps, , drop = FALSE]
GWAS_EAS_ld_sub <- GWAS_EAS_ld[common_snps, common_snps, drop = FALSE]


# 4. Run SuSiE fine-mapping with comprehensive error handling
d1 <- list(
    beta = gwas_loci_data_sub$beta,
    varbeta = gwas_loci_data_sub$varbeta,
    MAF = gwas_loci_data_sub$maf,
    pval = gwas_loci_data_sub$P,
    snp = rownames(gwas_loci_data_sub),
    LD = as.matrix(GWAS_EAS_ld_sub),
    N = gwas_loci_data_sub$N[1],
    type = gwas_qtype,
    chrom = chr,
    pos = position
)

fm_prefix <- sprintf(
    "input/finemapping_GWAS_1KGP_EAS/%s_%s",
    GWASNAME,
    SENTINELSNP
)
kriging_prefix <- sprintf(
    "input/finemapping_GWAS_1KGP_EAS/kriging_rss/%s_%s",
    GWASNAME,
    SENTINELSNP
)
if (!dir.exists(dirname(kriging_prefix))) {
    dir.create(dirname(kriging_prefix), recursive = TRUE, showWarnings = FALSE)
}

s1 <- tryCatch(
    {
        runsusie(
            d1,
            maxit = 1000,
            repeat_until_convergence = FALSE,
            R_finite = 504,
            R_mismatch = "eb"
        )
    },
    error = function(e) {
        message("SuSiE failed to converge: ", conditionMessage(e))
        NULL
    }
)

condz_in <- kriging_rss(d1$beta / sqrt(d1$varbeta), d1$LD, n = d1$N)
condz_in_df <- condz_in$conditional_dist
rownames(condz_in_df) <- d1$snp
fwrite(
    condz_in_df,
    file = sprintf("%s.txt.gz", kriging_prefix),
    sep = "\t",
    row.names = TRUE
)

d1$LD <- NULL

# Save result based on outcome
if (is.null(s1)) {
    # Save input data when SuSiE fails entirely
    save(d1, file = paste0(fm_prefix, ".noConverged.RData"))
} else {
    susie_summary <- tryCatch(summary(s1), error = function(e) NULL)

    if (is.null(susie_summary) || is.null(susie_summary$cs)) {
        # Save input data when no credible sets are returned
        save(d1, file = paste0(fm_prefix, ".noCS.RData"))
    } else {
        # Save successful SuSiE result
        save(s1, d1, file = paste0(fm_prefix, ".finemapping.RData"))
    }
}

# cs_corr <- get_cs_correlation(s1, Xcorr = d1$LD)
# round(cs_corr, 3)

# condz_in = kriging_rss(d1$beta / sqrt(d1$varbeta), d1$LD, n = d1$N)
# save(condz_in, file = paste0(kriging_prefix, ".RData"))

# check_alignment(d1, thr = 0.2, do_plot = FALSE)

# message("Processing completed for: ", GWASNAME, " at ", SENTINELSNP)

# gwas_fm_type <- case_when(
#     file.exists(paste0(fm_prefix, ".noConverged.RData")) ~ "noConverged",
#     file.exists(paste0(fm_prefix, ".noCS.RData")) ~ "noCS",
#     file.exists(paste0(fm_prefix, ".finemapping.RData")) ~ "finemapping",
#     file.exists(paste0(fm_prefix, ".susie_ser.RData")) ~ "susie_ser"
# )

# condz_in_res <- fread(sprintf("input/finemapping_GWAS_1KGP_EAS/kriging_rss/%s_%s.txt.gz", GWASNAME, SENTINELSNP))

# if(any(condz_in_res$logLR>2 & abs(condz_in_res$z>2))){
#   system(sprintf("mv %s.%s.RData %s.noConfidentLD.RData", fm_prefix, gwas_fm_type, fm_prefix))
# }
