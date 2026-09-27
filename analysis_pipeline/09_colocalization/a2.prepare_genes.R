#!/usr/bin/env Rscript

# Load required libraries
library(data.table)

# ------------------------------------------------------------
# Function: validate and parse command-line arguments
# Expected input order:
#   1. QTLTYPE (e.g., "SNV_eQTL")
#   2. TISSUENAME (e.g., "Pancreas_Body")
# ------------------------------------------------------------
argvs <- commandArgs(trailingOnly = TRUE)

# Check number of arguments
if (length(argvs) != 2) {
    stop(
        "Error: Exactly 2 command-line arguments are required: QTLTYPE TISSUENAME"
    )
}

QTLTYPE <- argvs[1]
TISSUENAME <- argvs[2]

output_file <- sprintf(
    "input/GWAS_QTL_overlap/%s/%s.txt",
    QTLTYPE,
    TISSUENAME
)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)


# if (!file.exists("input/EAS_GWAS_sig_variants.txt")) {
#     stop("xQTL file not found: input/EAS_GWAS_sig_variants.txt")
# }
# GWAS_sig_variants <- fread("input/EAS_GWAS_sig_variants.txt")
# GWAS_sig_variants <- GWAS_sig_variants[
#     !(GWAS_sig_variants$V2 == "chr6" &
#         GWAS_sig_variants$V3 > 28510120 &
#         GWAS_sig_variants$V3 < 33480577),
# ]

# sig_xqtl_path <- sprintf(
#     "/lustre/home/tzhang/2026-05-07-gtop_xqtl-Project/xqtl_atlas/gtop/%s/06_nominal_slim/nom_thresh/%s.nom_thresh.txt.gz",
#     QTLTYPE,
#     TISSUENAME
# )
# if (!file.exists(sig_xqtl_path)) {
#     stop("xQTL file not found: ", sig_xqtl_path)
# }
# sig_xQTL <- fread(sig_xqtl_path)
# sig_GWAS_xQTL <- sig_xQTL[variant_id %in% GWAS_sig_variants$V5]
# selected_gwas <- unique(GWAS_sig_variants$V1[GWAS_sig_variants$V5%in%sig_GWAS_xQTL$variant_id])

# Load GWAS loci once
message("Loading GWAS loci...")

GWAS_loci <- fread("input/EAS_GWAS_noMHC_all_sentinal_loci.txt")
# GWAS_loci <- GWAS_loci[GWAS_name%in%selected_gwas]

# Load xQTL lead variants
xqtl_path <- sprintf(
    "/lustre/home/tzhang/2026-05-07-gtop_xqtl-Project/xqtl_atlas/gtop/%s/04_xgene_info/%s.txt",
    QTLTYPE,
    TISSUENAME
)
if (!file.exists(xqtl_path)) {
    stop("xQTL file not found: ", xqtl_path)
}
lead_xQTL <- fread(xqtl_path)
# lead_xQTL <- lead_xQTL[phenotype_id %in% sig_GWAS_xQTL$phenotype_id]

# Parse chr and pos from chr_pos_ref_alt in one go (vectorized)
lead_xQTL[,
    c("lead_var_chr", "lead_var_pos") := tstrsplit(
        chr_pos_ref_alt,
        "_",
        keep = c(1, 2)
    )
]
lead_xQTL$lead_var_pos <- as.numeric(lead_xQTL$lead_var_pos)

GWAS_loci$tmp_start <- GWAS_loci$position - 250000
GWAS_loci$tmp_end <- GWAS_loci$position + 250000

# Use non-equi join instead of row-by-row loop (much faster!)
GWAS_QTL_pairdf <- GWAS_loci[
    lead_xQTL,
    on = .(
        chr == lead_var_chr,
        tmp_start < lead_var_pos,
        tmp_end > lead_var_pos
    ),
    .(
    xqtl_type = QTLTYPE,
    tissue = TISSUENAME,
    qtl_phenotype_id = phenotype_id,
    qtl_variant_id = variant_id,
    qtl_chr_pos_ref_alt = chr_pos_ref_alt,
    gwas_name = GWAS_name,
    chr=x.chr,
    position,
    rsid,
    loci_start = x.loci_start,
    loci_end = x.loci_end
    ),
    nomatch = 0L,
    allow.cartesian = TRUE
]

# Add metadata columns
# GWAS_QTL_pairdf[, `:=`(
#     xQTL_type = QTLTYPE,
#     tissue = TISSUENAME,
#     QTL_phenotype_id = phenotype_id,
#     QTL_variant_id = variant_id,
#     QTL_chr_pos_ref_alt = chr_pos_ref_alt
# )]

# # Select and reorder output columns
# result_cols <- c(
#     "xQTL_type",
#     "tissue",
#     "QTL_phenotype_id",
#     "QTL_variant_id",
#     "QTL_chr_pos_ref_alt",
#     "GWAS_name",
#     "gwas_chr",
#     "position",
#     "rsid",
#     "loci_start",
#     "loci_end"
# )
# GWAS_QTL_pairdf <- GWAS_QTL_pairdf[, ..result_cols]
# setnames(GWAS_QTL_pairdf, "gwas_chr", "chr")
GWAS_QTL_pairdf <- unique(GWAS_QTL_pairdf)

# Write output
fwrite(GWAS_QTL_pairdf, output_file, sep = "\t")
message("Output written to: ", output_file)
