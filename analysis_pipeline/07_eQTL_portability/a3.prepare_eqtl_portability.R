#!/usr/bin/env Rscript

CURRENT_DIR <- "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/2026-05-12-compare_with_gtex_revision"
setwd(CURRENT_DIR)

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
})

GTOP_lead_info <- fread(
    "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/xqtl_atlas/gtop/snv_eqtl/07_summary/lead_qtl.txt.gz"
)
GTOP_sig_info <- fread(
    "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/xqtl_atlas/gtop/snv_eqtl/07_summary/sig_qtl.txt.gz"
)

GTEx_lead_info <- fread(
    "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/xqtl_atlas/gtexv8_eqtl/07_summary/lead_qtl.txt.gz"
)

tissue_link <- fread("input/gtop_gtex_tissues.txt")

argvs <- commandArgs(trailingOnly = TRUE)
row_index <- as.numeric(argvs[1])

gtop_tissue <- tissue_link$GTOP_tissue[row_index]
gtex_tissue <- tissue_link$GTEx_tissue[row_index]

## QTL information of GTOP significant variants in GTEx
gtop_gtex_info <- fread(
    sprintf(
        "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/xqtl_atlas/gtexv8_eqtl/08_gtop_overlap/%s.txt.gz",
        gtop_tissue
    )
)
length(unique(gtop_gtex_info$stable_ensg))
summary(gtop_gtex_info)
gtop_gtex_info <- gtop_gtex_info[!is.na(gtop_gtex_info$gtex_se), ]

GTEx_threshold <- fread(
    sprintf(
        "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/xqtl_atlas/gtexv8_eqtl/03_permutations/%s.txt.gz",
        gtex_tissue
    )
)
GTEx_threshold$stable_ensg <- gsub("\\..+", "", GTEx_threshold$gene_id)
gtop_gtex_info <- merge(
    gtop_gtex_info,
    GTEx_threshold[, .(stable_ensg, gtex_threshold = pval_nominal_threshold)],
    by = "stable_ensg"
)
length(unique(gtop_gtex_info$stable_ensg))

GTOP_lead_variant <- GTOP_lead_info[tissue == gtop_tissue]
GTOP_sig_variant <- GTOP_sig_info[tissue == gtop_tissue]

gtop_gtex_info <- merge(
    gtop_gtex_info,
    GTOP_lead_variant[, .(
        stable_ensg,
        gtop_threshold = pval_nominal_threshold
    )],
    by = "stable_ensg"
)

gtop_gtex_info <- merge(
    gtop_gtex_info,
    GTOP_lead_variant[, .(stable_ensg, chr_pos_ref_alt, lead_gtop = "Yes")],
    by = c("stable_ensg", "chr_pos_ref_alt"),
    all.x = TRUE
)
gtop_gtex_info <- merge(
    gtop_gtex_info,
    GTOP_sig_variant[, .(stable_ensg, chr_pos_ref_alt, sig_gtop = "Yes")],
    by = c("stable_ensg", "chr_pos_ref_alt"),
    all.x = TRUE
)

gtop_gtex_info$lead_gtop[is.na(gtop_gtex_info$lead_gtop)] <- "No"
gtop_gtex_info$sig_gtop[is.na(gtop_gtex_info$sig_gtop)] <- "No"

## Lead eVariant
merged_df1 <- gtop_gtex_info[lead_gtop == "Yes"]
merged_df1$type1 <- dplyr::case_when(
    merged_df1$gtop_beta / merged_df1$gtex_beta > 0.5 &
        merged_df1$gtop_beta / merged_df1$gtex_beta < 2 ~ "consistent",
    TRUE ~ "inconsistent"
)

## mash results
mash_info <- readRDS(
    "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/2026-05-23-he_QTL_revision/output/data/mash/noct_eqtl/top_pairs_lead/m.s_zval.RDS"
)
mash_info <- mash_info$result$PosteriorMean
mash_info <- as.data.frame(mash_info)

beta_se_info <- read.table(
    "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/2026-05-23-he_QTL_revision/output/data/mash/noct_eqtl/Strong/strong_beta_se.MashR_input.txt.gz",
    row.names = 1,
    header = T,
    check.names = FALSE
)
all(rownames(mash_info) %in% rownames(beta_se_info))

beta_se_info <- beta_se_info[rownames(mash_info), ]

gtop_name1 <- sprintf("GTOP_%s_slope_se", gtop_tissue)
gtop_name2 <- sprintf("GTOP_%s", gtop_tissue)
gtex_name1 <- sprintf("GTEx_%s_slope_se", gtex_tissue)
gtex_name2 <- sprintf("GTEx_%s", gtex_tissue)

GTOP_beta_value <- beta_se_info[[gtop_name1]] * mash_info[[gtop_name2]]
GTEx_beta_value <- beta_se_info[[gtex_name1]] *
    mash_info[[gtex_name2]]

mash_beta_info <- data.table(
    "stable_ensg" = gsub(".+,", "", rownames(beta_se_info)),
    "chr_pos_ref_alt" = gsub(",.+", "", rownames(beta_se_info)),
    "GTOP_mash_beta" = GTOP_beta_value,
    "GTEx_mash_beta" = GTEx_beta_value
)
merged_df1 <- merge(
    merged_df1,
    mash_beta_info,
    by = c("stable_ensg", "chr_pos_ref_alt")
)

merged_df1$type2 <- dplyr::case_when(
    merged_df1$GTOP_mash_beta / merged_df1$GTEx_mash_beta > 0.5 &
        merged_df1$GTOP_mash_beta / merged_df1$GTEx_mash_beta <
            2 ~ "consistent",
    TRUE ~ "inconsistent"
)

## uncertainty-aware
z <- qnorm(0.975)
merged_df1$ratio_uncertainty <- merged_df1$gtex_beta / merged_df1$gtop_beta
merged_df1$ratio_se <- sqrt(
    merged_df1$gtex_se^2 /
        merged_df1$gtop_beta^2 +
        (merged_df1$gtex_beta^2 * merged_df1$gtop_se^2) / merged_df1$gtop_beta^4
)
merged_df1$lowCI <- merged_df1$ratio_uncertainty - z * merged_df1$ratio_se
merged_df1$upCI <- merged_df1$ratio_uncertainty + z * merged_df1$ratio_se

merged_df1$type3 <- dplyr::case_when(
    merged_df1$lowCI > 0.5 &
        merged_df1$upCI < 2 ~ "consistent",
    TRUE ~ "inconsistent"
)

# uncertainty-aware
merged_df1$type4 <- dplyr::case_when(
    merged_df1$gtex_p_value < merged_df1$gtex_threshold ~ "consistent",
    TRUE ~ "inconsistent"
)

# ## all eVariants
merged_df2 <- gtop_gtex_info[sig_gtop == "Yes"]
merged_df2$type1 <- dplyr::case_when(
    merged_df2$gtop_beta / merged_df2$gtex_beta > 0.5 &
        merged_df2$gtop_beta / merged_df2$gtex_beta < 2 ~ "consistent",
    TRUE ~ "inconsistent"
)
merged_df2$type2 <- dplyr::case_when(
    merged_df2$gtex_p_value < merged_df2$gtex_threshold ~ "consistent",
    TRUE ~ "inconsistent"
)

merged_df3 <- merged_df2 |>
    group_by(stable_ensg) |>
    summarise(sum_count = sum(type2 == "consistent"))

gtop_gtex_info <- gtop_gtex_info[!is.na(gtop_gtex_info$gtop_se), ]
GTEx_info_list <- split(gtop_gtex_info, gtop_gtex_info$stable_ensg)

library(future)
library(future.apply)
library(progressr)
plan(multisession, workers = 20)
handlers(global = TRUE)

with_progress({
    p <- progressor(along = GTEx_info_list)
    GTEx_coloc_info <- future_lapply(
        GTEx_info_list,
        function(tmp_data) {
            D1 <- list(
                snp = tmp_data$chr_pos_ref_alt,
                beta = tmp_data$gtex_beta,
                varbeta = tmp_data$gtex_se^2,
                MAF = tmp_data$gtex_maf,
                N = tmp_data$gtex_n[1],
                type = "quant"
            )
            D2 <- list(
                snp = tmp_data$chr_pos_ref_alt,
                beta = tmp_data$gtop_beta,
                varbeta = tmp_data$gtop_se^2,
                MAF = tmp_data$gtop_maf,
                N = tmp_data$gtop_n[1],
                type = "quant"
            )
            tmp_coloc <- coloc::coloc.abf(
                dataset1 = D1,
                dataset2 = D2,
                p1 = 1e-4,
                p2 = 1e-4,
                p12 = 1e-5
            )
            p()
            tmp_coloc$summary[6]
        },
        future.seed = TRUE
    )
})
plan(sequential)

res_list <- list(
    "portability_df" = gtop_gtex_info,
    "summary_count" = data.frame(
        "tissue" = gtop_tissue,
        "type1" = rep(
            c(
                "consistent",
                "inconsistent"
            ),
            6
        ),
        "type2" = c(
            rep("ratio_nominal", 2),
            rep("ratio_mash", 2),
            rep("ratio_uncertainty-aware", 2),
            rep("FDR", 2),
            rep("gene", 4)
        ),
        "type3" = c(
            rep("Lead eVariant", 2),
            rep("Lead eVariant", 2),
            rep("Lead eVariant", 2),
            # rep("All eVariants", 2),
            rep("Lead eVariant", 2),
            # rep("All eVariants", 2),
            rep("All eGenes", 2),
            rep("coloc genes", 2)
        ),
        "count" = c(
            sum(merged_df1$type1 == "consistent"),
            sum(merged_df1$type1 == "inconsistent"),
            # sum(merged_df2$type1 == "consistent"),
            # sum(merged_df2$type1 == "inconsistent"),
            sum(merged_df1$type2 == "consistent"),
            sum(merged_df1$type2 == "inconsistent"),
            sum(merged_df1$type3 == "consistent"),
            sum(merged_df1$type3 == "inconsistent"),
            sum(merged_df1$type4 == "consistent"),
            sum(merged_df1$type4 == "inconsistent"),
            # sum(merged_df2$type2 == "consistent"),
            # sum(merged_df2$type2 == "inconsistent"),
            sum(merged_df3$sum_count > 0),
            sum(merged_df3$sum_count == 0),
            sum(GTEx_coloc_info > 0.5),
            sum(GTEx_coloc_info <= 0.5)
        )
    )
)
saveRDS(
    res_list,
    sprintf("output/data/snv_eqtl/portability/%s.rds", gtop_tissue)
)
