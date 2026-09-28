#!/usr/bin/env Rscript

# ==============================================================================
# Meta-analysed summary statistics
# ==============================================================================

#%% ------------------------ 0. prepare files (packages, input files, output files)
CURRENT_DIR <- "/path/to/dir"
setwd(CURRENT_DIR)

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(patchwork)
    library(ggrastr)
})

## beta se values
strong_beta_se <- vroom::vroom(
    "/path/to/output/data/mash/noct_eqtl/Strong/strong_beta_se.MashR_input.txt.gz",
) %>%
    as.data.frame() %>%
    column_to_rownames("pair_id")

raw_z.stat <- strong_beta_se[, seq(3, ncol(strong_beta_se), 3)]
colnames(raw_z.stat) <- gsub("_zval", "", colnames(raw_z.stat))

raw_beta <- strong_beta_se[, seq(1, ncol(strong_beta_se), 3)]
colnames(raw_beta) <- gsub("_slope", "", colnames(raw_beta))

## lfsr
posterior_list <- readRDS(
    "/path/to/output/data/mash/noct_eqtl/top_pairs_lead/m.s_zval.RDS",
)
lfsr_mash <- posterior_list$result$lfsr
lfsr_mash <- as.data.frame(lfsr_mash)

lfsr_mash$stable_gene_id <- gsub(".+,", "", rownames(lfsr_mash))
lfsr_mash$chr_pos_ref_alt <- gsub(",.+", "", rownames(lfsr_mash))


#%% ------------------------ 1. analysis
tissue_change <- c(
    "Adipose" = "Adipose_Visceral_Omentum",
    "Adrenal_Gland" = "Adrenal_Gland",
    "Liver" = "Liver",
    "Muscle" = "Muscle_Skeletal",
    "Pancreas_Body" = "Pancreas",
    "Pancreas_Head" = "Pancreas",
    "Pancreas_Tail" = "Pancreas",
    "Skin" = "Skin_Not_Sun_Exposed_Suprapubic",
    "Spleen" = "Spleen",
    "Whole_Blood" = "Whole_Blood"
)
new_tissue_change <- paste0("GTEx_", tissue_change)
names(new_tissue_change) <- paste0("GTOP_", names(tissue_change))

sig_num_list <- rbindlist(lapply(names(new_tissue_change), function(x) {
    ## tissue information
    col1 <- x
    col2 <- new_tissue_change[x]

    portability_info <- readRDS(sprintf(
        "/path/to/output/data/snv_eqtl/portability/%s.rds",
        gsub("GTOP_", "", x)
    ))
    portability_info <- portability_info$portability_df
    nominal_gene <- unique(portability_info$stable_ensg[
        portability_info$sig_gtop == "Yes"
    ])

    ## Top eSNPs (raw portable) top SNPs
    portable_lead_df <- portability_info[lead_gtop == "Yes"]
    portable_lead_df$type1 <- dplyr::case_when(
        portable_lead_df$gtex_p_value <
            portable_lead_df$gtex_threshold ~ "consistent",
        TRUE ~ "inconsistent"
    )

    ## Portable Gene
    portable_gene_df <- portability_info[sig_gtop == "Yes"]
    portable_gene_df$type1 <- dplyr::case_when(
        portable_gene_df$gtex_p_value <
            portable_gene_df$gtex_threshold ~ "consistent",
        TRUE ~ "inconsistent"
    )
    portable_gene_df <- portable_gene_df |>
        group_by(stable_ensg) |>
        summarise(sum_count = sum(type1 == "consistent"))

    ## mash
    tmp_lfsr <- as.data.frame(lfsr_mash[, c(
        "stable_gene_id",
        "chr_pos_ref_alt",
        col1,
        col2
    )])

    tmp_lfsr_lead <- tmp_lfsr |>
        group_by(stable_gene_id) |>
        filter(!!sym(col1) == min(!!sym(col1)))

    tmp_lfsr_lead <- tmp_lfsr_lead[
        !duplicated(tmp_lfsr_lead$stable_gene_id),
    ]

    tmp_lfsr_lead <- tmp_lfsr_lead[tmp_lfsr_lead[[col1]] < 0.05, ]
    tmp_lfsr_lead$type1 <- dplyr::case_when(
        tmp_lfsr_lead[[col2]] < 0.05 ~ "consistent",
        TRUE ~ "inconsistent"
    )
    mashr_egene <- unique(tmp_lfsr_lead$stable_gene_id)

    ## eGenes
    tmp_lfsr_gene <- tmp_lfsr |>
        group_by(stable_gene_id) |>
        summarise(GTOP_min = min(!!sym(col1)), GTEx_min = min(!!sym(col2)))
    tmp_lfsr_gene <- tmp_lfsr_gene[tmp_lfsr_gene$GTOP_min < 0.05, ]
    tmp_lfsr_gene$type1 <- dplyr::case_when(
        tmp_lfsr_gene$GTEx_min < 0.05 ~ "consistent",
        TRUE ~ "inconsistent"
    )

    data.frame(
        "tissue" = gsub("GTOP_", "", x),
        "type1" = rep(
            c(
                "consistent",
                "inconsistent"
            ),
            4
        ),
        "type2" = c(
            rep("nominal", 2),
            rep("mash", 2),
            rep("nominal", 2),
            rep("mash", 2)
        ),
        "type3" = c(
            rep("Lead eVariant", 4),
            rep("All eGenes", 4)
        ),
        "value" = c(
            sum(portable_lead_df$type1 == "consistent"),
            sum(portable_lead_df$type1 == "inconsistent"),
            sum(tmp_lfsr_lead$type1 == "consistent"),
            sum(tmp_lfsr_lead$type1 == "inconsistent"),
            sum(portable_gene_df$sum_count > 0),
            sum(portable_gene_df$sum_count == 0),
            sum(tmp_lfsr_gene$type1 == "consistent"),
            sum(tmp_lfsr_gene$type1 == "inconsistent")
        )
    )
}))

sig_num_list$tissue <- factor(
    sig_num_list$tissue,
    levels = gsub("GTOP_", "", names(new_tissue_change))
)

sig_num_list$type2 <- factor(sig_num_list$type2, levels = c("nominal", "mash"))

mash_summary_df <- sig_num_list |>
    group_by(tissue, type2, type3) |>
    mutate(all_count = sum(value)) |>
    mutate(ratio = value / all_count)

mash_summary_df[mash_summary_df$type1 == "consistent", ] %>%
    group_by(type2, type3) %>%
    summarise(mean_ratio = mean(ratio, na.rm = T))
