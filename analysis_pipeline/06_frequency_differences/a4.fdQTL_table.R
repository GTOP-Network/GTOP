#!/usr/bin/env Rscript

# ==============================================================================
# fine-mapping QTL to independent loci
# ==============================================================================

#%% ------------------------ 0. prepare files (packages, input files, output files)
PROJECT_DIR <- "/path/to/dir"
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(ggrastr)
    library(igraph)
})

variant_freq <- fread(
    "../path/to/output/data/integrated_freq.txt"
)


fdQTL_merged_data <- lapply(
    c("snv_eqtl", "snv_juqtl", "snv_tuqtl"),
    function(tmp_xqtl_type) {
        fdQTL_res <- fread(
            sprintf("output/data/fm_freq/%s_fm_rare.txt", tmp_xqtl_type)
        )

        fdQTL_res$ensg <- gsub(";.+", "", fdQTL_res$gene_locus)
        fdQTL_res$indep_loci <- gsub(".+;", "", fdQTL_res$gene_locus)
        fdQTL_res$cs_id <- fdQTL_res$indep_loci
        fdQTL_res$qtl_type <- tmp_xqtl_type
        table(table(fdQTL_res$ensg))

        fdQTL_res <- fdQTL_res |> separate_rows(cs_id, sep = "\\|")
        fdQTL_res <- fdQTL_res |> separate_rows(cs_id, sep = "\\+")
        fdQTL_res <- fdQTL_res |> separate_rows(cs_id, sep = ",")
        fdQTL_res$tissue <- gsub("_L\\d+$", "", fdQTL_res$cs_id)
        fdQTL_res$cs_id <- gsub(".+_", "", fdQTL_res$cs_id)

        fm_info <- fread(
            sprintf(
                "output/data/cs_combining/%s/collapsed_cs_results.txt",
                tmp_xqtl_type
            )
        )
        all(fdQTL_res$indep_loci %in% fm_info$indep_loci)

        setDT(fdQTL_res)
        fdQTL_info <- merge(
            fdQTL_res[, .(
                gene_locus,
                tissue,
                ensg,
                indep_loci,
                cs_id,
                afr,
                eur,
                sas,
                amr
            )],
            fm_info[, .(
                xqtl_type,
                tissue,
                phenotype_id,
                ensg,
                symbol,
                variant_id,
                pip,
                indep_loci
            )],
            by = c("tissue", "ensg", "indep_loci")
        )
        length(unique(fdQTL_info$gene_locus))

        fdQTL_info <- merge(
            fdQTL_info,
            variant_freq[, .(
                variant_id,
                af_afr,
                af_amr,
                af_eas,
                af_eur,
                af_sas,
                af_neas
            )],
            by = "variant_id"
        )

        fdQTL_output <- fdQTL_info[,
            .(
                xqtl_type,
                tissue,
                phenotype_id,
                ensg,
                symbol,
                variant_id,
                cs_id,
                pip,
                indep_loci,
                maf_eas = signif(0.5 - abs(0.5 - af_eas), 3),
                maf_afr = signif(0.5 - abs(0.5 - af_afr), 3),
                maf_eur = signif(0.5 - abs(0.5 - af_eur), 3),
                maf_sas = signif(0.5 - abs(0.5 - af_sas), 3),
                maf_amr = signif(0.5 - abs(0.5 - af_amr), 3),
                afr,
                eur,
                sas,
                amr
            )
        ]

        cols <- c("afr", "eur", "sas", "amr")

        fdQTL_output[,
            (cols) := lapply(.SD, function(x) {
                ifelse(
                    x == 0,
                    "U",
                    ifelse(x == 0.01, "R", ifelse(x == 0.1, "C", x))
                )
            }),
            .SDcols = cols
        ]
        fdQTL_output
    }
)

juqtl_genes <- unique(fdQTL_merged_data[[2]]$ensg)

fdQTL_merged_df <- rbind(
    fdQTL_merged_data[[1]],
    fdQTL_merged_data[[2]],
    fdQTL_merged_data[[3]][!fdQTL_merged_data[[3]]$ensg %in% juqtl_genes]
)

tmp_names <- paste0(fdQTL_merged_df$ensg, "_", fdQTL_merged_df$indep_loci)
length(unique(tmp_names[fdQTL_merged_df$xqtl_type == "snv_eqtl"]))
length(unique(tmp_names[fdQTL_merged_df$xqtl_type != "snv_eqtl"]))

fdQTL_merged_df$xqtl_type <- case_when(
    fdQTL_merged_df$xqtl_type == "snv_eqtl" ~ "eQTL",
    fdQTL_merged_df$xqtl_type == "snv_juqtl" ~ "sQTL (junction usage)",
    fdQTL_merged_df$xqtl_type == "snv_tuqtl" ~ "sQTL (transcript usage)",
)

cols <- grep("^maf_", names(fdQTL_merged_df), value = TRUE)

fdQTL_merged_df[, (cols) := lapply(.SD, round, digits = 3), .SDcols = cols]

fdQTL_merged_df <- fdQTL_merged_df |>
    arrange(xqtl_type, phenotype_id, indep_loci, tissue)

fwrite(fdQTL_merged_df, "output/data/fdQTL_all_info.txt", sep = "\t")


uniqe_loci <- unique(fdQTL_merged_df[, .(
    xqtl_type,
    ensg,
    indep_loci,
    afr,
    eur,
    sas,
    amr
)])
table(uniqe_loci$xqtl_type == "eQTL")

uniqe_loci1 <- uniqe_loci[uniqe_loci$xqtl_type == "eQTL", ]
table(
    uniqe_loci1$afr != "C" &
        uniqe_loci1$eur != "C" &
        uniqe_loci1$sas != "C" &
        uniqe_loci1$amr != "C"
) # 76

uniqe_loci2 <- uniqe_loci[uniqe_loci$xqtl_type != "eQTL", ]
table(
    uniqe_loci2$afr != "C" &
        uniqe_loci2$eur != "C" &
        uniqe_loci2$sas != "C" &
        uniqe_loci2$amr != "C"
) # 184
