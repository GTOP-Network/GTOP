#!/usr/bin/env Rscript

# ==============================================================================
# Fine-mapping summary (Xudong Zou, GTOP xQTL fine-mapping)
# ==============================================================================

#%% ------------------------ 0. prepare files (packages, input files, output files)
PROJECT_DIR <- "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/2026-05-11-fine_mapping_revision"
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(ggrastr)
    library(pbmcapply)
    library(igraph)
})

# ------------ input
IN_EQTL_FILE <- "/media/london_B/zouxudong/2024-10-21-aGTEx-main/Revision_Nature/2026-05-25-joint-finemap/input/SuSiE_finemapping_summary/GTOP_finemapping.eQTL.snv_all_tissues.txt"
IN_JUQTL_FILE <- "/media/london_B/zouxudong/2024-10-21-aGTEx-main/Revision_Nature/2026-05-25-joint-finemap/input/SuSiE_finemapping_summary/GTOP_finemapping.ju_sQTL.snv_all_tissues.txt"
IN_TUQTL_FILE <- "/media/london_B/zouxudong/2024-10-21-aGTEx-main/Revision_Nature/2026-05-25-joint-finemap/input/SuSiE_finemapping_summary/GTOP_finemapping.tu_sQTL.snv_all_tissues.txt"

PLINK_BIN <- "/media/bora_A/zhangt/src/bin/plink2"
GTOP_BFILE <- "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/xqtl_atlas/data/genotype/gtop/gtop_snv.maf05"
KGP_BFILE <- "/media/bora_A/zhangt/src/data/1000G/five_ancestry_groups/EAS/1000G.EAS.maf01"

# ------------ output
OUT_FM_FILE <- "output/data/fine_mapped_gtop_xqtl.txt"

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


#%% ------------------------ 1. Load fine-mapping data (Add gene name, symbol)
eqtl_susie <- fread(IN_EQTL_FILE)[, .(
    xqtl_type = "snv_eqtl",
    tissue = Tissue,
    phenotype_id = locus_id,
    cs_id = cs,
    cs_size,
    variant_id,
    pip
)]
length(unique(paste0(
    eqtl_susie$tissue,
    "_",
    eqtl_susie$phenotype_id,
    "_",
    eqtl_susie$cs_id
))) # 29114

juqtl_susie <- fread(IN_JUQTL_FILE)[, .(
    xqtl_type = "snv_juqtl",
    tissue = Tissue,
    phenotype_id = locus_id,
    cs_id = cs,
    cs_size,
    variant_id,
    pip
)]

juqtl_susie$ensg <- gsub(".+:", "", juqtl_susie$phenotype_id)
length(unique(paste0(
    juqtl_susie$tissue,
    "_",
    juqtl_susie$ensg,
    "_",
    juqtl_susie$cs_id
))) # 30075

tuqtl_susie <- fread(IN_TUQTL_FILE)[, .(
    xqtl_type = "snv_tuqtl",
    tissue = Tissue,
    phenotype_id = locus_id,
    cs_id = cs,
    cs_size,
    variant_id,
    pip
)]
tuqtl_susie$ensg <- gsub(".+_", "", tuqtl_susie$phenotype_id)
length(unique(paste0(
    tuqtl_susie$tissue,
    "_",
    tuqtl_susie$phenotype_id,
    "_",
    tuqtl_susie$cs_id
))) # 16757

length(unique(c(
    paste0(
        juqtl_susie$tissue,
        "_",
        juqtl_susie$ensg,
        "_",
        juqtl_susie$cs_id
    ),
    paste0(
        tuqtl_susie$tissue,
        "_",
        tuqtl_susie$phenotype_id,
        "_",
        tuqtl_susie$cs_id
    )
))) # 46832


length(unique(eqtl_susie$phenotype_id)) # 8376

eqtl_lead_qtl <- fread(
    "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/xqtl_atlas/gtop/snv_eqtl/07_summary/lead_qtl.txt.gz"
)
juqtl_lead_qtl <- fread(
    "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/xqtl_atlas/gtop/snv_juqtl/07_summary/lead_qtl.txt.gz"
)
tuqtl_lead_qtl <- fread(
    "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/xqtl_atlas/gtop/snv_tuqtl/07_summary/lead_qtl.txt.gz"
)

eqtl_susie <- merge(
    eqtl_susie,
    eqtl_lead_qtl[, .(tissue, phenotype_id, ensg, stable_ensg, symbol)],
    by = c("tissue", "phenotype_id")
)

juqtl_susie <- merge(
    juqtl_susie,
    juqtl_lead_qtl[, .(tissue, phenotype_id, ensg, stable_ensg, symbol)],
    by = c("tissue", "phenotype_id")
)

tuqtl_susie <- merge(
    tuqtl_susie,
    tuqtl_lead_qtl[, .(tissue, phenotype_id, ensg, stable_ensg, symbol)],
    by = c("tissue", "phenotype_id")
)

susie_res <- rbind(eqtl_susie, juqtl_susie, tuqtl_susie)
print(table(susie_res$xqtl_type))
# eQTL: 1567091 | juQTL: 1212854 | tuQTL: 803532

cat(sprintf(
    "Unique genes — eQTL: %d | other xQTL: %d\n",
    uniqueN(susie_res$ensg[susie_res$xqtl_type == "snv_eqtl"]),
    uniqueN(susie_res$ensg[susie_res$xqtl_type != "snv_eqtl"])
))

dir.create("output/data", recursive = T)
fwrite(susie_res, OUT_FM_FILE, sep = "\t")


#%% ------------------------ 2. CS distribution and inter/intra-CS LD analysis
raw_susie_res <- fread(OUT_FM_FILE)

lapply(c("snv_eqtl", "snv_juqtl", "snv_tuqtl"), function(tmp_qtl_type) {
    # tmp_qtl_type <- "snv_eqtl"

    ld_tmp_dir <- sprintf("output/data/cs_combining/%s/ld_tmp", tmp_qtl_type)
    dir.create(ld_tmp_dir, recursive = TRUE, showWarnings = FALSE)

    out_figure_dir <- sprintf(
        "output/figures/cs_combining/%s/",
        tmp_qtl_type
    )
    dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)

    ## Load dat
    qtl_susie_res <- raw_susie_res[xqtl_type == tmp_qtl_type]

    ## CS count per eGene and tissue
    tissue_gene_order <- qtl_susie_res |>
        group_by(tissue) |>
        summarise(n_genes = uniqueN(ensg), .groups = "drop") |>
        arrange(n_genes)

    cs_number_df <- qtl_susie_res |>
        group_by(tissue, ensg) |>
        summarise(cs_number = uniqueN(cs_id), .groups = "drop") |>
        group_by(tissue, cs_number) |>
        summarise(gene_number = uniqueN(ensg), .groups = "drop") |>
        mutate(tissue = factor(tissue, levels = rev(tissue_gene_order$tissue)))

    multi_cs_pct <- cs_number_df |>
        group_by(tissue) |>
        summarise(
            pct_multi = sum(gene_number[cs_number > 1]) /
                sum(gene_number) *
                100,
            .groups = "drop"
        )
    range(100 - multi_cs_pct$pct_multi)
    cat(sprintf(
        "Mean proportion of eGenes with > 1 CS: %.4f%%\n",
        mean(multi_cs_pct$pct_multi)
    ))
    # 5.7624%

    cs_number_ratio <- cs_number_df |>
        group_by(tissue) |>
        mutate(total_gene = sum(gene_number)) |>
        mutate(ratio = gene_number / total_gene)

    ## CS count distribution per tissue
    ggplot(
        cs_number_df,
        aes(x = tissue, y = gene_number, fill = as.character(cs_number))
    ) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(
            values = c("#81cee8", "#449abd", "#ba4d8f", "#ffc748", "#896aa7")
        ) +
        theme_pub() +
        theme(
            axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                vjust = 1,
                color = "black"
            ),
            axis.text.y = element_text(color = "black"),
            axis.ticks = element_line(color = "black")
        ) +
        labs(x = "", y = "Proportion of eGene", fill = "CS number") +
        geom_text(
            data = tissue_gene_order,
            aes(x = tissue, y = 0.97, label = n_genes),
            size = 3.5,
            inherit.aes = FALSE
        )

    ggsave(paste0(out_figure_dir, "number_of_cs.pdf"), width = 6, height = 5)

    ## Compute pairwise LD for variants in multi-CS loci
    qtl_susie_res[,
        cs_n := uniqueN(cs_id),
        by = .(tissue, phenotype_id, xqtl_type)
    ]

    multi_leads <- qtl_susie_res[
        cs_n > 1L,
        {
            i <- which.max(pip)
            j <- which.min(pip)
            .(
                lead_var = variant_id[i],
                lead_pip = pip[i],
                min_var = variant_id[j],
                min_pip = pip[j]
            )
        },
        by = .(tissue, phenotype_id, ensg, xqtl_type, cs_id)
    ]
    cat(sprintf("Rows in multi-CS variant set: %d\n", nrow(multi_leads))) # 3622

    run_ld <- function(bfile, snp_file, out_pref) {
        ld_file <- paste0(out_pref, ".vcor")

        cmd <- paste(
            PLINK_BIN,
            "--bfile",
            bfile,
            "--extract",
            snp_file,
            "--r2-unphased",
            "--ld-window 999999",
            "--ld-window-kb 99999",
            "--ld-window-r2 0",
            "--out",
            out_pref,
            "--silent"
        )

        system(
            cmd,
            ignore.stdout = TRUE,
            ignore.stderr = TRUE
        )

        if (!file.exists(ld_file)) {
            return(
                data.table(
                    SNP_A = character(),
                    SNP_B = character(),
                    R2 = numeric()
                )
            )
        }

        fread(ld_file)[,
            .(
                SNP_A = ID_A,
                SNP_B = ID_B,
                R2 = UNPHASED_R2
            )
        ]
    }

    ## Run 1KGP-EAS first, then GTOP
    ## GTOP only fills LD pairs absent from 1KGP
    tissue_ld_list <- lapply(
        unique(multi_leads$tissue),
        function(tmp_tissue) {
            tissue_vars <- unique(c(
                multi_leads[
                    tissue == tmp_tissue,
                    lead_var
                ],
                multi_leads[
                    tissue == tmp_tissue,
                    min_var
                ]
            ))

            tag <- sprintf(
                "%s_%d",
                gsub("[^A-Za-z0-9._-]", "_", tmp_tissue),
                length(tissue_vars)
            )

            snp_file <- file.path(
                ld_tmp_dir,
                paste0(tag, ".snplist")
            )

            writeLines(
                tissue_vars,
                snp_file
            )

            ## 1. 1KGP-EAS
            kgp_pref <- file.path(
                ld_tmp_dir,
                paste0(tag, "_1KGP_EAS")
            )

            kgp_ld <- run_ld(
                bfile = KGP_BFILE,
                snp_file = snp_file,
                out_pref = kgp_pref
            )

            kgp_ld[,
                ld_source := "1KGP-EAS"
            ]

            ## 2. GTOP
            gtop_pref <- file.path(
                ld_tmp_dir,
                paste0(tag, "_GTOP")
            )

            gtop_ld <- run_ld(
                bfile = GTOP_BFILE,
                snp_file = snp_file,
                out_pref = gtop_pref
            )

            gtop_ld[,
                ld_source := "GTOP"
            ]

            ## 3. Define orientation-independent pair rs1-rs2 == rs2-rs1
            if (nrow(kgp_ld)) {
                kgp_ld[,
                    pair_id := paste(
                        pmin(SNP_A, SNP_B),
                        pmax(SNP_A, SNP_B),
                        sep = "||"
                    )
                ]
            }

            if (nrow(gtop_ld)) {
                gtop_ld[,
                    pair_id := paste(
                        pmin(SNP_A, SNP_B),
                        pmax(SNP_A, SNP_B),
                        sep = "||"
                    )
                ]
            }

            ## 4. Remove GTOP pairs already available  in 1KGP-EAS
            if (nrow(kgp_ld) && nrow(gtop_ld)) {
                gtop_ld <- gtop_ld[
                    !pair_id %chin% kgp_ld$pair_id
                ]
            }

            ## 5. Merge
            res <- rbindlist(
                list(
                    kgp_ld,
                    gtop_ld
                ),
                use.names = TRUE,
                fill = TRUE
            )

            res[,
                pair_id := NULL
            ]

            message(
                sprintf(
                    "%s: 1KGP=%d pairs; GTOP supplement=%d pairs; total=%d",
                    tmp_tissue,
                    nrow(kgp_ld),
                    nrow(gtop_ld),
                    nrow(res)
                )
            )
            res
        }
    )

    names(tissue_ld_list) <- unique(multi_leads$tissue)

    ##. Classify SNP pairs as intra- or inter-CS
    susie_cs_ld <- rbindlist(lapply(
        unique(qtl_susie_res$tissue),
        function(tmp_tissue) {
            multi_cs_tissue <- multi_leads[
                tissue == tmp_tissue
            ]

            ## SNP -> phenotype / CS membership
            cs_membership <- unique(rbind(
                multi_cs_tissue[, .(
                    phenotype_id,
                    cs_id,
                    variant_id = lead_var
                )],
                multi_cs_tissue[, .(
                    phenotype_id,
                    cs_id,
                    variant_id = min_var
                )]
            ))

            ## LD pairs
            ld_pairs <- copy(
                tissue_ld_list[[tmp_tissue]]
            )

            ##---------------------------------------
            ## Make LD pair orientation-independent
            ##---------------------------------------
            ld_pairs[, `:=`(
                var1 = pmin(SNP_A, SNP_B),
                var2 = pmax(SNP_A, SNP_B)
            )]

            ld_pairs <- unique(
                ld_pairs,
                by = c("var1", "var2")
            )

            ##---------------------------------------
            ## Construct SNP pairs within phenotype
            ##---------------------------------------
            pair_info <- merge(
                cs_membership[,
                    .(
                        phenotype_id,
                        cs1 = cs_id,
                        SNP_A = variant_id
                    )
                ],
                cs_membership[,
                    .(
                        phenotype_id,
                        cs2 = cs_id,
                        SNP_B = variant_id
                    )
                ],
                by = "phenotype_id",
                allow.cartesian = TRUE
            )

            ## remove self-pairs
            pair_info <- pair_info[
                SNP_A != SNP_B
            ]

            ## orientation-independent key
            pair_info[, `:=`(
                var1 = pmin(SNP_A, SNP_B),
                var2 = pmax(SNP_A, SNP_B)
            )]

            ## remove duplicated A-B / B-A
            pair_info <- unique(
                pair_info,
                by = c(
                    "phenotype_id",
                    "cs1",
                    "cs2",
                    "var1",
                    "var2"
                )
            )

            ##---------------------------------------
            ## Add LD
            ##---------------------------------------
            tissue_res <- merge(
                pair_info,
                ld_pairs[,
                    .(
                        var1,
                        var2,
                        R2,
                        ld_source
                    )
                ],
                by = c("var1", "var2")
            )

            ##---------------------------------------
            ## Classify
            ##---------------------------------------
            tissue_res[,
                pair_type := fifelse(
                    cs1 == cs2,
                    "Intra-credible set",
                    "Inter-credible set"
                )
            ]

            tissue_res[,
                tissue := tmp_tissue
            ]

            tissue_res
        }
    ))

    ## Plot: intra- vs. inter-CS LD
    ggplot(susie_cs_ld, aes(x = pair_type, y = R2)) +
        rasterise(geom_boxplot(outlier.size = 1), dpi = 600) +
        labs(x = "", y = "Unphased r2") +
        theme_pub()

    ggsave(
        paste0(out_figure_dir, "ld_r2_among_cs.pdf"),
        width = 3,
        height = 4
    )

    susie_cs_ld |>
        group_by(pair_type) |>
        summarise(mean_r2 = mean(R2), .groups = "drop") ## Inter-credible set  0.0808 Intra-credible set  0.839

    ## 3. PLINK-based pairwise r² (called ONLY for multi-CS genes)

    # tmp_qtl_type <- "snv_eqtl"

    R2_THRESH <- 0.01

    plink_r2 <- function(
        vars,
        tissue_name
    ) {
        all_v <- vars
        nv <- length(all_v)
        mat <- diag(nv)
        dimnames(mat) <- list(all_v, all_v)

        sub_multi_ld <- tissue_ld_list[[tissue_name]][
            SNP_A %in% all_v & SNP_B %in% all_v
        ]

        if (nrow(sub_multi_ld) > 0L) {
            ia <- match(sub_multi_ld$SNP_A, all_v)
            ib <- match(sub_multi_ld$SNP_B, all_v)
            mat[cbind(ia, ib)] <- sub_multi_ld$R2
            mat[cbind(ib, ia)] <- sub_multi_ld$R2
        }

        avail <- intersect(vars, all_v)
        mat[avail, avail, drop = FALSE]
    }

    # Graph-based merging (connected components)
    count_loci <- function(cs_dt, r2_mat, thresh = R2_THRESH) {
        cs_dt <- as.data.frame(cs_dt)
        rownames(cs_dt) <- cs_dt$lead_var

        n_cs <- nrow(cs_dt)
        if (n_cs == 1L) {
            return(1L)
        }

        if (is.null(r2_mat)) {
            warning("r² matrix unavailable — treating all CS as independent")
            return(n_cs)
        }

        vars <- cs_dt$lead_var
        avail <- intersect(vars, rownames(r2_mat))
        missing <- setdiff(vars, avail)

        if (length(missing)) {
            message(sprintf(
                "    [!] %d lead variant(s) absent from LD reference — kept as independent",
                length(missing)
            ))
        }

        if (length(avail) < 2L) {
            return(n_cs)
        }

        sub <- r2_mat[avail, avail, drop = FALSE]
        adj <- sub > thresh
        diag(adj) <- FALSE

        rownames(adj) <- cs_dt[rownames(adj), "cs_id"]
        colnames(adj) <- cs_dt[colnames(adj), "cs_id"]

        g <- graph_from_adjacency_matrix(adj, mode = "undirected")
        comp <- components(g)

        out_cs <- vector()
        for (i in 1:comp$no) {
            out_cs <- c(
                out_cs,
                paste0(
                    unique(names(comp$membership)[comp$membership == i]),
                    collapse = ","
                )
            )
        }
        out_cs <- paste0(out_cs, collapse = ";")

        out_cs
    }

    # LD merging loop: only multi-CS tissue-gene pairs
    ld_multi <- multi_leads[,
        {
            message(sprintf(
                "  %-12s | %-22s | %d CS",
                tissue[1],
                ensg[1],
                length(unique(cs_id))
            ))

            r2_mat <- tryCatch(
                plink_r2(lead_var, tissue[1]),
                error = function(e) {
                    warning(e$message)
                    NULL
                }
            )

            indep_loci <- count_loci(.SD, r2_mat)
            n_indep <- length(strsplit(indep_loci, split = ";")[[1]])

            cat(r2_mat)
            message(sprintf("    => %s independent locus/loci", n_indep))

            .(indep_loci = indep_loci, n_after = n_indep)
        },
        by = .(tissue, phenotype_id, ensg, xqtl_type)
    ]

    ## Combine: multi-CS (LD-merged)
    original <- unique(
        qtl_susie_res[, .(
            tissue,
            phenotype_id,
            ensg,
            xqtl_type,
            n_before = cs_n
        )]
    )

    results <- merge(
        original[n_before > 1L],
        ld_multi,
        by = c("tissue", "phenotype_id", "ensg", "xqtl_type")
    )

    results[, reduced := n_before - n_after]

    cat(sprintf(
        " Before (non-overlapping CS)    : %d loci\n",
        sum(results$n_before)
    ))
    cat(sprintf(
        " After  (LD-merged, r\u00b2 > %.2f)  : %d loci\n",
        R2_THRESH,
        sum(results$n_after)
    ))
    cat(sprintf(
        " Collapsed                       : %d loci (%.1f%% reduction)\n",
        sum(results$reduced),
        100 * sum(results$reduced) / sum(results$n_before)
    ))

    fwrite(
        results,
        paste0(dirname(ld_tmp_dir), "/ld_adjusted_loci_results.tsv"),
        sep = "\t"
    )

    ## Per-tissue grouped bar
    dt_total <- data.table(
        Method = factor(
            c(
                "Before\n(Non-overlapping CS)",
                sprintf("After\n(r\u00b2 > %.2f merged)", R2_THRESH)
            ),
            levels = c(
                "Before\n(Non-overlapping CS)",
                sprintf("After\n(r\u00b2 > %.2f merged)", R2_THRESH)
            )
        ),
        n = c(sum(results$n_before), sum(results$n_after))
    )
    dt_total
    # (3622-2260)/3622

    tissue_tbl <- results[,
        .(
            before = sum(n_before),
            after = sum(n_after),
            collapsed = sum(reduced),
            pct_reduc = round(100 * sum(reduced) / sum(n_before), 1)
        ),
        keyby = .(xqtl_type, tissue)
    ]

    ts_long <- melt(
        tissue_tbl,
        id.vars = c("tissue", "xqtl_type"),
        measure.vars = c("before", "after"),
        variable.name = "Method",
        value.name = "n"
    )
    ts_long[,
        Method := fifelse(
            Method == "before",
            "Before (non-overlapping CS)",
            sprintf("After (r\u00b2 > %.2f merged)", R2_THRESH)
        )
    ]

    ts_long$Method <- factor(
        ts_long$Method,
        levels = c(
            "Before (non-overlapping CS)",
            sprintf("After (r\u00b2 > %.2f merged)", R2_THRESH)
        )
    )

    ggplot(ts_long, aes(reorder(tissue, -n), n, fill = Method)) +
        geom_col(position = position_dodge(0.82), width = 0.78) +
        # geom_text(
        #     aes(label = n),
        #     position = position_dodge(0.82),
        #     vjust = -0.35,
        #     size = 3
        # ) +
        scale_fill_manual(
            values = c(
                "Before (non-overlapping CS)" = "#46b1c9",
                setNames(
                    "#d94b4f",
                    sprintf("After (r\u00b2 > %.2f merged)", R2_THRESH)
                )
            )
        ) +
        # scale_y_continuous(expand = expansion(c(0, .20))) +
        labs(
            title = "Per-tissue comparison",
            x = "",
            y = "Number of independent loci",
            fill = NULL
        ) +
        theme_pub() +
        theme(
            axis.text.x = element_text(angle = 40, hjust = 1),
            legend.position = "top",
            legend.direction = "horizontal"
        )

    ggsave(
        paste0(out_figure_dir, "ld_adjusted_loci_comparison.pdf"),
        width = 7,
        height = 5,
    )

    ts_long_value <- reshape2::dcast(ts_long, tissue ~ Method, value.var = "n")
    100 -
        mean(
            ts_long_value$`After (r² > 0.01 merged)` /
                ts_long_value$`Before (non-overlapping CS)` *
                100
        ) # 36.43562
})
