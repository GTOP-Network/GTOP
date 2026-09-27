#!/usr/bin/env Rscript

# ==============================================================================
# fine-mapping QTL to independent loci
# ==============================================================================

#%% ------------------------ 0. prepare files (packages, input files, output files)
PROJECT_DIR <- "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/2026-05-11-fine_mapping_revision"
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(ggrastr)
    library(igraph)
})

# ------------ input
IN_SUSIE_FILE <- "output/data/fine_mapped_gtop_xqtl.txt"
IN_AF_FILE <- "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/2026-05-09-gtop_gnomad_af_revision/output/data/integrated_freq.txt"

PLINK_BIN <- "/media/bora_A/zhangt/src/bin/plink2"
GTOP_BFILE <- "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/xqtl_atlas/data/genotype/gtop/gtop_snv.maf05"
KGP_BFILE <- "/media/bora_A/zhangt/src/data/1000G/five_ancestry_groups/EAS/1000G.EAS.maf01"


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
snv_freq_res <- fread(IN_AF_FILE)
raw_susie_res <- fread(IN_SUSIE_FILE)

lapply(c("snv_eqtl", "snv_juqtl", "snv_tuqtl"), function(tmp_qtl_type) {
    # tmp_qtl_type <- "snv_eqtl"
    PATH_COLLAPSE_RES <- sprintf(
        "output/data/cs_combining/%s/ld_adjusted_loci_results.tsv",
        tmp_qtl_type
    )

    ## result
    susie_collapse <- fread(PATH_COLLAPSE_RES)

    susie_collapse <- susie_collapse |> separate_rows(indep_loci, sep = ";")
    susie_collapse$cs_id <- susie_collapse$indep_loci
    susie_collapse <- susie_collapse |> separate_rows(cs_id, sep = ",")
    setDT(susie_collapse)

    susie_res <- merge(
        raw_susie_res[xqtl_type == tmp_qtl_type],
        susie_collapse[, .(tissue, phenotype_id, xqtl_type, cs_id, indep_loci)],
        by = c("tissue", "phenotype_id", "xqtl_type", "cs_id"),
        all.x = TRUE
    )

    susie_res$indep_loci[is.na(susie_res$indep_loci)] <- susie_res$cs_id[is.na(
        susie_res$indep_loci
    )]

    # if (tmp_qtl_type == "snv_juqtl") {
    #     susie_res$phenotype_id <- sapply(
    #         strsplit(susie_res$phenotype_id, ":"),
    #         function(x) {
    #             paste0(x[c(1, 2, 3, 6)], collapse = ":")
    #         }
    #     )
    # }

    ## gene count
    susie_res[, gene_n := uniqueN(tissue), by = .(ensg, xqtl_type)]
    susie_res$indep_loci <- paste0(susie_res$tissue, "_", susie_res$indep_loci)

    used_df1 <- susie_res[gene_n == 1]
    used_df1 <- used_df1[, .(
        xqtl_type,
        tissue,
        phenotype_id,
        variant_id,
        indep_loci,
        pip,
        ensg,
        symbol
    )]

    gene_multi_leads <- susie_res[
        gene_n > 1L,
        {
            i <- which.max(pip)
            .(lead_var = variant_id[i], lead_pip = pip[i])
        },
        by = .(tissue, ensg, symbol, xqtl_type, indep_loci)
    ]

    gene_multi_leads <- gene_multi_leads |>
        group_by(ensg) |>
        mutate(lead_unique = length(unique(lead_var)))
    setDT(gene_multi_leads)
    gene_multi_leads1 <- gene_multi_leads[lead_unique > 1]

    used_df2 <- gene_multi_leads[lead_unique == 1]
    used_df2 <- used_df2 |>
        group_by(ensg, xqtl_type) |>
        summarise(loci_id = paste0(sort(unique(indep_loci)), collapse = ","))
    used_df2 <- merge(
        susie_res,
        used_df2,
        by = c(
            "ensg",
            "xqtl_type"
        )
    )
    used_df2 <- used_df2[, .(
        xqtl_type,
        tissue,
        phenotype_id,
        variant_id,
        indep_loci = loci_id,
        pip,
        ensg,
        symbol
    )]

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
            "--ld-window-kb 1000",
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

    run_1KGP_GTOP <- function(tmp_snps) {
        ld_tmp_dir <- sprintf(
            "output/data/cs_combining/%s/ld_tmp",
            tmp_qtl_type
        )

        tag <- sprintf(
            "%s_%d",
            "test_gene",
            length(tmp_snps)
        )
        snp_file <- file.path(
            ld_tmp_dir,
            paste0(tag, ".snplist")
        )

        writeLines(
            tmp_snps,
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
                "1KGP=%d pairs; GTOP supplement=%d pairs; total=%d",
                nrow(kgp_ld),
                nrow(gtop_ld),
                nrow(res)
            )
        )
        return(res)
    }

    multi_ld <- run_1KGP_GTOP(tmp_snps = unique(gene_multi_leads1$lead_var))

    ## R2 filteration
    R2_THRESH <- 0.01

    tissue_plink_r2 <- function(
        vars
    ) {
        all_v <- vars
        nv <- length(all_v)
        mat <- diag(nv)
        dimnames(mat) <- list(all_v, all_v)

        sub_multi_ld <- multi_ld[
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

    ## STEP 4 — Graph-based merging (connected components)
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

        rownames(adj) <- cs_dt[rownames(adj), "indep_loci"]
        colnames(adj) <- cs_dt[colnames(adj), "indep_loci"]

        g <- graph_from_adjacency_matrix(adj, mode = "undirected")
        comp <- components(g)

        out_cs <- vector()
        for (i in 1:comp$no) {
            out_cs <- c(
                out_cs,
                paste0(
                    unique(names(comp$membership)[comp$membership == i]),
                    collapse = "|"
                )
            )
        }
        out_cs <- paste0(out_cs, collapse = ";")

        out_cs
    }

    # ==============================================================================
    # [F] STEP 5 — LD merging loop: only multi-CS tissue-gene pairs
    # ==============================================================================
    gene_multi_leads1 <- gene_multi_leads1 |>
        group_by(ensg, symbol, xqtl_type, lead_var) |>
        summarise(indep_loci = paste0(sort(indep_loci), collapse = "+"))
    setDT(gene_multi_leads1)

    ld_multi <- gene_multi_leads1[,
        {
            message(sprintf(
                " %-22s | %d CS",
                symbol[1],
                length(unique(indep_loci))
            ))

            r2_mat <- tryCatch(
                tissue_plink_r2(lead_var),
                error = function(e) {
                    warning(e$message)
                    NULL
                }
            )

            tissue_indep_loci <- count_loci(.SD, r2_mat)
            n_indep <- length(strsplit(tissue_indep_loci, split = ";")[[1]])
            message(sprintf("    => %s independent locus/loci", n_indep))

            .(tissue_indep_loci = tissue_indep_loci, n_after = n_indep)
        },
        by = .(ensg, symbol, xqtl_type)
    ]

    res_multi <- ld_multi |> separate_rows(tissue_indep_loci, sep = ";")
    res_multi$loci_id <- res_multi$tissue_indep_loci
    res_multi <- res_multi |> separate_rows(tissue_indep_loci, sep = "\\|")
    res_multi <- res_multi |> separate_rows(tissue_indep_loci, sep = "\\+")
    setDT(res_multi)

    used_df3 <- merge(
        susie_res,
        res_multi[, .(
            ensg,
            symbol,
            xqtl_type,
            indep_loci = tissue_indep_loci,
            n_after,
            loci_id
        )],
        by = c(
            "symbol",
            "ensg",
            "xqtl_type",
            "indep_loci"
        )
    )
    used_df3 <- used_df3[, .(
        xqtl_type,
        tissue,
        phenotype_id,
        variant_id,
        indep_loci = loci_id,
        pip,
        ensg,
        symbol
    )]

    used_df <- rbind(used_df1, used_df2, used_df3)

    fwrite(
        used_df,
        sprintf(
            "output/data/cs_combining/%s/collapsed_cs_results.txt",
            tmp_qtl_type
        ),
        sep = "\t"
    )

    ## freq
    fm_freq_res <- merge(
        used_df,
        snv_freq_res[, .(
            variant_id,
            type_afr,
            type_eur,
            type_sas,
            type_amr,
            type_neas
        )],
        by = "variant_id"
    )

    fm_freq_res <- fm_freq_res |>
        group_by(ensg, indep_loci) |>
        mutate(sum_pip = sum(pip))

    setDT(fm_freq_res)

    onehot_reshape <- function(idata, indep_name, icol_name) {
        #
        col_indep <- enquo(indep_name)
        col_sym <- enquo(icol_name)
        col_name_str <- rlang::as_name(col_sym)
        new_col_name <- paste0("new_", col_name_str)
        # 1. one-hot
        onehot_data <- idata %>%
            dplyr::select(!!col_sym) %>%
            model.matrix(~ . - 1, data = .)

        oh_cols <- colnames(onehot_data)
        if (length(oh_cols) == 0) {
            stop("No one-hot encoded columns generated.")
        }

        # 2. weight
        pip_weights <- idata$pip
        onehot_weighted <- sweep(onehot_data, 1, pip_weights, `*`)

        # 3. merge meta data
        combined <- idata %>%
            dplyr::select(!!col_indep, !!col_sym) %>%
            as.data.frame() %>%
            cbind(onehot_weighted)

        # 4. indep_loci
        result <- combined %>%
            dplyr::group_by(.data[[col_indep]]) %>%
            dplyr::summarise(
                across(all_of(oh_cols), sum, .names = "{.col}")
            )

        # 5. max weighted value
        value_matrix <- as.matrix(result[oh_cols])
        if (!any(grepl("U$", colnames(value_matrix)))) {
            f_af_type <- apply(value_matrix, 1, function(x) {
                if (x[2] / sum(x[1:2]) > 0.9) {
                    return("R")
                } else {
                    return("C")
                }
            })
        } else {
            f_af_type <- apply(value_matrix, 1, function(x) {
                if (sum(x[2:3]) / sum(x) > 0.9) {
                    if (x[3] / sum(x[2:3]) > 0.9) {
                        return("U")
                    } else {
                        return("R")
                    }
                } else {
                    return("C")
                }
            })
        }

        table(f_af_type)
        setDT(result)
        result[[new_col_name]] <- f_af_type

        return(result)
    }

    fm_freq_res$gene_locus <- paste0(
        fm_freq_res$ensg,
        ";",
        fm_freq_res$indep_loci
    )

    fm_type1 <- onehot_reshape(
        idata = fm_freq_res,
        indep_name = "gene_locus",
        icol_name = "type_afr"
    )
    fm_type2 <- onehot_reshape(
        idata = fm_freq_res,
        indep_name = "gene_locus",
        icol_name = "type_eur"
    )
    fm_type3 <- onehot_reshape(
        idata = fm_freq_res,
        indep_name = "gene_locus",
        icol_name = "type_sas"
    )
    fm_type4 <- onehot_reshape(
        idata = fm_freq_res,
        indep_name = "gene_locus",
        icol_name = "type_amr"
    )
    fm_type5 <- onehot_reshape(
        idata = fm_freq_res,
        indep_name = "gene_locus",
        icol_name = "type_neas"
    )

    fm_type <- cbind(
        fm_type1[, .SD, .SDcols = c(1, ncol(fm_type1))],
        fm_type2[, .SD, .SDcols = ncol(fm_type2)],
        fm_type3[, .SD, .SDcols = ncol(fm_type3)],
        fm_type4[, .SD, .SDcols = ncol(fm_type4)],
        fm_type5[, .SD, .SDcols = ncol(fm_type5)]
    )
    dir.create("output/data/fm_freq")
    fwrite(
        fm_type,
        sprintf("output/data/fm_freq/%s_res.txt", tmp_qtl_type),
        sep = "\t"
    )

    fm_type_df <- unique(fm_type[,
        .(
            gene_locus,
            afr = new_type_afr,
            eur = new_type_eur,
            sas = new_type_sas,
            amr = new_type_amr,
            neas = new_type_neas
        )
    ])

    fm_type_df[,
        (2:6) := lapply(.SD, function(x) {
            ifelse(
                x == "U",
                0,
                ifelse(x == "R", 0.01, ifelse(x == "C", 0.1, x))
            )
        }),
        .SDcols = 2:6
    ]

    table(rowSums(fm_type_df[, 2:6] < 0.1) > 0)
    table(rowSums(fm_type_df[, 2:5] < 0.1) > 0)

    fm_rare_type_df <- fm_type_df[rowSums(fm_type_df[, 2:5] < 0.1) > 0, ]

    fwrite(
        fm_rare_type_df,
        file = sprintf("output/data/fm_freq/%s_fm_rare.txt", tmp_qtl_type),
        sep = "\t"
    )
})


fd_eQTL <- fread(
    "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/2026-05-11-fine_mapping_revision/output/data/fm_freq/snv_eqtl_fm_rare.txt"
)
table(rowSums(fd_eQTL[, 2:6] < 0.1) > 0)
table(rowSums(fd_eQTL[, 2:5] < 0.1) > 0)
table(rowSums(fd_eQTL[, 2:5] < 0.1) == 4)
table(rowSums(fd_eQTL[, 2:5] < 0.01) == 4)

tmp_data <- fd_eQTL[, .(
    SNP = gene_locus,
    AFR = afr,
    EUR = eur,
    SAS = sas,
    AMR = amr
)]

fwrite(
    tmp_data,
    "/media/london_A/mengxin/GTOP_code/fig-4/input/Fig4d.txt",
    sep = "\t"
)

library(reticulate)

py_require("matplotlib==3.6.3")
py_require("git+https://github.com/aabiddanda/geovar")

py_config()

py_run_string(
    '
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ---- compatibility patches for geovar ----

# New NumPy removed np.row_stack
if not hasattr(np, "row_stack"):
    np.row_stack = np.vstack

# New Matplotlib removed matplotlib.cm.get_cmap
if not hasattr(cm, "get_cmap"):
    def get_cmap(name=None, lut=None):
        cmap = plt.colormaps[name if name is not None else "viridis"]
        if lut is not None:
            cmap = cmap.resampled(lut)
        return cmap
    cm.get_cmap = get_cmap

# Import geovar after compatibility patches
from geovar import *

plt.rcParams["pdf.fonttype"] = 42

geovar_test = GeoVar()
geovar_test.add_freq_mat("/media/london_A/mengxin/GTOP_code/fig-4/input/Fig4d.txt")
geovar_test.geovar_binning()

geovar_plot = GeoVarPlot()
geovar_plot.add_data_geovar(geovar_test)
geovar_plot.filter_data()
geovar_plot.add_cmap()

fig, ax = plt.subplots(1, 1, figsize=(3.5, 6))
geovar_plot.plot_geovar(ax)

ax.set_xticklabels(geovar_plot.poplist)

plt.savefig(
    "output/fd_eQTL_freq_data.pdf",
    dpi=600,
    bbox_inches="tight"
)

plt.close()
'
)
## tissue gene
tissue_gene_info <- reshape2::melt(tmp_data)
tissue_gene_info <- tissue_gene_info[tissue_gene_info$value < 0.1, ]
tissue_gene_info$gene <- gsub(";.+", "", tissue_gene_info$SNP)
tissue_gene_info$tissue <- gsub(".+;", "", tissue_gene_info$SNP)
tissue_gene_info <- tissue_gene_info |> separate_rows(tissue, sep = "\\+|\\|")
tissue_gene_info <- tissue_gene_info |> separate_rows(tissue, sep = ",")
tissue_gene_info <- tissue_gene_info[grepl("_", tissue_gene_info$tissue), ]
tissue_gene_info$tissue <- gsub("_L.+", "", tissue_gene_info$tissue)

tissue_gene_summary <- tissue_gene_info |>
    group_by(variable, tissue) |>
    summarise(count = length(unique(gene)))

tissue_gene_summary |> group_by(variable) |> summarise(mean = mean(count))

color_df <- fread(
    "/media/pacific/share/Datasets/Asian_GTEx/Metainfo/input/GMTiP_tissue_code_and_colors.csv"
)
color_vec <- paste0("#", color_df$Tissue_Color_Code)
names(color_vec) <- color_df$Tissue
tissue_gene_summary$variable

ggplot(
    tissue_gene_summary,
    aes(x = variable, y = count)
) +
    geom_boxplot(outlier.color = "NA", fill = "#e5e5e4") +
    geom_point(
        aes(group = tissue, color = tissue),
        position = position_dodge(width = 0.4)
    ) +
    scale_color_manual(values = color_vec) +
    theme_pub() +
    labs(x = "", y = "Gene number") +
    scale_fill_manual(values = c("#6874b4", "#4bb9b9"))

ggsave("output/figures/fm_freq/gene_count.pdf", width = 6, height = 4)

fwrite(
    tissue_gene_summary,
    "/media/london_A/mengxin/GTOP_code/fig-4/input/Fig4e.txt",
    sep = "\t"
)


## fd_sQTLs
fd_juQTL <- fread(
    "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/2026-05-11-fine_mapping_revision/output/data/fm_freq/snv_juqtl_fm_rare.txt"
)
fd_tuQTL <- fread(
    "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/2026-05-11-fine_mapping_revision/output/data/fm_freq/snv_tuqtl_fm_rare.txt"
)
fd_juQTL_gene <- gsub(";.+", "", fd_juQTL$gene_locus)
fd_tuQTL_gene <- gsub(";.+", "", fd_tuQTL$gene_locus)
fd_tuQTL_gene <- !(fd_tuQTL_gene %in% fd_juQTL_gene)
table(fd_tuQTL_gene)

fd_sQTL <- rbind(fd_juQTL, fd_tuQTL[fd_tuQTL_gene, ])

table(rowSums(fd_sQTL[, 2:5] < 0.1) > 0)
table(rowSums(fd_sQTL[, 2:5] < 0.1) == 4)
table(rowSums(fd_sQTL[, 2:5] < 0.01) == 4)

table(fd_sQTL$neas)
282 / (282 + 9277)

tmp_fd_sQTL <- fd_sQTL[, .(
    SNP = gene_locus,
    AFR = afr,
    EUR = eur,
    SAS = sas,
    AMR = amr
)]
fwrite(
    tmp_fd_sQTL,
    "/media/london_A/mengxin/GTOP_code/extend/extend_7/input/ext_Fig7a.txt",
    sep = "\t"
)


# library(tidyverse)

# # 输入示例：test.freq.csv
# # 至少应包含：AFR, AMR, EAS, EUR, SAS 五列，值为 0–1 的等位基因频率
# freq <- tmp_fd_sQTL #readr::read_csv("test.freq.csv", show_col_types = FALSE)

# # 列的显示顺序；与网页示例一致
# pops <- c("AFR", "EUR", "SAS", "AMR")

# # 0 -> U（未观察到）
# # (0, 0.05] -> R（rare）
# # (0.05, 1] -> C（common）
# to_geovar <- function(x, rare_cutoff = 0.05) {
#     case_when(
#         is.na(x) | x == 0 ~ "U",
#         x <= rare_cutoff ~ "R",
#         TRUE ~ "C"
#     )
# }

# # 统计每一种跨群体频率模式的比例
# plot_df <- freq |>
#     transmute(across(all_of(pops), to_geovar)) |>
#     unite("pattern", all_of(pops), sep = "") |>
#     count(pattern, name = "n") |>
#     mutate(
#         pct = n / sum(n),
#         # 按比例从大到小排列；如要复制不同排序规则，可改这里
#         pattern = fct_reorder(pattern, pct, .desc = TRUE),
#         ymin = cumsum(lag(pct, default = 0)),
#         ymax = cumsum(pct),
#         ymid = (ymin + ymax) / 2,
#         row = row_number()
#     ) |>
#     separate(pattern, into = pops, sep = c(1, 2, 3, 4), remove = FALSE) |>
#     pivot_longer(all_of(pops), names_to = "population", values_to = "class") |>
#     mutate(
#         population = factor(population, levels = pops),
#         x = as.numeric(population)
#     )

# # 可选：仅保留占比至少 0.1% 的模式，避免图过长
# # plot_df <- plot_df |> group_by(pattern) |> filter(first(pct) >= 0.001) |> ungroup()

# ggplot(plot_df) +
#     geom_tile(
#         aes(x = x, y = ymid, width = 1, height = pct, fill = class),
#         colour = "white",
#         linewidth = 0.35
#     ) +
#     geom_text(
#         aes(x = x, y = ymid, label = class, colour = class == "C"),
#         size = 3.2,
#         fontface = "bold"
#     ) +
#     scale_x_continuous(
#         breaks = seq_along(pops),
#         labels = pops,
#         expand = c(0, 0)
#     ) +
#     scale_y_continuous(
#         labels = scales::percent_format(accuracy = 1),
#         expand = c(0, 0)
#     ) +
#     scale_fill_manual(
#         values = c(U = "white", R = "#8CC5E8", C = "#0072B2"),
#         name = "Allele-frequency class"
#     ) +
#     scale_colour_manual(
#         values = c(`TRUE` = "white", `FALSE` = "black"),
#         guide = "none"
#     ) +
#     coord_cartesian(xlim = c(0.5, length(pops) + 0.5), clip = "off") +
#     labs(x = NULL, y = "Variants (%)") +
#     theme_minimal(base_size = 12) +
#     theme(
#         panel.grid = element_blank(),
#         axis.text.x = element_text(face = "bold"),
#         legend.position = "none",
#         plot.margin = margin(5.5, 5.5, 5.5, 5.5)
#     )

library(reticulate)
py_require("matplotlib==3.6.3")
py_require("git+https://github.com/aabiddanda/geovar")
py_config()

py_run_string(
    '
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# New NumPy removed np.row_stack
if not hasattr(np, "row_stack"):
    np.row_stack = np.vstack

# New Matplotlib removed matplotlib.cm.get_cmap
if not hasattr(cm, "get_cmap"):
    def get_cmap(name=None, lut=None):
        cmap = plt.colormaps[name if name is not None else "viridis"]
        if lut is not None:
            cmap = cmap.resampled(lut)
        return cmap
    cm.get_cmap = get_cmap

# Import geovar after compatibility patches
from geovar import *

plt.rcParams["pdf.fonttype"] = 42

geovar_test = GeoVar()
geovar_test.add_freq_mat("/media/london_A/mengxin/GTOP_code/extend/extend_7/input/ext_Fig7a.txt")
geovar_test.geovar_binning()

geovar_plot = GeoVarPlot()
geovar_plot.add_data_geovar(geovar_test)
geovar_plot.filter_data()
geovar_plot.add_cmap()

fig, ax = plt.subplots(1, 1, figsize=(3, 6))
geovar_plot.plot_geovar(ax)

ax.set_xticklabels(geovar_plot.poplist)

plt.savefig(
    "output/fd_sQTL_freq_data.pdf",
    dpi=600,
    bbox_inches="tight"
)

plt.close()
'
)

# CHRPOS_RSID <- fread(
#     "2025-06-11-specific_xQTL/input/GTOP/SNV/GTOP_chrpos_rsid.txt",
#     header = FALSE,
#     col.names = c("chr_pos_ref_alt", "variant_id")
# )

# ### plot
# tmp_data <- fread(
#     "/media/bora_A/zhangt/2026-05-07-gtop_xqtl-Project/2026-05-11-fine_mapping_revision/output/data/fm_freq/snv_eqtl_fm_rare.txt"
# )

# ## fm_res
# eQTL_fm_res <- fread(sprintf(
#     "2025-06-11-specific_xQTL/output/data/%s/finemapping_merge_variants.txt",
#     "eQTL"
# ))

# eQTL_fm_proxy <- eQTL_fm_res %>%
#     group_by(ensg) %>%
#     summarise(cs_num = length(unique(variant_proxy)))
# eQTL_fm_proxy$cs_num[eQTL_fm_proxy$cs_num >= 7] <- 7
# eQTL_fm_proxy <- eQTL_fm_proxy %>%
#     group_by(cs_num) %>%
#     summarise(phenotype_count = length(unique(ensg)))

# sQTL_fm_res <- fread(sprintf(
#     "2025-06-11-specific_xQTL/output/data/%s/finemapping_merge_variants.txt",
#     "sQTL"
# ))

# sQTL_fm_proxy <- sQTL_fm_res %>%
#     group_by(ensg) %>%
#     summarise(cs_num = length(unique(variant_proxy)))
# sQTL_fm_proxy$cs_num[sQTL_fm_proxy$cs_num >= 7] <- 7
# sQTL_fm_proxy <- sQTL_fm_proxy %>%
#     group_by(cs_num) %>%
#     summarise(phenotype_count = length(unique(ensg)))

# fm_proxy_count <- rbind(
#     eQTL_fm_proxy %>% mutate("xqtl_type" = "eQTL"),
#     sQTL_fm_proxy %>% mutate("xqtl_type" = "sQTL")
# )

# ggplot(fm_proxy_count, aes(x = cs_num, y = phenotype_count, fill = xqtl_type)) +
#     geom_col(position = "dodge") +
#     scale_x_continuous(breaks = 1:7, labels = c(1:6, ">=7")) +
#     theme_classic() +
#     labs(x = "Number of independent QTLs", y = "Gene number") +
#     scale_fill_manual(values = c("eQTL" = "#9ac294", "sQTL" = "#7b86a7"))

# ggsave(
#     "2025-06-11-specific_xQTL/output/result/Gene_number_of_independent_QTL.pdf",
#     width = 4,
#     height = 3
# )

# fwrite(
#     fm_proxy_count,
#     "/media/london_A/mengxin/GTOP_code/extend/extend_6/input/Extended_fig6a_data.txt",
#     sep = "\t"
# )

# ## eQTL for EAS-specific SNVs
# SNV_eQTL_fm_freq <- fread(
#     "2025-06-11-specific_xQTL/output/data/eQTL/finemapping_SNV_add_four_freq.txt"
# )
# length(unique(SNV_eQTL_fm_freq$variant_proxy)) # 11,121

# SNV_eQTL_fm_freq_lead <- unique(SNV_eQTL_fm_freq[
#     type_afr == new_type_afr &
#         type_nfe == new_type_nfe &
#         type_sas == new_type_sas &
#         type_amr == new_type_amr
# ])
# length(unique(SNV_eQTL_fm_freq_lead$variant_proxy)) # 11,105

# SNV_eQTL_fm_freq_lead <- SNV_eQTL_fm_freq_lead %>%
#     group_by(variant_proxy, ensg) %>%
#     mutate(
#         is_same = variant_id == variant_proxy
#     ) %>%
#     dplyr::slice(which.max(ifelse(is_same, 1, -pip))) %>%
#     dplyr::select(-is_same)
# setDT(SNV_eQTL_fm_freq_lead)

# SNV_eQTL_fm_freq_lead <- unique(SNV_eQTL_fm_freq_lead[,
#     .(
#         SNP = variant_id,
#         afr = new_type_afr,
#         nfe = new_type_nfe,
#         sas = new_type_sas,
#         amr = new_type_amr
#     )
# ])
# length(unique(SNV_eQTL_fm_freq_lead$SNP))
# nrow(SNV_eQTL_fm_freq_lead) # 11,057

# SNV_eQTL_fm_freq_lead[,
#     (2:5) := lapply(.SD, function(x) {
#         ifelse(x == "U", 0, ifelse(x == "R", 0.01, ifelse(x == "C", 0.1, x)))
#     }),
#     .SDcols = 2:5
# ]

# fwrite(
#     SNV_eQTL_fm_freq_lead,
#     file = "/media/london_A/mengxin/GTOP_code/extend/extend_6/input/Extended_fig6b_data.txt",
#     sep = "\t"
# )

# ### plot
# library(reticulate)

# if (!py_module_available("matplotlib")) {
#     py_install("matplotlib")
# }

# if (!py_module_available("geovar")) {
#     py_install("git+https://github.com/aabiddanda/geovar")
# }

# py_run_string(
#     '
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# # import pkg_resources
# from geovar import *

# plt.rcParams[\'pdf.fonttype\'] = 42

# geovar_test = GeoVar()
# geovar_test.add_freq_mat("/media/london_A/mengxin/GTOP_code/extend/extend_6/input/Extended_fig6b_data.txt")
# geovar_test.geovar_binning()

# geovar_plot = GeoVarPlot()
# geovar_plot.add_data_geovar(geovar_test)
# geovar_plot.filter_data()
# geovar_plot.add_cmap()

# fig, ax = plt.subplots(1,1,figsize=(3,6))
# geovar_plot.plot_geovar(ax)
# ax.set_xticklabels(geovar_plot.poplist)
# plt.savefig("./ExtendedFig6b.pdf", dpi=300, bbox_inches=\'tight\')
# '
# )

# ## finemapping_add_freq_sQTL
# SNV_sQTL_fm_freq <- fread(
#     "2025-06-11-specific_xQTL/output/data/sQTL/finemapping_SNV_add_four_freq.txt"
# )
# length(unique(SNV_sQTL_fm_freq$variant_proxy)) # 11760

# SNV_sQTL_fm_freq_lead <- unique(SNV_sQTL_fm_freq[
#     type_afr == new_type_afr &
#         type_nfe == new_type_nfe &
#         type_sas == new_type_sas &
#         type_amr == new_type_amr
# ])
# length(unique(SNV_sQTL_fm_freq_lead$variant_proxy))

# SNV_sQTL_fm_freq_lead <- SNV_sQTL_fm_freq_lead %>%
#     group_by(variant_proxy, ensg) %>%
#     mutate(
#         is_same = variant_id == variant_proxy
#     ) %>%
#     dplyr::slice(which.max(ifelse(is_same, 1, -pip))) %>%
#     dplyr::select(-is_same)
# setDT(SNV_sQTL_fm_freq_lead)

# SNV_sQTL_fm_freq_lead <- unique(SNV_sQTL_fm_freq_lead[,
#     .(
#         SNP = variant_id,
#         afr = new_type_afr,
#         nfe = new_type_nfe,
#         sas = new_type_sas,
#         amr = new_type_amr
#     )
# ])
# length(unique(SNV_sQTL_fm_freq_lead$SNP))
# nrow(SNV_sQTL_fm_freq_lead) # 11694

# SNV_sQTL_fm_freq_lead[,
#     (2:5) := lapply(.SD, function(x) {
#         ifelse(x == "U", 0, ifelse(x == "R", 0.01, ifelse(x == "C", 0.1, x)))
#     }),
#     .SDcols = 2:5
# ]

# fwrite(
#     SNV_sQTL_fm_freq_lead,
#     file = "/media/london_A/mengxin/GTOP_code/extend/extend_6/input/Extended_fig6c_data.txt",
#     sep = "\t"
# )

# py_run_string(
#     '
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# # import pkg_resources
# from geovar import *

# plt.rcParams[\'pdf.fonttype\'] = 42

# geovar_test = GeoVar()
# geovar_test.add_freq_mat("/media/london_A/mengxin/GTOP_code/extend/extend_6/input/Extended_fig6c_data.txt")
# geovar_test.geovar_binning()

# geovar_plot = GeoVarPlot()
# geovar_plot.add_data_geovar(geovar_test)
# geovar_plot.filter_data()
# geovar_plot.add_cmap()

# fig, ax = plt.subplots(1,1,figsize=(3,6))
# geovar_plot.plot_geovar(ax)
# ax.set_xticklabels(geovar_plot.poplist)
# plt.savefig("./ExtendedFig6c.pdf", dpi=300, bbox_inches=\'tight\')
# '
# )

# ## only EAS specific fine-mapping eQTLs
# SNV_eQTL_fm_freq_lead1 <- unique(SNV_eQTL_fm_freq[
#     type_afr == new_type_afr &
#         type_nfe == new_type_nfe &
#         type_sas == new_type_sas &
#         type_amr == new_type_amr
# ])
# length(unique(SNV_eQTL_fm_freq_lead1$variant_proxy)) # 11105

# SNV_eQTL_fm_freq_lead1 <- SNV_eQTL_fm_freq_lead1 %>%
#     group_by(variant_proxy, ensg) %>%
#     mutate(
#         is_same = variant_id == variant_proxy
#     ) %>%
#     dplyr::slice(which.max(ifelse(is_same, 1, -pip))) %>%
#     dplyr::select(-is_same)
# setDT(SNV_eQTL_fm_freq_lead1)

# SNV_eQTL_fm_freq_lead1 <- unique(SNV_eQTL_fm_freq_lead1[,
#     .(
#         SNP = variant_id,
#         afr = new_type_afr,
#         nfe = new_type_nfe,
#         sas = new_type_sas,
#         amr = new_type_amr,
#         neas = new_type_neas
#     )
# ])
# length(unique(SNV_eQTL_fm_freq_lead1$SNP))
# dim(SNV_eQTL_fm_freq_lead1) # 11057

# SNV_eQTL_fm_freq_lead2 <- SNV_eQTL_fm_freq_lead1[
#     apply(SNV_eQTL_fm_freq_lead1, 1, function(x) {
#         any(x[2:5] != "C")
#     }),
# ]
# dim(SNV_eQTL_fm_freq_lead2)

# 1150 / 11057 * 100
# table(apply(SNV_eQTL_fm_freq_lead1, 1, function(x) {
#     all(x[2:5] %in% c("U", "R"))
# }))
# 198 / 11057 * 100

# SNV_eQTL_fm_freq_lead3 <- SNV_eQTL_fm_freq_lead1[
#     apply(SNV_eQTL_fm_freq_lead1, 1, function(x) {
#         any(x[6] != "C")
#     }),
# ]
# dim(SNV_eQTL_fm_freq_lead3)

# 306 / 11057 * 100

# SNV_eQTL_fm_freq_lead1[,
#     (2:5) := lapply(.SD, function(x) {
#         ifelse(x == "U", 0, ifelse(x == "R", 0.01, ifelse(x == "C", 0.1, x)))
#     }),
#     .SDcols = 2:5
# ]

# fwrite(
#     SNV_eQTL_fm_freq_lead1,
#     file = "/media/london_A/mengxin/GTOP_code/fig-4/input/Fig4a.txt",
#     sep = "\t"
# )

# py_run_string(
#     '
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# # import pkg_resources
# from geovar import *

# plt.rcParams[\'pdf.fonttype\'] = 42

# geovar_test = GeoVar()
# geovar_test.add_freq_mat("/media/london_A/mengxin/GTOP_code/fig-4/input/Fig4a.txt")
# geovar_test.geovar_binning()

# geovar_plot = GeoVarPlot()
# geovar_plot.add_data_geovar(geovar_test)
# geovar_plot.filter_data()
# geovar_plot.add_cmap()

# fig, ax = plt.subplots(1,1,figsize=(3,6))
# geovar_plot.plot_geovar(ax)
# ax.set_xticklabels(geovar_plot.poplist)
# plt.savefig("./Fig4a.pdf", dpi=300, bbox_inches=\'tight\')
# '
# )

# SNV_sQTL_fm_freq_lead1 <- unique(SNV_sQTL_fm_freq[
#     type_afr == new_type_afr &
#         type_nfe == new_type_nfe &
#         type_sas == new_type_sas &
#         type_amr == new_type_amr
# ])
# length(unique(SNV_sQTL_fm_freq_lead1$variant_proxy))

# SNV_sQTL_fm_freq_lead1 <- SNV_sQTL_fm_freq_lead1 %>%
#     group_by(variant_proxy, ensg) %>%
#     mutate(
#         is_same = variant_id == variant_proxy
#     ) %>%
#     dplyr::slice(which.max(ifelse(is_same, 1, -pip))) %>%
#     dplyr::select(-is_same)
# setDT(SNV_sQTL_fm_freq_lead1)

# SNV_sQTL_fm_freq_lead1 <- unique(SNV_sQTL_fm_freq_lead1[,
#     .(
#         SNP = variant_id,
#         afr = new_type_afr,
#         nfe = new_type_nfe,
#         sas = new_type_sas,
#         amr = new_type_amr,
#         neas = new_type_neas
#     )
# ])
# length(unique(SNV_sQTL_fm_freq_lead1$SNP))
# dim(SNV_sQTL_fm_freq_lead1)

# SNV_sQTL_fm_freq_lead2 <- SNV_sQTL_fm_freq_lead1[
#     apply(SNV_sQTL_fm_freq_lead1, 1, function(x) {
#         any(x[2:5] != "C")
#     }),
# ]
# dim(SNV_sQTL_fm_freq_lead2)

# 1616 / 11694 * 100
# table(apply(SNV_sQTL_fm_freq_lead1, 1, function(x) {
#     all(x[2:5] %in% c("U", "R"))
# }))
# 377 / 11694 * 100

# SNV_sQTL_fm_freq_lead3 <- SNV_sQTL_fm_freq_lead1[
#     apply(SNV_sQTL_fm_freq_lead1, 1, function(x) {
#         any(x[6] != "C")
#     }),
# ]
# dim(SNV_sQTL_fm_freq_lead3)

# 534 / 11694 * 100

# apply(SNV_eQTL_fm_freq_lead1[, 2:5], 2, function(x) {
#     sum(x %in% c("R", "U"))
# })
# # afr nfe sas amr
# # 815 703 380 382
# apply(SNV_eQTL_fm_freq_lead1[, 2:5], 2, function(x) {
#     sum(x %in% c("R", "U")) / length(x)
# })
# apply(SNV_sQTL_fm_freq_lead1[, 2:5], 2, function(x) {
#     sum(x %in% c("R", "U")) / length(x)
# })

# ## eQTL information ------------------
# eQTL_fm_res <- fread(
#     sprintf(
#         "2025-06-11-specific_xQTL/output/data/%s/finemapping_merge_variants.txt",
#         "eQTL"
#     ),
#     sep = "\t"
# )

# eQTL_fm_tmp <- merge(
#     eQTL_fm_res,
#     GTOP_freq[, .(variant_id, chr_pos_ref_alt, compare1KG)],
#     by = "variant_id"
# )
# eQTL_fm_res1 <- onehot_reshape(idata = eQTL_fm_tmp, icol_name = "compare1KG")

# table(eQTL_fm_res1$compare1KG)
# table(eQTL_fm_res1$new_compare1KG)

# eQTL_fm_res1 <- unique(eQTL_fm_res1[compare1KG == new_compare1KG])
# length(unique(eQTL_fm_res1$variant_proxy)) # 11120

# eQTL_fm_res1 <- eQTL_fm_res1 %>%
#     group_by(variant_proxy, ensg) %>%
#     mutate(
#         is_same = variant_id == variant_proxy
#     ) %>%
#     dplyr::slice(which.max(ifelse(is_same, 1, -pip))) %>%
#     dplyr::select(-is_same)
# setDT(eQTL_fm_res1)

# eQTL_fm_res1 <- unique(eQTL_fm_res1[, .(
#     SNP = variant_id,
#     compare = new_compare1KG
# )])
# length(unique(SNV_eQTL_fm_freq_lead$SNP))
# nrow(SNV_eQTL_fm_freq_lead) # 11,057

# table(eQTL_fm_res1$compare)
# # No   Yes
# # 1107 10457

# a <- fread("")
