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

# ------------ input
IN_SUSIE_FILE <- "output/data/fine_mapped_gtop_xqtl.txt"
IN_AF_FILE <- "path/to/output/data/integrated_freq.txt"

PLINK_BIN <- "/path/to/src/bin/plink2"
GTOP_BFILE <- "/path/to/data/genotype/gtop/gtop_snv.maf05"
KGP_BFILE <- "/path/to/src/data/1000G/five_ancestry_groups/EAS/1000G.EAS.maf01"


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
    "/path/to/output/data/fm_freq/snv_eqtl_fm_rare.txt"
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
geovar_test.add_freq_mat("/path/to/GTOP_code/fig-4/input/Fig4d.txt")
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
    "/path/to/GMTiP_tissue_code_and_colors.csv"
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


## fd_sQTLs
fd_juQTL <- fread(
    "/path/to/output/data/fm_freq/snv_juqtl_fm_rare.txt"
)
fd_tuQTL <- fread(
    "/path/to/output/data/fm_freq/snv_tuqtl_fm_rare.txt"
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
geovar_test.add_freq_mat("/path/to/GTOP_code/extend/extend_7/input/ext_Fig7a.txt")
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


