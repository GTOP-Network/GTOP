#!/usr/bin/env Rscript

# ==============================================================================
# eQTL protablility
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

#%% ------------------------ 1. Load QTL
tissue_link <- fread("input/gtop_gtex_tissues.txt")

portability_info_list <- lapply(
    1:nrow(tissue_link),
    function(row_index) {
        gtop_tissue <- tissue_link$GTOP_tissue[row_index]
        tmp_data <- readRDS(sprintf(
            "output/data/snv_eqtl/portability/%s.rds",
            gtop_tissue
        ))
        return(tmp_data)
    }
)
names(portability_info_list) <- tissue_link$GTOP_tissue

portable_summary <- rbindlist(lapply(portability_info_list, function(x) {
    x[[2]]
}))

summary_df1 <- portable_summary |>
    group_by(tissue, type2, type3) |>
    mutate(all_count = sum(count)) |>
    mutate(ratio = count / all_count)

summary_df1$type2 <- factor(
    summary_df1$type2,
    levels = c(
        "ratio_nominal",
        "ratio_mash",
        "ratio_uncertainty-aware",
        "FDR",
        "gene"
    )
)
summary_df1$type3 <- factor(
    summary_df1$type3,
    levels = c("Lead eVariant", "All eVariants", "All eGenes", "coloc genes")
)
summary_df1$tissue <- factor(
    summary_df1$tissue,
    levels = tissue_link$GTOP_tissue
)

p1 <- ggplot(
    summary_df1[summary_df1$type3 == "Lead eVariant", ],
    aes(type1, count / 1000, fill = type1)
) +
    geom_col(width = 0.8) +
    facet_grid(
        tissue ~ type2,
        scales = "free_y",
        axes = "all_y",
        axis.labels = "margins"
    ) +
    theme_pub() +
    scale_fill_manual(values = c("#ddbf93", "#af5543")) +
    labs(x = "", y = "Number of eVariants")

p2 <- ggplot(
    summary_df1[summary_df1$type3 != "Lead eVariant", ],
    aes(type1, count / 1000, fill = type1)
) +
    geom_col(width = 0.8) +
    facet_grid(
        tissue ~ type3,
        scales = "free_y",
        axes = "all_y",
        axis.labels = "margins"
    ) +
    theme_pub() +
    scale_fill_manual(values = c("#ddbf93", "#af5543")) +
    labs(x = "", y = "Number of eGenes")

p1 + p2 + plot_layout(guides = "collect", widths = c(2, 1))

ggsave(
    "output/figures/eqtl_portability/portability_three_metrics.pdf",
    width = 9,
    height = 10
)

saveRDS(
    summary_df1,
    "/path/to/GTOP_code/supp/supp_fig34.Portability_three_metrics/input/Fig34_data.rds"
)

summary_df1[summary_df1$type1 == "consistent", ] %>%
    group_by(type2, type3) %>%
    summarise(mean_ratio = mean(ratio, na.rm = T))


p3 <- ggplot(
    summary_df1[
        summary_df1$type3 == "Lead eVariant",
    ],
    aes(type2, ratio, fill = type1)
) +
    stat_summary(
        fun = mean,
        geom = "bar",
        position = position_dodge(width = 0.7),
        width = 0.6
    ) +
    stat_summary(
        fun.data = mean_se,
        geom = "errorbar",
        position = position_dodge(width = 0.7),
        width = 0.2,
        linewidth = 0.6
    ) +
    lims(y = c(0, 1)) +
    theme_pub() +
    scale_fill_manual(
        values = c("#e9bd76", "#3e71a6")
    ) +
    labs(x = "", y = "Number of eVariants")


p4 <- ggplot(
    summary_df1[
        summary_df1$type3 != "Lead eVariant",
    ],
    aes(type3, ratio, fill = type1)
) +
    stat_summary(
        fun = mean,
        geom = "bar",
        position = position_dodge(width = 0.7),
        width = 0.6
    ) +
    stat_summary(
        fun.data = mean_se,
        geom = "errorbar",
        position = position_dodge(width = 0.7),
        width = 0.2,
        linewidth = 0.6
    ) +
    lims(y = c(0, 1)) +
    theme_pub() +
    scale_fill_manual(
        values = c("#e9bd76", "#3e71a6")
    ) +
    labs(x = "", y = "Number of eGenes")

p3 +
    p4 +
    plot_layout(
        guides = "collect",
        widths = c(2, 1)
    )

ggsave(
    "output/figures/eqtl_portability/portability_three_metrics_boxplot1.pdf",
    width = 10,
    height = 3
)
fwrite(summary_df1, "/media/london_A/mengxin/GTOP_code/fig-4/input/Fig4g.txt")


#%% ------------------------ 2. eQTL portability (GTOP)
calc_var_x <- function(maf) {
    2 * maf * (1 - maf)
}

calc_minus_log10p <- function(beta, se) {
    z <- beta / se
    -log10(2 * pnorm(-abs(z)))
}

add_classification <- function(discovery_sig, replication_sig) {
    case_when(
        discovery_sig & replication_sig ~ "Both Significant",
        discovery_sig & !replication_sig ~ "Discovery Only",
        !discovery_sig & replication_sig ~ "Prediction Only",
        TRUE ~ "Neither Significant"
    )
}

calc_accuracy <- function(pred_sig, replication_sig) {
    round(mean(pred_sig == replication_sig), 4)
}

count_quadrants <- function(class_vec) {
    tibble(class = class_vec) %>%
        count(class, name = "n") %>%
        complete(
            class = c(
                "Both Significant",
                "Discovery Only",
                "Prediction Only",
                "Neither Significant"
            ),
            fill = list(n = 0)
        )
}


make_annotation_df <- function(class_vec, accuracy) {
    cnt <- count_quadrants(class_vec)

    both <- cnt$n[cnt$class == "Both Significant"]
    disconly <- cnt$n[cnt$class == "Discovery Only"]
    predonly <- cnt$n[cnt$class == "Prediction Only"]
    neither <- cnt$n[cnt$class == "Neither Significant"]

    tibble(
        label = c(
            "Exp. Port.",
            "Exp. Non-Port.",
            sprintf("Accuracy: %.2f", accuracy)
        ),
        non_port = c(disconly, neither, NA),
        port = c(both, predonly, NA)
    )
}

plot_portability <- function(
    df,
    xvar,
    yvar,
    class_var,
    tissue,
    xlab,
    ylab,
    accuracy,
    xlim = c(0, 100),
    ylim = c(0, 100)
) {
    ann <- make_annotation_df(df[[class_var]], accuracy)

 
    x_left <- xlim[2] * 0.01
    x_col1 <- xlim[2] * 0.30
    x_col2 <- xlim[2] * 0.42

    y_header <- ylim[2] * 0.98
    y_row1 <- ylim[2] * 0.94
    y_row2 <- ylim[2] * 0.89
    y_acc <- ylim[2] * 0.84

    ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
        rasterise(
            geom_point(aes(color = .data[[class_var]]), size = 1, alpha = 0.7),
            dpi = 600
        ) +
        scale_color_manual(
            values = c(
                "Both Significant" = "#3b52a7",
                "Discovery Only" = "#bd5a45",
                "Prediction Only" = "#e7c798",
                "Neither Significant" = "#18c4c2"
            )
        ) +
        coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
        theme_classic() +
        theme(
            axis.text = element_text(color = "black"),
            axis.ticks = element_line(color = "black"),
            legend.position = "none",
            plot.title = element_text(face = "bold")
        ) +
        labs(
            title = tissue,
            x = xlab,
            y = ylab
        ) +
        # 表头
        # annotate(
        #     "text",
        #     x = x_col1,
        #     y = y_header,
        #     label = "Non-Port.",
        #     hjust = 0.5,
        #     size = 3
        # ) +
        # annotate(
        #     "text",
        #     x = x_col2,
        #     y = y_header,
        #     label = "Port.",
        #     hjust = 0.5,
        #     size = 3
        # ) +
        # # 第1行
        # annotate(
        #     "text",
        #     x = x_left,
        #     y = y_row1,
        #     label = ann$label[1],
        #     hjust = 0,
        #     size = 3
        # ) +
        # annotate(
        #     "text",
        #     x = x_col1,
        #     y = y_row1,
        #     label = ann$non_port[1],
        #     color = "#bd5a45",
        #     size = 3
        # ) +
        # annotate(
        #     "text",
        #     x = x_col2,
        #     y = y_row1,
        #     label = ann$port[1],
        #     color = "#3b52a7",
        #     size = 3
        # ) +
        # # 第2行
        # annotate(
        #     "text",
        #     x = x_left,
        #     y = y_row2,
        #     label = ann$label[2],
        #     hjust = 0,
        #     size = 3
        # ) +
        # annotate(
        #     "text",
        #     x = x_col1,
        #     y = y_row2,
        #     label = ann$non_port[2],
        #     color = "#18c4c2",
        #     size = 3
        # ) +
        # annotate(
        #     "text",
        #     x = x_col2,
        #     y = y_row2,
        #     label = ann$port[2],
        #     color = "#e7c798",
        #     size = 3
        # ) +
        # Accuracy
        annotate(
            "text",
            x = x_left,
            y = y_acc,
            label = ann$label[3],
            hjust = 0,
            size = 3
        )
}

LDscore_df <- fread(
    "/path/to/src/data/gnomAD/ld_scores/gnomad.genomes.r2.1.1.nfe.adj.ld_scores.hg38.ldscore"
)
LDscore_df$chr_pos_ref_alt <- paste0(
    "chr",
    LDscore_df$CHR,
    "_",
    LDscore_df$BP,
    "_",
    LDscore_df$ref,
    "_",
    LDscore_df$alt
)

portable_info_list <- lapply(seq_len(nrow(tissue_link)), function(i) {
    tissue <- tissue_link$GTOP_tissue[i]

    df <- portability_info_list[[tissue]]$portability_df[lead_gtop == "Yes"] %>%
        dplyr::filter(
            !is.na(gtex_beta) &
                !is.na(gtex_se)
        ) %>%
        mutate(
            var_x_d = calc_var_x(gtop_maf),
            var_x_r = calc_var_x(gtex_maf)
        )

    N_GTOP <- df$gtop_n[1]
    N_GTEx <- df$gtex_n[1]

    df <- df %>%
        mutate(
            minus_log10p_raw = -log10(gtop_p_value),
            minus_log10p_gtex = -log10(gtex_p_value),

            # sample size only
            se_pred1 = gtop_se * sqrt(N_GTOP / N_GTEx),
            minus_log10p_pred1 = calc_minus_log10p(gtop_beta, se_pred1),

            # MAF only
            se_pred2 = gtop_se * sqrt(var_x_d / var_x_r),
            minus_log10p_pred2 = calc_minus_log10p(gtop_beta, se_pred2),

            # sample size + MAF
            se_pred3 = gtop_se * sqrt((N_GTOP / N_GTEx) * (var_x_d / var_x_r)),
            minus_log10p_pred3 = calc_minus_log10p(gtop_beta, se_pred3)
        ) %>%
        na.omit()

    replication_sig <- df$gtex_p_value < df$gtex_threshold
    raw_sig <- df$gtop_p_value < df$gtop_threshold
    maf_sig <- df$minus_log10p_pred2 > -log10(df$gtop_threshold)
    n_maf_sig <- df$minus_log10p_pred3 > -log10(df$gtop_threshold)

    df <- df %>%
        mutate(
            raw_type = add_classification(raw_sig, replication_sig),
            maf_type = add_classification(maf_sig, replication_sig),
            n_maf_type = add_classification(n_maf_sig, replication_sig)
        )

    acc_raw <- calc_accuracy(raw_sig, replication_sig)
    acc_maf <- calc_accuracy(maf_sig, replication_sig)
    acc_n_maf <- calc_accuracy(n_maf_sig, replication_sig)

    summary_df <- tibble(
        tissue = tissue,
        type = c("raw", "MAF", "sample_size+MAF", "gtop_n", "gtex_n"),
        value = c(
            round(c(acc_raw, acc_maf, acc_n_maf) * 100, 2),
            df$gtop_n[1],
            df$gtex_n[1]
        )
    )

    df <- merge(
        df,
        LDscore_df[, .(chr_pos_ref_alt, L2)],
        by = "chr_pos_ref_alt",
        all.x = TRUE
    )

    list(
        summary = summary_df,
        data = df
    )
})

names(portable_info_list) <- tissue_link$GTOP_tissue

portable_info_df <- rbindlist(lapply(portable_info_list, function(x) {
    x[[1]]
}))
portable_info_df1 <- portable_info_df[
    portable_info_df$type %in% c("raw", "MAF", "sample_size+MAF"),
]
portable_info_df1 %>%
    group_by(type) %>%
    summarise(mean_ratio = mean(value, na.rm = T))

#   type            mean_ratio
#   <chr>                <dbl>
# 1 MAF                   74.8
# 2 raw                   72.3
# 3 sample_size+MAF       74.6

portable_info_df1$tissue <- factor(
    portable_info_df1$tissue,
    levels = tissue_link$GTOP_tissue
)
portable_info_df1$type <- factor(
    portable_info_df1$type,
    levels = c("raw", "MAF", "sample_size+MAF")
)

color_df <- fread(
    "/path/to/GMTiP_tissue_code_and_colors.csv"
)
color_vec <- paste0("#", color_df$Tissue_Color_Code)
names(color_vec) <- color_df$Tissue

count_df <- portable_info_df1[portable_info_df1$type != "MAF", ]

overview_plot <- ggplot(
    count_df,
    aes(type, value, fill = type)
) +
    geom_boxplot() +
    geom_point(
        aes(group = tissue, color = tissue),
        position = position_dodge(width = 0.4)
    ) +
    scale_color_manual(values = color_vec) +
    theme_pub() +
    geom_line(
        aes(group = tissue, color = tissue),
        position = position_dodge(width = 0.4),
        alpha = 0.5
    ) +
    labs(x = "", y = "Proportion of portable eQTLs") +
    ggpubr::stat_compare_means(
        comparisons = list(c("raw", "sample_size+MAF")),
        method = "t.test",
        paired = T
    ) +
    scale_fill_manual(values = c("#6874b4", "#e9bd76"))

overview_plot

ggsave(
    "output/figures/eqtl_portability/portability_correction.pdf",
    width = 4.2,
    height = 4
)

fwrite(count_df, "/path/to/GTOP_code/fig-4/input/Fig4h.txt")


p1_list <- lapply(names(portable_info_list), function(x) {
    tmp_data <- portable_info_list[[x]][[2]]
    tmp_summary <- portable_info_list[[x]][[1]]
    plot_portability(
        df = tmp_data,
        xvar = "minus_log10p_gtex",
        yvar = "minus_log10p_raw",
        class_var = "raw_type",
        tissue = x,
        xlab = sprintf(
            "-log(p) GTEx (n=%s)",
            tmp_summary$value[tmp_summary$type == "gtex_n"]
        ),
        ylab = sprintf(
            "-log(p) GTOP (n=%s)",
            tmp_summary$value[tmp_summary$type == "gtop_n"]
        ),
        accuracy = tmp_summary$value[tmp_summary$type == "raw"]
    )
})

p2_list <- lapply(names(portable_info_list), function(x) {
    tmp_data <- portable_info_list[[x]][[2]]
    tmp_summary <- portable_info_list[[x]][[1]]
    plot_portability(
        df = tmp_data,
        xvar = "minus_log10p_gtex",
        yvar = "minus_log10p_pred2",
        class_var = "maf_type",
        tissue = x,
        xlab = sprintf(
            "-log(p) GTEx (n=%s)",
            tmp_summary$value[tmp_summary$type == "gtex_n"]
        ),
        ylab = "-log(p) GTOP (adjusted for MAF)",
        accuracy = tmp_summary$value[tmp_summary$type == "MAF"]
    )
})

p3_list <- lapply(names(portable_info_list), function(x) {
    tmp_data <- portable_info_list[[x]][[2]]
    tmp_summary <- portable_info_list[[x]][[1]]
    plot_portability(
        df = tmp_data,
        xvar = "minus_log10p_gtex",
        yvar = "minus_log10p_pred3",
        class_var = "n_maf_type",
        tissue = x,
        xlab = sprintf(
            "-log(p) GTEx (n=%s)",
            tmp_summary$value[tmp_summary$type == "gtex_n"]
        ),
        ylab = "-log(p) GTOP\n(adjusted for MAF and study size)",
        accuracy = tmp_summary$value[tmp_summary$type == "sample_size+MAF"]
    )
})
names(p1_list) <- names(p2_list) <- names(p3_list) <- names(portable_info_list)

overview_plot +
    p1_list$Adrenal_Gland +
    p2_list$Adrenal_Gland +
    p3_list$Adrenal_Gland +
    plot_layout(nrow = 1, widths = c(0.5, 1, 1, 1))

ggsave(
    "output/figures/eqtl_portability/Adrenal_Gland_portability.pdf",
    width = 15,
    height = 4
)

saveRDS(
    portable_info_list,
    file = "/path/to/input/Fig35_data.rds"
)

portable_data <- rbindlist(lapply(names(portable_info_list), function(x) {
    tmp_data <- portable_info_list[[x]][[2]]
    tmp_data$tissue <- x
    tmp_data
}))

# bp1 <- ggplot(portable_data, aes(tissue, L2, fill = raw_type)) +
#     geom_boxplot() +
#     theme_pub() +
#     scale_fill_manual(
#         values = c(
#             "Both Significant" = "#3b52a7",
#             "Discovery Only" = "#bd5a45",
#             "Prediction Only" = "#e7c798",
#             "Neither Significant" = "#18c4c2"
#         )
#     )

bp2 <- ggplot(
    portable_data[
        !is.na(portable_data$L2) &
            !portable_data$tissue %in% c("Muscle", "Skin"),
    ],
    aes(n_maf_type, L2, fill = n_maf_type)
) +
    geom_boxplot(outlier.size = 0.8) +
    theme_pub() +
    scale_fill_manual(
        values = c(
            "Both Significant" = "#3b52a7",
            "Discovery Only" = "#bd5a45",
            "Prediction Only" = "#e7c798",
            "Neither Significant" = "#18c4c2"
        )
    ) +
    ggpubr::stat_compare_means(
        comparisons = list(
            c("Both Significant", "Prediction Only"),
            c("Discovery Only", "Neither Significant")
        ),
        method = "t.test"
    ) +
    labs(x = "", y = "LD score")

bp2

ggsave(
    "path/to/portability_correction_LDscore.pdf",
    width = 5,
    height = 4
)

LD_value <- portable_data |>
    group_by(tissue, n_maf_type) |>
    summarise(mean_LD = mean(L2, na.rm = T))
LD_value <- reshape2::dcast(
    LD_value,
    tissue ~ n_maf_type,
    value.var = "mean_LD"
)

mean(LD_value$`Both Significant`) # 179.2248
mean(LD_value$`Discovery Only`) # 157.6043
mean(LD_value$`Neither Significant`, na.rm = T) # 56.13948
mean(LD_value$`Prediction Only`, na.rm = T) # 54.52136
