#===============================================#
# Assessment of eQTL sharing between GTOP and GTEx #
# Supp-Figure-35#
#===============================================#

library(tidyverse)
library(data.table)
library(patchwork)
library(ggrastr)

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
  
  #
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

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig35.Portability_correction//")

portable_info_list <- readRDS("input/Fig35_data.rds")

portable_info_df <- rbindlist(lapply(portable_info_list, function(x) {
  x[[1]]
}))
portable_info_df1 <- portable_info_df[
  portable_info_df$type %in% c("raw", "MAF", "sample_size+MAF"),
]
portable_info_df1 %>%
  group_by(type) %>%
  summarise(mean_ratio = mean(value, na.rm = T))

portable_info_df1$type <- factor(
  portable_info_df1$type,
  levels = c("raw", "MAF", "sample_size+MAF")
)

color_vec <- readRDS("../../fig-4/input/tissue_color.RDS")

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
  theme_classic() +
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

p_combined <- overview_plot +
  p1_list$Adrenal_Gland +
  p2_list$Adrenal_Gland +
  p3_list$Adrenal_Gland +
  plot_layout(nrow = 1, widths = c(0.5, 1, 1, 1))

print(p_combined)
