#===============================================#
# Assessment of eQTL sharing between GTOP and GTEx #
# Supp-Figure-34#
#===============================================#

library(tidyverse)
library(data.table)
library(patchwork)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig34.Portability_three_metrics/")

fig_data <- readRDS("input/Fig34_data.rds")

fig_data$type2 <- factor(
  fig_data$type2,
  levels = c(
    "ratio_nominal",
    "ratio_mash",
    "ratio_uncertainty-aware",
    "FDR",
    "gene"
  )
)
fig_data$type3 <- factor(
  fig_data$type3,
  levels = c("Lead eVariant", "All eVariants", "All eGenes", "coloc genes")
)

p1 <- ggplot(
  fig_data[fig_data$type3 == "Lead eVariant", ],
  aes(type1, count / 1000, fill = type1)
) +
  geom_col(width = 0.8) +
  facet_grid(
    tissue ~ type2,
    scales = "free_y",
    axes = "all_y",
    axis.labels = "margins"
  ) +
  theme_classic() +
  scale_fill_manual(values = c("#ddbf93", "#af5543")) +
  labs(x = "", y = "Number of eVariants") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

p2 <- ggplot(
  fig_data[fig_data$type3 != "Lead eVariant", ],
  aes(type1, count / 1000, fill = type1)
) +
  geom_col(width = 0.8) +
  facet_grid(
    tissue ~ type3,
    scales = "free_y",
    axes = "all_y",
    axis.labels = "margins"
  ) +
  theme_classic() +
  scale_fill_manual(values = c("#ddbf93", "#af5543")) +
  labs(x = "", y = "Number of eGenes") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

p1 + p2 + plot_layout(guides = "collect", widths = c(2, 1))
