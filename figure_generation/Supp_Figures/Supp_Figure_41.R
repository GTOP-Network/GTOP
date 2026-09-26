#==============================================#
# Coloc_information #
# Supp-Figure-41#
#==============================================#
library(data.table)
library(ggplot2)


setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig41.Coloc_information/input")

# Supp.Fig.41a:  ---------------------------------------------------
locus_lvl <- readRDS("supp_fig41a_data.rds")

tmp_plot_df <- locus_lvl[, .N, by = .(phecode_abbr, best_status)] |>
  arrange(best_status, desc(N))

tmp_plot_df$phecode_abbr <- factor(
  tmp_plot_df$phecode_abbr,
  levels = unique(tmp_plot_df$phecode_abbr)
)
tmp_plot_df$best_status <- factor(
  tmp_plot_df$best_status,
  levels = c("no_eQTL", "eQTL_no_coloc", "coloc")
)
ggplot(tmp_plot_df, aes(phecode_abbr, N, fill = best_status)) +
  geom_col() +
  labs(y = "Number of GWAS loci", x = NULL, fill = "Locus class") +
  theme_classic() +
  scale_fill_viridis_d() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))


# Supp.Fig.41b:  ---------------------------------------------------
or_results <- readRDS("supp_fig41b_data.rds")

ggplot(or_results, aes(x = OR, y = feature)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.7) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.16,
    linewidth = 0.8
  ) +
  geom_point(size = 3.4) +
  scale_x_log10() +
  labs(
    title = "Locus features associated with colocalization",
    # subtitle = "Odds ratio for colocalization per 1-SD increase in each feature",
    x = "Odds ratio (95% CI; log scale)",
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.line.y = element_blank()
  )


# Supp.Fig.41c:  ---------------------------------------------------
gene_assignment_df <- data.frame(
  "type" = rep(
    c(
      "All fine-mapped QTL signals",
      "All colocalized fine-mapped GWAS loci",
      "All colocalized GWAS loci"
    ),
    each = 2
  ),
  "type2" = rep(c("Single eGene", "Multiple eGenes"), 3),
  "value" = c(16880, 21454 - 16880, 295, 369 - 295, 533, 635 - 533)
)
ggplot(gene_assignment_df, aes(value, type, fill = type2)) +
  geom_col(position = "fill") +
  labs(y = "", x = "Gene assianment", fill = "Locus class") +
  scale_fill_manual(values = c("#e38785", "#95bce4")) +
  theme_classic()

