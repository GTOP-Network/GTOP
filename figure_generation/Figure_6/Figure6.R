#==================================#
# Disease association #
# Figure-6 #
#==================================#

setwd("/path/to/GTOP_code/fig-6")

library(data.table)
library(tidyverse)
library(ggradar)
library(ggpubr)
library(ComplexHeatmap)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(ggvenn)
library(ggradar)

## Fig.6a, SLDSC results -------------------------------------------------------

plotdf <- fread("./input/Fig6a.txt")
plotdf %>%
    mutate(
        method = factor(method, levels = c("allsig", "finemap", "finemappip")),
        qtltype = factor(qtltype, levels = c("eQTL", "juQTL", "tuQTL"))
    ) %>%
    ggplot(.) +
    geom_pointrange(
        aes(
            x = method,
            y = logen,
            ymax = logCI_upper,
            ymin = logCI_lower,
            color = qtltype
        ),
        position = position_dodge(width = .8),
        size = .3
    ) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "gray60",
        alpha = 0.5
    ) +
    scale_color_manual(
        values = c("eQTL" = "#a2c398", "juQTL" = "#5181b1", "tuQTL" = "#b97975")
    ) +
    ylab("log2(Heritability enrichment)") +
    xlab("") +
    theme_pubr()


## Fig.6b, coloc and SMR results -----------------------------------------------

plot_data <- readRDS("./input/Fig6b.RDS")

ggvenn::ggvenn(
    plot_data,
    fill_color = c("#9ac294", "#c2c48e", "#7b86a7", "#8c7c9c"),
    show_percentage = F
)

plot_GWAS_loci <- fread("./input/Fig6b.txt")
rect_df <- unique(plot_GWAS_loci[, c("Trait_class", "Phecode")])
rect_df <- as.data.frame(table(rect_df$Trait_class)[unique(
    rect_df$Trait_class
)])
rect_df$data2 <- cumsum(rect_df$Freq) + 0.5
rect_df$data1 <- rect_df$data2 - (rect_df$Freq)
rect_df$Phecode <- rect_df$Var1

plot_GWAS_loci$Phecode <- factor(
    plot_GWAS_loci$Phecode,
    levels = unique(plot_GWAS_loci$Phecode)
)
plot_GWAS_loci$type <- factor(
    plot_GWAS_loci$type,
    levels = c("eQTL", "sQTL", "eQTL,sQTL")
)

ggplot(plot_GWAS_loci, aes(x = Phecode)) +
    geom_col(aes(y = qtl_loci_count, fill = type), width = 1) +
    labs(
        title = "",
        x = "Trait category",
        fill = "Category",
        y = "Number of GWAS loci"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    geom_rect(
        data = rect_df,
        aes(xmin = data1, xmax = data2, ymin = -5, ymax = -1, fill = Phecode),
        color = "white",
        size = 0
    ) +
    scale_fill_manual(
        values = c(
            "Quantitative_trait" = "#2F4F4F",
            "Circulatory_system" = "#487c51",
            "Dermatologic" = "#7690a4",
            "Digestive" = "#f0d795",
            "Endocrine_metabolic" = "#4682B4",
            "Genitourinary" = "#B0C4DE",
            "Hematopoietic" = "#00bfc4",
            "Infectious_diseases" = "#9e4832",
            "Mental_disorders" = "#7c776b",
            "Musculoskeletal" = "#DB7093",
            "Neoplasms" = "#7371a3",
            "Neurological" = "#00008B",
            "Respiratory" = "#5192c1",
            "Sense_organs" = "#c4598a",
            "Symptoms" = "#2F4F4F",
            "eQTL" = "#9ac294",
            "sQTL" = "#7b86a7",
            "eQTL,sQTL" = "#bc9cc1"
        )
    )


## Fig.6c, Compared with GTEx, JCTF, and MAGE ----------------------------------

compare_ratio_df <- fread("./input/Fig6c.txt")

ggradar::ggradar(
    compare_ratio_df,
    values.radar = c("0", "30", "60"),
    grid.min = 0,
    grid.mid = 30,
    grid.max = 60,
    group.line.width = 1,
    group.point.size = 3,
    group.colours = c("#5baab8", "#dbb938", "#e15f28"),
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    legend.position = "bottom"
)


## Fig.6d, Example of GTOP specific compared with GTOP -------------------------

plot_df <- fread("./input/Fig6d.txt")

ggplot(plot_df, aes(-log10(GWAS_p), -log10(QTL_p))) +
    geom_point(aes(color = type), size = 2, alpha = 0.9) +
    theme_classic(base_size = 14) +
    scale_color_manual(values = c("GTEx" = "#e7c798", "GTOP" = "#bd5a45")) +
    theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")
    )


## Fig.6e, enrichment in GWAS and TPMI -----------------------------------------
data <- fread("./input/Fig6e.txt")

data$type <- factor(
    data$type,
    levels = unique(data$type)
)
data$variable <- factor(
    data$variable,
    levels = unique(data$variable)
)
ggplot(data, aes(x = value, y = type, fill = variable)) +
    geom_bar(stat = "identity", position = "fill", alpha = 0.8) +
    labs(y = "", x = "Percentage of significant variants (%)") +
    theme_classic() +
    scale_fill_manual(values = c("tested" = "#e7c798", "signif" = "#bd5a45")) +
    theme(axis.text = element_text(colour = "black"), legend.position = "none")


## Fig.6f, GWAS rare variants ------------------------------------------------
plot_df <- fread("./input/Fig6f.txt", sep = "\t")

plot_df$type1 <- factor(
    plot_df$type1,
    levels = c(
        "At least \n one ancestry \n MAF <0.01",
        "All ancestries \n MAF <0.01"
    )
)
plot_df$type2 <- factor(
    plot_df$type2,
    levels = c("eQTL", "sQTL", "eQTL + sQTL")
)

ggplot(plot_df, aes(type1, value, fill = type2)) +
    geom_col(position = position_dodge(width = 0.9), width = 0.8) +
    theme_classic() +
    labs(x = "", y = "Number of colocalized GWAS loci") +
    scale_fill_manual(
        values = c(
            "eQTL" = "#95c595",
            "sQTL" = "#7b86a7",
            "eQTL + sQTL" = "#eab096"
        )
    ) +
    theme(
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")
    )


## Fig.6g, Number of disease-related SVs and TRs -------------------------------
df_coloc_stat <- readRDS("./input/Fig6g.RDS")

upset_list_snv <- df_coloc_stat$gene_tissue[df_coloc_stat$have_SNV]
upset_list_SV <- df_coloc_stat$gene_tissue[df_coloc_stat$have_SV]
upset_list_TR <- df_coloc_stat$gene_tissue[df_coloc_stat$have_TR]

gene_list <- list(SNV = upset_list_snv, SV = upset_list_SV, TR = upset_list_TR)

library(ComplexHeatmap)
mat <- make_comb_mat(gene_list, mode = "distinct")
UpSet(mat, comb_order = order(-comb_size(mat)))


## Fig.6h, Example of SVs ------------------------------------------------------

load("./input/Fig6h.RData")

## loci information
SNP_name <- "rs77303550"
SV_name <- "chr16_72056410_DEL_CM400_1724"
loci_start <- SNV_eQTL_data$pos[which.min(SNV_eQTL_data$p)] - 120000
loci_end <- SNV_eQTL_data$pos[which.min(SNV_eQTL_data$p)] + 200000
GWAS_name <- "Hyperlipidemia"
chr_name <- GWAS_data$chrom[1]

## GWAS locuszoom
GWAS_locus <- locus(
    data = as.data.frame(GWAS_data),
    xrange = c(loci_start, loci_end),
    seqname = chr_name,
    index_snp = SNP_name,
    ens_db = "EnsDb.Hsapiens.v86"
)
GWAS_locus$data$ld <- LD_info[GWAS_locus$data$rsid, "R2"]
GWAS_locus$data$col <- "transparent"
GWAS_plot <- gg_scatter(
    GWAS_locus,
    pcutoff = FALSE,
    yzero = T,
    labels = SNP_name,
    color = "black",
    legend_pos = "right",
    LD_scheme = c(
        "#e5e5e5",
        "#e5e5e5",
        "#3e70b4",
        "#3f7d1d",
        "orange",
        "red",
        "red"
    )
) +
    annotate(
        "text",
        x = loci_end / 10^6,
        y = max(-log10(GWAS_data$p)),
        label = GWAS_name,
        hjust = "right"
    ) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())

## eQTL locuszoom
SNV_eQTL_locus <- locus(
    data = as.data.frame(SNV_eQTL_data),
    xrange = c(loci_start, loci_end),
    seqname = chr_name,
    index_snp = SNP_name,
    ens_db = "EnsDb.Hsapiens.v86"
)
SNV_eQTL_locus$data$ld <- LD_info[SNV_eQTL_locus$data$rsid, "R2"]
SNV_eQTL_locus$data$col <- "transparent"
SNV_eQTL_plot <- gg_scatter(
    SNV_eQTL_locus,
    pcutoff = FALSE,
    yzero = T,
    labels = SNP_name,
    color = "black",
    legend_pos = "right",
    LD_scheme = c(
        "#e5e5e5",
        "#e5e5e5",
        "#3e70b4",
        "#3f7d1d",
        "orange",
        "red",
        "red"
    )
) +
    annotate(
        "text",
        x = loci_end / 10^6,
        y = max(-log10(SNV_eQTL_data$p)),
        label = "eQTL_Whole Blood",
        hjust = "right"
    ) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())

## SV-eQTL locuszoom
SV_eQTL_locus <- locus(
    data = as.data.frame(SV_eQTL_data),
    xrange = c(loci_start, loci_end),
    seqname = chr_name,
    index_snp = SV_name,
    ens_db = "EnsDb.Hsapiens.v86"
)
SV_eQTL_locus$data$col <- "transparent"
SV_eQTL_plot <- gg_scatter(
    SV_eQTL_locus,
    pcutoff = FALSE,
    labels = SV_name,
    color = "black",
    legend_pos = "right",
    LD_scheme = c(
        "#e5e5e5",
        "#e5e5e5",
        "#3e70b4",
        "#3f7d1d",
        "orange",
        "red",
        "red"
    )
) +
    annotate(
        "text",
        x = loci_end / 10^6,
        y = max(-log10(SV_eQTL_data$p)),
        label = "eQTL_Whole Blood",
        hjust = "right"
    ) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())

## joint-fine-mapping
joint_fm_locus <- locus(
    data = as.data.frame(joint_fm_data),
    yvar = "PIP",
    xrange = c(loci_start, loci_end),
    seqname = chr_name, #index_snp = SNP_name,
    ens_db = "EnsDb.Hsapiens.v86"
)
joint_fm_locus$data$col <- "transparent"
joint_fm_plot <- gg_scatter(
    joint_fm_locus,
    pcutoff = FALSE,
    yzero = T,
    labels = SV_name,
    color = "black",
    legend_pos = "right",
    LD_scheme = c(
        "#e5e5e5",
        "#e5e5e5",
        "#3e70b4",
        "#3f7d1d",
        "orange",
        "red",
        "red"
    )
) +
    annotate(
        "text",
        x = loci_end / 10^6,
        y = max(joint_fm_data$PIP),
        label = "Joint fine-mapping",
        hjust = "right"
    ) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())

## gene structure
gene_plot <- gg_genetracks(
    SNV_eQTL_locus,
    highlight = "HP",
    filter_gene_name = c("HP"),
    filter_gene_biotype = c("protein_coding")
)

## plot
wrap_plots(
    list(GWAS_plot, SNV_eQTL_plot, SV_eQTL_plot, joint_fm_plot, gene_plot),
    ncol = 1,
    heights = c(3, 3, 3, 3, 1)
)
