#==============================================#
# evoluation #
# Supp-Figure-27#
#==============================================#

library(ggplot2)
library(tidyverse)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig27")


# Supp.Fig.27a evolution -----------------------------------------
selection_res <- fread("./input/supp_fig27a_data.txt")

ggplot(selection_res) +
  geom_col(aes(y=reorder(tissue, qtl_number), x=qtl_number)) +
  geom_col(aes(y=tissue, x=-gene_number)) +
  facet_grid(.~Type, scales = "free") +
  theme_classic() +
  labs(x="Number of variants/eGenes", y="")


gene_list <- readRDS("./input/supp_fig27b_data.RDS")
ggvenn::ggvenn(gene_list)


# Example of fd-QTL under evoluation --------------------------------------------

gene_qtl_subdf1 <- fread("./input/supp_fig27c_data.txt")

p1 <- ggplot(gene_qtl_subdf1) +
  geom_point(aes(x=pos, maf_GTOP), color="#bd5a45", alpha=0.3) +
  geom_segment(aes(x=pos, xend=pos, y=0, yend=maf_GTOP), color="#bd5a45", alpha=0.3) +
  geom_point(aes(x=pos, maf_GTEx), color="#e7c798", alpha=0.3) +
  geom_segment(aes(x=pos, xend=pos, y=0, yend=maf_GTEx), color="#e7c798", alpha=0.3) +
  theme_classic() +
  labs(x="", y="MAF") 

p2 <- ggplot(gene_qtl_subdf1) +
  geom_point(aes(x=pos, -log10(gene_qtl_subdf1$p_value)), color="grey") +
  theme_classic() +
  labs(x="", y="-log10(P)")

loci_start <- min(as.numeric(gene_qtl_subdf1$pos))
loci_end <- max(as.numeric(gene_qtl_subdf1$pos))
SNV_eQTL_locus <- locus(data = as.data.frame(gene_qtl_subdf1[, .(chrom=chr, pos=as.numeric(pos), 
                                                                 rsid=variant_id, p=p_value)]),
                        xrange = c(loci_start, loci_end),
                        seqname = gene_qtl_subdf1$chr[1],
                        ens_db = "EnsDb.Hsapiens.v86")
gene_plot <- gg_genetracks(SNV_eQTL_locus, highlight = example_gene, 
                           filter_gene_name = c(example_gene),
                           filter_gene_biotype = c("protein_coding"))

p1+p2+gene_plot+plot_layout(ncol = 1, heights = c(3,3,3,0.8))
