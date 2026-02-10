

## popSpecific SV and popSpecific SNV
## popSpecific SV and non-popSpecific SNV


#!/usr/bin/env Rscript

setwd("/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project")

# Load required libraries
library(data.table)
library(tidyverse)

## Tissue color: color_vec
color_df <- fread("/media/pacific/share/Datasets/Asian_GTEx/Metainfo/input/GMTiP_tissue_code_and_colors.csv")
color_vec <- paste0("#", color_df$Tissue_Color_Code)
names(color_vec) <- color_df$Tissue
rm(color_df)

## Allele frequency of SNV and SV: freq_data
freq_data <- fread("2025-06-11-specific_xQTL/output/data_plot/data/SNP_SV_addgroup_information.txt")

## Finemapping result of SNV-xQTL
finemapping_SNV_SVs <- fread("2025-06-11-specific_xQTL/output/data_plot/data/SNV_SV_variant_xQTL_finemapping_add_type.txt")
finemapping_SNV_SVs <- finemapping_SNV_SVs[!pip_Type%in%c("only_SNV")]
length(unique(finemapping_SNV_SVs$variant_id))  # 489726 --> 158517
ggvenn::ggvenn(list("freq"=unique(freq_data$variant_id[freq_data$Type=="SV"]), 
                    "fine"=unique(finemapping_SNV_SVs$variant_id[finemapping_SNV_SVs$Type=="SV"])))

## Allele frequency + finemapping
freq_finemapping <- merge(freq_data[,.(variant_id, Type, af_GMTiP, af_EUR, two_group)],
                          finemapping_SNV_SVs[, .(variant_id, phenotype_id, xGene, symbol, pip, tissue, xQTL_type, cs_id, cs_size, gene_type, pip_Type)],
                          by="variant_id")

ff_SV <- freq_finemapping[Type=="SV"]
length(unique(ff_SV$variant_id))
table(sapply(strsplit(unique(ff_SV$variant_id), split = "_" ), function(x){x[3]}))
# chr19_21647331_INV_3BM5C9_415127

# chr2_97202193_DEL_37M6F_5609  ANKRD36
# chr19_12053911_DEL_0M5BC_12871 ZNF844
# chr14_52785748_DEL_1DM437_2613

write.table(unique(ff_SV$xGene), "2025-06-11-specific_xQTL/output/data_plot/data/specific_SV_genes.txt", quote = F, row.names = F, col.names = F)

## --
tmp_A <- ff_SV[two_group=="EUR_rare"]
length(unique(tmp_A$variant_id))

tmp_A <- finemapping_SNV_SVs %>% group_by(phenotype_id, xGene, symbol, tissue, xQTL_type, cs_id, cs_size, gene_type, pip_Type) %>% 
  summarise(SV_list=variant_id[Type=="SV"&pip==max(pip[Type=="SV"])][1],
            SV_pip=pip[Type=="SV"&pip==max(pip[Type=="SV"])][1],
            SNV_list=variant_id[Type=="SNV"&pip==max(pip[Type=="SNV"])][1],
            SNV_pip=pip[Type=="SNV"&pip==max(pip[Type=="SNV"])][1]) %>% 
  as.data.table()
tmp_A <- tmp_A[pip_Type!="only_SV"]
tmp_A$SV_type <- sapply(strsplit(tmp_A$SV_list, split = "_" ), function(x){x[3]})

tmp_A <- merge(tmp_A[, .(SV_type, phenotype_id, xGene, symbol, tissue, cs_id, gene_type, pip_Type, SV_list, SNV_list)], 
               unique(freq_finemapping[, .(SV_list=variant_id, SV_group=two_group)]), 
               by="SV_list")
tmp_A <- merge(tmp_A[, .(SV_type, phenotype_id, xGene, symbol, tissue, cs_id, gene_type, pip_Type, SV_list, SV_group, SNV_list)], 
               unique(freq_finemapping[, .(SNV_list=variant_id, SNV_group=two_group)]), 
               by="SNV_list")
tmp_A <- tmp_A[SV_type != "INV"]

tres <- vector()
for(i in unique(tmp_A$SV_type)){
  for(j in unique(tmp_A$pip_Type)){
    tmp_df <- tmp_A[SV_type==i&pip_Type==j]
    val1 <- nrow(tmp_df[SV_group=="EUR_rare"&SNV_group=="EUR_rare"])
    val2 <- nrow(tmp_df[SV_group=="EUR_common"&SNV_group=="EUR_rare"])
    val3 <- nrow(tmp_df[SV_group=="EUR_rare"&SNV_group=="EUR_common"])
    val4 <- nrow(tmp_df[SV_group=="EUR_common"&SNV_group=="EUR_common"])
    tdata <- fisher.test(matrix(c(val1, val2, val3, val4), byrow = T, nrow = 2))
    tres <- c(tres, i, j, val1, val2, val3, val4, 
              tdata$p.value, tdata$estimate, tdata$conf.int[1], tdata$conf.int[2])
  }
}
tres_df <- as.data.frame(matrix(tres, ncol = 10, byrow = T))
colnames(tres_df) <- c("Type", "Type2", "v1", "v2", "v3", "v4",
                       "pvalue", "OR", "lower_CI", "upper_CI")
for(i in 3:10){tres_df[,i] <- as.numeric(tres_df[,i])}
tres_df$Psig <- case_when(
  tres_df$pvalue < 0.05 ~ "Sig",
  TRUE ~ "NS"
)
# tres_df <- tres_df[tres_df$v1 >=10, ]
tres_df$se_value <- (tres_df$upper_CI - tres_df$lower_CI) / (2 * 1.96)
unique(tres_df$Type2)
tres_df$Type2 <- factor(tres_df$Type2, levels = rev(c(
  "lead_SNV", "lead_SV", "equal_SNV_SV"
)))
tres_df <- tres_df[order(tres_df$Type, decreasing = F), ]

ggforestplot::forestplot(
  df = tres_df,
  name = Type,
  estimate = OR,
  se = se_value,
  pvalue = pvalue,
  psignif = 0.05,
  logodds = FALSE,
  xlab = "Odds ratio",
  title = "",
  colour = Type2,
  size = 0.5
)

ggsave("2025-06-11-specific_xQTL/output/data_plot/plot/specific_common_SV_enrichment.pdf", width = 4, height = 3)


tmp_B <- freq_finemapping %>% group_by(tissue, pip_Type, two_group) %>% 
  summarise(count=length(unique(xGene[Type=="SV"])))
tmp_B <- reshape2::dcast(tmp_B, tissue + pip_Type ~ two_group, value.var = "count")

ggplot(data = tmp_B, aes(x=EUR_common, y=EUR_rare, color=tissue, shape=pip_Type)) +
  geom_point(size=3) +
  theme_classic() +
  scale_color_manual(values = color_vec) +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  facet_wrap(pip_Type~., scales = "free")

ggsave("2025-06-11-specific_xQTL/output/data_plot/plot/specific_common_SV_gene_correlation.pdf", width = 7, height = 5)

