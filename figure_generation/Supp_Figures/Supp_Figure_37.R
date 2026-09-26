#===============================================#
# Fine-mapping-filter-poor-coverage #
# Supp-Figure-37                         #
#===============================================#
library(data.table)
library(dplyr)
library(stringi)
library(stringr)
library(pbapply)
library(ggpubr)
library(purrr)
library(forcats)
library(tidyr)
library(ggplot2)
library(magrittr)
library(reshape2)


setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig37.Fine-mapping-filter-poor-coverage/")


# Supp.Fig.37a: Number of cs within SV/TR group by Tissue ------------------------------------------------------------

tissue_order <- c("Whole_Blood","Pancreas_Tail","Pancreas_Head","Pancreas_Body",
                  "Adrenal_Gland","Liver","Adipose","Gallbladder","Spleen","Muscle","Skin" )

fread("input/Supp_Fig37a.number_of_CS.txt") %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order),
         class = factor(class)) %>%
  ggplot(aes(x = Tissue, y = n, fill = class)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = rep("#8090B4", 3), guide = "none") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Tissues", y = "# of credible sets")

fread("input/Supp_Fig37a.proportion_of_CS.txt") %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order)) %>%
  ggplot(aes(x = Tissue, y = prop, fill = type)) +
  geom_col(position = "stack") +
  scale_fill_manual(breaks = c("PC_snv", "PC_sv_tr", "PC_sv_trLead"),
                    values = c("#bdbdbd", "#fc9272", "#de2d26")) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Tissues", y = "Proportion of credible sets")



# Supp.Fig.37b: enrichment by tissue and by SV or TR -----------------------

dat <- readRDS("input/Supp_Fig37b.SVTR_lead_vs_SNV_lead_enrichment.RDS")
color <- read.csv("/media/london_A/mengxin/GTOP_code/supp/supp_fig40.R4.Heritability-long-reads-coloc/input/GTOP_tissue_coloc_code", header = T) %>%
  filter(Tissue %in% unique(dat$Tissue))
dat$PIP <- factor(dat$PIP)
model_p <- summary(lm(OR ~ PIP, data=dat[dat$Type=="SVTR" & dat$Pval<0.05,]))$fstatistic
pval <- pf( model_p[1], model_p[2],model_p[3],lower.tail = FALSE)
adjr2 <- summary(lm(OR ~ PIP, data=dat[dat$Type=="SVTR" & dat$Pval<0.05,]))$adj.r.squared

ggplot(dat[dat$Pval<0.05,], aes(x = PIP, y = OR, fill = as.numeric(as.character(PIP)))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  geom_jitter(aes(color = Tissue), width = 0.2, size = 1.2) +
  scale_fill_gradient(low = "#fff7ec", high = "#b30000", name = "PIP threshold") +
  scale_color_manual(values = setNames(paste0("#", color$Tissue_Color_Code), color$Tissue),guide="none") +
  theme_classic() +
  labs(x = "PIP threshold", y = "Odds Ratio (lead SV or TR vs lead SNV)") +
  theme(axis.text.x = element_text(colour = "black"),legend.position = "top")+
  annotate("text",x=1,y=max(dat$OR[dat$Pval<0.05],na.rm=T)*1.15,
    label=paste0("P = ",format.pval(pval,digits=3),"\nAdj. R² = ", round(adjr2,3)),size=4 )

  








