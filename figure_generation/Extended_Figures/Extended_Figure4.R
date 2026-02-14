#==============================================#
# Extended Fig.4 #
#==============================================#

library(data.table)
library(ggplot2)
library(stringi)
library(stringr)
library(dplyr)
library(ggsci)
library(tidyverse)
library(ggpubr)
library(magrittr)
library(ggbreak)
library(ggrastr)
setwd("/media/london_A/mengxin/GTOP_code/extend/extend_145/input")

# Extended.Fig.4a ---------------------------------------------------------

qq_plot_list <- readRDS("qq_plot.all_tissue.RDS")
p <- ggplot(qq_plot_list,aes(x=ExpP,y=ObservedP,color=Tissue)) + 
  geom_point_rast(raster.dpi = getOption("ggrastr.default.dpi", 300)) + 
  geom_abline(intercept = 0,slope = 1,color="grey") +
  theme_pubr() + scale_color_manual(values=tissue_colors)
ggsave(
  filename = "eQTL_genomic_inflation_lambda.pdf",
  plot = p,              
  width = 6,
  height = 6,
  units = "in"
)


# Extended Fig.4b:TSS supp ------------------------------------------------------------------

dat <- readRDS("ExtendFig4b.RDS")

ggplot(dat,aes(x=distance,color=group)) + geom_density(size=1.5) + 
  theme_pubr() + 
  scale_color_manual(values = c("#7B88A8","#609561","#983927"))

# Extended.Fig.4c:  torus enrichment ------------------------------------------------------------------

order <- rev(c("enhancer","promoter","open chromatin region","CTCF binding site","TF binding site","3 prime UTR","5 prime UTR","frameshift","intron","missense","NC transcript","splice acceptor","splice donor","splice region","stop gained", "synonymous"))

plotdf <- fread("Exfig 4c.left.txt")
p1 <- plotdf %>%
  dplyr::mutate(Ann=factor(Ann, levels=order)) %>% 
  ggplot(.) +
  geom_pointrange(aes(x=Ann, y = logmFC, ymin=lFC, ymax=hFC, color = QTL, shape=QTL), 
                  position=position_dodge(width=0.8), size = 0.5)+
  scale_color_manual(values=c("eQTL"="#a2bf98" , "juQTL"="#7d8bad",
                              "tuQTL"="#b47973"))+
  scale_shape_manual(values = c("eQTL"=16, "juQTL"=15, "tuQTL"=17))+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  # ylim(-8,20)+
  ylab(expression("Log"[2]*"(Fold Enrichment)"))+
  xlab("")+
  coord_flip()+
  theme_pubr()+
  theme(
    axis.text.x = element_text(color="black"),
    axis.text.y = element_text(color="black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
p1

plotdf <- fread("Exfig 4c.right.txt")
p2 <- plotdf %>% 
  dplyr::mutate(feature=factor(feature, levels=order)) %>% 
  ggplot(.)+
  geom_bar(aes(x=meanratio, y=feature, fill=type), stat = "identity", position=position_dodge(width=0.9))+
  geom_errorbar(aes(xmax=meanratio+sdratiio, xmin=meanratio-sdratiio, y=feature, group=type), position=position_dodge(width=0.9),width=0)+
  scale_fill_manual(values=c("eQTL"="#a2bf98" , "juQTL"="#7d8bad",
                             "tuQTL"="#b47973"))+
  # scale_x_break(c(0.12,0.65),
  #               space = 0.3,
  #               scales = .5)+
  theme_pubr(legend = "top")+
  xlab("Proportion of variants")+
  ylab("")+
  theme(axis.text.y = element_blank())
p2

library(patchwork)
all <- wrap_elements(p1) + wrap_elements(p2) +
  plot_layout(widths = c(2.5, 1))

all
# Extended Fig.4d:  eVariant distince from TSS ------------------------------------------------------------------

effect_data <- fread("ExtendFig4d-e.txt")
ggplot(effect_data)+ geom_violin(aes(y=QTL, x=slope, fill=VarSubType), 
                                       position = position_dodge())+
  xlab("Effect Size")+
  ylab("Molecular Trait")+
  theme_classic()+theme(legend.position = "none")+
  scale_fill_manual(values = c("SNV"="#c9c9c9","STR"="#0f3c7a","VNTR"="#ab889a","INS"="#b55f60","DEL"="#577b95"))


# Extended.Fig.4e: larger sv length with higher effect size  --------------------------------------------------------
sv_qtl <- effect_data[effect_data$VarType=="SV" & effect_data$QTL=="eQTL",] %>%
  mutate(length=as.integer(word(variant_id,start=5,end=5,sep="\\_")),
         group=case_when(length < 100 ~ "<100bp",
                         length < 200 & length >= 100 ~ "<200bp",
                         length < 500 & length >= 200 ~ "<500bp",
                         length < 10000 & length >= 500 ~ "<10kb",
                         length >= 10000 ~ ">10kb"))

sv_qtl$group <- factor(sv_qtl$group,levels = c("<100bp","<200bp","<500bp","<10kb",">10kb"))
groups2 <- levels(sv_qtl$group)
comparisons2 <- combn(groups2, 2, simplify = FALSE)
ggplot(sv_qtl, aes(x = group, y = abs(slope), fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.4, color = "black") +
  scale_fill_manual(values = c("#f1eef6","#bdc9e1","#74a9cf","#2b8cbe","#045a8d")) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "|Effect Size|") +
  stat_compare_means( comparisons = comparisons2,
                      method = "wilcox.test", label = "p.format",step.increase = 0.08 )



# Extended.Fig.4f: eGene constraint score ---------------------------------------------------

fet_plot <- fread("ExtendFig4f.eGene_constraint.txt") %>%
  mutate(log2OR = log2(OR),
         p_star = case_when(FDR < 0.001 ~ "***",FDR < 0.01  ~ "**",FDR < 0.05  ~ "*",TRUE ~ "ns"),
         log2OR=ifelse(p_star=="ns",0,log2OR)) %>% 
  mutate(VarSubType=factor(VarSubType,levels=rev(c("VNTR","STR","DEL","INS","BND","INV","DUP","SNV"))))

ggplot(fet_plot, aes(y = VarSubType, x = log2OR, fill = VarSubType)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +  
  geom_text(aes(label = p_star),
            position = position_dodge(width = 0.8),
            hjust = ifelse(fet_plot$log2OR >= 0, -0.1, 1.1), 
            vjust = 0.5, size = 5, color = "black") +
  facet_wrap(~score_type, scales = "free_x") +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none") +
  labs(y = "Variant Type", x = "log2(OR) compared to non-eGenes")+
  scale_fill_manual(values = c("SNV"="#c9c9c9","STR"="#0f3c7a","VNTR"="#ab889a","INS"="#b55f60","DEL"="#577b95"))


# Extended.Fig.4g: eGene constraint score correlated with effect size---------------------------------------------
library(broom)
bin_stats <- readRDS("ExtendFig4g.RDS")

fit_lines <- bin_stats %>% group_by(score_type) %>%
  do({fit <- lm(mean_score_adj ~ mean_effect, data = .)
  newdat <- data.frame(mean_effect = seq(min(.$mean_effect), max(.$mean_effect), length.out = 100))
  preds <- predict(fit, newdata = newdat, interval = "confidence", level = 0.95)
  data.frame(score_type = unique(.$score_type), mean_effect = newdat$mean_effect,
             fitted = preds[,"fit"],ci_low = preds[,"lwr"],ci_high = preds[,"upr"])}) %>%ungroup()

lm_results <- bin_stats %>%
  group_by(score_type) %>%
  do(tidy(lm(mean_score_adj ~ mean_effect, data = .))) %>%
  filter(term == "mean_effect") %>%
  ungroup() %>% mutate(FDR = p.adjust(p.value, method = "BH"))

label_df <- bin_stats %>%
  group_by(score_type) %>%
  summarise(x = max(mean_effect) * 0.95,y = max(mean_score) * 0.95) %>%
  left_join(lm_results, by = "score_type") %>%
  mutate(label = paste0("β = ", round(estimate, 3),"\nP = ", signif(p.value, 3)))

ggplot(bin_stats, aes(x = mean_effect, y = mean_score, color = score_type)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0) +
  geom_ribbon(data = fit_lines, aes(x = mean_effect, ymin = ci_low, ymax = ci_high, fill = score_type),
              alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = fit_lines, aes(x = mean_effect, y = fitted, color = score_type), size = 1) +
  geom_text(data = label_df,
            aes(x = x, y = y, label = label, color = score_type),
            hjust = 1, vjust = 1, size = 5, fontface = "bold") +
  scale_color_manual(values = c("pLI"="#bb0021","pRec"="darkgreen","pNull"="#3b4992")) +
  scale_fill_manual(values = c("pLI"="#bb0021","pRec"="darkgreen","pNull"="#3b4992")) +
  theme_classic(base_size = 14) +
  labs(x = "Absolute effect size", y = "Mean(Exac Score)") +
  theme(legend.position = "top")

# Extended.Fig.4h functional enrichment -----------------------------------

group_order <- c("Splice Site","5UTR","Coding Exon",  "3UTR", "Non-coding Exon",
                 "Intron", "Upstream","Downstream", "Intergenic")

group_labels <- c("Splice Site","5UTR","Coding Exon",  "3UTR", "Non-coding Exon",
                  "Intron", "Upstream","Downstream", "Intergenic")
results <- fread("./ExtendFig4h.txt")
results$Group <- factor(results$Group, levels = rev(group_order), labels = rev(group_labels))

ggplot(results, aes(x = OR, y = Group, color = QTL_type)) +
  geom_point(size = 3.5, position = position_dodge(0.7)) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high),
                 position = position_dodge(0.7), height = 0.2, size = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray30") +
  scale_color_manual(values = c("sQTL" = "#7c8bad", "eQTL" = "#9fc999")) +
  scale_x_log10() +
  labs(title = "Enrichment of QTL TRs in Genomic Regions",
       x = "Odds Ratio (95% CI)", 
       y = "",
       color = "QTL Type") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11),
    legend.position = "top")


# Extended.Fig.4i motif enrichment ----------------------------------------


motif_enrichment_results <- fread("./ExtendFig4i.txt")
significant_motifs <- motif_enrichment_results %>%
  filter(P_adjusted < 0.05, Enrichment_Fold > 1) %>%
  arrange(P_adjusted, desc(Enrichment_Fold))%>%
  mutate(sig_label = case_when(
    P_adjusted < 0.001 ~ "***",
    P_adjusted < 0.01  ~ "**",
    P_adjusted < 0.05  ~ "*",
    TRUE               ~ ""),
    Enrichment_Fold_log2 = log2(Enrichment_Fold),
    CI_lower_log2 = log2(CI_lower),
    CI_upper_log2 = log2(CI_upper))

ggplot(significant_motifs[1:20,], aes(x = reorder(Motif, Enrichment_Fold_log2), y = Enrichment_Fold_log2)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = CI_lower_log2, ymax = CI_upper_log2), width = 0.3) +
  geom_text(aes(y = Enrichment_Fold_log2 + 0.1, label = sig_label),) +
  coord_flip() +
  theme_classic() +
  labs(y = "log2(OR)")


# Extended.Fig.4j RBP enrichment ------------------------------------------


dat <- fread("./ExtendFig4j.txt")
colnames(dat) <- c("RBP","Tissue","peak","qtl_peak","shuffle_peak","qtl","pvalue","odds_ratio","ci_down","ci_up","stats")

dat <- dat %>%
  mutate(cell=word(RBP,start=1,end=1,sep="\\-"),
         RBP=word(RBP,start=2,end=2,sep="\\-"),
         odds_ratio=ifelse(pvalue < 0.05,odds_ratio,1),
         odds_ratio=ifelse(is.infinite(odds_ratio),1,odds_ratio))

dat[, logOR := log2(odds_ratio)]

cell_types <- unique(dat$cell)
cell_type="HepG2"
dat_cell <- dat[cell == cell_type]
motif_order <- dat_cell %>%
  group_by(RBP) %>%
  summarise(mean_logOR = mean(logOR, na.rm = TRUE),tissue_count = n_distinct(Tissue)  ) %>%
  arrange(desc(mean_logOR), desc(tissue_count)) %>% head(15) %>% 
  pull(RBP)

tissue_order<-c("Liver","Pancreas_Tail","Pancreas_Head","Spleen","Adrenal_Gland","Pancreas_Body","Adipose","Whole_Blood","Gallbladder","Muscle","Skin")

plot_data <- dat_cell %>%filter(RBP %in% motif_order) %>% 
  mutate(RBP = factor(RBP, levels = rev(motif_order)),
         Tissue = factor(Tissue, levels = rev(tissue_order)),
         signif_label = case_when(
           pvalue < 0.001 ~ "***",
           pvalue < 0.01 ~ "**",
           pvalue < 0.05 ~ "*",
           TRUE ~ ""), fill_value = logOR)

ggplot(plot_data, aes(x = Tissue, y = RBP, fill = fill_value)) +
  geom_tile(color = "white", linewidth = 0.1) +
  #geom_text(aes(label = signif_label), size = 3, color = "black",vjust = 0.8) +
  #scale_fill_gradient2(low = "#D73027", mid = "white",high = "#4575B4",midpoint = 0,
  #  name = "log2(OR)",na.value = "gray90") +
  scale_fill_steps2(low = "#D73027",mid = "white",
                    high = "#0a3a6e",midpoint = 0,
                    n.breaks = 20, name = "log2(OR)",na.value = "gray90" ) +
  labs() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.text.y = element_text(size = 8),legend.position = "right")







