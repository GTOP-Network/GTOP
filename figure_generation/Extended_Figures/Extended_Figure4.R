#==============================================#
#  Characteristics of molecular QTLs across variant classes#
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

#input data include 3 column: ExpP Tissue ObservedP
p <- ggplot(df_plot,aes(x=ExpP,y=ObservedP,color=Tissue)) + 
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

dat <- fread("ExtendFig4b.txt")
dat$VarType<-"SNV"
dat[dat$VarSubType %in% c("STR","VNTR"),]$VarType<-"TR"
dat[dat$VarSubType %in% c("DEL","INS"),]$VarType<-"SV"

ggplot(dat,aes(x=start_distance,color=VarType)) + geom_density(size=1.5) + theme_pubr() +
  scale_color_manual(values = c("#7B88A8","#609561","#983927"))+
  labs(x="Relative distance from TSS(Mb)")

# Extended.Fig.4c:  torus enrichment ------------------------------------------------------------------

order <- rev(c("enhancer","promoter","open chromatin region","CTCF binding site","TF binding site","3 prime UTR","5 prime UTR","frameshift","intron","missense","NC transcript","splice acceptor","splice donor","splice region","stop gained", "synonymous"))

plotdf <- fread("Exfig 4c.left.txt")
plotdf$QTL<-factor(plotdf$QTL,levels=c("tuQTL","juQTL","eQTL"))
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
plotdf$type<-factor(plotdf$type,levels = c("tuQTL","juQTL","eQTL"))
p2 <- plotdf %>% 
  dplyr::mutate(feature=factor(feature, levels=order)) %>% 
  ggplot(.)+
  geom_bar(aes(x=meanratio, y=feature, fill=type), stat = "identity", position=position_dodge(width=0.9))+
  geom_errorbar(aes(xmax=meanratio+sdratiio, xmin=meanratio-sdratiio, y=feature, group=type), position=position_dodge(width=0.9),width=0)+
  scale_fill_manual(values=c("eQTL"="#a2bf98" , "juQTL"="#7d8bad",
                             "tuQTL"="#b47973"))+
  scale_x_break(c(0.125,0.7),
                space = 0.3,
                scales = .5)+
  theme_pubr(legend = "top")+
  xlab("Proportion of variants")+
  ylab("")+
  theme(axis.text.y = element_blank())
p2

# Extended Fig.4d:  effect size distrbution ------------------------------------------------------------------

effect_data <- fread("ExtendFig4d.txt.gz")
effect_data$VarSubType<-factor(effect_data$VarSubType,levels = rev(c("INS","DEL","VNTR","STR","SNV")))
effect_data$QTL<-factor(effect_data$QTL,levels = rev(c("eQTL","juQTL","tuQTL")))
p1 <- ggplot(effect_data)+ geom_violin(aes(y=QTL, x=slope, fill=VarSubType), 
                                       position = position_dodge())+
  xlab("Effect Size")+
  ylab("Molecular Trait")+
  theme_classic()+theme(legend.position = "none")+
  scale_fill_manual(values = c("SNV"="#7e8daf","STR"="#0f3c7a","VNTR"="#ab889a","INS"="#b55f60","DEL"="#577b95"))
num <- effect_data %>% group_by(QTL,VarSubType) %>% dplyr::summarise(num=n()) %>% arrange(num)
p2 <- num %>% mutate(QTL=factor(QTL,levels=c("tuQTL","juQTL","eQTL")),
                     VarSubType=factor(VarSubType,levels=rev(c("INS","DEL","VNTR","STR","SNV")))) %>% 
  ggplot(aes(y = QTL, x = num, fill=VarSubType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("SNV"="#7e8daf","STR"="#0f3c7a","VNTR"="#ab889a","INS"="#b55f60","DEL"="#577b95"))+
  scale_x_log10()+
  labs(x = "QTL count", y = "") +
  theme_classic()
cowplot::plot_grid(p1,p2,ncol = 2)

# Extended.Fig.4e: eGene constraint score ---------------------------------------------------

fet_plot <- fread("ExtendFig4e.eGene_constraint.txt") %>%
  dplyr::filter(score_type %in% c("pLI"))%>%
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
  scale_fill_manual(values = c("SNV"="#7e8daf","STR"="#0f3c7a","VNTR"="#ab889a","INS"="#b55f60","DEL"="#577b95"))


# Extended.Fig.4f: eGene constraint score correlated with effect size---------------------------------------------
library(broom)
bin_stats <- readRDS("ExtendFig4f.RDS")

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

# Extended.Fig.4g  TR-sQTL functional enrichment -----------------------------------

dat <- fread("ExtendFig4g.txt")
group_order <-  c("5UTR","Splice_site","3UTR","Coding_exon", "Intron", "Upstream","Downstream", "Intergenic")
group_labels <- c("5UTR","Splice_site","3UTR","Coding_exon", "Intron", "Upstream","Downstream", "Intergenic")
dat$Group <- factor(dat$Group, levels = rev(group_order), labels = rev(group_labels))
p1 <- ggplot(dat, aes(x = OR, y = Group, color = QTL_type)) +
  geom_point(size = 3.5, position = position_dodge(0.7)) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high),
                 position = position_dodge(0.7), height = 0.2, size = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray30") +
  scale_color_manual(values = c("sQTL" = "#7c8bad", "eQTL" = "#9fc999")) +
  scale_x_log10() +
  labs(x = "Odds Ratio (95% CI)",  y = "",color = "QTL Type") +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11),
    legend.position = "top");p1

plot <- dat %>% dplyr::select(Group, QTL_type, Proportion_in_QTL,Count_in_QTL,Total_in_genome)%>%
  mutate(Proportion = Proportion_in_QTL,Ratio = Count_in_QTL / Total_in_genome )
plot$Group <- factor(plot$Group,levels = c("Intergenic","Downstream","Upstream","Intron", "Coding_exon","3UTR","Splice_site","5UTR"))
scale_factor <- max(plot$Proportion) / max(plot$Ratio)

p2 <- ggplot(plot, aes(y = Group)) +
  geom_bar(aes(x = Proportion*100, fill = QTL_type),
           stat = "identity",
           position = position_dodge(width = 0.75),
           width = 0.7 ) +
  labs(x = NULL,y = "Proportion in QTLs (%)" ) +
  scale_x_break(c(15, 60),scales = "free") +
  scale_x_continuous(breaks = seq(0, 70, by = 5))+
  theme_classic(base_size = 13) +
  scale_fill_manual(values = c("eQTL"="#9bc396","sQTL"="#7a88a8"))+
  scale_color_manual(values = c("eQTL"="#9bc396","sQTL"="#7a88a8"))+
  theme(axis.text = element_text(colour = "black")) +
  theme(legend.position = "top")+
  theme(
    axis.text.x.top = element_blank(),   
    axis.ticks.x.top = element_blank(), 
    axis.line.x.top = element_blank()   
  );p2

library(patchwork)
p1 + p2
# Extended.Fig.4h motif enrichment ----------------------------------------

motif_enrichment_results <- fread("./ExtendFig4h.txt")
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
  geom_bar(stat = "identity",fill="#7c8aaa") +
  geom_errorbar(aes(ymin = CI_lower_log2, ymax = CI_upper_log2), width = 0.3) +
  geom_text(aes(y = Enrichment_Fold_log2 + 0.1, label = sig_label),) +
  coord_flip() +
  theme_classic() +
  labs(y = "log2(OR)")


# Extended.Fig.4i RBP enrichment ------------------------------------------


plot_data <- fread("./ExtendFig4i.txt")
motif_order <- plot_data %>%
  group_by(RBP) %>%
  summarise(mean_logOR = mean(logOR, na.rm=TRUE),tissue_count = n_distinct(Tissue)) %>%
  arrange(desc(mean_logOR), desc(tissue_count)) %>%
  pull(RBP)

tissue_order <- c("Pancreas_Body","Pancreas_Tail","Pancreas_Head","Whole_Blood","Adipose",
                  "Adrenal_Gland","Gallbladder","Liver","Muscle","Skin","Spleen")

plot_data <- plot_data %>%
  mutate(RBP = factor(RBP, levels = rev(motif_order)), Tissue = factor(Tissue, levels = tissue_order))

ggplot(plot_data,aes(x=Tissue,y=RBP,fill=logOR))+
  geom_tile(color="white",linewidth=0.1)+
  scale_fill_steps2(low="#D73027",mid="white", high="#0a3a6e",
                    midpoint=0,n.breaks=6, name="log2(OR)",na.value="gray90")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90,hjust=1,size=10),
        axis.text.y=element_text(size=8),legend.position="right" )







