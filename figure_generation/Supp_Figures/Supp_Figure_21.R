
#==============================================#
# tuQTL & tuVariants #
# Supp-Figure-2#
#==============================================#

library(patchwork)
library(ggplot2)
library(ggpubr)
library(ggbreak)
setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig21")

# Supp.Fig.21a left torus result of known/novel tuQTL ---------------------
order <- rev(c("enhancer","promoter","open chromatin region","CTCF binding site","TF binding site","3 prime UTR","5 prime UTR","frameshift","intron","missense","NC transcript","splice acceptor","splice donor","splice region","stop gained", "synonymous"))
plotdf <- fread("./input/Figure S21a.left.txt")

p1 <- plotdf %>%
  mutate(Ann=factor(Ann, levels=order)) %>% 
  ggplot(.) +
  geom_pointrange(aes(x=Ann, y = logmFC, ymin=lFC, ymax=hFC, color = QTL, shape=QTL), 
                  position=position_dodge(width=0.8), size = 0.5)+
  scale_color_manual(values=rev(c("tuQTL novel"="#20407b" , "tuQTL known"="#8892ae")))+
  scale_shape_manual(values = c("tuQTL novel"=16, "tuQTL known"=15))+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  # ylim(-8,20)+
  ylab(expression("Log"[2]*"(Fold Enrichment)"))+
  coord_flip()+
  theme_pubr(legend = "right")+
  theme(
    axis.text.x = element_text(color="black"),
    axis.text.y = element_text(color="black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
p1


# Supp.Fig.21a right known or novel tuVariants proportion -----------------


plotdf <- fread("./input/Figure S21a.right.txt")

p2 <- plotdf %>% 
  mutate(feature=factor(feature, levels=order)) %>% 
  ggplot(.)+
  geom_bar(aes(x=meanratio, y=feature, fill=type), stat = "identity", position=position_dodge(width=0.9))+
  geom_errorbar(aes(xmax=meanratio+sdratiio, xmin=meanratio-sdratiio, y=feature, group=type), position=position_dodge(width=0.9),width=0)+
  scale_fill_manual(values=rev(c("novel"="#20407b" , "known"="#8892ae")))+
  # scale_x_break(c(0.1,0.70),
  #               space = 0.3,
  #               scales = 1.5)+
  theme_pubr(legend = "right")+
  xlab("Proportion of variants")+
  ylab("")
p2


# Supp.Fig.21b distribution of credible set variants relative to s --------

density_df <- fread("./input/Figure S21b.txt")

density_df %>% 
  ggplot(., aes(x = pos_index, color=qtltype)) +
  geom_density( alpha = 0.5, adjust = 0.5) +
  labs(x = "Position", y = "Density of sQTL") +
  scale_color_manual(values = c("SV" ="#477980" , "TR" = "#9d3928", "Small_variant" = "#7d8bad"))+
  geom_vline(xintercept = 50, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = 71, linetype = "dashed", color = "black", alpha = 0.5) +
  theme_pubr()
