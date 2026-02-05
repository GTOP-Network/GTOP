#===============================================#
# MS peptide  #
# Supp Figure 11 #
#===============================================#
library(ggplot2)
library(ggpubr)
library(tidyverse)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig11/input")

# Supp.Fig.11a ------------------------------------------------------------

plot_df <- fread("supp11a.MS_peptide_length_dist.txt")

ggplot(plot_df)+
  geom_bar(aes(x=peptide_length, y=count), stat = "identity")+
  xlab("Peptide length (aa)")+
  ylab("Number of peptides")+
  theme_pubr()


# Supp.Fig.b --------------------------------------------------------------

df <- fread("supp11b.MS_unique_peptide_support_novel_transcript.txt")
plot_df <- df %>% 
  pivot_longer(cols = -c("tissue"), names_to = "class",
               values_to = "number")
plot_df %>% 
  mutate(class=factor(class, levels=c("> 5", "2 - 5", "1"))) %>% 
  ggplot(.)+
  geom_bar(aes(x=tissue, y=number, fill=class), stat = "identity")+
  scale_fill_manual(values = c("1"="#cedbeb", "2 - 5"="#82aece", "> 5"="#4174ad"))+
  xlab("")+
  ylab("Number of novel transcripts")+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
