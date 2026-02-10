#==============================================#
# Sites #
# Supp-Figure-7#
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

setwd("/path/to/GTOP_code/supp/supp_fig7/input")

color <- c(
  "intron"="#5574A6","intergenic"="#cbcbcb",
  "3UTR"="#329262","Coding"= "#02599C",
  "Other"="#cccccc","downstream"="#AFC8E2",
  "non_coding_transcript_exon"="#DEECF6","non_coding_transcript"="#DEECF6",
  "splice_acceptor"="#F9D5D5","splice_donor"="#FF7F50",
  "splice_region"="#8A4B43",
  "upstream"="#A2B5CD","5UTR"="#651067",
  "NMD_transcript_variant"="#003264",
  "coding_sequence"="#25a7e1",
  "incomplete_terminal_codon"="#C8A096",
  "stop_lost"="#B482FF",
  "start_lost"="#FFAAC8",
  "synonymous_variant"="#A0B4AA",
  "missense_variant"="#465A6E",
  "intron_variant"="#82966E",
  "stop_gained"="#CBE5C7",
  "frameshift"="#ff69b4","inframe_insertion"="#ff69b4","inframe_deletion"="#ff69b4",
  "protein_altering"="#66aa00",
  "synonymous"="#76ee00",
  "stop_retained"="#ff0000",
  "transcript_amplification"="#CBE5C7",
  "transcript_ablation"="#CBE5C7",
  "start_retained"="#CBE5C7")

# Supp.Fig.7a:  VEP annotation ------------------------------------------------------------------

vep <- fread("Supp_Fig7a.txt")
consequence_order <- vep %>%
  group_by(consequence_severe) %>%
  summarise(total_percentage = sum(percentage)) %>%
  arrange(total_percentage) %>%  
  pull(consequence_severe)

vep$consequence_severe <- factor(vep$consequence_severe,levels = consequence_order)
vep$type <- factor(vep$type,levels = c("SNV","SV","TR"))

ggplot(vep, aes(x = type, y = percentage, fill = consequence_severe)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "All Variants", y = "Percentage of Sites", fill = "Region") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))+
  theme(legend.position = "right")+
  scale_fill_manual(values = color)


# Supp.Fig.7b: SV subtype feature  --------------------------------------------------------

type_pro <- fread("Supp_Fig7b.txt")
type_pro$region <- factor(type_pro$consequence_severe,
                          levels = rev(c("Intron","Intergenic","Upstream",
                                         "Downstream","3UTR","Coding","5UTR")))
type_pro$type <- factor(type_pro$type,levels = c("INS","DEL","DUP","INV","BND"))

p1 <- ggplot(type_pro, aes(x = type, y = count, fill = region)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "SV Type", y = "Number of Sites", fill = "Region") +
  theme_classic() +
  theme(legend.position = "top")+
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))+
  coord_flip()

p2 <- ggplot(type_pro, aes(x = type, y = pro, fill = region)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "SV Type", y = "Percentage of Sites", fill = "Region") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))+
  theme(legend.position = "top")+
  coord_flip()
cowplot::plot_grid(p2,p1,ncol = 2)


# Supp.Fig.7c: TR subtype feature  --------------------------------------------------------

subtype_pro <- fread("Supp_Fig7c.txt")

subtype_pro$region <- factor(subtype_pro$consequence_severe,
                             levels = rev(c("Intron","Intergenic","Upstream",
                                            "Downstream","3UTR","Coding","5UTR")))
subtype_pro$type <- factor(subtype_pro$type,levels = c("2","3","4","5","6","VNTR"))

p1 <- ggplot(subtype_pro, aes(x = type, y = count, fill = region)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "TR Type", y = "Number of Sites", fill = "Region") +
  theme_classic() +
  theme(legend.position = "top")+
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))+
  coord_flip()

p2 <- ggplot(subtype_pro, aes(x = type, y = pro, fill = region)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "TR Type", y = "Percentage of Sites", fill = "Region") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))+
  theme(legend.position = "top")+
  coord_flip()
cowplot::plot_grid(p2,p1,ncol = 2)

