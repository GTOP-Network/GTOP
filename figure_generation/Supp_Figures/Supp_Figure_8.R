#==============================================#
# VEP annotation #
# Supp-Figure-8#
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
setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig2567/input")

color <- c(
  "intron"="#5574A6","intergenic"="#cbcbcb",
  "3UTR"="#329262","Coding"= "#02599C",
  "Other"="#cccccc","downstream"="#AFC8E2",
  "non_coding_transcript_exon"="#DEECF6","non_coding_transcript"="#DEECF6",
  "splice_acceptor"="#F9D5D5","splice_donor"="#FF7F50",
  "splice_region"="#8A4B43",
  "upstream"="#A2B5CD","5UTR"="#651067",
  "NMD_transcript_variant"="#003264",
  "Coding"="#25a7e1","coding_sequence"="#25a7e1","coding_transcript"="#25a7e1",
  "incomplete_terminal_codon"="#C8A096",
  "stop_lost"="#B482FF",
  "start_lost"="#FFAAC8",
  "synonymous"="#A0B4AA",
  "missense"="#465A6E",
  #"intron_variant"="#82966E",
  "stop_gained"="#CBE5C7",
  "frameshift"="#ff69b4","inframe_insertion"="#ff69b4","inframe_deletion"="#ff69b4",
  "protein_altering"="#66aa00",
  "synonymous"="#76ee00",
  "stop_retained"="#ff0000",
  "transcript_amplification"="#cccccc",
  "transcript_ablation"="#cccccc",
  "start_retained"="#CBE5C7")

# Supp.Fig.8a:  VEP annotation ------------------------------------------------------------------
vep <- fread("Supp_Fig8a.txt")
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


# Supp.Fig.8b: SV subtype feature  --------------------------------------------------------

subtype_pro <- fread("Supp_Fig8b.txt") %>%
  dplyr::filter(type %in% c("INS","DEL","DUP","INV"))%>%
  mutate(Consequence=str_replace_all(Consequence,"_variant",""),
         Consequence = recode(Consequence,
                              "5_prime_UTR" = "5UTR",
                              "3_prime_UTR" = "3UTR",
                              "downstream_gene" = "downstream",
                              "upstream_gene" = "upstream",
                              # splice related
                              "splice_donor_region" = "splice_donor",
                              "splice_donor_5th_base" = "splice_donor",
                              "splice_polypyrimidine_tract" = "splice_region",
                              "inframe_insertion"="frameshift","inframe_deletion"="frameshift",
                              "transcript_amplification"="Other",
                              "transcript_ablation"="Other","incomplete_terminal_codon"= "Other",
                              "mature_miRNA"="Other",
                              "coding_sequence"="Coding","coding_transcript"="Coding",
                              # keep unchanged
                              .default = Consequence))

region_order <- subtype_pro %>% group_by(Consequence) %>%summarise(mean_pro = mean(pro)) %>%arrange(mean_pro) %>%pull(Consequence)
subtype_pro$Consequence <- factor(subtype_pro$Consequence,levels = region_order)
subtype_pro$type <- factor(subtype_pro$type,levels = c("INS","DEL","DUP","INV"))

p1 <- ggplot(subtype_pro, aes(x = type, y = count, fill = Consequence)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "SV Type", y = "Number of Sites", fill = "Region") +
  theme_classic() +
  theme(legend.position = "top")+
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))+
  scale_fill_manual(values = color)+
  coord_flip()

p2 <- ggplot(subtype_pro, aes(x = type, y = pro, fill = Consequence)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "SV Type", y = "Percentage of Sites", fill = "Region") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))+
  theme(legend.position = "top")+
  scale_fill_manual(values = color)+
  coord_flip()

cowplot::plot_grid(p2,p1,ncol = 2)


# Supp.Fig.8c: TR subtype feature  --------------------------------------------------------

subtype_pro <- fread("Supp_Fig8c.txt") %>%
  mutate(Consequence=str_replace_all(Consequence,"_variant",""),
         Consequence = recode(Consequence,"5_prime_UTR" = "5UTR",
                              "3_prime_UTR" = "3UTR",
                              "downstream_gene" = "downstream",
                              "upstream_gene" = "upstream",
                              # splice related
                              "splice_donor_region" = "splice_donor",
                              "splice_donor_5th_base" = "splice_donor",
                              "splice_polypyrimidine_tract" = "splice_region",
                              "inframe_insertion"="frameshift","inframe_deletion"="frameshift",
                              "transcript_amplification"="Other",
                              "transcript_ablation"="Other","incomplete_terminal_codon"= "Other",
                              "mature_miRNA"="Other",
                              # keep unchanged
                              .default = Consequence))

region_order <- subtype_pro %>% group_by(Consequence) %>%summarise(mean_pro = mean(pro)) %>%arrange(mean_pro) %>%pull(Consequence)
subtype_pro$Consequence <- factor(subtype_pro$Consequence,levels = region_order)
subtype_pro$type <- factor(subtype_pro$type,levels = c("2","3","4","5","6","VNTR"))

p1 <- ggplot(subtype_pro, aes(x = type, y = count, fill = Consequence)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "TR Type", y = "Number of Sites", fill = "Region") +
  theme_classic() +
  theme(legend.position = "top")+
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))+
  scale_fill_manual(values = color)+
  coord_flip()

p2 <- ggplot(subtype_pro, aes(x = type, y = pro, fill = Consequence)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "TR Type", y = "Percentage of Sites", fill = "Region") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))+
  theme(legend.position = "top")+
  scale_fill_manual(values = color)+
  coord_flip()
cowplot::plot_grid(p2,p1,ncol = 2)

