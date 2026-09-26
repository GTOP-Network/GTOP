#==============================================#
# LRS specific QTL tagged by nearby samll variant #
# Supp-Figure-27#
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
library(scales)
library(ggrastr)
library(ggupset)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig27. R1.7-LRS-specific.QTL/input/")

# Supp.Fig.27a: SV eQTL --------------------------------------------------------

dat <- fread("Supp_Fig27a.sv_eqtl.txt") %>% mutate(R2=ifelse(is.na(R2),0,R2))%>%
  mutate(combo_group = case_when( Class == "Overlap"  ~ "Overlap (all)",Class == "Specific" & R2 >= 0.8 ~ "Specific_HighLD",
                                  Class == "Specific" & R2 < 0.8  ~ "Specific_LowLD", TRUE    ~ NA_character_ )) %>%
  distinct(tissue, phenotype_id, variant_id, combo_group) %>%
  group_by(tissue, combo_group) %>%
  summarise(n_pairs = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(total_pairs_tissue = sum(n_pairs),proportion = n_pairs / total_pairs_tissue) %>%
  ungroup()

ggplot(dat %>% mutate(tissue = factor(tissue, levels = dat %>%
                                                  group_by(tissue) %>%
                                                  summarise(total_variants = sum(n_pairs), .groups = "drop") %>%
                                                  arrange(desc(total_variants)) %>%
                                                  pull(tissue)),
                                combo_group  = factor(combo_group, levels = c( "Specific_LowLD","Specific_HighLD","Overlap (all)"))), 
       aes(x = tissue, y = n_pairs, fill = combo_group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "",y = "Number of eQTLs",fill = "Class" ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "right")+
  scale_fill_manual(values = c("Specific_HighLD"="#cf928f","Specific_LowLD"="#933628","Overlap (all)"= "#8b9dc5"))

# Supp.Fig.27b: TR eQTL --------------------------------------------------------

dat <- fread("Supp_Fig27b.tr_eqtl.txt") %>% mutate(R2=ifelse(is.na(R2),0,R2))%>%
  mutate(combo_group = case_when( Class == "Overlap"  ~ "Overlap (all)",Class == "Specific" & R2 >= 0.8 ~ "Specific_HighLD",
                                  Class == "Specific" & R2 < 0.8  ~ "Specific_LowLD", TRUE    ~ NA_character_ )) %>%
  distinct(tissue, phenotype_id, variant_id, combo_group) %>%
  group_by(tissue, combo_group) %>%
  summarise(n_pairs = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(total_pairs_tissue = sum(n_pairs),proportion = n_pairs / total_pairs_tissue) %>%
  ungroup()


ggplot(dat %>% mutate(tissue = factor(tissue, levels = dat %>%
                                                             group_by(tissue) %>%
                                                             summarise(total_variants = sum(n_pairs), .groups = "drop") %>%
                                                             arrange(desc(total_variants)) %>%
                                                             pull(tissue)),
                                           combo_group  = factor(combo_group, levels = c( "Specific_LowLD","Specific_HighLD","Overlap (all)"))), 
                  aes(x = tissue, y = n_pairs, fill = combo_group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "",y = "Number of eQTLs",fill = "Class" ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "right")+
  scale_fill_manual(values = c("Specific_HighLD"="#cf928f","Specific_LowLD"="#933628","Overlap (all)"= "#8b9dc5"))


