#===============================================#
# ASE longcallR #
# Supp-Figure-20                         #
#===============================================#
library(ggplot2)
library(ggpubr)
library(ComplexUpset)
library(tidyverse)
library(data.table)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig20.ASE-longcallR/input")


# Supp.Fig.20a --the overlap of genes with significant ASE or ASTS events----------------------------
library(ggVennDiagram)
df_list<-readRDS("supp20a.RDS")
p4 <- ggVennDiagram(df_list,force_upset = T,order.set.by = "name");p4


# Supp.Fig.20b --comparison of allelic fold change for overlapped ASE events--------------------------

m12.nna<-readRDS("supp20b.RDS")
p5 <- ggplot(m12.nna,aes(x=logafc,y=logFC_trans,color=Significant)) + geom_point(size=2.5) + theme_classic() + 
  xlab("Loras aFC") + ylab("LongcallR aFC") + geom_hline(yintercept = 0,linetype="dashed") + geom_vline(xintercept = 0,linetype="dashed") + 
  scale_color_manual(values = c("grey","darkblue"));p5


# Supp.Fig.20c --the proportion of single-donor ASE events--------------------
df_stat.ase3<-readRDS("supp20c.RDS")
df_g1 <- df_stat.ase3 %>% filter(count == 1)


group_2 <- df_g1 %>% filter(count_total == 1)   
group_3 <- df_g1 %>% filter(count_total > 1)    

group_2_4 <- group_2 %>% filter(count_sample == 1)  
group_2_5 <- group_2 %>% filter(count_sample > 1)   

group_3_4 <- group_3 %>% filter(count_sample == 1)
group_3_5 <- group_3 %>% filter(count_sample > 1)

cat("group_1 -> group_2:", nrow(group_2), "\n")   
cat("group_1 -> group_3:", nrow(group_3), "\n")   
cat("group_2 -> group_4:", nrow(group_2_4), "\n") 
cat("group_2 -> group_5:", nrow(group_2_5), "\n") 
cat("group_3 -> group_4:", nrow(group_3_4), "\n") 
cat("group_3 -> group_5:", nrow(group_3_5), "\n") #

# group_1: ASE sig in one donor; group_2: tested in one donor; group_3:tested in >1 donors; group_4: test in one tissue; group_5: test in > 1 tissues 
links <- data.frame(
  source=c("group_1","group_1", "group_2", "group_2", "group_3", "group_3"), 
  target=c("group_2","group_3", "group_4", "group_5", "group_4", "group_5"), 
  value  = c(
    nrow(group_2),
    nrow(group_3),
    nrow(group_2_4),
    nrow(group_2_5),
    nrow(group_3_4),
    nrow(group_3_5)
  )
)

nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
library(networkD3)
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE)

p


# Supp.Fig.20d -----Enrichment of genes with ASE or ASTS events-----------------------------------

load("supp20d.RData")

p<-ggplot(df_plot, aes(x = group, y = OR)) +
  geom_pointrange(
    aes(ymin = CI_low, ymax = CI_high),
    size = 0.8
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Odds Ratio (95% CI)",
    title = "Enrichment of ASE / ASTS signals in QTLs"
  );p
