#==============================================#
# GTOP-eQTL-correlation-SNV #
# Supp-Figure-22#
#==============================================#
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(data.table)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig22")
# Supp.Fig.22a correlation of eQTL effects between pancreas -----------



# Supp.Fig.22b rb of eQTL between tissues -----------

df <- fread("./input/Figure S22b.txt")

ggplot(df, aes(x = tis1, y = tis2, fill = r_b)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), 
            size = 3, 
            color = "black") +
  scale_fill_gradient2(
    high = "#2171B5",
    mid = "#9ECAE1",
    low = "white",
    midpoint = median(df$r_b, na.rm = TRUE),
    name = "Rb"
  )+
  labs(
    x = "Discovery",
    y = "Replication",
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  coord_fixed(ratio = 1)

