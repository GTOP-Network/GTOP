#==============================================#
# evoluation #
# Supp-Figure-27#
#==============================================#

library(ggplot2)
library(tidyverse)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig26")


# Supp.Fig.26a AF difference -----------------------------------------
diff_freq_group_df <- fread("./input/supp_fig26a_data.txt")

diff_freq_group_df$proxy_group <- factor(diff_freq_group_df$proxy_group, 
                                         levels = c("EAS > EUR", "EAS < EUR"))
ggplot(diff_freq_group_df, aes(x=proxy_group, y=count, fill=proxy_delta_group)) +
  geom_col(position = "fill") +
  scale_fill_manual(values =  c("#f8e4da", "#f1cf9c", "#e6986d", 
                                "#cf5e47", "#a43331", "#5c95d6")) +
  theme_classic() +
  labs(x="", y="Percentage of eQTLs")


# Supp.Fig.26b Number of variants --------------------------------------------
diff_freq_group_df <- fread("./input/supp_fig26b_data.txt")

ggplot(freq_group_count_df, aes(x=Type_GTEx, y=allcount, 
                                fill=proxy_delta_group)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values =  c("#f8e4da", "#f1cf9c", "#e6986d", 
                                "#cf5e47", "#a43331", "#5c95d6")) +
  theme_classic() +
  labs(x="", y="Number of eQTLs")

