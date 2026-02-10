#==============================================#
# evoluation #
# Supp-Figure-25#
#==============================================#

library(ggplot2)
library(tidyverse)

setwd("/path/to/GTOP_code/supp/supp_fig25")


# Supp.Fig.25a AF difference -----------------------------------------
diff_freq_group_df <- fread("./input/supp_fig25a_data.txt")

diff_freq_group_df$proxy_group <- factor(diff_freq_group_df$proxy_group, 
                                         levels = c("EAS > EUR", "EAS < EUR"))
ggplot(diff_freq_group_df, aes(x=proxy_group, y=count, fill=proxy_delta_group)) +
  geom_col(position = "fill") +
  scale_fill_manual(values =  c("#f8e4da", "#f1cf9c", "#e6986d", 
                                "#cf5e47", "#a43331", "#5c95d6")) +
  theme_classic() +
  labs(x="", y="Percentage of eQTLs")


# Supp.Fig.25b Number of variants --------------------------------------------
diff_freq_group_df <- fread("./input/supp_fig25b_data.txt")

ggplot(diff_freq_group_df, aes(x=new_Type_GTEx, y=allcount, 
                                fill=proxy_delta_group)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values =  c("#f8e4da", "#f1cf9c", "#e6986d", 
                                "#cf5e47", "#a43331", "#5c95d6")) +
  theme_classic() +
  labs(x="", y="Number of eQTLs")

