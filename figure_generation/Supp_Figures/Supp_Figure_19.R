#==============================================#
# Proportion of expression variance captured by the PCs  #
# Supp-Figure-19#
#==============================================#

library(data.table)
library(dplyr)
library(ggpubr)
library(cowplot)

setwd("/path/to/GTOP_code/supp/supp_fig19")


# Supp.Fig.19 Cumulative PVE-------------------------------------------------


df <- fread("./input/Figure S19.txt")
figlist <- list()
type <- unique(df$qtltype)[c(2, 1, 3)]

for (qtltype_t in type) {
  plotdf <- df %>% dplyr::filter(qtltype==qtltype_t)
  
  plotdf$tissue <- factor(plotdf$tissue, levels = plotdf$tissue)
  
  
  scale_factor <- max(plotdf$num_PCs) / max(plotdf$Cumulative_PVE)
  
  figlist[[qtltype_t]] <- ggplot(plotdf, aes(x = tissue)) +
    geom_bar(aes(y = Cumulative_PVE, fill = "pro_explain"), 
             stat = "identity") +
    geom_point(aes(y = num_PCs / scale_factor, color = "PCnumber"),
               shape = 19) +
    scale_y_continuous(
      name = "Proportion of expression variance captured by the PCs",
      limits = c(0, 0.8),
      sec.axis = sec_axis(~. * scale_factor, 
                          name = "Number of PCs in xQTL mapping")
    ) +
    scale_fill_manual(values = c("pro_explain" = "#c9c9c9")) +
    scale_color_manual(values = c("PCnumber" = "#9d3929")) +
    theme_pubr()+
    theme(axis.text.x = element_blank(),
          legend.position = "none")+
    ggtitle(qtltype_t)
}

cowplot::plot_grid(plotlist = figlist[type], ncol = 1)




