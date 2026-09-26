#==============================================#
# GTOP-eQTL-correlation-SNV #
# Supp-Figure-29#
#==============================================#
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(data.table)
library(corrplot)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig29.Figure.S21-internal-replication-eQTL")
# Supp.Fig.29a correlation of eQTL effects between pancreas -----------

dat.w <- fread("./input/Figure S29a.cov.txt")
rownames(dat.w) <- dat.w$replication_tissue
dat.w$replication_tissue <- NULL
corrplot(as.matrix(dat.w),method = "number",is.corr = F,col = COL1('Blues', 200)[50:200],
         type = "upper",diag = F)


min_val <- min(as.matrix(dat.w) ,na.rm = TRUE)
max_val <- max(as.matrix(dat.w), na.rm = TRUE)

scatterdf <- fread("./input/Figure S29a.scatter.txt")
#df_tissue_pairs <- scatterdf %>% distinct(rep, dis)

scatter_list <- list()

for(tissuePair in unique(scatterdf$tissuepair)){
  
  print(tissuePair)
  dat <- scatterdf %>% dplyr::filter(tissuepair==tissuePair)
  rho <- dat$rho[1]

  prop <- (rho - min_val) / (max_val - min_val)
  color_index <- round(1 + prop * (151 - 1))
  color <- COL1("Blues", 200)[50:200][color_index]

  p <- ggplot(dat,aes(x=x_var,y=y_var)) + geom_point(color = color, size = 0.5) + theme_pubr() + 
  ggtitle(label = tissuePair)
  scatter_list[[tissuePair]] <- p
}

cowplot::plot_grid(plotlist = scatter_list,ncol = 4,align = "hv")


# Supp.Fig.29b rb of eQTL between tissues -----------

df <- fread("./input/Figure S29b.txt")
df[, label := gsub("\\\\n", "\n", label)]

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

