#==============================================#
# LRS WGS Population#
# Supp-Figure-5#
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
library(ggrepel)
library(scales)

setwd("/path/to/GTOP_code/supp/supp_fig5/input")


# Supp.Fig.5a: Gtop Genotype PCA with 1KGP and GTEx in EAS population ------------------------------------------------------------------

df_pca <- fread("Supp_Fig5a.txt",
                     header=TRUE,check.names=FALSE)
title <- colnames(df_pca)
df_pca$pop <- ifelse(df_pca$superpop=="GTEX","GTEX",df_pca$pop)
df_pca <- rbind(df_pca[df_pca$pop=="GTOP",],
                df_pca[df_pca$pop %in% c("CDX","CHB","CHS","JPT","KHV","GTEX"),])

df_pca$pop <- ifelse(df_pca$pop %in% c("CDX","CHB","CHS","JPT","KHV","GTEX"),df_pca$pop,"GTOP")
colors_0 <- c("CDX"="#f39b7f","GTOP"="#E64B35FF","CHB"="#00a087","GTEX"="#4DBBD5FF",
              "CHS"="#3d5488","KHV"="#8391b4","JPT"="#91d1c2")

p <- ggplot(df_pca,aes(x=PC1,y=PC2)) +
  geom_point(aes(color=pop)) +
  xlab('PC1 (0.8%)') + ylab('PC2 (0.4%)') +
  #scale_color_npg() +
  scale_color_manual(values=colors_0) +
  theme_classic() + 
  theme(axis.line=element_line(color='black'),
        axis.text=element_text(color='black'),
        axis.title=element_text(face='bold',size=rel(1.25)),
        legend.position=c(0.9,0.3),
        legend.text=element_text(face='bold',size=rel(1.1)),
        legend.key=element_blank() )

# Supp.Fig.5b: GTOP Genotype PCA  --------------------------------------------------------

df_pca <- read.table("Supp_Fig5b.pca_stats.txt",header=FALSE,row.names=1)
datap <- data.frame(t(df_pca))
datap$PC <- seq(1,dim(df_pca)[2])
scale_factor <- max(datap$prop_var, na.rm = TRUE) / max(datap$cumm_prop, na.rm = TRUE)
p <- ggplot(datap, aes(x = PC)) +
  geom_col(aes(y = prop_var), fill = "grey") +
  geom_line(aes(y = cumm_prop * scale_factor), color = "black", size = 1) +
  geom_point(aes(y = cumm_prop * scale_factor), color = "black", size = 1) +
  scale_y_continuous(
    name = "PVE",
    sec.axis = sec_axis(~ . / scale_factor,name = "Cumulative proportion explained",
                        labels = percent_format(accuracy = 1))) +
  labs(x = "Genotype PC index",) +
  theme_classic() +
  theme(axis.title.y.left  = element_text(color = "black"),
        axis.title.y.right = element_text(color = "black"))
p
# Supp.Fig.5c: ADMIXTURE cv error ---------------------------------------------------

cv_df <- fread("Supp_Fig5c.cv_error.txt")
p <- ggplot(cv_df, aes(x = K, y = CV)) +
  geom_line(color = "#3d5488", size = 1) +
  geom_point(color = "#3d5488", size = 3) +
  theme_classic() +
  scale_x_continuous(breaks = seq(min(cv_df$K), max(cv_df$K), by = 1))+
  labs(x = "Number of Ancestral Populations (K)",
       y = "Cross-Validation Error") +
  theme(
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold", size = 13))
p
# Supp.Fig.5d: variant number and length ----------------------------------------

k <- 5
datap <- read.table(
  "structure.K5.out"
)
col_popu  <- 3
col_start <- 4
popus   <- unique(datap[, col_popu])
popu_n  <- table(datap[, col_popu])
n_popu  <- length(popus)
poses_border <- seq_len(n_popu)
datap_q <- datap[, col_start:ncol(datap)]

colors_n <- c(652,83,115,312,258,257,28,632,103,116,595,598,114,131,497,24,376,511,141,137)
colors_0 <- colors()[colors_n[1:k]]
par(mar = c(4, 0, 0, 5), bty = "n", xpd = NA)
plot(c(0, n_popu), c(0, 1),type = "n",axes = FALSE,xlab = "", ylab = "")

x1 <- 0
for (i in seq_len(n_popu)) {
  idx <- datap[, col_popu] == popus[i]
  q <- as.matrix(datap_q[idx, , drop = FALSE])
  storage.mode(q) <- "numeric"
  n_i <- nrow(q)
  w   <- 1 / n_i
  for (j in seq_len(n_i)) {
    y0 <- 0
    for (kk in seq_len(k)) {
      rect( x1, y0,
        x1 + w, y0 + q[j, kk],
        col    = colors_0[kk],
        border = NA)
      y0 <- y0 + q[j, kk]
    }
    x1 <- x1 + w
  }
}
rect(0, 0, n_popu, 1)
abline(v = poses_border[-n_popu])
mtext("K = 5", side = 4, cex = 2, line = -1)
for (i in seq_len(n_popu)) {
  text(i - 0.5, -0.1, labels = popus[i],
       srt = 45, adj = c(1, NA), cex = 1)
}

axis(1, at = seq(0, n_popu), labels = FALSE, tcl = -0.8)










