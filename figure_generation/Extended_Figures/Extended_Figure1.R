#==============================================#
# Extended Fig.1 #
#==============================================#

setwd("/media/london_A/mengxin/GTOP_code/extend/extend_145/input")
library(data.table)
library(ggplot2)
library(stringi)
library(stringr)
library(dplyr)
library(ggsci)
library(tidyverse)
library(ggpubr)
library(magrittr)

# Extended.Fig.1a: variant number and length ----------------------------------------
dat <- fread("ExtendFig1a.PathogenicTR.txt")
gene_order <- dat %>%
  group_by(gene) %>%summarise(mean_CNV = mean(CNV_LRS, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_CNV)) %>%pull(gene)

dat <- dat %>%mutate(gene = factor(gene, levels = gene_order))

shared_bg <- dat %>%
  filter(class == "Shared") %>%
  group_by(gene) %>%
  summarise(xmin = 0, xmax = max(CNV_LRS), .groups = "drop") %>%
  mutate(gene = factor(gene, levels = gene_order))

max_CNV <- max(dat$CNV_LRS, na.rm = TRUE)

lrs_bg <- dat %>%
  filter(class != "Shared") %>%
  group_by(gene) %>%
  summarise(xmin = min(CNV_LRS), xmax = max_CNV, .groups = "drop") %>%
  mutate(gene = factor(gene, levels = gene_order))

map_CNV <- function(x) ifelse(x <= 50, x, 50 + (x - 50) / 10)

dat <- dat %>%mutate(CNV_plot = map_CNV(CNV_LRS))

shared_bg <- shared_bg %>%mutate(xmin_plot = map_CNV(xmin), xmax_plot = map_CNV(xmax))

lrs_bg <- lrs_bg %>%mutate(xmin_plot = map_CNV(xmin), xmax_plot = map_CNV(xmax))

breaks1 <- c(1,10,20,30,40,50)
labels1 <- as.character(breaks1)

max_val <- ceiling(max(dat$CNV_LRS, na.rm = TRUE) / 50) * 50
breaks2_orig <- seq(100, max_val, by = 50)
breaks2_mapped <- 50 + (breaks2_orig - 50) / 5
labels2 <- as.character(breaks2_orig)

breaks_all <- c(breaks1, breaks2_mapped)
labels_all <- c(labels1, labels2)
dat$class <- factor(dat$class, levels = c("LRS higher", "LRS specific", "Shared"))
p1 <- ggplot() +
  geom_rect(data = lrs_bg,
            aes( xmin = xmin_plot,xmax = xmax_plot,
                 ymin = as.numeric(gene) - 0.4,ymax = as.numeric(gene) + 0.4),
            fill = "grey90",alpha = 0.7) +
  geom_point(data = dat,aes(x = CNV_plot, y = gene, color = class),size = 1) +
  scale_color_manual(
    values = c("LRS specific" = "#994999",
               "LRS higher" = "red",
               "Shared" = "grey60")) +
  scale_x_continuous(breaks = breaks_all, labels = labels_all) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "CNV_LRS", y = "Gene", color = "Class") +
  coord_flip()
p1



