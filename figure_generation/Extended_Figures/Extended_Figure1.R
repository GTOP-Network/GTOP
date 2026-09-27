#==============================================#
# Extended Fig.1 #
#==============================================#

setwd("/path/to/GTOP_code/extend/extend_145/input")
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
lrs_bg <-dat %>%
  filter(class != "Shared") %>%
  group_by(gene) %>%
  summarise(xmin = min(CNV_LRS), xmax = max(dat$CNV_LRS, na.rm = TRUE), .groups = "drop") %>%
  mutate(gene = factor(gene, levels = gene_order))


map_CNV <- function(x){
  ifelse(x <= 50, x, 50 + (x - 50) / 10)
}

valid_genes <- dat %>%
  mutate( CNV_plot = map_CNV(CNV_LRS), class = factor(class,levels = c("LRS higher", "LRS specific", "Shared") ) )
shared_bg <- shared_bg %>% mutate( xmin_plot = map_CNV(xmin), xmax_plot = map_CNV(xmax))
lrs_bg <- lrs_bg %>% mutate( xmin_plot = map_CNV(xmin), xmax_plot = map_CNV(xmax))
breaks1 <- c(1,10,20,30,40,50)
max_val <- ceiling(max(valid_genes$CNV_LRS, na.rm = TRUE)/50)*50
breaks2_orig <- seq(100,max_val,by=50)
breaks2_mapped <- map_CNV(breaks2_orig)
breaks_all <- c(breaks1, breaks2_mapped)
labels_all <- c( as.character(breaks1),as.character(breaks2_orig))

p1 <- ggplot() +
  geom_rect(data = shared_bg,aes( xmin = xmin_plot,xmax = xmax_plot,
                                  ymin = as.numeric(gene)-0.4, ymax = as.numeric(gene)+0.4 ), fill="grey90", alpha=0.7) +
  geom_rect( data = lrs_bg,aes(xmin=xmin_plot,xmax=xmax_plot,ymin=as.numeric(gene)-0.4,ymax=as.numeric(gene)+0.4), fill="grey90",alpha=0.7) +
  geom_point( data=valid_genes,aes( x=CNV_plot,y=gene, color=class ),size=1 ) +
  scale_color_manual( values=c("LRS specific"="#994999", "LRS higher"="red", "Shared"="grey60") ) +
  scale_x_continuous( breaks=breaks_all, labels=labels_all, limits=c(0,max(breaks_all))) +
  theme_classic() +
  theme(axis.text=element_text(color="black"),axis.text.x=element_text(angle=90, hjust=1)) +
  labs(x="CNV_LRS",y="Gene", color="Class" ) +
  
  coord_flip()


p1



