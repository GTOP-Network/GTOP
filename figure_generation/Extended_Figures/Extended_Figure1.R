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

# Extended.Fig.1a:  SV length distrbution ------------------------------------------------------------------

df_sv <- readRDS("ExtendFig1a.SV_length.profile.LRS_vs_SRS.RDS")
df_sv.del <- df_sv %>% filter(SVtype=="DEL")
df_sv.ins <- df_sv %>% filter(SVtype=="INS")
df_sv.del$abs_len <- abs(as.numeric(df_sv.del$length))
df_sv.ins$length <- as.numeric(df_sv.ins$length)

df_sv.del$lenBin <- cut(df_sv.del$abs_len,breaks = c(49,60,70,80,90,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
                                                     11000,12000,13000,14000,15000,
                                                     20000,25000,30000,35000,40000,45000,50000,100000,200000,300000,400000,500000,1000000,1936387),
                        labels = c(60,70,80,90,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
                                   11000,12000,13000,14000,
                                   15000,20000,25000,30000,35000,40000,45000,50000,100000,200000,300000,400000,500000,1000000,2000000))
df_sv.ins$lenBin <- cut(df_sv.ins$length,breaks = c(49,60,70,80,90,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
                                                    11000,12000,13000,14000,15000,
                                                    20000,25000,30000,35000),
                        labels = c(60,70,80,90,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
                                   11000,12000,13000,14000,15000,20000,25000,
                                   30000,35000))

df_sv.del$lenBin <- as.numeric(as.character(df_sv.del$lenBin))
df_sv.ins$lenBin <- as.numeric(as.character(df_sv.ins$lenBin))

df_plot.del <- as.data.frame(table(df_sv.del$lenBin,df_sv.del$tech))
df_plot.ins <- as.data.frame(table(df_sv.ins$lenBin,df_sv.ins$tech))

df_plot.del$Var1 <- as.numeric(as.character(df_plot.del$Var1))
df_plot.ins$Var1 <- as.numeric(as.character(df_plot.ins$Var1))
p1 <- ggplot(df_plot.del,aes(x=Var1,y=Freq,color=Var2)) + geom_line(size=1.5) + theme_pubr() + 
  scale_y_log10(breaks=c(1,10,100,1000,10000),
                labels=c("1","10","100","1K","10K")) + 
  scale_x_log10(breaks=c(100,1000,10000,100000,500000,100000000),
                labels=c("100","1K","10K","100K","500K",">=1M")) + 
  scale_color_manual(breaks = c("LRS","SRS"),values = c("#9c3929","#1c9099"));p1

p2 <- ggplot(df_plot.ins,aes(x=Var1,y=Freq,color=Var2)) + geom_line(size=1.5) + theme_pubr() + 
  scale_y_log10(breaks=c(1,10,100,1000,10000),
                labels=c("1","10","100","1K","10K")) + 
  scale_x_log10(breaks=c(100,1000,10000,100000),
                labels=c("100","1K","10K","100K")) + 
  scale_color_manual(breaks = c("LRS","SRS"),values = c("#9c3929","#1c9099"));p2

cowplot::plot_grid(p1,p2,ncol = 2,align = "h")



# Extended.Fig.1b-c: SV/TR number compare with GTEx  --------------------------------------------------------
var <- data.frame(type = c("SV", "STR", "VNTR"),Count = c(95957,1056776,38615))
var$label <- paste0(var$Count,"\n", " (", scales::percent(var$Count/sum(var$Count)), ")")
var$dataset <- "GTOP"
var2 <- data.frame(type = c("SV", "STR", "VNTR"),Count = c( 23602, 175226, 10264))
var2$label <- paste0(var2$Count,"\n", " (", scales::percent(var2$Count/sum(var2$Count)), ")")
var2$dataset <- "GTEx"
dat <- rbind(var,var2)
colors <- c("GTOP" = "#a5382a","GTEx" = "#517f84")
dat$dataset <- factor(dat$dataset,levels = c("GTOP","GTEx"))
dat$type <- factor(dat$type,levels = c("SNP","Indel","STR","VNTR","SV"))
p1 <- ggplot(dat[dat$type=="SV",], aes(x = type, y = Count, fill = dataset)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  scale_fill_manual(values = colors) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3) +
  labs(x = "Variant Type",y = "# of Variants") +
  #scale_y_log10(labels = scales::comma) +  
  theme_classic() +
  theme(legend.position = "none",axis.text = element_text(color = "black", size = 10), 
        axis.ticks = element_line(color = "black"))

p2 <- ggplot(dat[dat$type=="STR" | dat$type=="VNTR",], 
             aes(x = type, y = Count, fill = dataset)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  scale_fill_manual(values = colors) +
  facet_wrap(~type ,scales = "free_y", nrow = 1) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3) +
  labs(x = "Variant Type",y = "# of Variants") +
  #scale_y_log10(labels = scales::comma) +  
  theme_classic() +
  theme(legend.position = "none",axis.text = element_text(color = "black", size = 10), 
        axis.ticks = element_line(color = "black"))

p1
p2

# Extended.Fig.1d: LRS specific allele diversity ---------------------------------------------------

ac_df <- fread("ExtendFig1d.txt")
p4 <- ggplot(ac_df, aes(x = factor(threshold), y = OR)) +
  geom_point(size = 3, color = "#993828") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "#993828") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Allele count threshold",
    y = "Odds Ratio (LRS-specific vs Shared)") +
  theme_classic(base_size = 14)

p4

# Extended.Fig.1e: variant number and length ----------------------------------------
info <- fread("ExtendFig1e_PathogenicTR.txt")
info <- info %>% mutate(gene=word(word(V4,start=1,end=1,sep=":"),start=-1L,end=-1L,sep="\\="),
                        TRID=word(word(V4,start=1,end=1,sep=";"),start=-1L,end=-1L,sep="\\:"))

datt <- fread("ExtendFig1e_PathogenicTR_CNV.txt") %>%mutate(class = case_when(
  is.na(Sequencing_SRS) ~ "LRS specific",(CNV_LRS - CNV_SRS) > 2 ~ "LRS higher",
  TRUE ~ "Shared")) %>%
  select(-gene) %>%left_join(info %>% select(gene, TRID), by = "TRID") %>%
  mutate(motif = word(TRID, start = -1L, end = -1L, sep = "_"),
         RU = paste0(nchar(motif), "bp", ifelse(nchar(motif) > 6, "VNTR", "STR")))

valid_genes <- datt %>%filter(!is.na(CNV_LRS))

gene_order <- valid_genes %>%
  group_by(gene) %>%summarise(mean_CNV = mean(CNV_LRS, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_CNV)) %>%pull(gene)

valid_genes <- valid_genes %>%mutate(gene = factor(gene, levels = gene_order))

shared_bg <- valid_genes %>%
  filter(class == "Shared") %>%
  group_by(gene) %>%
  summarise(xmin = 0, xmax = max(CNV_LRS), .groups = "drop") %>%
  mutate(gene = factor(gene, levels = gene_order))

max_CNV <- max(valid_genes$CNV_LRS, na.rm = TRUE)

lrs_bg <- valid_genes %>%
  filter(class != "Shared") %>%
  group_by(gene) %>%
  summarise(xmin = min(CNV_LRS), xmax = max_CNV, .groups = "drop") %>%
  mutate(gene = factor(gene, levels = gene_order))

map_CNV <- function(x) ifelse(x <= 50, x, 50 + (x - 50) / 10)

valid_genes <- valid_genes %>%mutate(CNV_plot = map_CNV(CNV_LRS))

shared_bg <- shared_bg %>%mutate(xmin_plot = map_CNV(xmin), xmax_plot = map_CNV(xmax))

lrs_bg <- lrs_bg %>%mutate(xmin_plot = map_CNV(xmin), xmax_plot = map_CNV(xmax))

breaks1 <- 1:50
labels1 <- as.character(breaks1)

max_val <- ceiling(max(valid_genes$CNV_LRS, na.rm = TRUE) / 50) * 50
breaks2_orig <- seq(100, max_val, by = 50)
breaks2_mapped <- 50 + (breaks2_orig - 50) / 5
labels2 <- as.character(breaks2_orig)

breaks_all <- c(breaks1, breaks2_mapped)
labels_all <- c(labels1, labels2)

p1 <- ggplot() +
  geom_rect(data = lrs_bg,
            aes( xmin = xmin_plot,xmax = xmax_plot,
                 ymin = as.numeric(gene) - 0.4,ymax = as.numeric(gene) + 0.4),
            fill = "grey90",alpha = 0.7) +
  geom_point(data = valid_genes,aes(x = CNV_plot, y = gene, color = class),size = 1) +
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




