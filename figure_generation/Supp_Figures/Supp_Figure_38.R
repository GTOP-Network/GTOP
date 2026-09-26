#==============================================#
# Joint_Finemap #
# Supp-Figure-38#
#==============================================#
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(data.table)
library(dplyr)
library(magrittr)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig38.Joint_Finemap/input")

# Supp.Fig.38a ------------------------------------------------------------

df_cs.filter <- readRDS("./SNV_SVTR_mixed_CS.RDS")
Tissue_List <- unique(df_cs.filter$Tissue)
fit_df <- data.frame()
for(i in 1:length(Tissue_List)){
  df_tissue <- df_cs.filter %>% filter(Tissue == Tissue_List[i])
  fit <- glm(
    svtr_lead ~ offset(offset),
    data = df_tissue,
    family = binomial()
  )
  
  beta <- coef(fit)["(Intercept)"]
  
  se <- summary(fit)$coefficients[
    "(Intercept)", 
    "Std. Error"
  ]
  
  OR <- exp(beta)
  
  CI_low <- exp(beta - 1.96 * se)
  CI_high <- exp(beta + 1.96 * se)
  
  fit_df <- rbind(fit_df,list(oddsratio=OR,CI95_low=CI_low,CI95_high=CI_high))
}

fit_df$Tissue <- Tissue_List

fit_df <- fit_df[order(fit_df$oddsratio),]

fit_df$Tissue <- factor(fit_df$Tissue,levels = fit_df$Tissue)

df_cs.stratified <- readRDS("./SNV_SVTR_mixed_CS.R2_stratified.RDS")
df_cs.ld_high <- df_cs.stratified %>% filter(R2>=0.8)
df_cs.ld_low <- df_cs.stratified %>% filter(R2<0.8)

# fit model and plot
fit_df.high <- data.frame()
for(i in 1:length(Tissue_List)){
  df_tissue <- df_cs.ld_high %>% filter(Tissue == Tissue_List[i])
  fit <- glm(
    svtr_lead ~ offset(offset),
    data = df_tissue,
    family = binomial()
  )
  
  beta <- coef(fit)["(Intercept)"]
  
  se <- summary(fit)$coefficients[
    "(Intercept)", 
    "Std. Error"
  ]
  
  OR <- exp(beta)
  
  CI_low <- exp(beta - 1.96 * se)
  CI_high <- exp(beta + 1.96 * se)
  
  fit_df.high <- rbind(fit_df.high,list(oddsratio=OR,CI95_low=CI_low,CI95_high=CI_high))
}

fit_df.high$Tissue <- Tissue_List

fit_df.low <- data.frame()
for(i in 1:length(Tissue_List)){
  df_tissue <- df_cs.ld_low %>% filter(Tissue == Tissue_List[i])
  fit <- glm(
    svtr_lead ~ offset(offset),
    data = df_tissue,
    family = binomial()
  )
  
  beta <- coef(fit)["(Intercept)"]
  
  se <- summary(fit)$coefficients[
    "(Intercept)", 
    "Std. Error"
  ]
  
  OR <- exp(beta)
  
  CI_low <- exp(beta - 1.96 * se)
  CI_high <- exp(beta + 1.96 * se)
  
  fit_df.low <- rbind(fit_df.low,list(oddsratio=OR,CI95_low=CI_low,CI95_high=CI_high))
}

fit_df.low$Tissue <- Tissue_List

fit_df.high$group <- "High_LD"
fit_df.low$group <- "Low_LD"

fit_df$group <- "unstratificated"

fit_all <- rbind(fit_df,fit_df.high,fit_df.low)
fit_all$group <- factor(fit_all$group,levels = c("unstratificated","High_LD","Low_LD"))

# plot
library(ggplot2)
p <- ggplot(fit_all,aes(x=Tissue,y=oddsratio)) + geom_pointrange(aes(ymin = CI95_low,ymax = CI95_high)) + theme_classic() + 
  geom_hline(yintercept = 0,linetype="dashed") + coord_flip() + facet_wrap(~group);p

# Supp.Fig.38b ------------------------------------------------------------

df_m <- readRDS("./Simulation.enrichment.RDS")

df_real <- readRDS("./observed_enrichment.RDS")

df_real$PIP <- factor(df_real$PIP,levels = c(0,0.2,0.4,0.6,0.8))

library(ggplot2)
df_m$PIP <- factor(df_m$PIP,levels = c(0,0.2,0.4,0.6,0.8))
p1 <- ggplot(df_m,aes(x=PIP,y=OR)) + geom_boxplot(width = 0.35, alpha = 1, outlier.shape = NA,color="#D6D6D6") + 
  geom_jitter(aes(x = PIP, y = OR),color = "#555555",size = 1,width = 0.1,alpha = 1
  ) + geom_point(data = df_real,aes(x = PIP, y = OR),color = "#C76B5A",size = 5,shape = 18) + theme_classic();p1


# Supp.Fig.38c ------------------------------------------------------------

extract_tissue <- function(x){
  return(strsplit(x,split = ":",fixed = T)[[1]][1])
}

extract_run <- function(x){
  return(strsplit(x,split = ":",fixed = T)[[1]][2])
}

df_pip0 <- readRDS("downsample_SVTR_enrich_causal_CS.compare_to_SNV.0.RDS")
df_pip2 <- readRDS("downsample_SVTR_enrich_causal_CS.compare_to_SNV.0.2.RDS")
df_pip4 <- readRDS("downsample_SVTR_enrich_causal_CS.compare_to_SNV.0.4.RDS")
df_pip6 <- readRDS("downsample_SVTR_enrich_causal_CS.compare_to_SNV.0.6.RDS")
df_pip8 <- readRDS("downsample_SVTR_enrich_causal_CS.compare_to_SNV.RDS")

df_pip0$TissueName <- sapply(as.character(df_pip0$Tissue),extract_tissue)
df_pip0$RUN <- sapply(as.character(df_pip0$Tissue),extract_run)

df_pip2$TissueName <- sapply(as.character(df_pip2$Tissue),extract_tissue)
df_pip2$RUN <- sapply(as.character(df_pip2$Tissue),extract_run)

df_pip4$TissueName <- sapply(as.character(df_pip4$Tissue),extract_tissue)
df_pip4$RUN <- sapply(as.character(df_pip4$Tissue),extract_run)

df_pip6$TissueName <- sapply(as.character(df_pip6$Tissue),extract_tissue)
df_pip6$RUN <- sapply(as.character(df_pip6$Tissue),extract_run)

df_pip8$TissueName <- sapply(as.character(df_pip8$Tissue),extract_tissue)
df_pip8$RUN <- sapply(as.character(df_pip8$Tissue),extract_run)

df_pip0$Tissue <- NULL
df_pip2$Tissue <- NULL
df_pip4$Tissue <- NULL
df_pip6$Tissue <- NULL
df_pip8$Tissue <- NULL


df_pip0  <- df_pip0  %>% mutate(PIP_cutoff = 0)
df_pip02 <- df_pip2 %>% mutate(PIP_cutoff = 0.2)
df_pip04 <- df_pip4 %>% mutate(PIP_cutoff = 0.4)
df_pip06 <- df_pip6 %>% mutate(PIP_cutoff = 0.6)
df_pip08 <- df_pip8 %>% mutate(PIP_cutoff = 0.8)

df_all <- bind_rows(
  df_pip0,
  df_pip02,
  df_pip04,
  df_pip06,
  df_pip08
)

df_heatmap <- df_all %>%
  mutate(
    log2OR = log2(OR),
    significant = Pval < 0.05
  ) %>%
  group_by(TissueName, PIP_cutoff) %>%
  summarise(
    median_log2OR = median(log2OR, na.rm = TRUE),
    median_OR = median(OR, na.rm = TRUE),
    n_sig = sum(significant, na.rm = TRUE),
    n_run = n(),
    .groups = "drop"
  )


df_heatmap <- df_heatmap %>%
  mutate(
    PIP_cutoff = factor(
      PIP_cutoff,
      levels = c(0, 0.2, 0.4, 0.6, 0.8),
      labels = c(
        "PIP ≥ 0",
        "PIP ≥ 0.2",
        "PIP ≥ 0.4",
        "PIP ≥ 0.6",
        "PIP ≥ 0.8"
      )
    )
  )


tissue_order <- df_heatmap %>%
  filter(PIP_cutoff == "PIP ≥ 0") %>%
  arrange(median_OR) %>%
  pull(TissueName)

df_heatmap <- df_heatmap %>%
  mutate(
    TissueName = factor(
      TissueName,
      levels = tissue_order
    )
  )


df_heatmap <- df_heatmap %>%
  mutate(
    label = paste0(
      sprintf("%.2f", median_OR),
      "\n",
      n_sig, "/", n_run
    )
  )

# plot
p <- ggplot(
  df_heatmap,
  aes(
    x = PIP_cutoff,
    y = TissueName,
    fill = median_log2OR
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.7
  ) +
  
  geom_text(
    aes(label = label),
    size = 3.2,
    lineheight = 0.9
  ) +
  
  scale_fill_gradient2(
    name = "Median\nlog2(OR)",
    midpoint = 0
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    
    axis.text.y = element_text(
      color = "black"
    ),
    
    legend.title = element_text(
      size = 10
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    axis.ticks = element_blank(),
    
    plot.margin = margin(
      5, 5, 5, 5
    )
  )

p



# Supp.Fig.38d --SV-SNV----------------------------------------------------

df_condition_res <- readRDS("df_condition_res.sv.RDS")
df_condition_res %<>% filter(!is.na(Beta_0)) %>% filter(!is.na(Beta_1))
df_condition_res$Pval0_log <- -log10(df_condition_res$Pval_0)
df_condition_res$Pval1_log <- -log10(df_condition_res$Pval_1)
df_condition_res %<>% mutate(R2=SV_SNV^2)

df_condition_res.rev <- readRDS("df_condition_res.sv_rev.RDS")
df_condition_res.rev %<>% filter(!is.na(Beta_0)) %>% filter(!is.na(Beta_1))
df_condition_res.rev$Pval0_log <- -log10(df_condition_res.rev$Pval_0)
df_condition_res.rev$Pval1_log <- -log10(df_condition_res.rev$Pval_1)
df_condition_res.rev %<>% mutate(R2=SV_SNV^2)

df_plot <- data.frame(Count=c(
  dim(df_condition_res %>% filter(Pval_1<0.01))[1],
  dim(df_condition_res)[1]-dim(df_condition_res %>% filter(Pval_1<0.01))[1],
  dim(df_condition_res.rev %>% filter(Pval_1<0.01))[1],
  dim(df_condition_res.rev)[1]-dim(df_condition_res.rev %>% filter(Pval_1<0.01))[1]),
  Group=c("SV-SNV","SV-SNV","SNV-SV","SNV-SV"),
  Effect=c("Sig","Not_Sig","Sig","Not_Sig"))
df_plot$Group <- factor(df_plot$Group,levels = c("SV-SNV","SNV-SV"))
df_plot$Effect <- factor(df_plot$Effect,levels = c("Sig","Not_Sig"))

p3 <- ggplot(df_plot,aes(x=Group,y=Count,fill=Effect)) + geom_bar(stat="identity",width=.8) + theme_classic() + 
  scale_fill_manual(breaks = c("Not_Sig","Sig"),values = c("#878787","#bdb4ad"));p3


# Supp.Fig.38e TR-SNV ---------------------------------------------------

df_condition_res <- readRDS("df_condition_res.tr.RDS")
df_condition_res %<>% filter(!is.na(Beta_0)) %>% filter(!is.na(Beta_1))
df_condition_res$Pval0_log <- -log10(df_condition_res$Pval_0)
df_condition_res$Pval1_log <- -log10(df_condition_res$Pval_1)
df_condition_res %<>% mutate(R2=TR_SNV^2)
df_condition_res<-df_condition_res[!is.na(df_condition_res$TR_SNV),]

df_condition_res.rev <- readRDS("df_condition_res.tr.Rev.RDS")
df_condition_res.rev %<>% filter(!is.na(Beta_0)) %>% filter(!is.na(Beta_1))
df_condition_res.rev$Pval0_log <- -log10(df_condition_res.rev$Pval_0)
df_condition_res.rev$Pval1_log <- -log10(df_condition_res.rev$Pval_1)
df_condition_res.rev %<>% mutate(R2=TR_SNV^2)

df_plot <- data.frame(Count=c(
  dim(df_condition_res %>% filter(Pval_1<0.01))[1],
  dim(df_condition_res)[1]-dim(df_condition_res %>% filter(Pval_1<0.01))[1],
  dim(df_condition_res.rev %>% filter(Pval_1<0.01))[1],
  dim(df_condition_res.rev)[1]-dim(df_condition_res.rev %>% filter(Pval_1<0.01))[1]),
  Group=c("TR-SNV","TR-SNV","SNV-TR","SNV-TR"),
  Effect=c("Sig","Not_Sig","Sig","Not_Sig"))
df_plot$Group <- factor(df_plot$Group,levels = c("TR-SNV","SNV-TR"))
df_plot$Effect <- factor(df_plot$Effect,levels = c("Sig","Not_Sig"))

p3 <- ggplot(df_plot,aes(x=Group,y=Count,fill=Effect)) + geom_bar(stat="identity",width=.8) + theme_classic() + 
  scale_fill_manual(breaks = c("Not_Sig","Sig"),values = c("#878787","#bdb4ad"));p3

