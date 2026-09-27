#!/bin/Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(R.utils)
library(mashr)
library(doParallel)
library(foreach)
library(tidyverse)

"%&%" = function(a, b) { paste0(a, b) }
readFile = function(f){
    if(grepl(".RDS", f, ignore.case=T)){
        DF = readRDS(f)
    } else{
        DF = fread(f)
    }
    return(DF)
}

strong_rmna <- function(tmp_df){
    gene_group <- stringr::word(rownames(tmp_df), start = 2, end = 2, sep = ",")
    gene_nanum <- apply(tmp_df, 1, function(x){sum(!is.na(x))})
    gene_df <- data.frame("gene"=gene_group, "nonNAnum"=gene_nanum, "QTL"=rownames(tmp_df))
    gene_filter <- gene_df %>% group_by(gene) %>% filter(nonNAnum==max(nonNAnum))
    tmp_output <- tmp_df[gene_filter$QTL, ]
    tmp_output
}

strong_zscore <- function(tmp_df){
    na_count <- rowSums(is.na(tmp_df))
    tmp_df <- tmp_df[na_count!=ncol(tmp_df),]
    gene_group <- stringr::word(rownames(tmp_df), start = 2, end = 2, sep = ",")
    source_group <- stringr::word(colnames(tmp_df), start = 1, end = 1, sep = "_")
    # length(unique(gene_group))
    gene_zcore <- apply(tmp_df, 1, function(x){max(abs(na.omit(x)))})
    source_support <- apply(tmp_df, 1, function(x){length(unique(source_group[!is.na(x)]))})
    gene_df <- data.frame("gene"=gene_group, "zcore"=gene_zcore, "support"=source_support, "QTL"=rownames(tmp_df))
    gene_df$score <- gene_df$zcore*gene_df$support
    gene_filter <- gene_df %>% group_by(gene) %>% filter(score==max(score))
    gene_filter <- gene_filter %>% group_by(gene) %>% filter(QTL==QTL[1])
    tmp_output <- tmp_df[gene_filter$QTL, ]
    tmp_output
}

strong_unique <- function(tmp_df){
    QTL_group <- stringr::word(rownames(tmp_df), start = 2, end = 2, sep = ",")
    QTL_nanum <- apply(tmp_df, 1, function(x){median(abs(x), na.rm = T)})
    QTL_df <- data.frame("gene"=QTL_group, "nonNAnum"=QTL_nanum, "QTL"=rownames(tmp_df))
    QTL_filter <- QTL_df %>% group_by(gene) %>% filter(nonNAnum==max(nonNAnum))
    table(table(QTL_filter$gene))
    QTL_filter <- QTL_filter %>% group_by(gene) %>% filter(QTL==QTL[1])
    table(table(QTL_filter$gene))
    tmp_output <- tmp_df[QTL_filter$QTL, ]
    tmp_output
}

#----------------------------------------------------------------------------
ARGS <- commandArgs(trailingOnly = TRUE)
file_strong = ARGS[1] 
file_random = ARGS[2]
dropna = as.logical(as.numeric(ARGS[3]))
dir_output = ARGS[4]

# file_strong = "2023-11-04-xQTL_characteristics/2024-03-20-mashr/output/CReQTL_subtype/Strong/strong_pairs.MashR_input.txt.gz"
# file_random = "2023-11-04-xQTL_characteristics/2024-03-20-mashr/output/CReQTL_subtype/Random/MashR.random_subset_1000000.RDS"
# dropna = as.logical(as.numeric("0"))
# dir_output = "2023-11-04-xQTL_characteristics/2024-03-20-mashr/output/CReQTL_subtype/top_pairs/"
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

if (length(dir_output) == 0){
    dir_output = "./"
} else {
    dir_output = dir_output %&% "/"
}
dir.create(dir_output, recursive = T)

set.seed(9823)

############ strong pairs
raw.zval_data.strong = as.data.frame(readFile(file_strong))
rownames(raw.zval_data.strong) = raw.zval_data.strong$pair_id
raw.zval_data.strong$pair_id = NULL

zval_data.strong <- strong_zscore(tmp_df = raw.zval_data.strong)
na_count <- rowSums(!is.na(zval_data.strong))
ncol(zval_data.strong)
hist(na_count)
table(na_count)

if(grepl("subtype", file_strong)){
    zval_data.strong <- zval_data.strong[na_count>=2, ]
}

# zval_data.strong1 <- strong_rmna(tmp_df = raw.zval_data.strong)
# ggvenn::ggvenn(list("zscore"=rownames(zval_data.strong), "rmna"=rownames(zval_data.strong1)))
# setdiff(rownames(zval_data.strong), rownames(zval_data.strong1))[11:22]
# rownames(zval_data.strong1)[grepl("CCNL2:PDIUI", rownames(zval_data.strong1))]
# rownames(zval_data.strong)[grepl("CCNL2:PDIUI", rownames(zval_data.strong))]
# barplot(as.numeric(raw.zval_data.strong["1:821532:A:G,CCNL2:PDIUI", ]))
# barplot(as.numeric(raw.zval_data.strong["1:901516:T:C,CCNL2:PDIUI", ]))

if(dropna){
    zval_data.strong = zval_data.strong[rowSums(is.na(zval_data.strong))==0,]
} else {
    zval_data.strong[is.na(zval_data.strong)] = 0
}
zval_data.strong = as.matrix(zval_data.strong)

############ all pairs of a random subset of 1000000 tests
zval_data.random.subset = as.data.frame(readFile(file_random))
rownames(zval_data.random.subset) = zval_data.random.subset$pair_id
zval_data.random.subset$pair_id = NULL
zval_data.random.subset = as.matrix(zval_data.random.subset)

zval_data.random.subset <- zval_data.random.subset[, colnames(zval_data.strong)]
if (sum(colnames(zval_data.random.subset) != colnames(zval_data.strong))>0){
    stop("Different colnames between random subset and strong subset.")
}

## Correlation structure
data.temp = mash_set_data(Bhat=zval_data.random.subset,alpha=1,zero_Bhat_Shat_reset=1e6)
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

data.random = mash_set_data(Bhat=zval_data.random.subset,alpha=1,V=Vhat,zero_Bhat_Shat_reset=1e6)
data.strong = mash_set_data(Bhat=zval_data.strong,alpha=1,V=Vhat,zero_Bhat_Shat_reset=1e6)

# U.flash <- cov_flash(data.strong, factors="nonneg", tag="non_neg")

## Data driven covariances
# if (file.exists(dir_output %&% "U.ed.RDS")){
#     U.ed = readRDS(dir_output %&% "U.ed.RDS")
# } else {
    U.pca = cov_pca(data.strong, 5)
    U.ed = cov_ed(data.strong, U.pca, algorithm = "teem")
    ##
    attr_names <- attr(U.ed,"names")
    attr_names <- attr_names[!is.na(attr_names)]
    U.ed_new <- list()
    for(i in 1:length(attr_names)){
        U.ed_new[[attr_names[i]]] <- U.ed[,,i]
    }
    saveRDS(U.ed_new, dir_output %&% "U.ed.RDS")
    U.ed = readRDS(dir_output %&% "U.ed.RDS")
# }

## Fit mash model (estimate mixture proportions)
U.c = cov_canonical(data.random)

## create study covariance matrices
# name_source <- colnames(data.random$Bhat)
# name_source <- sapply(strsplit(name_source, "_"), function(x){x[1]})
# U.study <- list()
# for(i in unique(name_source)){
#     U.tmp = matrix(0, nrow = length(name_source), ncol = length(name_source))
#     tmp_index <- name_source%in%i
#     U.tmp[tmp_index, tmp_index] <- 1
#     U.study[[i]] <- U.tmp
# }

if(grepl("edQTL_subtype", file_strong)){
    for(i in names(U.ed)){
        U.ed[[i]] <- round(U.ed[[i]], digits = 13)
    }
}

m.r = mash(data.random, Ulist = c(U.ed, U.c), outputlevel = 1)
# m.r <- readRDS(dir_output %&% "m.r.RDS")
saveRDS(m.r, dir_output %&% "m.r.RDS")

## Compute posterior summaries
m.s = mash(data.strong, g=get_fitted_g(m.r), fixg=TRUE)

saveRDS(m.s, dir_output %&% "m.s_zval.RDS")
saveRDS(get_lfsr(m.s), dir_output %&% "lfsr_m.s.RDS")
# saveRDS(get_pm(m.s), dir_output %&% "pm_m.s.RDS")
# saveRDS(get_psd(m.s), dir_output %&% "psd_m.s.RDS")
# saveRDS(get_significant_results(m.s), dir_output %&% "significant_results_m.s.RDS")
# saveRDS(get_pairwise_sharing(m.s), dir_output %&% "pairwise_sharing_m.s.RDS")
# saveRDS(get_pairwise_sharing(m.s, factor=0), dir_output %&% "pairwise_sharing_factor0_m.s.RDS")
# saveRDS(get_loglik(m.s), dir_output %&% "loglik_m.s.RDS")
# saveRDS(get_estimated_pi(m.s), dir_output %&% "estimated_pi_m.s.RDS")
