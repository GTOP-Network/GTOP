library(ggplot2)
library(scales)
library(ggrastr)
library(data.table)
library(stringi)
library(ggplot2)
library(patchwork)
library(viridis)
library(tidyr)
library(dplyr)

setwd("/path/to/dir/")
tissue_info <- fread("/path/to/input/GTOP_tissue_info.csv")

# load data and filter low expression by Tissue -------------------------------------------------

read_mat <- function(file){
  dt <- fread(file, data.table = FALSE)
  rownames(dt) <- dt[,1]
  mat <- as.matrix(dt[,-1,drop=FALSE])
  mode(mat) <- "numeric"
  return(mat)
}

### ========== 1. input file ==========
salmon_tpm_file <- "./input/enhanced.transcript.tpm.salmon.tsv"
rsem_tpm_file <- "./input/enhanced.transcript.tpm.rsem.tsv"

salmon_count_file <- "./input/enhanced.transcript.count.salmon.tsv"
rsem_count_file <- "./input/enhanced.transcript.count.rsem.tsv"

### ========== 2. load data ==========
salmon_tpm <- read_mat(salmon_tpm_file)
rsem_tpm <- read_mat(rsem_tpm_file)

salmon_count <- read_mat(salmon_count_file)
rsem_count <- read_mat(rsem_count_file)

### ========== 3. transcript & sample intersect ==========
common_tx <- Reduce(intersect, list(rownames(salmon_tpm), rownames(rsem_tpm)))

salmon_tpm <- salmon_tpm[common_tx, ]
rsem_tpm <- rsem_tpm[common_tx, ]

salmon_count <- salmon_count[common_tx, ]
rsem_count <- rsem_count[common_tx, ]

common_samples <- Reduce(intersect, list(colnames(salmon_tpm),colnames(rsem_tpm),
                                         colnames(salmon_count), colnames(rsem_count)))

salmon_tpm <- salmon_tpm[, common_samples]
rsem_tpm <- rsem_tpm[, common_samples]

salmon_count <- salmon_count[, common_samples]
rsem_count <- rsem_count[, common_samples]

### ========== 4. mean TPM with Salmon and RSEM ==========
salmon_rsem_tpm <- (salmon_tpm + rsem_tpm) / 2
salmon_rsem_tpm_df <- data.frame(
  transcript_id = rownames(salmon_rsem_tpm),
  salmon_rsem_tpm,
  check.names = FALSE)
dir.create("./output",showWarnings = F)
fwrite(salmon_rsem_tpm_df,
       "./input/enhanced.transcript.tpm.salmon_rsem_average.tsv",sep = "\t")

### ========== 5. define tissue(sample ID split) ==========
tissue <- sapply(strsplit(colnames(salmon_tpm), "-"), `[`, 3)

### ========== 6. filter function  ==========
filter_tx <- function(tpm, count, tpm_cutoff = 0.1, count_cutoff = 6, frac = 0.2){
  min_n <- ncol(tpm) * frac
  keep <- rowSums(tpm > tpm_cutoff & count >= count_cutoff) >= min_n
  return(keep)
}

### ========== 7. correlation function ==========
calc_cor <- function(m1, m2){
  if(nrow(m1) == 0) return(numeric(0))
  out <- numeric(nrow(m1))
  for(i in seq_len(nrow(m1))){
    out[i] <- suppressWarnings(
      cor(as.numeric(m1[i,]),as.numeric(m2[i,]), method="spearman", use="pairwise.complete.obs"
      ) )
  }
  
  return(out)
}

# Correlation by Tissue -------------------------------------------------

res_list <- list()

for(t in unique(tissue)){
  cat("\nProcessing tissue:", t, "\n")
  idx <- which(as.character(tissue) == as.character(t))
  s1 <- salmon_tpm[,idx,drop=FALSE]
  s2 <- rsem_tpm[,idx,drop=FALSE]
  
  c1 <- salmon_count[,idx,drop=FALSE]
  c2 <- rsem_count[,idx,drop=FALSE]
  
  keep1 <- filter_tx(s1, c1)
  keep2 <- filter_tx(s2, c2)
  
  keep <- keep1 & keep2
  num_keep <- sum(keep)  
  print(paste0("Salmon and RSEM all suported transcript number：", num_keep))
  if(sum(keep) == 0){cat("Skip tissue:", t, "- no transcripts passed filter\n");next}
  tpms1 <- s1[keep,,drop=FALSE]
  if(ncol(tpms1) <= 2){cat("Skip tissue:", t, "- not enough samples\n");next}
  tpms2 <- s2[keep,,drop=FALSE]
  
  log1 <- log2(tpms1 + 1)
  log2 <- log2(tpms2 + 1)
  
  cor_df <- data.frame(
    Transcript = rownames(log1),
    Tissue = t,
    MeanExpr = rowMeans(tpms1),
    Salmon_RSEM = calc_cor(log1, log2) )
  
  res_list[[t]] <- cor_df
}

final_cor_df <- rbindlist(res_list)

# MARD by Tissue -------------------------------------------------

calc_mard <- function(m1, m2, eps = 0){
  m1[m1 < eps] <- 0
  m2[m2 < eps] <- 0
  mean1 <- rowMeans(m1)
  mean2 <- rowMeans(m2)
  denom <- mean1 + mean2
  mard <- ifelse(denom == 0, 0, abs(mean1 - mean2) / denom)
  return(mard)
}

res_list <- list()

for(t in unique(tissue)){
  cat("Processing tissue:", t, "\n")
  idx <- which(as.character(tissue) == as.character(t))
  s1 <- salmon_tpm[, idx, drop = FALSE]
  s2 <- rsem_tpm[, idx, drop = FALSE]
  
  c1 <- salmon_count[, idx, drop = FALSE]
  c2 <- rsem_count[, idx, drop = FALSE]
  keep1 <- filter_tx(s1, c1)
  keep2 <- filter_tx(s2, c2)
  keep <- keep1 & keep2
  
  if(sum(keep) == 0){
    cat("Skip tissue (no transcripts passed):", t, "\n")
    next
  }
  
  s1 <- s1[keep,,drop = FALSE]
  s2 <- s2[keep,,drop = FALSE]
  
  mard_12 <- calc_mard(s1, s2)
  
  df <- data.frame(Transcript = rownames(s1),Tissue = t,
                   MeanExpr = rowMeans(s1), MARD_Salmon_RSEM = mard_12)
  res_list[[t]] <- df
}

final_mard <- rbindlist(res_list)

# MAE by Tissue ------------------------------------------------------------

calc_pairwise_abs_diff <- function(log1, log2) {
  mean_log1 <- rowMeans(log1, na.rm = TRUE)
  mean_log2 <- rowMeans(log2, na.rm = TRUE)
  diff_SR <- abs(mean_log1 - mean_log2)
  return(list(SR = diff_SR))
}

res_list_pairwise <- list()
for(t in unique(tissue)){
  cat("\nProcessing tissue:", t, "\n")
  idx <- which(as.character(tissue) == as.character(t))
  s1 <- salmon_tpm[,idx,drop=FALSE]
  s2 <- rsem_tpm[,idx,drop=FALSE]
  c1 <- salmon_count[,idx,drop=FALSE]
  c2 <- rsem_count[,idx,drop=FALSE]
  
  keep1 <- filter_tx(s1, c1)
  keep2 <- filter_tx(s2, c2)
  keep <- keep1 & keep2
  
  if(sum(keep) == 0){
    cat("Skip tissue:", t, "- no transcripts passed filter\n")
    next
  }
  
  tpms1 <- s1[keep,,drop=FALSE]
  tpms2 <- s2[keep,,drop=FALSE]
  log1 <- log2(tpms1 + 1)
  log2 <- log2(tpms2 + 1)
  
  diffs <- calc_pairwise_abs_diff(log1, log2)
  
  mean_tpm <- rowMeans(tpms1)  
  pair_df <- data.frame( Transcript = rownames(log1),Tissue = t,
                         MeanExpr = mean_tpm,diff_Salmon_RSEM = diffs$SR )
  
  res_list_pairwise[[t]] <- pair_df
}

pairwise_all_mae <- do.call(rbind, res_list_pairwise)


# combine three metrics -------------------------------------------------

merged <- final_mard %>% 
  dplyr::left_join(final_cor_df, by = c("Transcript", "Tissue", "MeanExpr")) %>%
  dplyr::left_join(pairwise_all_mae, by = c("Transcript", "Tissue", "MeanExpr"))

setnames(merged, old = c("diff_Salmon_RSEM", "Salmon_RSEM"), new = c("MAE_Salmon_RSEM", "Cor_Salmon_RSEM"))

merged <- merged %>% mutate(
  Tier1 = "Yes",
  Tier2 = ifelse(Cor_Salmon_RSEM > 0.5 | is.na(Cor_Salmon_RSEM), "Yes", "No"),
  Tier3 = ifelse(Tier2 == "Yes" & (MeanExpr < 5 | MARD_Salmon_RSEM < 0.33 | MAE_Salmon_RSEM < 0.33), "Yes", "No"))%>% 
  mutate(Tissue = as.integer(Tissue)) %>%
  dplyr::left_join(tissue_info %>% dplyr::select(Tissue, Tissue_Code) %>% setnames("TissueName", "Tissue"), by = "Tissue")

# Filtering transcripts with highly inconsistent -----------------------------------------

ave_tpm <- as.data.table(salmon_rsem_tpm_df)
setnames(ave_tpm, 1,"Transcript")

sample_info <- data.table( Sample = colnames(ave_tpm)[-1])
sample_info[, Tissue := sapply( strsplit(Sample, "-"),function(x) x[3])]
sample_info$Tissue <- as.integer(sample_info$Tissue)

outdir <- "./output/Salmon_RSEM_Tier3_average_TPM"
dir.create(outdir,showWarnings = FALSE)

for(t in unique(sample_info$Tissue)){
  
  cat("Processing tissue:",t,"\n")
  samples <- sample_info[Tissue==t, Sample]
  keep_tx <- merged[ Tissue==t &Tier3=="Yes", unique(Transcript)]
  
  tpm <- ave_tpm[Transcript %in% keep_tx, c("Transcript",samples),with=FALSE ]
  tissue_name <- unique(merged[Tissue == t, TissueName] )
  
  outfile <- file.path( outdir, paste0(tissue_name, ".filter.average_Salmon_RSEM_TPM.bed" ) )
  fwrite(tpm, outfile, sep="\t" )
  
}






