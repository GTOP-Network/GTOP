#!/usr/bin/env Rscript

# Load required libraries
library(data.table)
library(dplyr)

# ------------------------------------------------------------
# Function: validate and parse command-line arguments
# Expected input order:
#   1. INPUTDIR (e.g., "slurm/coloc/task1/output")
#   2. OUTPUTDIR (e.g., "output")
# ------------------------------------------------------------
argvs <- commandArgs(trailingOnly = TRUE)

INPUTDIR <- argvs[1]
OUTPUTDIR <- argvs[2]
if(!file.exists(OUTPUTDIR)){dir.create(OUTPUTDIR, recursive = TRUE)}

coloc_files <- list.files(INPUTDIR)
coloc_res <- rbindlist(lapply(coloc_files, function(x){
  f_info <-  fread(file.path(INPUTDIR, x))
  f_info
}))
colnames(coloc_res) <- tolower(colnames(coloc_res))
colnames(coloc_res)[colnames(coloc_res)=="loci"] <- "locus_id"

GWAS_info <- fread("input/collected_GWAS_names.txt")
colnames(GWAS_info) <- tolower(colnames(GWAS_info))

coloc_res <- merge(GWAS_info[, .(gwas_name, phecode_abbr, disease_trait, trait_class=category)], coloc_res, by="gwas_name")

for(i in unique(coloc_res$xqtl_type)){
  sub_data <- coloc_res[xqtl_type==i]
  sub_data <- sub_data %>% arrange(gwas_name)
  fwrite(sub_data, file = paste0(OUTPUTDIR, "/", i, "_coloc_raw.txt"), sep = "\t")
}

# length(unique(paste0(coloc_res$gwas_name, "_", coloc_res$gwas_sentinel))) # 4240
# length(unique(paste0(coloc_res$phecode_abbr, "_", coloc_res$loci))) # 3782

# a <- fread("/lustre/home/tzhang/2026-05-07-gtop_xqtl-Project/2025-09-28-coloc/input/run_data/coloc_for_gtop.txt")
# length(unique(paste0(a$V6, "_", a$V9))) # 4240
# length(unique(paste0(coloc_res$phecode_abbr, "_", coloc_res$loci))) # 3782
