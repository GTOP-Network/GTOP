#!/usr/bin/env Rscript

setwd("/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project")

# Load required libraries
library(data.table)
library(tidyverse)

CHRPOS_RSID <- fread("2025-06-11-specific_xQTL/input/GMTiP/SNP/GMTiP_SNP_chrpos_rsid.txt", 
                     header = FALSE, col.names = c("SNP", "rsid"))

## Allele frequency of SNV and SV: freq_data ----------------------------------
freq_data <- fread("2025-06-11-specific_xQTL/output/data_plot/data/SNP_SV_addgroup_information.txt")

## Finemapping result of SNV-xQTL
finemapping_SNVs <- fread("2025-06-11-specific_xQTL/output/data_plot/data/SNV_variant_xQTL_finemapping.txt")
table(finemapping_SNVs$xQTL_type)

selected_xQTL <- "SNV_eQTL"
finemapping_SNVs <- finemapping_SNVs[xQTL_type==selected_xQTL]
finemapping_SNVs_control <- finemapping_SNVs[sig_type!="xGene"]
finemapping_SNVs <- finemapping_SNVs[sig_type=="xGene"]

# tmp_df <- finemapping_SNVs %>% group_by(gene_name, cs_id, tissue) %>% summarise(max_pip=max(pip))
# hist(tmp_df$max_pip, breaks = 100)
# length(unique(finemapping_SNVs[sig_type=="xGene"&pip>0]$variant_id)) # 365925
# length(unique(finemapping_SNVs[sig_type=="xGene"]$variant_id)) # 372324

finemapping_SNVs <- finemapping_SNVs %>% 
  group_by(phenotype_id, gene_name, gene_symbol, gene_type, 
           cs_id, cs_size, tissue, xQTL_type) %>% 
  mutate(variant_proxy=paste0(variant_id[pip==max(pip)], collapse = ","))
setDT(finemapping_SNVs)
finemapping_SNVs$variant_proxy <-  gsub(",.+", "", finemapping_SNVs$variant_proxy)

tmp_SNVs <- unique(finemapping_SNVs[variant_proxy==variant_id, .(variant_proxy, pip)])
tmp_SNVs <- tmp_SNVs %>% group_by(variant_proxy) %>% summarise(max_pip=max(pip)) %>% as.data.frame()
setDT(tmp_SNVs)
tmp_SNVs <- merge(tmp_SNVs[, .(rsid=variant_proxy, P=max_pip)], CHRPOS_RSID, by="rsid")
tmp_SNVs[, c("chr", "pos", "ref", "alt") := tstrsplit(SNP, "_")]
tmp_SNVs <- tmp_SNVs[, .(CHR=chr,SNP=rsid,BP=pos,FRQ="-",
                         BETA="-",SE="-",P)]
tmp_SNVs$P <- 1-tmp_SNVs$P
tmp_SNVs$P[tmp_SNVs$P>=1] <- 1
fwrite(tmp_SNVs, "2025-06-11-specific_xQTL/tmp/indepedent_SNVs.input", sep = "\t")
system(command = sprintf("%s --bfile %s --clump 2025-06-11-specific_xQTL/tmp/indepedent_SNVs.input %s --out %s",
                         "/media/bora_A/zhangt/src/bin/plink2",
                         "/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/genotype/GMTiP_SNV",
                         "--clump-r2 0.8 --clump-p1 1 --clump-p2 1 --clump-kb 1000",
                         "2025-06-11-specific_xQTL/tmp/indepedent_SNVs.res"))
tmp_SNVs_res <- fread("2025-06-11-specific_xQTL/tmp/indepedent_SNVs.res.clumps")
tmp_SNVs_res <- tmp_SNVs_res[SP2!="."] %>% separate_rows(SP2, sep=",") %>% as.data.frame() %>% column_to_rownames("SP2")

finemapping_SNVs$variant_proxy[finemapping_SNVs$variant_proxy%in%rownames(tmp_SNVs_res)] <- 
  tmp_SNVs_res[finemapping_SNVs$variant_proxy[finemapping_SNVs$variant_proxy%in%rownames(tmp_SNVs_res)], "ID"]
length(unique(finemapping_SNVs$variant_proxy)) # 8470


## freq_finemapping ------------------------------------------------------------
freq_finemapping <- merge(finemapping_SNVs, 
                          freq_data[,.(variant_proxy=variant_id, 
                                       proxy_af_GMTiP=af_GMTiP, 
                                       proxy_af_EUR=af_EUR, 
                                       proxy_delta_group=delta_group)],
                          by="variant_proxy")

ff_control <- merge(finemapping_SNVs_control, 
                    freq_data[,.(variant_id, 
                                 af_GMTiP, 
                                 af_EUR, 
                                 delta_group,
                                 two_group)],
                    by="variant_id")

length(unique(freq_finemapping$gene_name)) # 7739
freq_finemapping <- merge(freq_finemapping, 
                          freq_data[, .(variant_id, af_GMTiP, af_EUR, two_group)], 
                          by="variant_id")
freq_finemapping$two_group <- case_when(
  freq_finemapping$af_EUR == 0 | freq_finemapping$af_EUR == 1 ~ "EUR_unobserved",
  TRUE ~ freq_finemapping$two_group
)

data_onehot <- freq_finemapping[, .(variant_id, two_group)] 
data_onehot <- model.matrix(~ two_group - 1, data = data_onehot)
data_onehot <- data_onehot*abs(freq_finemapping$pip)
data_onehot <- cbind(freq_finemapping[, .(variant_proxy, phenotype_id)], data_onehot)
data_onehot <- data_onehot %>% group_by(variant_proxy, phenotype_id) %>% 
  mutate(EUR_common=sum(two_groupEUR_common), 
         EUR_rare=sum(two_groupEUR_rare),
         EUR_unobserved=sum(two_groupEUR_unobserved))
freq_finemapping$new_two_group <- apply(
  data_onehot, 1, function(x){
    if(as.numeric(x[7])==max(as.numeric(x[6:8]))){
      return("EUR_rare")
    }else if(as.numeric(x[8])==max(as.numeric(x[6:8]))){
      return("EUR_unobserved")
    }else{
      return("-")
    }
  }
)


## mash 
# eQTL_mash_model <- readRDS("/media/london_B/zouxudong/2024-10-21-aGTEx-main/2025-05-28-mashR-eQTL/output/mashr_flashr_workflow_output/eQTL_file_list.snv.mash.EZ.V_identity.mash_model.rds")
# mash_res <- readRDS("/media/london_B/zouxudong/2024-10-21-aGTEx-main/2025-05-28-mashR-eQTL/output/mashr_flashr_workflow_output/eQTL.nominal.mash.final.results.snv.rds")
# mash_for_ffdata = mash(data.strong, g=get_fitted_g(m.r), fixg=TRUE)


## Add other annotations
GTOP_LD_extend <- function(variant_list1, variant_list2, population){
  variant_list2 <- setdiff(variant_list2, variant_list1)
  tmp_var_list <- unique(c(variant_list1, variant_list2))
  write.table(tmp_var_list, "2025-06-11-specific_xQTL/tmp/tmp_var.list", 
              quote = F, row.names = F, col.names = F)
  write.table(variant_list1, "2025-06-11-specific_xQTL/tmp/tmp_query.list", 
              quote = F, row.names = F, col.names = F)
  system(sprintf("%s --bfile %s --extract %s --ld-snp-list %s --r2-unphased --ld-window-kb 1000 --ld-window-r2 0.8 --out %s",
                 "/media/bora_A/zhangt/src/bin/plink2",
                 "/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/data/genotype/GMTiP_SNV", 
                 "2025-06-11-specific_xQTL/tmp/tmp_var.list", "2025-06-11-specific_xQTL/tmp/tmp_query.list", 
                 "2025-06-11-specific_xQTL/tmp/tmp_var"))
  tmp_var_LD <- fread("2025-06-11-specific_xQTL/tmp/tmp_var.vcor")
  tmp_var_LD <- tmp_var_LD[(tmp_var_LD$ID_A%in%variant_list1 &
                              tmp_var_LD$ID_B%in%variant_list2) | 
                             (tmp_var_LD$ID_B%in%variant_list1 &
                                tmp_var_LD$ID_A%in%variant_list2) ]
  out_df <- tmp_var_LD[, .(ID_A, ID_B, UNPHASED_R2)]
  colnames(out_df) <- c( "IDA_rsid", "IDB_rsid", "UNPHASED_R2")
  return(out_df)
}


PS_SNVs <- fread("2025-06-11-specific_xQTL/output/data_plot/data/positive_selection_SNV.txt")
heQTL_SNVs <- fread("2025-09-30-mash/output/data_plot/EAS_biased_eQTL.txt")
heQTL_SNVs$tissue <- gsub("GMTiP_", "", heQTL_SNVs$tissue)
archaic_SNVs <- fread("2025-06-13-evolution/input/SNP_archaic_orgin.txt")
CHRPOS_RSID$tmp_id <- gsub("_[ATCG].+", "", CHRPOS_RSID$SNP)
CHRPOS_RSID$tmp_id <- gsub("_", ":", CHRPOS_RSID$tmp_id)
CHRPOS_RSID$tmp_id <- gsub("chr", "", CHRPOS_RSID$tmp_id)
archaic_SNVs <- merge(archaic_SNVs, CHRPOS_RSID[,.(ID=tmp_id, rsid)], by="ID")

## calculate high LD variants 
ps_SNV_LD <- GTOP_LD_extend(variant_list1 = PS_SNVs$variant_id, variant_list2 = unique(freq_finemapping$variant_id))
length(unique(PS_SNVs$variant_id))
length(unique(c(PS_SNVs$variant_id, ps_SNV_LD$IDB_rsid)))

he_SNV_LD <- GTOP_LD_extend(variant_list1 = heQTL_SNVs$rsid, variant_list2 = unique(freq_finemapping$variant_id))
heQTL_SNVs1 <- inner_join(heQTL_SNVs[, .(tissue, rsid, Gene)], unique(he_SNV_LD[, .(rsid=IDA_rsid, IDB_rsid)]), by="rsid")
heQTL_SNVs1 <- rbind(heQTL_SNVs[, .(tissue, rsid, Gene)], heQTL_SNVs1[, .(tissue, rsid=IDB_rsid, Gene)])
heQTL_SNVs1 <- paste0(heQTL_SNVs1$tissue, "_", heQTL_SNVs1$rsid, "_", heQTL_SNVs1$Gene)

length(unique(heQTL_SNVs$rsid))
length(unique(c(he_SNV_LD$rsid, he_SNV_LD$IDB_rsid)))

arc_SNV_LD <- GTOP_LD_extend(variant_list1 = archaic_SNVs$rsid, variant_list2 = unique(freq_finemapping$variant_id))
length(unique(archaic_SNVs$rsid))
length(unique(c(archaic_SNVs$rsid, arc_SNV_LD$IDB_rsid)))

freq_finemapping$posSelection <- "-"
freq_finemapping$posSelection[freq_finemapping$variant_id%in%unique(c(PS_SNVs$variant_id, ps_SNV_LD$IDB_rsid))] <- "Yes"

freq_finemapping$popBias <- "-"
tmp_id <- paste0(freq_finemapping$tissue, "_", freq_finemapping$variant_id, "_", freq_finemapping$gene_name)
table(tmp_id%in%heQTL_SNVs1)
freq_finemapping$popBias[tmp_id%in%heQTL_SNVs1] <- "Yes"

freq_finemapping$arcSNP <- "-"
freq_finemapping$arcSNP[freq_finemapping$variant_id%in%unique(c(archaic_SNVs$rsid, arc_SNV_LD$IDB_rsid))] <- "Yes"

data_onehot <- freq_finemapping[, .(variant_id, posSelection)] 
data_onehot <- model.matrix(~ posSelection - 1, data = data_onehot)
data_onehot <- data_onehot*abs(freq_finemapping$pip)
data_onehot <- cbind(freq_finemapping[, .(variant_proxy, phenotype_id)], data_onehot)
data_onehot <- data_onehot %>% group_by(variant_proxy, phenotype_id) %>% 
  mutate(SumNo=sum(`posSelection-`), 
         SumYes=sum(posSelectionYes))
freq_finemapping$new_posSelection <- apply(
  data_onehot, 1, function(x){
    if(as.numeric(x[6])==max(as.numeric(x[5:6]))){
      return("Yes")
    }else{
      return("-")
    }
  }
)

rm(data_onehot)
data_onehot <- freq_finemapping[, .(variant_id, popBias)] 
data_onehot <- model.matrix(~ popBias - 1, data = data_onehot)
data_onehot <- data_onehot*abs(freq_finemapping$pip)
data_onehot <- cbind(freq_finemapping[, .(variant_proxy, phenotype_id)], data_onehot)
data_onehot <- data_onehot %>% group_by(variant_proxy, phenotype_id) %>% 
  mutate(SumNo=sum(`popBias-`), 
         SumYes=sum(popBiasYes))
freq_finemapping$new_popBias <- apply(
  data_onehot, 1, function(x){
    if(as.numeric(x[6])==max(as.numeric(x[5:6]))){
      return("Yes")
    }else{
      return("-")
    }
  }
)


rm(data_onehot)
data_onehot <- freq_finemapping[, .(variant_id, arcSNP)] 
data_onehot <- model.matrix(~ arcSNP - 1, data = data_onehot)
data_onehot <- data_onehot*abs(freq_finemapping$pip)
data_onehot <- cbind(freq_finemapping[, .(variant_proxy, phenotype_id)], data_onehot)
data_onehot <- data_onehot %>% group_by(variant_proxy, phenotype_id) %>% 
  mutate(SumNo=sum(`arcSNP-`), 
         SumYes=sum(arcSNPYes))
freq_finemapping$new_arcSNP <- apply(
  data_onehot, 1, function(x){
    if(as.numeric(x[6])==max(as.numeric(x[5:6]))){
      return("Yes")
    }else{
      return("-")
    }
  }
)


## compare with GTEx
tissue_change <- c("Adipose"="Adipose_Subcutaneous",
                   "Adrenal_Gland"="Adrenal_Gland",
                   "Liver"="Liver",
                   "Muscle"="Muscle_Skeletal",
                   "Pancreas_Body"="Pancreas",
                   "Pancreas_Head"="Pancreas",
                   "Pancreas_Tail"="Pancreas",
                   "Skin"="Skin_Not_Sun_Exposed_Suprapubic",
                   "Spleen"="Spleen",
                   "Whole_Blood"="Whole_Blood")

GTEx_support_res <- pbmcapply::pbmclapply(seq_len(length(tissue_change)), function(f_index){
  GTEx_V8_egenes <- fread(sprintf("2025-06-11-specific_xQTL/data/GTEx_v8/GTEx_Analysis_v8_eQTL/%s.v8.egenes.txt.gz", 
                                  tissue_change[f_index]))
  GTEx_V8_egenes <- GTEx_V8_egenes[qval < 0.05]
  
  GTEx_V8_sig <- fread(sprintf("2025-06-11-specific_xQTL/data/GTEx_v8/GTEx_Analysis_v8_eQTL/%s.v8.signif_variant_gene_pairs.txt.gz", 
                               tissue_change[f_index]))
  GTEx_V8_sig <- GTEx_V8_sig[, c("gene_id", "variant_id")]
  
  GTEx_V8_sig <- inner_join(GTEx_V8_sig, GTEx_V8_egenes[, c("gene_id", "gene_name")], by="gene_id")
  GTEx_V8_sig$tissue <- names(tissue_change)[f_index]
  return(GTEx_V8_sig)
}, mc.cores = 10, mc.preschedule = FALSE)
GTEx_support_df <- rbindlist(GTEx_support_res)
GTEx_support_df$variantId <- gsub("_b38", "", GTEx_support_df$variant_id)
GTEx_support_df <- merge(GTEx_support_df, CHRPOS_RSID[, .(variantId=SNP, rsid)], by="variantId")
GTEx_support_id <- paste0(GTEx_support_df$tissue, "_", GTEx_support_df$rsid, "_", GTEx_support_df$gene_name)

freq_finemapping$GTEx_support <- "-"
tmp_id <- paste0(freq_finemapping$tissue, "_", freq_finemapping$variant_id, "_", freq_finemapping$gene_symbol)
table(tmp_id%in%GTEx_support_id)
freq_finemapping$GTEx_support[tmp_id%in%GTEx_support_id] <- "Yes"

GTEx_expGene_res <- pbmcapply::pbmclapply(seq_len(length(tissue_change)), function(f_index){
  GTEx_V8_egenes <- fread(sprintf("2025-06-11-specific_xQTL/data/GTEx_v8/GTEx_Analysis_v8_eQTL/%s.v8.egenes.txt.gz", 
                                  tissue_change[f_index]))
  GTEx_V8_egenes <- GTEx_V8_egenes
  
  GTEx_V8_egenes$tissue <- names(tissue_change)[f_index]
  return(GTEx_V8_egenes)
}, mc.cores = 10, mc.preschedule = FALSE)
GTEx_expGene_df <- rbindlist(GTEx_expGene_res)
GTEx_expGene_id <- paste0(GTEx_expGene_df$tissue, "_", GTEx_expGene_df$gene_name)

freq_finemapping$GTEx_expGene <- "-"
tmp_id <- paste0(freq_finemapping$tissue, "_", freq_finemapping$gene_symbol)
table(tmp_id%in%GTEx_expGene_id)
freq_finemapping$GTEx_expGene[tmp_id%in%GTEx_expGene_id] <- "Yes"


GTEx_eGene_res <- pbmcapply::pbmclapply(seq_len(length(tissue_change)), function(f_index){
  GTEx_V8_egenes <- fread(sprintf("2025-06-11-specific_xQTL/data/GTEx_v8/GTEx_Analysis_v8_eQTL/%s.v8.egenes.txt.gz", 
                                  tissue_change[f_index]))
  GTEx_V8_egenes <- GTEx_V8_egenes[qval < 0.05]
  
  GTEx_V8_egenes$tissue <- names(tissue_change)[f_index]
  return(GTEx_V8_egenes)
}, mc.cores = 10, mc.preschedule = FALSE)
GTEx_eGene_df <- rbindlist(GTEx_eGene_res)
GTEx_eGene_id <- paste0(GTEx_eGene_df$tissue, "_", GTEx_eGene_df$gene_name)

freq_finemapping$GTEx_eGene <- "-"
tmp_id <- paste0(freq_finemapping$tissue, "_", freq_finemapping$gene_symbol)
table(tmp_id%in%GTEx_eGene_id)
freq_finemapping$GTEx_eGene[tmp_id%in%GTEx_eGene_id] <- "Yes"

fwrite(freq_finemapping, file = "2025-06-11-specific_xQTL/output/data_plot/data/freq_finemapping.txt", sep = "\t")

