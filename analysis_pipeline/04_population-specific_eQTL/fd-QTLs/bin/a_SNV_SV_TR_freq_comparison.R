
# SNP Analysis Pipeline: Allele Frequency Comparisons
# Date: 2025-06-21

library(data.table)
library(dplyr)
library(Biostrings)
library(ggplot2)
library(ggpubr)

# SECTION 1: SNP FREQUENCY ANALYSIS (1000 Genomes & gnomAD) ---------------

#' Load and process GMTiP SNP frequency data
#' Note: Original dataset contains ~5.8 million SNPs
process_gtop_snps <- function() {
  # Read SNP rsID mapping
  CHRPOS2rsid <- fread("2025-06-11-specific_xQTL/input/GMTiP/SNP/GMTiP_SNP_chrpos_rsid.txt",
                       col.names = c("SNP", "rsid"))
  
  # Load allele frequencies from GMTiP dataset
  GTOP_freq <- fread("2025-06-11-specific_xQTL/input/GMTiP/SNP/SNP_EAS.afreq")
  GTOP_freq <- merge(GTOP_freq[, .(rsid = ID, af_GMTiP = ALT_FREQS)], 
                      CHRPOS2rsid,
                      by = "rsid")[, .(SNP, rsid, af_GMTiP)]
  
  return(GTOP_freq)
}


#' Load and process 1000 Genomes population frequencies
process_1kg_snps <- function() {
  populations <- c("AFR", "EAS", "EUR")
  freq_list <- rbindlist(lapply(populations, function(pop) {
    file_path <- sprintf("/media/bora_A/zhangt/single_cell_datasets/seven_genotypes_summary/population/1000G_%s.afreq", pop)
    fread(file_path)[, Population := pop]
  }))
  
  # Reshape to wide format (one row per SNP)
  OKG_freq <- dcast(freq_list, ID + REF + ALT ~ Population, value.var = "ALT_FREQS") %>%
    .[, SNP := paste0("chr", gsub(":", "_", ID))] %>%
    .[, c("chr", "pos") := tstrsplit(SNP, "_", keep = 1:2)] %>% 
    .[, .(SNP, chr, pos, af_AFR = AFR, af_EAS = EAS, af_EUR = EUR)]
    
  OKG_freq <- OKG_freq[order(chr, as.numeric(pos))]
  
  return(OKG_freq)
}


#' Load and process gnomAD population frequencies
process_gnomad_snps <- function() {
  # Read gnomAD allele frequencies
  gnomad_freq <- fread("2025-06-11-specific_xQTL/input/gnomAD/SNP/gnomAD_SNP_AF.txt", 
                       col.names = c("SNP", "rsid", "ac_AFR", "an_AFR", "af_AFR", "ac_EAS", "an_EAS", "af_EAS", "ac_EUR", "an_EUR", "af_EUR"))
  
  # Format SNP: convert "chr:pos" to "chr_pos"
  gnomad_freq[, SNP := gsub(":", "_", SNP)]

  # We'll keep rsid only for GMTiP matching, then remove it
  gnomad_freq <- gnomad_freq[, .(SNP, af_AFR, af_EAS, af_EUR, ac_AFR, ac_EAS, ac_EUR, an_AFR, an_EAS, an_EUR)]
  
  gnomad_freq <- gnomad_freq[, c("chr", "pos") := tstrsplit(SNP, "_", keep = 1:2)] %>% 
    .[, .(SNP, chr, pos, af_AFR, af_EAS, af_EUR, population, ac_AFR, ac_EAS, ac_EUR, an_AFR, an_EAS, an_EUR)]
  
  gnomad_freq <- gnomad_freq[order(chr, as.numeric(pos))]
  return(gnomad_freq)
}

# Main SNP processing pipeline -------------------------------------------------
SNP_freq_1KG <- process_1kg_snps()
fwrite(SNP_freq_1KG, "2025-06-11-specific_xQTL/output/data_plot/1KG_SNP_information.txt", sep = "\t")

freq_gnomAD <- process_gnomad_snps()
fwrite(freq_gnomAD, "2025-06-11-specific_xQTL/output/data_plot/gnomAD_SNP_information.txt", sep = "\t")

freq_1KG_EAS <- fread("2025-06-11-specific_xQTL/output/data_plot/1KG_SNP_information.txt") # 12867045
freq_1KG_EAS <- freq_1KG_EAS[af_EAS > 0 | af_EUR > 0 ] # 10044462

freq_gnomAD_EAS <- fread("2025-06-11-specific_xQTL/output/data_plot/gnomAD_SNP_information.txt") # 27568964
freq_gnomAD_EAS <- freq_gnomAD_EAS[af_EAS > 0 | af_EUR > 0 ] # 27122917
freq_gnomAD_EAS <- freq_gnomAD_EAS[an_EUR > 30000 & an_EAS > 5000 ] # 19012414  # an_AFR > 10000 & 

GTOP_freq <- process_gtop_snps() # 6037504
GTOP_freq[, c("chr", "pos", "ref", "alt") := tstrsplit(SNP, "_", keep = 1:4)]


ggvenn::ggvenn(list("GTOP"=GTOP_freq$SNP, "1KGP"=freq_1KG_EAS$SNP, "gnomAD"=freq_gnomAD_EAS$SNP), show_percentage = F)


# Process 1000 Genomes data
freq_1KG_res <- merge(GTOP_freq, freq_1KG_EAS, by = c("SNP", "chr", "pos")) %>%
  .[, source := "1KG"]  # Add source column

freq_gnomAD_EAS <- freq_gnomAD_EAS[!SNP %in% freq_1KG_res$SNP][
  , -c("ac_AFR", "ac_EAS", "ac_EUR", "an_AFR", "an_EAS", "an_EUR")]

# Process gnomAD data (with different frequency threshold)
freq_gnomAD_res <- merge(GTOP_freq, freq_gnomAD_EAS, by = c("SNP", "chr", "pos")) %>%
  .[, source := "gnomAD"]


# Merge results from both resources
combined_results <- list(
  freq_1KG_res,
  freq_gnomAD_res_clean  # Add only new SNPs
) %>% rbindlist()

# Mark SNPs found in both datasets
dual_source_ids <- intersect(freq_1KG_res$SNP, freq_gnomAD_res$SNP)
combined_results[SNP %in% dual_source_ids, source := "1KG&gnomAD"]

combined_results <- combined_results[order(chr, as.numeric(pos))]
combined_results[, var_type := fifelse(
  nchar(ref) == 1L & nchar(alt) == 1L, 
  "SNP", 
  "INDEL"
)]
table(combined_results$var_type)

combined_results <- combined_results[, .(SNP, rsid, chr, pos, var_type, af_GMTiP)]
combined_results$pos <- as.numeric(combined_results$pos)

table(combined_results$source)
fwrite(combined_results, "2025-06-11-specific_xQTL/output/data_plot/SNP_information.txt", sep = "\t")

# combined_results <- fread("2025-06-11-specific_xQTL/output/data_plot/SNP_information.txt")
length(unique(combined_results$SNP)) # 5630528
table(abs(combined_results$af_GMTiP-combined_results$af_EAS) < 0.2)

check_SNPs <- subset(GTOP_freq, !SNP%in%combined_results$SNP)


selected_SNPs <- combined_results[population!="-"]
table(selected_SNPs$var_type) # INDEL:67484 SNP:1039979

table(selected_SNPs$source)
table(selected_SNPs$population)

sub_SNPs <- selected_SNPs[var_type=="SNP"]
cor.test(sub_SNPs$af_GMTiP, sub_SNPs$af_EAS) # cor: 0.9935447; p-value < 2.2e-16
png("2025-06-11-specific_xQTL/output/data_plot/plot/AF_correlation_SNP_GMTiP_1KG_gnomAD.png", width = 3000, height = 3500, res = 600)
plot(sub_SNPs$af_GMTiP, sub_SNPs$af_EAS, 
     xlab="AF(GMTiP)", ylab="AF(1KG|gnomAD)", las=1, pch=21, cex = 0.9)
dev.off()

sub_INDELs <- selected_SNPs[var_type=="INDEL"]
cor.test(sub_INDELs$af_GMTiP, sub_INDELs$af_EAS) # cor: 0.9944385; p-value < 2.2e-16
png("2025-06-11-specific_xQTL/output/data_plot/plot/AF_correlation_INDEL_GMTiP_1KG_gnomAD.png", width = 3000, height = 3500, res = 600)
plot(sub_INDELs$af_GMTiP, sub_INDELs$af_EAS, 
     xlab="AF(GMTiP)", ylab="AF(1KG|gnomAD)", las=1, pch=21, cex = 0.9)
dev.off()

# Filter and save significant SNPs
table(abs(selected_SNPs$af_GMTiP-selected_SNPs$af_EAS) < 0.2)  # FASLE: 71 TRUE: 1107392

# selected_SNPs <- selected_SNPs[abs(af_GMTiP-af_EAS) < 0.2]
fwrite(selected_SNPs, "2025-06-11-specific_xQTL/output/data_plot/SNP_EAS_specific.txt", sep = "\t")

pdf("2025-06-11-specific_xQTL/output/data_plot/plot/Proportion_of_SNP_INDEL.pdf", width = 4, height = 4)
plot(eulerr::euler(list("Other SNPs"=combined_results$SNP[combined_results$var_type=="SNP"], 
                        "EAS_specific"=selected_SNPs$SNP[selected_SNPs$var_type=="SNP"]), 
                   shape = "circle"), quantities = TRUE, 
     main=paste0("Proportion: ", round(length(selected_SNPs$SNP[selected_SNPs$var_type=="SNP"])/
                   length(combined_results$SNP[combined_results$var_type=="SNP"])*100, 2), "%"))

plot(eulerr::euler(list("Other INDELs"=combined_results$SNP[combined_results$var_type=="INDEL"], 
                        "EAS_specific"=selected_SNPs$SNP[selected_SNPs$var_type=="INDEL"]), 
                   shape = "circle"), quantities = TRUE, 
     main=paste0("Proportion: ", round(length(selected_SNPs$SNP[selected_SNPs$var_type=="INDEL"])/
                                         length(combined_results$SNP[combined_results$var_type=="INDEL"])*100, 2), "%"))
dev.off()


# STRUCTURAL VARIANT (SV) ANALYSIS ----------------------------------------
# Execute SV processing
# GMTiP_SV_freq <- fread("2025-06-11-specific_xQTL/input/GMTiP/SV/SV_EAS.afreq") %>%
#   .[!is.na(ALT_FREQS) & !duplicated(ID), .(ID, af_GMTiP = ALT_FREQS)]
# 
# GMTiP_SV_freq$tmpID <- apply(GMTiP_SV_freq, 1, function(x){tmpx <- strsplit(x[1], split="_")[[1]];
# tmpy <- round(0.5-abs(0.5-as.numeric(x[2])), digits = 3)
# paste0(c(tmpx[1], tmpx[2], tmpx[3], tmpy), collapse = "_")})
# 
# 
# SV_ID_compare <- fread("/media/london_B/lixing/2024-08-29-TR-AsianGTEX/2025-01-16-LRS-SRS-SV-Compare/output/SV_bed_diff_study/GMTiP_LRS.inter.othreStudy_sameType.results")
# SV_ID_compare$tmpID <- sapply(strsplit(SV_ID_compare$ID, split=":|_"), function(x){
#   paste0(c(x[1], x[2], x[4], round(as.numeric(x[7]), digits = 3)), collapse = "_")
# })
# 
# ggvenn::ggvenn(list("GMTiP_z"=GMTiP_SV_freq$tmpID, "GMTiP_l"=SV_ID_compare$tmpID))
# 
# setdiff(GMTiP_SV_freq$tmpID, SV_ID_compare$tmpID)
# SV_ID_compare$tmpID[SV_ID_compare$tmpID=="chr3_90855061_INS_0.212"] <- "chr3_90855061_INS_0.213"
# 
# SV_1KG_compare <- SV_ID_compare[overlap_1KG=="Reported", ]
# SV_1KG_compare <- unique(SV_1KG_compare[, .(tmpID, `1KG-SVID`)])
# SV_1KG_compare$`1KG-SVID` <- sapply(strsplit(SV_1KG_compare$`1KG-SVID`, split=":"), function(x){x[2]})
# 
# GTOP_freq_res1 <- merge(GMTiP_SV_freq, 
#                          SV_1KG_compare[, .(tmpID, db_ID = `1KG-SVID`)], 
#                          by = "tmpID") %>%
#   merge(SV_freq_1KG, by.x = "db_ID", by.y = "CHRPOS") %>%
#   # classify_population(freq_diff = 0.5) %>%
#   .[, .(ID, db_ID, af_GMTiP, af_AFR, af_EAS, af_EUR)]
# GTOP_freq_res1[, source := "1KG"]
# 
# plot(GTOP_freq_res1$af_GMTiP, GTOP_freq_res1$af_EAS)
# 
# 
# 
# SV_gnomAD_compare <- SV_ID_compare[overlap_gnomAD=="Reported", ]
# SV_gnomAD_compare <- unique(SV_gnomAD_compare[, .(tmpID, `gnomAD-SVID`)])
# SV_gnomAD_compare$`gnomAD-SVID` <- sapply(strsplit(SV_gnomAD_compare$`gnomAD-SVID`, split=":"), function(x){x[1]})
# 
# gnomAD_clean <- fread("/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project/2025-06-11-specific_xQTL/input/gnomAD/SV/gnomAD_SV.bed")
# gnomAD_clean$tmpID <- apply(gnomAD_clean, 1, function(x){
#   gsub(" ", "", paste0(c(x[1], as.numeric(x[2])+1, x[3], strsplit(x[4], split="_")[[1]][3]), collapse = "_"))
# })
# 
# GTOP_freq_res2 <- merge(GMTiP_SV_freq, 
#                          SV_gnomAD_compare[, .(tmpID, dbx_ID = `gnomAD-SVID`)], 
#                          by = "tmpID") %>%
#   # merge(SV_freq_1KG, by.x = "db_ID", by.y = "CHRPOS") %>%
#   merge(gnomAD_clean[, .(dbx_ID=tmpID, db_ID = V4, af_AFR = V9, af_EAS = V12, af_EUR = V15,
#                          ac_AFR = V7, ac_EAS = V10, ac_EUR = V13, 
#                          an_AFR = V8, an_EAS = V11, an_EUR = V14)],
#   by = "dbx_ID"
# )
# 
# plot(GTOP_freq_res1$af_GMTiP, GTOP_freq_res1$af_EAS, cex=0.5)
# 
# length(unique(c(GTOP_freq_res1$ID, GTOP_freq_res2$ID)))


GMTiP_SV_freq <- fread("2025-06-11-specific_xQTL/input/GMTiP/SV/SV_EAS.afreq") %>%
  .[!is.na(ALT_FREQS) & !duplicated(ID), .(ID, af_GMTiP = ALT_FREQS)]

# Process 1000 Genomes SV data
freq_list <- rbindlist(lapply(c("AFR", "EUR", "EAS"), function(pop) {
  fread(sprintf("2025-06-11-specific_xQTL/input/1KG/SV/1000G_%s.afreq", pop))[, Population := pop]
}))

SV_freq_1KG <- dcast(freq_list, ID + REF + ALT ~ Population, 
                  value.var = "ALT_FREQS")[, .(CHRPOS = ID, af_AFR = AFR, af_EAS = EAS, af_EUR = EUR)]

# GTOP_freq_res1 <- merge(GMTiP_SV_freq, 
#                          GMTiP21KG_clean[, .(ID = V4, db_ID = V10)], 
#                          by = "ID") %>%
#   merge(SV_freq_1KG, by.x = "db_ID", by.y = "CHRPOS") %>%
#   classify_population(freq_diff = 0.5) %>%
#   .[, .(ID, db_ID, af_GMTiP, af_AFR, af_EAS, af_EUR, population)] %>%
#   .[, source := "1KG"]


# Load mapping between GMTiP and 1000 Genomes SVs
GMTiP21KG <- fread("2025-06-11-specific_xQTL/input/GMTiP/SV/GMTiP_LRS_inter_1KG_SV.txt")
type_list <- sapply(strsplit(GMTiP21KG$V4, split="_"), function(x){x[3]})
table(type_list)

GMTiP21KG_list <- split(GMTiP21KG, by = "V4")
GMTiP21KG_cleanlist <- pbmcapply::pbmclapply(seq_len(length(GMTiP21KG_list)), function(list_index) {
  tmp_data <- GMTiP21KG_list[[list_index]]
  tmp_type <- strsplit(tmp_data$V4[1], split="_")[[1]][3]
  tmp_name <- sapply(strsplit(tmp_data$V4, split="_"), function(x){paste0(x[c(1,2,3,5)], collapse = "_")})
  index1 <- tmp_name == tmp_data$V10
  if(any(index1)){
    return(tmp_data[which(index1)[1], ])
  }else{
    pos1 <- as.numeric(strsplit(tmp_data$V4, split="_")[[1]][2])
    pos2 <- as.numeric(sapply(strsplit(tmp_data$V10, split="_"), function(x){x[2]}))
    index2 <- pos1 == pos2
    index3 <- abs(nchar(tmp_data$V6) - nchar(tmp_data$V12)) <= 50 &
      abs(nchar(tmp_data$V5) - nchar(tmp_data$V11)) <= 50
    if(any( index2 & index3 )){
      sub_tmp_data <- tmp_data[index2 & index3, ]
      if(tmp_type == "INS"){
        index4 <- apply(sub_tmp_data, 1, function(row_info){
          alignment <- pairwiseAlignment(pattern = row_info[6], subject = row_info[12])
          nmatch(alignment)/min(nchar(row_info[6]), nchar(row_info[12])) > 0.8
        })
        if(any(index4)){
          return(sub_tmp_data[which(index4)[1], ])
        }
      }else if(tmp_type == "DEL"){
        index4 <- apply(sub_tmp_data, 1, function(row_info){
          alignment <- pairwiseAlignment(pattern = row_info[5], subject = row_info[11])
          nmatch(alignment)/min(nchar(row_info[5]), nchar(row_info[11])) > 0.8
        })
        if(any(index4)){
          return(sub_tmp_data[which(index4)[1], ])
        }
      }else{
        len1 <- as.numeric(strsplit(sub_tmp_data$V4, split="_")[[1]][5])
        len2 <- as.numeric(sapply(strsplit(sub_tmp_data$V10, split="_"), function(x){x[4]}))
        return(sub_tmp_data[which.min(abs(len1-len2))[1], ])
      }
    } else if(any(index3)) {
      sub_tmp_data <- tmp_data[index3, ]
      len1 <- as.numeric(strsplit(sub_tmp_data$V4, split="_")[[1]][5])
      len2 <- as.numeric(sapply(strsplit(sub_tmp_data$V10, split="_"), function(x){x[4]}))
      if(tmp_type == "INS"){
        index4 <- apply(sub_tmp_data, 1, function(row_info){
          alignment <- pairwiseAlignment(pattern = row_info[6], subject = row_info[12])
          nmatch(alignment)/min(nchar(row_info[6]), nchar(row_info[12])) > 0.8
        })
        if(any(index4)){
          return(sub_tmp_data[which.min(abs(len1-len2))[1], ])
        }
      }else if(tmp_type == "DEL"){
        index4 <- apply(sub_tmp_data, 1, function(row_info){
          alignment <- pairwiseAlignment(pattern = row_info[5], subject = row_info[11])
          nmatch(alignment)/min(nchar(row_info[5]), nchar(row_info[11])) > 0.8
        })
        if(any(index4)){
          return(sub_tmp_data[which(index4)[1], ])
        }
      }else{
        len1 <- as.numeric(strsplit(sub_tmp_data$V4, split="_")[[1]][5])
        len2 <- as.numeric(sapply(strsplit(sub_tmp_data$V10, split="_"), function(x){x[4]}))
        return(sub_tmp_data[which.min(abs(len1-len2))[1], ])
      }
    }
  }
}, mc.cores = 40)

GMTiP21KG_clean <- rbindlist(GMTiP21KG_cleanlist)

# Merge datasets with validation
GTOP_freq_res1 <- merge(GMTiP_SV_freq,
      GMTiP21KG_clean[, .(ID = V4, db_ID = V10)],
      by = "ID") %>%
  merge(SV_freq_1KG, by.x = "db_ID", by.y = "CHRPOS") %>%
  classify_population(freq_diff = 0.5) %>%
  .[, .(ID, db_ID, af_GMTiP, af_AFR, af_EAS, af_EUR, population)] %>%
  .[, source := "1KG"]


# Process gnomAD SV data
GMTiP2gnomAD <- fread("2025-06-11-specific_xQTL/input/GMTiP/SV/GMTiP_LRS_inter_gnomAD_SV.txt")

GMTiP2gnomAD_list <- split(GMTiP2gnomAD, by = "V4")
GMTiP2gnomAD_cleanlist <- pbmcapply::pbmclapply(seq_len(length(GMTiP2gnomAD_list)), function(list_index) {
  tmp_data <- GMTiP2gnomAD_list[[list_index]]
  tmp_name <- sapply(strsplit(tmp_data$V4, split="_"), function(x){paste0(x[c(1,2,3,5)], collapse = "_")})
  index1 <- tmp_name == tmp_data$V10
  if(any(index1)){
    return(tmp_data[which(index1)[1], ])
  }else{
    pos1 <- as.numeric(strsplit(tmp_data$V4, split="_")[[1]][2])
    pos2 <- as.numeric(sapply(strsplit(tmp_data$V10, split="_"), function(x){x[2]}))
    index2 <- pos1 == pos2
    len1 <- as.numeric(strsplit(tmp_data$V4, split="_")[[1]][5])
    len2 <- as.numeric(sapply(strsplit(tmp_data$V10, split="_"), function(x){x[4]}))
    index3 <- abs(len1 - len2) < 50
    if(any( index2 & index3 )){
        sub_tmp_data <- tmp_data[index2 & index3, ]
        len21 <- as.numeric(sapply(strsplit(sub_tmp_data$V10, split="_"), function(x){x[4]}))
        return(sub_tmp_data[which.min(abs(len1-len21))[1], ])
    } else if(any(index3)) {
        sub_tmp_data <- tmp_data[index3, ]
        len21 <- as.numeric(sapply(strsplit(sub_tmp_data$V10, split="_"), function(x){x[4]}))
        return(sub_tmp_data[which.min(abs(len1-len21))[1], ])
    }
  }
}, mc.cores = 40)

GMTiP2gnomAD_clean <- rbindlist(GMTiP2gnomAD_cleanlist)
freq_gnomAD_res <- merge(
  GMTiP_SV_freq, 
  GMTiP2gnomAD_clean[, .(ID = V4, db_ID = V10, af_AFR = V15, af_EAS = V18, af_EUR = V21,
                         ac_AFR = V13, ac_EAS = V16, ac_EUR = V19, 
                         an_AFR = V14, an_EAS = V17, an_EUR = V20)],
  by = "ID"
)

plot(freq_gnomAD_res$af_GMTiP, freq_gnomAD_res$af_EAS)
test_data <- freq_gnomAD_res[an_AFR > 10000 & an_EUR > 30000 & an_EAS > 2000 & ac_EAS > 100 ]
plot(test_data$af_GMTiP, test_data$af_EAS)
freq_gnomAD_res <- test_data[, -c("ac_AFR", "ac_EAS", "ac_EUR", "an_AFR", "an_EAS", "an_EUR")]

GTOP_freq_res2 <- classify_population(freq_gnomAD_res, freq_diff = 0.5)[
  , .(ID, db_ID, af_GMTiP, af_AFR, af_EAS, af_EUR, population)
][, source := "gnomAD"]

GTOP_freq_res1[ID %in% GTOP_freq_res2$ID, source := "1KG&gnomAD"]
GMTiP_SV_freq_res <- rbind(
  GTOP_freq_res1,
  GTOP_freq_res2[!ID %in% GTOP_freq_res1$ID]
)

# Final SV processing and output
fwrite(GMTiP_SV_freq_res, "2025-06-11-specific_xQTL/output/data_plot/SV_information.txt", sep = "\t")
plot(GMTiP_SV_freq_res$af_GMTiP, GMTiP_SV_freq_res$af_EAS)


# GMTiP_SV_freq_res <- fread("2025-06-11-specific_xQTL/output/data_plot/SV_information.txt")
selected_SVs <- GMTiP_SV_freq_res[population != "-"]
cor.test(selected_SVs$af_GMTiP, selected_SVs$af_EAS) # cor: 0.9382915; p-value < 2.2e-16

png("2025-06-11-specific_xQTL/output/data_plot/plot/AF_correlation_SV_GMTiP_1KG_gnomAD.png", width = 3000, height = 3500, res = 600)
plot(selected_SVs$af_GMTiP, selected_SVs$af_EAS, 
     xlab="AF(GMTiP)", ylab="AF(1KG|gnomAD)", las=1, pch=21, cex = 0.9)
dev.off()

# Filter and save significant SNPs
table(abs(selected_SVs$af_GMTiP-selected_SVs$af_EAS) < 0.2)  # FASLE: 64 TRUE: 2143

# selected_SVs <- selected_SVs[abs(af_GMTiP-af_EAS) < 0.2]
fwrite(selected_SVs, "2025-06-11-specific_xQTL/output/data_plot/SV_EAS_specific.txt", sep = "\t")

pdf("2025-06-11-specific_xQTL/output/data_plot/plot/Proportion_of_SV.pdf", width = 4, height = 4)
plot(eulerr::euler(list("Other SVs"=GMTiP_SV_freq_res$ID, 
                        "EAS_specific"=selected_SVs$ID), 
                   shape = "circle"), quantities = TRUE, 
     main=paste0("Proportion: ", round(length(selected_SVs$ID)/
                                         length(GMTiP_SV_freq_res$ID)*100, 2), "%"))
dev.off()


## Merge SNV and SV
freq_SNP_data <- fread("2025-06-11-specific_xQTL/output/data_plot/data/SNP_information.txt")[, Type := "SNP"]
freq_SV_data <- fread("2025-06-11-specific_xQTL/output/data_plot/data/SV_information.txt")[, Type := "SV"]
freq_SNP_data$len <- sapply(strsplit(freq_SNP_data$SNP, split = "_"), function(x){max(nchar(x[3:4]))})
freq_SV_data$len <- sapply(strsplit(freq_SV_data$ID, split = "_"), function(x){x[length(x)]})

freq_data <- rbind(freq_SNP_data[, .(variant_id=rsid, af_GMTiP, af_EAS, af_EUR, Type, len)],
                   freq_SV_data[, .(variant_id=ID, af_GMTiP, af_EAS, af_EUR, Type, len)]) %>% 
  mutate(
    EAS_af = af_GMTiP,
    EUR_af = af_EUR
  ) %>% 
  mutate(
    delta_af = abs(EAS_af - EUR_af),
    diff_af = EAS_af - EUR_af
  ) %>%
  mutate(
    diff_group = case_when(
      # EAS_af > 0.05 & EUR_af < 0.01 ~ "EAS_specific (low EUR)",
      diff_af >= 0.50 ~ "AF > 0.5 (low EUR)",
      diff_af >= 0.40 ~ "AF > 0.4 (low EUR)",
      diff_af >= 0.30 ~ "AF > 0.3 (low EUR)",
      diff_af >= 0.20 ~ "AF > 0.2 (low EUR)",
      diff_af >= 0.10 ~ "AF > 0.1 (low EUR)",
      diff_af > 0 ~ "AF > 0 (low EUR)",
      # EAS_af < 0.95 & EUR_af > 0.99 ~ "EAS_specific (high EUR)",
      diff_af <= -0.50 ~ "AF > 0.5 (high EUR)",
      diff_af <= -0.40 ~ "AF > 0.4 (high EUR)",
      diff_af <= -0.30 ~ "AF > 0.3 (high EUR)",
      diff_af <= -0.20 ~ "AF > 0.2 (high EUR)",
      diff_af <= -0.10 ~ "AF > 0.1 (high EUR)",
      diff_af <= 0 ~ "AF > 0 (high EUR)",
    )
  ) %>%
  mutate(
    delta_group = case_when(
      # EAS_af > 0.05 & EUR_af < 0.01 ~ "EAS_specific",
      # EAS_af < 0.95 & EUR_af > 0.99 ~ "EAS_specific",
      delta_af > 0.50 ~ "|ΔAF| > 0.5",
      delta_af > 0.40 ~ "|ΔAF| > 0.4",
      delta_af > 0.30 ~ "|ΔAF| > 0.3",
      delta_af > 0.20 ~ "|ΔAF| > 0.2",
      delta_af > 0.10 ~ "|ΔAF| > 0.1",
      TRUE ~ "|ΔAF| > 0"
    )
  ) %>% 
  mutate(
    two_group = case_when(
      EAS_af >= 0.05 & EUR_af < 0.01 ~ "EUR_rare",
      EAS_af <= 0.95 & EUR_af > 0.99 ~ "EUR_rare",
      TRUE ~ "EUR_common"
    )
  )
fwrite(freq_data, "2025-06-11-specific_xQTL/output/data_plot/data/SNP_SV_addgroup_information.txt", sep = "\t")



# TANDEM REPEAT (TR) ANALYSIS --------------------------------------------------
# Load GMTiP TR frequencies
GMTiP_TR_info <- fread("/media/london_B/lixing/2024-08-29-TR-AsianGTEX/2024-12-11-TR-LongReads/output/QTL_Filter/GMTip_LRS_TR_131INDs.Miss85.AF95.AC1.dosage.impute.simpleTR.filter.txt")

GMTiP_TR_info$chr <- gsub("_.+", "", GMTiP_TR_info$TRID)
GMTiP_TR_info_list <- split(GMTiP_TR_info, by="chr") 

GMTiP_TR_all_freq <- pbmcapply::pbmclapply(GMTiP_TR_info_list, function(f_sub_df){
  f_sub_df[, as.data.table(table(as.numeric(.SD))), by = TRID, .SDcols = -c("TRID", "chr")]
}, mc.cores = 22, mc.preschedule = FALSE)
GMTiP_TR_all_df <- rbindlist(GMTiP_TR_all_freq)
colnames(GMTiP_TR_all_df) <- c("TRID", "Copies", "af_GMTiP")
range(GMTiP_TR_all_df$af_GMTiP)

GMTiP_samplesize <- 130
GMTiP_TR_all_df$af_GMTiP <- GMTiP_TR_all_df$af_GMTiP/GMTiP_samplesize

unique_motif <- names(table(GMTiP_TR_all_df$TRID)[table(GMTiP_TR_all_df$TRID) == 1])
GMTiP_TR_all_df <- GMTiP_TR_all_df[!TRID%in%unique_motif]

fwrite(GMTiP_TR_all_df, "2025-06-11-specific_xQTL/input/GMTiP/TR/all_length_TR_freq.txt", sep = "\t")

# GMTiP_TR_all_df <- fread("2025-06-11-specific_xQTL/input/GMTiP/TR/all_length_TR_freq.txt")

GMTiP_TR_names <- unique(GMTiP_TR_all_df$TRID)
GMTiP_TR_names <- data.frame("ID"=GMTiP_TR_names)
GMTiP_TR_names <- tidyr::separate(GMTiP_TR_names, col = "ID", into = c("chrom", "start", "end", "motif"), sep = "_", remove = FALSE)
GMTiP_TR_names <- GMTiP_TR_names[, c("chrom", "start", "end", "ID", "motif")]

## 1KG TR data
dist_1KG_res <- fread("2025-06-11-specific_xQTL/input/1KG/TR/all_TR_dist.txt", sep="\t")
dist_1KG_res$End <- sapply(strsplit(dist_1KG_res$ID, split = "_"), function(x){x[3]})
# library(tidyverse)
# dist_1KG_res <- dist_1KG_res %>%
#   dplyr::mutate(
#     # Primary classification
#     population = dplyr::case_when(
#       (EAS_AFR > 2 & EAS_EUR > 2) ~ "S1",
#       (EAS_AFR > 2 | EAS_EUR > 2) ~ "S2",
#       TRUE ~ "-"
#     )
#   ) %>%
#   as.data.table()

# selected_1KG <- dist_1KG_res[population!="-"][, .(ID, Motif, population, EAS_AFR, EAS_EUR)]

TRID_1KG_res <- dist_1KG_res[, .(chrom=Chrom, start=Start, end=End, CHRSE=ID, motif=Motif)]

freq_1KG_res <- fread("2025-06-11-specific_xQTL/input/1KG/TR/all_length_TR_freq.txt")
only_one_size <- names(table(freq_1KG_res$ID)[table(freq_1KG_res$ID) == 1])

##
library(Biostrings)
GMTiP_1KG_overlap <- bedtoolsr::bt.intersect(GMTiP_TR_names, TRID_1KG_res, wo = TRUE, f = 0.8, F = 0.8)

motif_list <- unique(GMTiP_1KG_overlap$V5)
motif_rc_list <- sapply(motif_list, function(x){as.character(Biostrings::reverseComplement(DNAString(x)))})
GMTiP_1KG_overlap$GMTiP_motif_rc <- motif_rc_list[GMTiP_1KG_overlap$V5]
index1 <- GMTiP_1KG_overlap$V5 == GMTiP_1KG_overlap$V10 | GMTiP_1KG_overlap$GMTiP_motif_rc == GMTiP_1KG_overlap$V10
table(index1)
setDT(GMTiP_1KG_overlap)
GMTiP_1KG_overlap <- GMTiP_1KG_overlap[index1][, .(TRID=V4, KG1ID=V9, KG1motif=V10)]
GMTiP_1KG_overlap <- GMTiP_1KG_overlap[!KG1ID%in%only_one_size]

GMTiP_1KG_overlap1 <- merge(GMTiP_1KG_overlap, GMTiP_TR_all_df[, .(TRID, Copies, af_GMTiP)], by=c("TRID"))
freq_1KG_res$Copies <- as.numeric(freq_1KG_res$Copies)
GMTiP_1KG_overlap1$Copies <- as.numeric(GMTiP_1KG_overlap1$Copies)
GMTiP_1KG_overlap1 <- merge(GMTiP_1KG_overlap1, 
                            freq_1KG_res[, .(KG1ID=ID, KG1motif=Motif, 
                                             Copies, AFR, EAS, EUR)], by = c("KG1ID", "KG1motif", "Copies"), all = TRUE)
# GMTiP_1KG_overlap1 <- GMTiP_1KG_overlap1[!KG1ID%in%only_one_size]
GMTiP_1KG_overlap1$Chrom <- gsub("_.+", "", GMTiP_1KG_overlap1$KG1ID)
GMTiP_1KG_overlap1_list <- split(GMTiP_1KG_overlap1, by = "Chrom")

library(waddR)

safe_wasserstein <- function(x, y) {
  if (length(x) == 0 | length(y) == 0) return(-99)
  wasserstein_metric(x, y, p = 2)
  # return("1")
}

# AFR: 893; AMR: 490; EAS: 585, SAS: 601; EUR: 633; H3Africa: 348
dist_G1_list <- pbmcapply::pbmclapply(GMTiP_1KG_overlap1_list, function(f_sub_df){
  f_output <- f_sub_df[, {
    GMTiP_ID <- unique(na.omit(TRID))
    GMTiP_vec = rep(Copies[!is.na(af_GMTiP)], round(af_GMTiP[!is.na(af_GMTiP)] * GMTiP_samplesize))
    eas_vec = rep(Copies[!is.na(EAS)], round(EAS[!is.na(EAS)] * 585))
    afr_vec = rep(Copies[!is.na(AFR)], round(AFR[!is.na(AFR)] * 893))
    eur_vec = rep(Copies[!is.na(EUR)], round(EUR[!is.na(EUR)] * 633))
    
    max_GMTiP <- Copies[which.max(af_GMTiP)]
    af_GMTiP <- af_GMTiP[which.max(af_GMTiP)]
    max_EAS <- Copies[which.max(EAS)]
    af_EAS <- EAS[which.max(EAS)]
    
    .(
      TRID=GMTiP_ID,
      Nsize = length(Copies),
      max_GMTiP = max_GMTiP,
      af_GMTiP = af_GMTiP,
      max_EAS = max_EAS,
      af_EAS = af_EAS,
      GMTiP_AFR = safe_wasserstein(GMTiP_vec, afr_vec),
      GMTiP_EAS = safe_wasserstein(GMTiP_vec, eas_vec),
      GMTiP_EUR = safe_wasserstein(GMTiP_vec, eur_vec),
      KG1_EAS_EUR = safe_wasserstein(eas_vec, eur_vec)
    )
  }, by = .(KG1ID, KG1motif)]
  f_output
}, mc.cores = 22, mc.preschedule = FALSE)

dist_G1_res <- rbindlist(dist_G1_list)
dist_G1_res$source <- "1KG"

fwrite(dist_G1_res, "2025-06-11-specific_xQTL/input/GMTiP/TR/GMTiP_1KG_all_TR_dist.txt", sep = "\t")

# dist_G1_res <- fread("2025-06-11-specific_xQTL/input/GMTiP/TR/GMTiP_1KG_all_TR_dist.txt")
dist_G1_res <- dist_G1_res[TRID!=""]
dist_G1_res <- dist_G1_res[!is.na(TRID)]

fwrite(dist_G1_res, "2025-06-11-specific_xQTL/output/data_plot/TR_information.txt", sep = "\t")


sub_dist_G1_res <- merge(dist_G1_res, selected_1KG[, .(KG1ID=ID, KG1motif=Motif, population, EAS_AFR, EAS_EUR)], by=c("KG1ID", "KG1motif"))
fwrite(sub_dist_G1_res, "2025-06-11-specific_xQTL/output/data_plot/TR_EAS_specific.txt", sep = "\t")

boxplot(sub_dist_G1_res[, c("GMTiP_AFR", "GMTiP_EAS", "GMTiP_EUR")])


library(ggplot2)
library(data.table)
library(ggpubr)

#
plot_data <- melt(sub_dist_G1_res[, .(GMTiP_AFR, GMTiP_EAS, GMTiP_EUR)], 
                  measure.vars = c("GMTiP_AFR", "GMTiP_EAS", "GMTiP_EUR"),
                  variable.name = "Group",
                  value.name = "Value")

ggplot(plot_data, aes(x = Group, y = log2(Value+1), fill = Group)) +
  geom_jitter(width = 0.2, alpha = 0.4, color = "grey80") +
  geom_boxplot(width = 0.5, alpha = 0.7, outliers = FALSE) +
  labs(x = "",
       y = "log2(2-Wasserstein distance+1)") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  theme(text = element_text(family = "Arial"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  stat_compare_means(
    comparisons = list(c("GMTiP_AFR", "GMTiP_EAS"),
                       c("GMTiP_EAS", "GMTiP_EUR"),
                       c("GMTiP_AFR", "GMTiP_EUR")),
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.01
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

ggsave("2025-06-11-specific_xQTL/output/data_plot/plot/AF_dist_TR_GMTiP_1KG.png", width = 4, height = 4)


pdf("2025-06-11-specific_xQTL/output/data_plot/plot/Proportion_of_TR.pdf", width = 4, height = 4)
plot(eulerr::euler(list("Other TRs"=unique(dist_G1_res$TRID), 
                        "EAS_specific"=unique(sub_dist_G1_res$TRID)), 
                   shape = "circle"), quantities = TRUE, 
     main=paste0("Proportion: ", round(length(unique(sub_dist_G1_res$TRID))/
                                         length(unique(dist_G1_res$TRID))*100, 2), "%"))
dev.off()



## -----------------------------------------------------
library(ggplot2)
library(ggalluvial)
library(data.table)

df_SNP <- fread("2025-06-11-specific_xQTL/output/data_plot/SNP_EAS_specific.txt")
table(df_SNP$source);table(df_SNP$population)
df_SV <- fread("2025-06-11-specific_xQTL/output/data_plot/SV_EAS_specific.txt")
table(df_SNP$source);table(df_SNP$population)
df_TR <- fread("2025-06-11-specific_xQTL/output/data_plot/TR_EAS_specific.txt")
table(df_TR$population)


df_merge <- rbind(df_SNP[var_type=="SNP"][, type := "SNP"][,.(ID=rsid, source, type, population)],
                  df_SNP[var_type=="INDEL"][, type := "INDEL"][,.(ID=rsid, source, type, population)],
                  df_TR[, type := "TR"][,.(ID=TRID, source, type, population)],
                  df_SV[, type := "SV"][,.(ID, source, type, population)])

df_merge$new_type <- dplyr::case_when(
  df_merge$population %in% c("S1") ~ "EAS_specific_to_both",
  df_merge$population %in% c("S2") ~ "EAS_specific_to_either",
  df_merge$population %in% c("H1", "L1") ~ "EAS_diff_to_both",
  df_merge$population %in% c("H2", "L2") ~ "EAS_diff_to_either"
)

data.table::fwrite(df_merge, "2025-06-11-specific_xQTL/output/data_plot/Variant_type.txt", sep = "\t")


df_count <- df_merge %>% group_by(type, new_type) %>% 
  summarise(count=length(unique(ID))) %>% ungroup() %>% 
  group_by(type) %>% mutate(all=sum(count))

df_count <- df_count[order(factor(df_count$type, levels = rev(c("SNP", "INDEL", "TR", "SV")))),]

df_count$Freq <- df_count$count/df_count$all
df_count$type <- paste0(df_count$type, " (", df_count$all, ")")
df_count$type <- factor(df_count$type, levels = unique(df_count$type))

df_count$new_type <- factor(df_count$new_type, 
                            levels = rev(c("EAS_specific_to_both", "EAS_specific_to_either",
                                       "EAS_diff_to_both", "EAS_diff_to_either")))

ggplot(df_count, aes(type, Freq))+
  geom_bar(aes(fill=new_type), stat = "identity", position = "fill")+
  coord_polar(theta = "y")+
  geom_text(hjust = 1, size = 3,
            aes(x = type, y = 0, label = type, fontface="italic")) +
  
  theme_minimal() +
  labs(x="", y="", fill="Type") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "right") +
  scale_fill_manual(values = c("#94c2e1", "#71a260","#e8c147", "#fb7f00"))

ggsave("2025-06-11-specific_xQTL/output/data_plot/plot/Count_and_frequency.pdf", width = 5, height = 3)



##
df_plot <- to_lodes_form(df_merge[,2:ncol(df_merge)], id = "value")
unique(df_plot$stratum)

df_plot$stratum <- factor(df_plot$stratum, 
                          levels = c("1KG&gnomAD", "1KG", "gnomAD", 
                                     "SNP", "INDEL", "TR", "SV", 
                                     "S1", "S2", "H1", "H2", "L1", "L2"))

col<- c('#b9b9b9', '#b9b9b9', '#b9b9b9', '#cd5535', '#4061a7','#ccda4b', '#7f7f7f',
        '#913386', '#a071a9', '#488c70', '#5aad5a', '#f1c545', '#ffc719')

ggplot(df_plot, aes(x = x, fill=stratum, label=stratum,
                    stratum = stratum, alluvium  = value))+
  geom_flow(width = 0.1,
            curve_type = "sine",
            alpha = 0.5,
            color = 'white',
            size = 0.1)+
  geom_stratum(width = 0.15, color="white", size=0.2)+
  geom_text(stat = 'stratum', size = 2, color = 'black')+
  scale_fill_manual(values = col)+
  theme_void()+
  theme(legend.position = 'none')

ggsave("2025-06-11-specific_xQTL/output/data_plot/plot/Variant_sankey.pdf", width = 4, height = 4)


## Comparison of variant length
df_SNP[, c("ref") := tstrsplit(SNP, "_", keep = 3)]
df_TR[, c("start", "end") := tstrsplit(TRID, "_", keep = 2:3)]
df_SV[, c("length") := tstrsplit(ID, "_", keep = 5)]

summ_data <- data.frame(
  "type" = c("SNP", "TR", "SV"),
  "Count" = c(nrow(df_SNP), nrow(df_TR), nrow(df_SV)),
  "Length" = c(sum(nchar(df_SNP$ref)), 
               sum(as.numeric(df_TR$end)-as.numeric(df_TR$start)),
               sum(as.numeric(df_SV$length)))
)

summ_data <- reshape2::melt(summ_data)
summ_data$type <- factor(summ_data$type, levels = c("SV", "TR", "SNP"))

ggplot(data=summ_data, aes(x=variable, y=value, fill=type)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_classic() +
  scale_fill_viridis_d(direction = -1) +
  theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  labs(x = "", y = "Proportion")

ggsave("2025-06-11-specific_xQTL/output/data_plot/plot/Variant_length_comparsion.pdf", width = 3, height = 4)

