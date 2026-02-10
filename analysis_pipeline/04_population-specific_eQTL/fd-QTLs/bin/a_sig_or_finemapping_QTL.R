#!/usr/bin/env Rscript

setwd("/media/bora_A/zhangt/2025-05-07-EAS_specific_xQTL-Project")

# Load required libraries
library(data.table)

clean_anno <- fread("/media/pacific/share/Datasets/Asian_GTEx/Analysis/annotation/gencode.v47.gene_info.txt") %>% as.data.frame()
rownames(clean_anno) <- clean_anno$V5

SNP_posrsid <- "2025-06-11-specific_xQTL/input/GMTiP/SNP/GMTiP_SNP_chrpos_rsid.txt"
CHRPOS_RSID <- fread(SNP_posrsid, header = FALSE)
colnames(CHRPOS_RSID) <- c("SNP", "rsid")

transcript_events <- fread("/media/bora_A/wangyn/gencode/novel_gencode47.with_novel_gene_all_exons.txt.gz", header = T)
transcript_events <- transcript_events[, .(chr, start, end, Name=transcript_id, gene_name=gene_id, strand, gene_symbol=gene_name)]
transcript_events <- transcript_events %>% 
  group_by(chr, Name, gene_name, strand, gene_symbol) %>% 
  summarise(start=min(start), end=max(end)) %>% as.data.frame()
rownames(transcript_events) <- transcript_events$Name


## lead and significant xQTLs ------------------------------------------------
QTL_info <- pbmcapply::pbmclapply(c("SNV_eQTL", "TR_eQTL", "SV_eQTL", 
                                    "SNV_juQTL", "TR_juQTL", "SV_juQTL", 
                                    "SNV_tuQTL", "TR_tuQTL", "SV_tuQTL"), function(xQTL_type){ 
  
  sig_variant_files <- list.files(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/%s/nominal_slim/nom_thresh/", xQTL_type), full.names = T)
  sig_variants <- rbindlist(lapply(sig_variant_files, function(x){
    tmp_df <- fread(x)
    tmp_df$tissue <- sub(".nom_thresh.txt.gz", "", basename(x))
    tmp_df
  }))
  sig_variants$sig_type <- "sig"
  sig_variants$variant_type <- gsub("_.+", "", xQTL_type)
  sig_variants$QTL_type <- gsub(".+_", "", xQTL_type)
  sig_variants <- sig_variants[,.(tissue, variant_id=SNP, phenotype_id=gene, gene_symbol=symbol, 
                                  beta, se, p_value, variant_type, QTL_type, sig_type)]
  
  lead_variant_files <- list.files(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/%s/clean_xGene", xQTL_type), full.names = T)
  lead_variants <- rbindlist(lapply(lead_variant_files, function(x){
    if(file.size(x) > 0){
      tmp_df <- fread(x)
      tmp_df$tissue <- sub(".txt", "", basename(x))
      # tmp_df
      tmp_df$sig_type <- "xGene"
      
      all_gene_info <- fread(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/%s/permutation/%s.gz", xQTL_type, basename(x)))
      all_gene_info <- all_gene_info[!all_gene_info$phenotype_id%in%tmp_df$xGene, 
                                     .(xGene=phenotype_id, lead_variant=variant_id, beta=slope, 
                                       p_value=pval_nominal, se=slope_se, num_var, pval_nominal_threshold, 
                                       gene_name=NA, gene_symbol=NA, tissue=sub(".txt", "", basename(x)), sig_type="Non-xGene")]
      
      rbind(tmp_df, all_gene_info)
    }
  }))
  lead_variants <- lead_variants[sig_type == "xGene"]
  tmp_id1 <- paste0(lead_variants$tissue, "_", lead_variants$lead_variant, "_", lead_variants$xGene)
  tmp_id2 <- paste0(sig_variants$tissue, "_", sig_variants$variant_id, "_", sig_variants$phenotype_id)
  
  sig_variants$sig_type[tmp_id2%in%tmp_id1] <- "lead"
  
  sig_variants
}, mc.cores = 6, mc.preschedule = FALSE)

QTL_info_df <- rbindlist(QTL_info)
fwrite(QTL_info_df, "2025-06-11-specific_xQTL/output/data_plot/data/tensorQTL_signif_QTL.txt", sep="\t")



## significant ------------------------------------------------
threshold_info <- pbmcapply::pbmclapply(c("SNV_eQTL", "TR_eQTL", "SV_eQTL", 
                                    "SNV_juQTL", "TR_juQTL", "SV_juQTL", 
                                    "SNV_tuQTL", "TR_tuQTL", "SV_tuQTL"), 
                                  function(xQTL_type){ 
                                      variant_files <- list.files(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/%s/permutation/", xQTL_type), full.names = T)
                                      sig_variants <- rbindlist(lapply(variant_files, function(x){
                                        if(file.size(x)>0){
                                          tmp_df <- fread(x)
                                          if(any(tmp_df$qval < 0.05)){
                                            tmp_df$tissue <- sub(".txt.gz", "", basename(x))
                                            tmp_df[qval < 0.05, .(tissue, phenotype_id, variant_id, pval_nominal, qval, pval_nominal_threshold)]
                                          }
                                        }
                                      }))
                                      sig_variants$variant_type <- gsub("_.+", "", xQTL_type)
                                      sig_variants$QTL_type <- gsub(".+_", "", xQTL_type)
                                      sig_variants
                                    }, mc.cores = 6, mc.preschedule = FALSE)

threshold_info_df <- rbindlist(threshold_info)
fwrite(threshold_info_df, "2025-06-11-specific_xQTL/output/data_plot/data/tensorQTL_QTL_threshold.txt", sep="\t")



QTL_info_df$gene_name[!QTL_info_df$xQTL_type%in%c("SNV_tuQTL", "SV_tuQTL", "TR_tuQTL")] <- 
  gsub(".+:", "", QTL_info_df$phenotype_id[!QTL_info_df$xQTL_type%in%c("SNV_tuQTL", "SV_tuQTL", "TR_tuQTL")])
QTL_info_df$gene_name[QTL_info_df$xQTL_type%in%c("SNV_tuQTL", "SV_tuQTL", "TR_tuQTL")] <- 
  transcript_events[QTL_info_df$phenotype_id[QTL_info_df$xQTL_type%in%c("SNV_tuQTL", "SV_tuQTL", "TR_tuQTL")], "gene_name"]
QTL_info_df$gene_symbol <- clean_anno[QTL_info_df$gene_name, "V6"]
QTL_info_df$gene_type <- clean_anno[QTL_info_df$gene_name, "V7"]

QTL_info_df <- QTL_info_df[, .(variant_id, phenotype_id, gene_name, gene_symbol, gene_type, 
                               beta, se, p_value, tissue, xQTL_type, variant_type, sig_type)]

fwrite(QTL_info_df, "2025-06-11-specific_xQTL/output/data_plot/data/SNP_lead_xQTL_summary.txt", 
       row.names = F, sep = "\t")



## fine-mapping -------------------------------------------------------------
## fine mapping SNV-xQTL results
susie_res <- lapply(c("SNV_eQTL", "SNV_juQTL", "SNV_tuQTL"), function(x){
  tissue_list <- list.files(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/%s/susie/", x), pattern = ".txt")
  t_df <- lapply(tissue_list, function(y){
    t_data <- fread(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/%s/susie/%s", x, y))
    t_data$tissue <- sub(".txt", "", y)
    t_data
  })
  t_df <- do.call('rbind', t_df)
  t_df$xQTL_type <- x
  t_df
})

susie_res <- do.call('rbind', susie_res)
susie_res <- susie_res %>% group_by(phenotype_id, tissue, xQTL_type, cs_id) %>% 
  mutate(cs_size=length(unique(variant_id))) %>% as.data.table()

susie_res$gene_name[susie_res$xQTL_type%in%c("SNV_eQTL", "SNV_juQTL")] <- gsub(".+:", "", susie_res$phenotype_id[susie_res$xQTL_type%in%c("SNV_eQTL", "SNV_juQTL")])
susie_res$gene_name[susie_res$xQTL_type%in%c("SNV_tuQTL")] <- transcript_events[susie_res$phenotype_id[susie_res$xQTL_type%in%c("SNV_tuQTL")], "gene_name"]
susie_res$variant_type <- "SNV"
susie_res$gene_symbol <- clean_anno[susie_res$gene_name, "V6"]
susie_res$gene_type <- clean_anno[susie_res$gene_name, "V7"]

ggvenn::ggvenn(list("A"=susie_res$phenotype_id[susie_res$xQTL_type=="SNV_eQTL"],
                    "B"=QTL_info_df$phenotype_id[QTL_info_df$xQTL_type=="SNV_eQTL"]))

combine_fm_tensorqtl <- function(susie_input, qtl_input){
  qtl_input <- qtl_input[qtl_input$xQTL_type%in%susie_input$xQTL_type, ]
  qtl_input_type <- paste0(qtl_input$phenotype_id, "_", qtl_input$tissue, "_", qtl_input$xQTL_type)
  susie_input_type <- paste0(susie_input$phenotype_id, "_", susie_input$tissue, "_", susie_input$xQTL_type)
  qtl_subdf <- qtl_input[!qtl_input_type%in%susie_input_type, ]
  
  rbind(susie_input[, .(variant_id, phenotype_id, gene_name, gene_symbol, gene_type, 
                        pip, cs_id, cs_size, tissue, xQTL_type, variant_type, sig_type="xGene")], 
        qtl_subdf[, .(variant_id, phenotype_id, gene_name, gene_symbol, gene_type, 
                      pip=-1, cs_id=NA, cs_size=NA, tissue, xQTL_type, variant_type, sig_type)])
}

merged_qtl <- combine_fm_tensorqtl(susie_input = susie_res, qtl_input = QTL_info_df)
length(unique(merged_qtl$phenotype_id)) # 83727

fwrite(merged_qtl, "2025-06-11-specific_xQTL/output/data_plot/data/SNV_variant_xQTL_finemapping.txt", sep = "\t")



## fine mapping joint SV-SNV-xQTL results
susie_res <- lapply(c("SNV_eQTL", "SNV_juQTL", "SNV_tuQTL"), function(x){
  tissue_list <- list.files(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/%s/combined_SNV_SV_susie/", x), pattern = ".txt")
  t_df <- lapply(tissue_list, function(y){
    t_data <- fread(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/%s/combined_SNV_SV_susie/%s", x, y))
    t_data$tissue <- sub(".txt", "", y)
    t_data
  })
  t_df <- do.call('rbind', t_df)
  t_df$xQTL_type <- gsub(".+_", "", x)
  t_df
})

susie_res <- do.call('rbind', susie_res)
susie_res <- susie_res %>% group_by(phenotype_id, tissue, xQTL_type, cs_id) %>% 
  mutate(cs_size=length(unique(variant_id))) %>% as.data.table()
susie_res$xGene[susie_res$xQTL_type%in%c("eQTL", "juQTL")] <- gsub(".+:", "", susie_res$phenotype_id[susie_res$xQTL_type%in%c("eQTL", "juQTL")])
susie_res$xGene[susie_res$xQTL_type%in%c("tuQTL")] <- transcript_events[susie_res$phenotype_id[susie_res$xQTL_type%in%c("tuQTL")], "gene_name"]
susie_clean_res <- merge(susie_res, clean_anno[, .(xGene=V5, symbol=V6, gene_type=V7)], by="xGene")
susie_clean_res$Type <- "SV"
susie_clean_res$Type[susie_clean_res$variant_id%in%CHRPOS_RSID$rsid] <- "SNV"

fwrite(susie_clean_res, "2025-06-11-specific_xQTL/output/data_plot/data/SNV_SV_variant_xQTL_finemapping.txt", sep = "\t")

## finemapping type
finemapping_SV_SNV_raw <- fread("2025-06-11-specific_xQTL/output/data_plot/data/SNV_SV_variant_xQTL_finemapping.txt")
finemapping_SV_SNV_raw <- finemapping_SV_SNV_raw[, .(phenotype_id, variant_id, xGene, symbol, gene_type, 
                                                     cs_id, cs_size, pip, Type, tissue, xQTL_type)]

finemapping_SV_index <- finemapping_SV_SNV_raw %>%
  group_by(phenotype_id, tissue, cs_id, xQTL_type) %>%
  summarise(SV_pip= max(pip[Type=="SV"]),
            SNV_pip = max(pip[Type=="SNV"]),
            SV_count= sum(Type=="SV"),
            SNV_count = sum(Type=="SNV"))
setDT(finemapping_SV_index)

finemapping_SV_index$pip_Type <- case_when(
  finemapping_SV_index$SV_count == 0 ~ "only_SNV",
  finemapping_SV_index$SNV_count == 0 ~ "only_SV",
  finemapping_SV_index$SV_pip > finemapping_SV_index$SNV_pip ~ "lead_SV",
  finemapping_SV_index$SV_pip < finemapping_SV_index$SNV_pip ~ "lead_SNV",
  finemapping_SV_index$SV_pip == finemapping_SV_index$SNV_pip ~ "equal_SNV_SV"
)
table(finemapping_SV_index$pip_Type)
# finemapping_SV_index <- finemapping_SV_index[count > 0.1]

finemapping_SV_SNV <- merge(finemapping_SV_SNV_raw, finemapping_SV_index[, .(phenotype_id, tissue, cs_id, xQTL_type, pip_Type)], 
                            by=c("phenotype_id", "tissue", "cs_id", "xQTL_type"))
# finemapping_SV <- merge(finemapping_SV, sigQTL_data[, .(phenotype_id=xGene, variant_id=SNP, tissue)], 
#                         by=c("phenotype_id", "variant_id", "tissue"))
# length(unique(finemapping_SV$phenotype_id))
# finemapping_SV$xQTL_type[finemapping_SV$xQTL_type%in%c("juQTL", "tuQTL")] <- "sQTL"

fwrite(finemapping_SV_SNV, "2025-06-11-specific_xQTL/output/data_plot/data/SNV_SV_variant_xQTL_finemapping_add_type.txt", sep = "\t")


## fine mapping TR-xQTL results  
susie_res <- lapply(c("SNV_eQTL", "SNV_juQTL", "SNV_tuQTL"), function(x){
  tissue_list <- list.files(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/%s/combined_SNV_TR_susie/", x), pattern = ".txt")
  t_df <- lapply(tissue_list, function(y){
    t_data <- fread(sprintf("2025-06-11-specific_xQTL/data/tensorqtl/%s/combined_SNV_TR_susie/%s", x, y))
    t_data$tissue <- sub(".txt", "", y)
    t_data
  })
  t_df <- do.call('rbind', t_df)
  t_df$xQTL_type <- gsub(".+_", "", x)
  t_df
})

susie_res <- do.call('rbind', susie_res)
susie_res <- susie_res %>% group_by(phenotype_id, tissue, xQTL_type) %>% 
  mutate(cs_size=length(unique(variant_id))) %>% as.data.table()
susie_res$xGene[susie_res$xQTL_type%in%c("eQTL", "juQTL")] <- gsub(".+:", "", susie_res$phenotype_id[susie_res$xQTL_type%in%c("eQTL", "juQTL")])
susie_res$xGene[susie_res$xQTL_type%in%c("tuQTL")] <- transcript_events[susie_res$phenotype_id[susie_res$xQTL_type%in%c("tuQTL")], "gene_name"]

susie_clean_res <- merge(susie_res, clean_anno[, .(xGene=V5, symbol=V6, gene_type=V7)], by="xGene")

fwrite(susie_clean_res, "2025-06-11-specific_xQTL/output/data_plot/data/SNV_TR_variant_xQTL_finemapping.txt", sep = "\t")
