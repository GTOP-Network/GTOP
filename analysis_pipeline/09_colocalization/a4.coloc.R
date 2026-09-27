#!/usr/bin/env Rscript

# Load required libraries
library(data.table)
library(coloc)
library(dplyr)
library(pbmcapply)

# ------------------------------------------------------------
# Function: validate and parse command-line arguments
# Expected input order:
#   1. QTLTYPE (e.g., "SNV_eQTL")
#   2. QTLTISSUE (e.g., "Adipose")
#   3. QTLGENE (char, e.g., "ENSG00000117620.15")
#   4. GWASNAME (char, e.g., "Lymphocyte@count-Sakaue-2021")
#   5. SENTINELSNP (char, e.g., "rs77046277")
#   6. CHRNAME (char, e.g., "chr1")
#   7. GWASSTART (numeric, e.g., 100279090)
#   8. GWASEND (numeric, e.g., 102279090)
# ------------------------------------------------------------
argvs <- commandArgs(trailingOnly = TRUE)

gwas_info <- fread("input/EAS_GWAS_noMHC_all_sentinal_loci.txt")

gwas_qtl_pair <- fread(argvs[1], header = FALSE)
GWAS_LD_source <- as.character(argvs[2])
QTL_LD_source <- as.character(argvs[3])
p1_setting <- as.numeric(argvs[4])
p2_setting <- as.numeric(argvs[5])
p12_setting <- as.numeric(argvs[6])
thread_number <- as.numeric(argvs[7])
output_file <- argvs[8]

message(c("p1: ", p1_setting, "\tp2: ", p2_setting, "\tp12: ", p12_setting))

gwas_qtl_pair <- merge(
    gwas_qtl_pair,
    gwas_info[, .(
        V6 = GWAS_name,
        V9 = rsid,
        adjusted_loci_start,
        adjusted_loci_end
    )],
    by = c("V6", "V9")
)

merged_coloc_result <- pbmclapply(
    1:nrow(gwas_qtl_pair),
    function(row_index) {
        QTLTYPE <- as.character(gwas_qtl_pair$V1[row_index])
        QTLTISSUE <- gwas_qtl_pair$V2[row_index]
        QTLGENE <- gwas_qtl_pair$V3[row_index]
        GWASNAME <- gwas_qtl_pair$V6[row_index]
        CHRNAME <- gwas_qtl_pair$V7[row_index]
        # SENTINELPOS <- as.numeric(gwas_qtl_pair$V8[row_index])
        SENTINELSNP <- gwas_qtl_pair$V9[row_index]

        LOCISTART <- gwas_qtl_pair$adjusted_loci_start[row_index]
        LOCIEND <- gwas_qtl_pair$adjusted_loci_end[row_index]

        LOCINAME <- sprintf("%s:%s-%s", CHRNAME, LOCISTART, LOCIEND)

        GWASpath <- sprintf(
            "input/%s/%s_%s",
            GWAS_LD_source,
            GWASNAME,
            SENTINELSNP
        )
        QTLpath <- sprintf(
            "input/%s/%s/%s/%s",
            QTL_LD_source,
            QTLTYPE,
            QTLTISSUE,
            QTLGENE
        )

        GWAS_fm_type <- case_when(
            file.exists(paste0(GWASpath, ".noConverged.RData")) ~ "noConverged",
            file.exists(paste0(GWASpath, ".noCS.RData")) ~ "noCS",
            file.exists(paste0(GWASpath, ".finemapping.RData")) ~ "finemapping",
            file.exists(paste0(GWASpath, ".susie_ser.RData")) ~ "susie_ser"
        )

        ## QTL information
        QTL_fm_type <- case_when(
            file.exists(paste0(QTLpath, ".noConverged.RData")) ~ "noConverged",
            file.exists(paste0(QTLpath, ".noCS.RData")) ~ "noCS",
            file.exists(paste0(QTLpath, ".finemapping.RData")) ~ "finemapping"
        )

        load(paste0(GWASpath, ".", GWAS_fm_type, ".RData"))

       qtl_load_success <-  tryCatch(
            {
                load(paste0(QTLpath, ".", QTL_fm_type, ".RData"))
            },
            error = function(e) {
                cat("QTL load failed:", QTLpath, "\n")
                return(FALSE)
            }
        )

        if ( !isFALSE(qtl_load_success) ) {
            d2$varbeta[is.na(d2$varbeta)] <- 1
            d1$MAF[d1$MAF == 0] <- 0.00001
            d1$beta[is.na(d1$beta)] <- 0

            if (length(intersect(d1$snp, d2$snp)) >= 10) {
                ##
                simple_coloc_index <- FALSE
                if (
                    GWAS_fm_type %in%
                        c("noConverged", "noCS") |
                        QTL_fm_type %in% c("noConverged", "noCS")
                ) {
                    simple_coloc_index <- TRUE
                }
                if (GWAS_fm_type == "finemapping") {
                    if (is.null(s1$sets$cs)) {
                        simple_coloc_index <- TRUE
                    }
                }
                ## coloc
                if (simple_coloc_index) {
                    simple_coloc <- coloc.abf(
                        dataset1 = d1,
                        dataset2 = d2,
                        p1 = p1_setting,
                        p2 = p2_setting,
                        p12 = p12_setting
                    )$summary %>%
                        as.data.frame() %>%
                        t() %>%
                        as.data.frame()

                    GWAS_high_snp <- d1$snp[d1$pval == min(d1$pval)][1]
                    QTL_high_sentinel <- d2$snp[d2$pval == min(d2$pval)][1]

                    coloc_res <- simple_coloc %>%
                        mutate(
                            GWAS_name = GWASNAME,
                            loci = LOCINAME,
                            xQTL_type = QTLTYPE,
                            tissue = QTLTISSUE,
                            phenotype_id = QTLGENE,
                            GWAS_fm_type = GWAS_fm_type,
                            QTL_fm_type = QTL_fm_type,
                            coloc_type = "simple_coloc",
                            GWAS_sentinel = SENTINELSNP,
                            GWAS_hit = GWAS_high_snp,
                            GWAS_beta = d1$beta[d1$snp == GWAS_high_snp],
                            GWAS_se = sqrt(d1$varbeta[d1$snp == GWAS_high_snp]),
                            GWAS_pvalue = d1$pval[d1$snp == GWAS_high_snp],
                            xQTL_hit = QTL_high_sentinel,
                            xQTL_beta = d2$beta[d2$snp == QTL_high_sentinel],
                            xQTL_se = sqrt(d2$varbeta[
                                d2$snp == QTL_high_sentinel
                            ]),
                            xQTL_pvalue = d2$pval[d2$snp == QTL_high_sentinel]
                        ) %>%
                        select(
                            GWAS_name,
                            loci,
                            xQTL_type,
                            tissue,
                            phenotype_id,
                            GWAS_fm_type,
                            QTL_fm_type,
                            coloc_type,
                            GWAS_sentinel,
                            GWAS_hit,
                            GWAS_beta,
                            GWAS_se,
                            GWAS_pvalue,
                            xQTL_hit,
                            xQTL_beta,
                            xQTL_se,
                            xQTL_pvalue,
                            nsnps,
                            PP.H0.abf,
                            PP.H1.abf,
                            PP.H2.abf,
                            PP.H3.abf,
                            PP.H4.abf
                        )

                    return(coloc_res)
                } else {
                    ##
                    ## coloc-susie
                    susie_coloc <- coloc.susie(
                        dataset1 = s1,
                        dataset2 = s2,
                        p1 = p1_setting,
                        p2 = p2_setting,
                        p12 = p12_setting
                    )
                    susie_coloc$summary

                    susie_coloc_res <- susie_coloc$summary

                    if (is.na(susie_coloc_res$PP.H4.abf[1])) {
                        simple_coloc <- coloc.abf(
                            dataset1 = d1,
                            dataset2 = d2,
                            p1 = p1_setting,
                            p2 = p2_setting,
                            p12 = p12_setting
                        )$summary %>%
                            as.data.frame() %>%
                            t() %>%
                            as.data.frame()

                        GWAS_high_snp <- d1$snp[d1$pval == min(d1$pval)][1]
                        QTL_high_sentinel <- d2$snp[d2$pval == min(d2$pval)][1]

                        coloc_res <- simple_coloc %>%
                            mutate(
                                GWAS_name = GWASNAME,
                                loci = LOCINAME,
                                xQTL_type = QTLTYPE,
                                tissue = QTLTISSUE,
                                phenotype_id = QTLGENE,
                                GWAS_fm_type = GWAS_fm_type,
                                QTL_fm_type = QTL_fm_type,
                                coloc_type = "simple_coloc",
                                GWAS_sentinel = SENTINELSNP,
                                GWAS_hit = GWAS_high_snp,
                                GWAS_beta = d1$beta[d1$snp == GWAS_high_snp],
                                GWAS_se = sqrt(d1$varbeta[
                                    d1$snp == GWAS_high_snp
                                ]),
                                GWAS_pvalue = d1$pval[d1$snp == GWAS_high_snp],
                                xQTL_hit = QTL_high_sentinel,
                                xQTL_beta = d2$beta[
                                    d2$snp == QTL_high_sentinel
                                ],
                                xQTL_se = sqrt(d2$varbeta[
                                    d2$snp == QTL_high_sentinel
                                ]),
                                xQTL_pvalue = d2$pval[
                                    d2$snp == QTL_high_sentinel
                                ]
                            ) %>%
                            select(
                                GWAS_name,
                                loci,
                                xQTL_type,
                                tissue,
                                phenotype_id,
                                GWAS_fm_type,
                                QTL_fm_type,
                                coloc_type,
                                GWAS_sentinel,
                                GWAS_hit,
                                GWAS_beta,
                                GWAS_se,
                                GWAS_pvalue,
                                xQTL_hit,
                                xQTL_beta,
                                xQTL_se,
                                xQTL_pvalue,
                                nsnps,
                                PP.H0.abf,
                                PP.H1.abf,
                                PP.H2.abf,
                                PP.H3.abf,
                                PP.H4.abf
                            )
                        return(coloc_res)
                    } else {
                        coloc_res <- susie_coloc_res %>%
                            mutate(
                                GWAS_name = GWASNAME,
                                loci = LOCINAME,
                                xQTL_type = QTLTYPE,
                                tissue = QTLTISSUE,
                                phenotype_id = QTLGENE,
                                GWAS_fm_type = GWAS_fm_type,
                                QTL_fm_type = QTL_fm_type,
                                coloc_type = "susie_coloc",
                                GWAS_sentinel = SENTINELSNP,
                                GWAS_hit = hit1,
                                GWAS_beta = d1$beta[match(hit1, d1$snp)],
                                GWAS_se = sqrt(d1$varbeta[match(hit1, d1$snp)]),
                                GWAS_pvalue = d1$pval[match(hit1, d1$snp)],
                                xQTL_hit = hit2,
                                xQTL_beta = d2$beta[match(hit2, d2$snp)],
                                xQTL_se = sqrt(d2$varbeta[match(hit2, d2$snp)]),
                                xQTL_pvalue = d2$pval[match(hit2, d2$snp)]
                            ) %>%
                            select(
                                GWAS_name,
                                loci,
                                xQTL_type,
                                tissue,
                                phenotype_id,
                                GWAS_fm_type,
                                QTL_fm_type,
                                coloc_type,
                                GWAS_sentinel,
                                GWAS_hit,
                                GWAS_beta,
                                GWAS_se,
                                GWAS_pvalue,
                                xQTL_hit,
                                xQTL_beta,
                                xQTL_se,
                                xQTL_pvalue,
                                nsnps,
                                PP.H0.abf,
                                PP.H1.abf,
                                PP.H2.abf,
                                PP.H3.abf,
                                PP.H4.abf
                            )
                        return(coloc_res)
                    }
                }
            }
        }
    },
    mc.cores = thread_number,
    mc.preschedule = FALSE
)
merged_coloc_result <- rbindlist(merged_coloc_result)

fwrite(merged_coloc_result, file = output_file, sep = "\t")
