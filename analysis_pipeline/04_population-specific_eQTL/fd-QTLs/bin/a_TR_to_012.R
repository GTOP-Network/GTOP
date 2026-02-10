library(dplyr)
library(matrixStats)


TR_info <- fread("/media/london_B/lixing/2024-08-29-TR-AsianGTEX/2024-12-11-TR-LongReads/output/QTL_Filter/GMTip_LRS_TR_131INDs.Miss85.AF95.AC1.dosage.impute.simpleTR.filter.txt")

mat <- as.matrix(TR_info %>% column_to_rownames("TRID"))

q1 <- rowQuantiles(mat, probs = 1/3, na.rm = TRUE)
q2 <- rowQuantiles(mat, probs = 2/3, na.rm = TRUE)

# 向量化分类
result <- matrix(1, nrow = nrow(mat), ncol = ncol(mat))
result[mat < q1] <- 0
result[mat > q2] <- 2

# 转换回数据框
df_categorized <- as.data.frame(result)
colnames(df_categorized) <- colnames(mat)
df_categorized$TRID <- rownames(mat)
df_categorized <- df_categorized %>% select(TRID, everything())

fwrite(df_categorized, "2025-06-11-specific_xQTL/input/GMTiP/TR/GMTip_LRS_TR_131INDs.Miss85.AF95.AC1.dosage.impute.simpleTR.filter.012.txt", sep = "\t")

rm(mat, q1, q2, result)
gc()
