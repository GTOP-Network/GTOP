library(data.table)
library(tidyverse)

files <- list.files("./SNP_ASE", pattern="*all*", full.names=T)
df <- data.frame(ase_readcount_file=files)
df2 <- df %>% mutate(sample_id=str_split( basename(ase_readcount_file) , "\\.", simplify=T)[, 1])

tissueinfo <- fread("GTOP_tissue_code_and_colors.csv", colClasses=rep("character", 4))

df2 <- df2 %>% mutate(Tissue_Code=str_split(sample_id, "-", simplify=T)[, 3])%>% left_join(tissueinfo[, c("Tissue", "Tissue_Code")])

final <- df2 %>% dplyr::select(sample_id, tissue_site_detail=Tissue, ase_readcount_file)

write_tsv(final, "./input/SNP_ASE.files.list")
