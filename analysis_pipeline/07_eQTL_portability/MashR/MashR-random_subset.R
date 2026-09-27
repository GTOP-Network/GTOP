#!/BIGDATA2/scau_hzhang_1/software/microsoft-r-open/microsoft-r-open-MRO4.0.2-X11/bin/Rscript

options(scipen = 999)

library(data.table)
"%&%" = function(a,b) paste0(a,b)

ARGS <- commandArgs(trailingOnly = TRUE)
dir_output = ARGS[1]
subset_size = as.numeric(ARGS[2])

DF_all_pairs = NULL

tmp = fread(dir_output %&% "/nominal_pairs.MashR_input.txt.gz",sep="\t",nThread=20)
colnames(tmp) = gsub(".nominal_pairs_zval","",colnames(tmp))
DF_all_pairs = tmp

if (nrow(DF_all_pairs) < subset_size){
    print("# of rows: " %&% nrow(DF_all_pairs))
    random_subset = DF_all_pairs
} else {
    set.seed(9823)
    random_subset = DF_all_pairs[floor(seq(1,nrow(DF_all_pairs),length.out=subset_size)),]
}
saveRDS(random_subset, dir_output %&% paste0("/MashR.random_subset_", subset_size, ".RDS"))

cat("Done\n")

