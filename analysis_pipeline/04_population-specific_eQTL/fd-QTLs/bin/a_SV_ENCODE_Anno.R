
library(data.table)

encode_data <- fread("/media/bora_A/zhangt/src/data/ENCODE_cCREs_regulation_UCSC") #42126344

encode_cleandata <- encode_data[, .(`#Chrom`=`#chrom`, Start=chromStart, End=chromEnd, cCRE=encodeLabel)]
encode_cleandata$`#Chrom` <- gsub("chr", "", encode_cleandata$`#Chrom`)
fwrite(encode_cleandata, "../src/bin/AnnotSV/share/AnnotSV/Annotations_Human/Users/GRCh38/AnyOverlap/ENCODE_cCREs.sorted.bed", sep = "\t")
