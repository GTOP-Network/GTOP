library(bambu)

args <- commandArgs(trailingOnly = TRUE)
threads <- args[1]
outdir <- args[2]
samples <- args[3:length(args)]  # args[-1]
cat(sprintf("####output directory: %s\n", outdir))
cat("####input files:\n")
print(samples)

workdir <- '/flashfs1/scratch.global/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/bambu'
refg <- file.path(workdir, 'ref/genome.fa')
refanno <- file.path(workdir, 'ref2/bambuAnnotations.rds')

if (!dir.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
}

bambuAnnotations <- readRDS(refanno)
# Transcript discovery only (no quantification)
se <- bambu(reads=samples, annotations=bambuAnnotations, genome=refg, verbose=TRUE, ncore=threads, lowMemory=TRUE, NDR=1, quant=FALSE)
se.filtered <- se[(!is.na(mcols(se)$NDR) & mcols(se)$NDR<0.5) | is.na(mcols(se)$NDR)]
writeToGTF(se.filtered, file.path(outdir, "discoveryOnly.ndr0.5.gtf"))
