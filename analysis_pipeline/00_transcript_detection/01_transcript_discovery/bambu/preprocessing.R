library(bambu)

args <- commandArgs(trailingOnly = TRUE)
threads <- args[1]
outdir <- args[2]
sample <- args[3]
cat(sprintf("####output directory: %s\n", outdir))
cat("####input file:\n")
print(sample)

workdir <- '/flashfs1/scratch.global/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/bambu'
refg <- file.path(workdir, 'ref/genome.fa')
refanno <- file.path(workdir, 'ref2/bambuAnnotations.rds')

if (!dir.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
}

bambuAnnotations <- readRDS(refanno)
se <- bambu(reads=sample, rcOutDir=outdir, annotations=bambuAnnotations, genome=refg, verbose=TRUE, ncore=threads)