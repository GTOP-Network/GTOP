library(bambu)

args <- commandArgs(trailingOnly = TRUE)
refgtf <- args[1]
refanno <- args[2]

cat(sprintf("####reference gtf: %s\n", refgtf))
cat(sprintf("####bambu annotations: %s\n", refanno))

if (!file.exists(refgtf)) {
    stop(sprintf("Reference GTF not found: %s", refgtf))
}

outdir <- dirname(refanno)
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
}

bambuAnnotations <- prepareAnnotations(refgtf)
saveRDS(bambuAnnotations, file=refanno)
