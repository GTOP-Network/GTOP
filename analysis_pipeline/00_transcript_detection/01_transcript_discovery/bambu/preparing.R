library(bambu)

workdir <- '/flashfs1/scratch.global/lhgong/longrw/myprojs/gtop/20251108-LR-RNAseq/bambu'
refgtf <- file.path(workdir, 'ref/gencode.v47.annotation.gtf')
outdir <- file.path(workdir, 'ref2')

if (!dir.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
}

bambuAnnotations <- prepareAnnotations(refgtf)
saveRDS(bambuAnnotations, file=file.path(outdir, 'bambuAnnotations.rds'))