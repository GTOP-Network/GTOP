library(bambu)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
    stop("Usage: bambu.R <threads> <outdir> <bam> [ref_genome_fa] [ref_gtf] [bambu_annotations_rds]")
}

threads <- as.integer(args[1])
outdir <- args[2]
sample <- args[3]
refg <- args[4]
refgtf <- args[5]
refanno <- args[6]

cat(sprintf("####output directory: %s\n", outdir))
cat("####input files:\n")
print(sample)
cat(sprintf("####reference genome: %s\n", refg))
cat(sprintf("####reference gtf: %s\n", refgtf))
cat(sprintf("####bambu annotations: %s\n", refanno))

if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}

if (!file.exists(sample)) {
    stop(sprintf("Input BAM not found: %s", sample))
}
if (!file.exists(refg)) {
    stop(sprintf("Reference genome FASTA not found: %s", refg))
}

if (file.exists(refanno)) {
    bambuAnnotations <- readRDS(refanno)
} else {
    if (!file.exists(refgtf)) {
        stop(sprintf("Reference GTF not found and annotation cache does not exist: %s", refgtf))
    }
    refanno_dir <- dirname(refanno)
    if (!dir.exists(refanno_dir)) {
        dir.create(refanno_dir, recursive = TRUE)
    }
    bambuAnnotations <- prepareAnnotations(refgtf)
    saveRDS(bambuAnnotations, file = refanno)
}

# General Usage
se <- bambu(reads=sample, annotations=bambuAnnotations, genome=refg, verbose=TRUE, ncore=threads, lowMemory=TRUE, quant=TRUE)
saveRDS(se, file=file.path(outdir, 'bambu.rds'))

se <- readRDS(file.path(outdir, 'bambu.rds'))

writeToGTF(rowRanges(se), file.path(outdir, "discoveryOnly.ndr.default.gtf"))

counts <- assays(se)$counts
write.csv(
  data.frame(
    transcript_id = rownames(counts),
    read_count = counts[, 1]
  ),
  file = file.path(outdir, "transcript_counts.csv"),
  row.names = FALSE,
  quote = FALSE
)

counts <- assays(se)$fullLengthCounts
write.csv(
  data.frame(
    transcript_id = rownames(counts),
    read_count = counts[, 1]
  ),
  file = file.path(outdir, "full_length_support.csv"),
  row.names = FALSE,
  quote = FALSE
)

