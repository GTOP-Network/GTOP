# GTOP Splicing Quantification

This directory contains code to replicate the splicing quantification procedure done as part of the initial GTOP publication. This procedure is broadly split into three sections:
1. Variant-aware read mapping with STAR and filtering reads
2. Annotation-agnostic estimation of intron excision with [Leafcutter2](https://github.com/cfbuenabadn/leafcutter2/)
3. Filtering and normalization for sQTL analsyses

*Note: these steps should be run in order*<br><br>

## Variant-aware read mapping with STAR

The `01_align_dev.sh` script uses VCF files to perform variant-aware alignment with STAR and keep WASP pass reads. The script assumes that `bcftools`, `htslib`, `vcftools`, and `STAR` are in the user's path. The script requires four arguments:

1. `workDir`: Working directory where outputs will be generated
2. `VCFdir`: Directory containing per doner VCF files
3. `starIndex`: Path to the STAR index directory
4. `fastqDir`: Path to directory containing FASTQ files, organized into subdirectories based on sequencing batch



## Splicing quantification with Leafcutter2

The `02_intron_usage.sh` script will use the `STAR` alignments generated in the previous section to estimate intron cluster usage per sample. Following, intron usage estimates are aggregated across all samples from a tissue typeand consolidated into a single file. The focal outputs from this script are `leafcutter2.junction_counts.gz` and `leafcutter2.cluster_ratios.gz` which contain ***counts*** of intron cluster usage and ***ratios***, respectively.<br>

The script takes six arguments:
1. `STAR_DIR`: The `STAR` alignments generated in the previous section
2. `Meta_data`: The second positional argument passed to the script, which specifies a metadata file containing two columns: Column 1: "tissue"; Column 2: "tissuecode", without header
3. `juncDir`: Directory containing junction files
4. `intronDir`: Directory to write output intron files
5. `GTFfile`: The GTF file
6. `reference`: The reference genome fasta file



*Before running this script, install `leafcutter2` following the instructions [here](https://github.com/cfbuenabadn/leafcutter2).*<br><br>

## Filtering and normalization for sQTL analyses

The `03_map_and_phenotype.sh` script will filter out lowly expressed and low-complexity splicing clusters. It will output three files:
1. Raw filtered intron excision ratios (no normalization)
2. Normalized filtered intron excision ratios for sQTL calling
3. A file mapping introns to genes

The script takes nine arguments:
1. `GTFfile`: The GTF file, used for mapping introns to genes
2. `leafcutterOutDir`: The directory containing leafcutter2 files from the previous step
3. `leafcutterOutPrefix`: Prefix of Leafcutter2 output files (i.e. everything before `.junction_counts.gz` and `.cluster_ratios.gz`)
4. `sampleLookupFile`: A file mapping columns of the leafcutter files to sampleIDs (to be included in the output file). Should be a two-column tab-separated file
5. `leafcutterScriptDir`: Path to leafcutter foldr, which is used for normalization (and some filtering)
6. `exonFile`: This gzipped TSV file contains annotation information for exons, with each row corresponding to a single exon entry. The file includes the following columns:chr: Chromosome name; start: 0-based genomic start position of the exon; end: 0-based exclusive genomic end position of the exon; strand: Strand of the exon on the chromosome, indicated by "+"; gene_id: Unique identifier for the gene associated with the exon (e.g., Ensembl gene ID: ENSG00000123456); gene_name: Common/gene symbol name of the associated gene (e.g., TP53, BRCA1); transcript_id: Unique identifier for the transcript that the exon belongs to (e.g., Ensembl transcript ID: ENST00000345678), linking the exon to its parent transcript.
7. `geneinfoFile`: This TSV file contains annotation information for genes, with each row corresponding to a single gene. The file includes the following columns:chr: gene_name: Common/gene symbol name of the associated gene (e.g., TP53, BRCA1);gene_type: biotype of gene; gene_id: Unique identifier for the gene associated with the exon (e.g., Ensembl gene ID: ENSG00000123456); Chromosome name; start: 1-based genomic start position of the exon; end: 1-based inclusive genomic end position of the exon; strand: Strand of the exon on the chromosome, indicated by "+".
7. `outDir`: Path to directory to write output files to
<br><br>

*Before running this script, install `leafcutter` following the instructions [here](https://davidaknowles.github.io/leafcutter/articles/Installation.html).*<br><br>