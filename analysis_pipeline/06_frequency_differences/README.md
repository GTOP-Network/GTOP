# Frequency-Differentiated QTL Analysis Pipeline

This directory contains a suite of R scripts for comprehensive **frequency-differentiated quantitative trait locus (fdQTL) analysis**. The pipeline covers allele frequency comparison across ancestral populations, QTL fine-mapping and credible-set processing, cross-tissue locus aggregation, and the identification of frequency-differentiated eQTLs and sQTLs based on fine-mapping posterior inclusion probability (PIP) values.

The scripts are organized into four sequential analytical stages (`a1`–`a4`), reflecting the logical workflow of the analysis.

## Script Overview

### `a1.gtop_gnomad_af_comparison.R`

Computes and compares allele frequencies of genetic variants across different ancestral populations, providing the population-frequency information required for downstream ancestry-aware QTL analyses.

### `a2.qtl_fine_mapping_summary.R`

Processes QTL fine-mapping results by merging credible-set variants with variants in linkage disequilibrium (LD; R² > 0.01), generating an expanded summary of fine-mapped QTL signals.

### `a3.qtl_mege_to_loci.R`

Aggregates and merges credible sets across different tissues to define QTL loci and consolidate shared or overlapping fine-mapping signals.

### `a4.fdQTL_table.R`

Identifies frequency-differentiated quantitative trait loci (fdQTLs) by integrating population allele frequency information with fine-mapped QTL signals.

## Usage Recommendations

Run the scripts sequentially to follow the standard analytical workflow:

`a1` → `a2` → `a3` → `a4`

Required R packages include, but are not limited to:

* `data.table`
* `dplyr`
* `ggplot2`
* `mashr`
* `qvalue`

Please ensure that all required dependencies are installed before running the pipeline.
