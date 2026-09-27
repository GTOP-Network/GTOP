# eQTL Portability Analysis Pipeline

This directory contains a suite of scripts for **eQTL portability and cross-tissue sharing analyses**. The pipeline covers data preprocessing, multivariate adaptive shrinkage (mash) analysis, and the evaluation and summarization of eQTL portability across populations and tissues.

The scripts are organized into five sequential analytical stages (`a1`–`a5`), reflecting the logical workflow of the analysis.

## Script Overview

### `a1.GTEx_data.R`

Identifies GTOP eGene–variant pairs that are also present in GTEx and prepares the overlapping QTL data for downstream portability analyses.

### `a2.mash_adjustment.sh`

Applies the **multivariate adaptive shrinkage (mash)** framework to model QTL effects across multiple tissues and calculates local false sign rates (lfsr) for downstream analyses.

### `a3.prepare_eqtl_portability.R`

Prepares the processed QTL data and evaluates eQTL portability using six complementary metrics.

### `a4.eqtl_portability.R`

Summarizes eQTL portability across the six metrics and generates the final portability results.

### `a5.mash_analysis.R`

Evaluates the portability and cross-tissue sharing patterns of sQTLs and eGenes based on mash-adjusted effect estimates and statistical significance.

## Usage Recommendations

Run the scripts sequentially to follow the standard analytical workflow:

`a1` → `a2` → `a3` → `a4` → `a5`

Required R packages include, but are not limited to:

* `data.table`
* `dplyr`
* `ggplot2`
* `mashr`
* `qvalue`

Please ensure that all required dependencies are installed before running the pipeline.

Scripts involving mash modeling and large-scale cross-tissue analyses can be computationally intensive. Running these steps on a high-performance computing (HPC) system is recommended when possible.
