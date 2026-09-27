## cis-QTL and GWAS Colocalization

This directory contains a series of R scripts for **fine-mapping and colocalization analyses of genome-wide association study (GWAS) and cis-QTL signals**. The pipeline integrates GWAS and QTL summary statistics to prioritize putative causal variants and evaluate whether GWAS and molecular QTL associations share the same underlying causal signals.

## Script Overview

### `a1.finemapping_for_GWAS.R`

Performs fine-mapping of GWAS summary statistics using **SuSiE**, with 1000 Genomes Project (1KGP) East Asian (EAS) genotypes as the linkage disequilibrium (LD) reference panel. This step identifies credible sets and estimates posterior inclusion probabilities (PIPs) for candidate causal variants within GWAS-associated loci.

### `a2.prepare_genes.R`

Identifies genes overlapping GWAS-associated loci and prepares the corresponding gene–locus pairs for downstream QTL fine-mapping and colocalization analyses.

### `a3.finemapping_for_QTL_GTOP_LD.R`

Performs fine-mapping of cis-QTL signals, including eQTLs and sQTLs, to identify credible sets and prioritize candidate causal variants underlying gene expression and splicing associations.

### `a4.coloc.R`

Performs colocalization analysis between GWAS and cis-QTL signals using **SuSiE-coloc**, evaluating whether the two association signals are consistent with sharing the same underlying causal variant.

### `a5.merge_result.R`

Integrates colocalization results across GWAS loci and cis-QTL datasets into a unified summary table, facilitating the interpretation and prioritization of shared genetic signals between complex traits and molecular phenotypes.

## Workflow Summary

`GWAS Fine-Mapping` → `Gene–Locus Preparation` → `QTL Fine-Mapping` → `Colocalization Analysis` → `Result Integration`

Corresponding scripts:

`a1.finemapping_for_GWAS.R` → `a2.prepare_genes.R` → `a3.finemapping_for_QTL_GTOP_LD.R` → `a4.coloc.R` → `a5.merge_result.R`
