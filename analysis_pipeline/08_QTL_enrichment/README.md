# GTOP QTL Enrichment Analysis

This directory contains code for the enrichment analyses:

1. **Genomic annotation enrichment**: Bayesian enrichment analysis using [torus](https://github.com/xqwen/torus).
    
2. **S-LDSC enrichment analysis**: Stratified LD score regression for heritability enrichment using [LDSC](https://github.com/bulik/ldsc).
  

---

## Part 1: Genomic annotation enrichment

The `run_torus.sh` script performs Bayesian enrichment analysis using [torus](https://github.com/xqwen/torus). The pipeline consists of two main steps, each implemented as a function within the master script.


|Function|Description|
|---|---|
|`prepare`|Processes tensorQTL nominal results to generate torus input format.|
|`torus`|Executes torus in estimation mode (`-est`) to calculate enrichment parameters for each annotation category.|

---

## Part 2: S-LDSC enrichment analysis

The `run_sldsc_enrichment.sh` script performs stratified LD score regression using  [LDSC](https://github.com/bulik/ldsc) to estimate heritability enrichment of genomic annotations. The pipeline consists of three main steps.


|Function|Description|
|---|---|
|`generate_snp_f`|Extracts sets of variants for LDSC annotation: (1) significant variants; (2) variants in credible sets from fine-mapping results; (3) variants in credible sets from fine-mapping results along with their maximum PIP posterior inclusion probabilities.|
|`ldsc_annot_f`|Formats the extracted variant sets for LD score calculation in the subsequent step.|
|`make_annot_f`|Calculates LD scores for each annotation using `ldsc.py --l2`.|
|`enrichment_f`|Performs S-LDSC enrichment analysis using `ldsc.py --h2` to estimate enrichment.|

