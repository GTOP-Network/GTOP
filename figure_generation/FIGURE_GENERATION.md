# GTOP Figure Generation — Detailed Documentation

## 1. Purpose

The `figure_generation/` directory contains the visualization layer for the GTOP Phase-I manuscript. It is separated from `analysis_pipeline/`, which contains the primary statistical and computational analyses.

## 2. Figure-to-script map

### Main figures

| Figure | Script |
|---|---|
| Figure 1 | `Figure_1/Figure1.R` |
| Figure 2 | `Figure_2/Figure2.R` |
| Figure 3 | `Figure_3/Figure3.R` |
| Figure 4 | `Figure_4/Figure4.R` |
| Figure 5 | `Figure_5/Figure5.R` |
| Figure 6 | `Figure_6/Figure6.R` |

### Extended figures

| Figure | Script |
|---|---|
| Extended Figure 1 | `Extended_Figures/Extended_Figure1.R` |
| Extended Figure 2 | `Extended_Figures/Extended_Figure2.R` |
| Extended Figure 3 | `Extended_Figures/Extended_Figure3.R` |
| Extended Figure 4 | `Extended_Figures/Extended_Figure4.R` |
| Extended Figure 5 | `Extended_Figures/Extended_Figure5.R` |
| Extended Figure 6 | `Extended_Figures/Extended_Figure6.R` |
| Extended Figure 7 | `Extended_Figures/Extended_Figure7.R` |
| Extended Figure 8 | `Extended_Figures/Extended_Figure8.R` |
| Extended Figure 9 | `Extended_Figures/Extended_Figure9.R` |

## 3. Main figure content

### Figure 1 — GTOP cohort and genetic variation

`Figure_1/Figure1.R` generates panels covering genotype PCA, RNA-seq multidimensional scaling, tissue clustering, variant number/length distributions, per-genome variant counts, SV comparisons, and LRS/SRS TR comparisons.

### Figure 2 — Long-read transcriptome

`Figure_2/Figure2.R` covers LRS sample structure, transcript discovery and novelty, SQANTI3 structural/coding annotation, alternative splicing, transcript tissue breadth, peptide support, transcript-level co-expression/GO analysis, and ASE/ASTS summaries.

### Figure 3 — Molecular QTLs

`Figure_3/Figure3.R` covers eQTL/sQTL summary statistics, TSS distance of fine-mapped QTLs, GTOP–GTEx effect-size comparison, eGene sharing across variant classes, SV effect size versus SV length, pathogenic TR-associated QTLs, representative splice-junction signals, and MASH-based cross-tissue analyses.

### Figure 4 — Population-specific QTLs

`Figure_4/Figure4.R` covers fine-mapping credible-set size and PIP, frequency-differentiated QTLs, representative frequency-differentiated loci, and GTOP–GTEx eQTL portability.

### Figure 5 — SV/TR fine-mapping

`Figure_5/Figure5.R` evaluates SV/TR tagging by small variants, poorly tagged signals, and the composition and enrichment of joint fine-mapping credible sets.

### Figure 6 — Disease association

`Figure_6/Figure6.R` integrates S-LDSC enrichment, GWAS/QTL colocalization, comparisons with external QTL resources, GTOP-specific signals, GWAS enrichment, rare-variant analyses, disease-associated SV/TR signals, and representative SV-associated loci.

## 4. Extended figures

The extended-figure scripts provide supporting analyses for the major conclusions:

- **Extended Figure 1:** genetic variation and pathogenic TR-related analyses.
- **Extended Figure 2:** additional long-read/short-read transcriptome characterization.
- **Extended Figure 3:** long-read ASE and ASTS analyses.
- **Extended Figure 4:** molecular QTL characteristics and genomic/functional enrichment.
- **Extended Figure 5:** QTL tissue sharing and tissue specificity.
- **Extended Figure 6:** cross-ancestry fine-mapping and MPRA-related analyses.
- **Extended Figure 7:** frequency-differentiated sQTLs and mash-based eQTL portability.
- **Extended Figure 8:** SV-sQTL and TR-sQTL signals poorly tagged by small variants.
- **Extended Figure 9:** SMR and GWAS/QTL colocalization analyses.

## 5. Supplementary figures

The repository currently distributes scripts for S2–S42, with the actual set of scripts listed in `figure_generation/README.md`. Major groups include sequencing/genetic QC (S2–S8), transcript discovery and characterization (S11–S17), allele-specific analyses (S18–S21), QTL phenotype/discovery analyses (S22–S31), fine-mapping/portability (S34–S39), and disease association (S40–S42).

## 6. Input organization

Inputs are stored within the relevant figure directories or in shared figure input directories:

```text
Figure_X/input/
Extended_Figures/input/
Supp_Figures/input/
```

The inputs are processed analysis outputs rather than raw sequencing data and include `.txt`, `.gz`, `.RDS`, `.RData`, and related formats.

## 7. Reproducibility workflow

```text
Upstream statistical analysis
            ↓
Processed result tables / objects
            ↓
Figure-specific input/
            ↓
Figure-specific R script
            ↓
Publication figure
```

For an individual figure:

1. Open the corresponding figure directory.
2. Inspect the input files loaded by the R script.
3. Update local paths if required.
4. Install the packages imported by the script.
5. Run the R script.

Example:

```bash
cd /path/to/GTOP/figure_generation/Figure_3
Rscript Figure3.R
```

## 8. Important distinction

The figure-generation scripts are **visualization workflows**. They do not replace the upstream analyses that generate QTLs, fine-mapping results, enrichment statistics, or colocalization results. Those analyses are documented under `analysis_pipeline/`.

## 9. Reuse

Individual plotting scripts can be adapted to related datasets when the expected input structures are retained. For substantial reuse, inspect the input-loading and transformation sections of the relevant script before substituting datasets.

## 10. Citation

Please cite the GTOP study when reusing these figure-generation workflows or derived results, together with the original software packages used for the corresponding analyses.
