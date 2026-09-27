# GTOP Figure Generation — Detailed Documentation

## 1. Purpose

The `figure_generation/` directory separates visualization from upstream statistical analysis. Each manuscript figure has a dedicated directory containing the main plotting script and the processed inputs needed for figure generation.

This structure allows a reader to trace a figure from the manuscript to its plotting code and then to the processed analysis results.

## 2. Figure-to-script map

| Manuscript figure | Script |
|---|---|
| Figure 1 | `Figure_1/Figure1.R` |
| Figure 2 | `Figure_2/Figure2.R` |
| Figure 3 | `Figure_3/Figure3.R` |
| Figure 4 | `Figure_4/Figure4.R` |
| Figure 5 | `Figure_5/Figure5.R` |
| Figure 6 | `Figure_6/Figure6.R` |
| Extended Figure 1 | `Extended_Figures/Extended_Figure1.R` |
| Extended Figure 2 | `Extended_Figures/Extended_Figure2.R` |
| Extended Figure 3 | `Extended_Figures/Extended_Figure3.R` |
| Extended Figure 4 | `Extended_Figures/Extended_Figure4.R` |
| Extended Figure 5 | `Extended_Figures/Extended_Figure5.R` |
| Extended Figure 6 | `Extended_Figures/Extended_Figure6.R` |
| Extended Figure 7 | `Extended_Figures/Extended_Figure7.R` |
| Extended Figure 8 | `Extended_Figures/Extended_Figure8.R` |
| Extended Figure 9 | `Extended_Figures/Extended_Figure9.R` |
| Supplementary Figures | `Supp_Figures/Supp_Figure_*.R` |

## 3. Main figure content

### Figure 1 — GTOP cohort and genetic variation

The script combines genotype PCA, RNA-expression MDS, tissue clustering, genome-wide variant summaries, per-sample variant counts, SV concordance, and TR detection comparisons.

### Figure 2 — Long-read transcriptome and allele-specific regulation

The script summarizes transcript-level sample structure, transcript discovery relative to GENCODE/GTEx references, SQANTI3 classifications, alternative splicing, tissue distribution, peptide support, WGCNA cardiac modules, ASE/ASTS detection, and a representative allele-specific locus.

### Figure 3 — Molecular QTL landscape

The script summarizes QTL counts by variant class, genomic distance to TSS, GTOP–GTEx effect-size concordance, eGene overlap among SNV/SV/TR QTLs, SV-length effects, disease-associated TRs, splice-junction regulation, and tissue sharing.

### Figure 4 — Population-specific regulatory effects

The script combines:

- comparison of credible-set size across GTOP, GTEx, and cross-ancestry fine-mapping;
- maximum-PIP summaries;
- a representative fine-mapped locus;
- geographic allele-frequency patterns;
- GWAS/QTL colocalization;
- GTOP–GTEx allele-frequency comparison;
- eQTL portability summaries.

### Figure 5 — Structural variation and tandem repeats in fine-mapping

The script evaluates:

- LD between SV/TR eQTLs and their best-tagging SNVs;
- poorly tagged SV/TR signals;
- variant composition of joint credible sets;
- tissue-level credible-set counts;
- credible-set classes containing SVs and/or TRs;
- enrichment of SV/TR variants among fine-mapped signals.

### Figure 6 — Disease genetics

The script integrates:

- S-LDSC enrichment;
- GWAS loci colocalized with eQTLs and sQTLs;
- comparison with external resources;
- representative locus-level colocalization;
- GWAS enrichment;
- QTL overlap with risk variants across frequency groups;
- SV/TR-containing colocalization signals;
- representative SV-sQTL/sQTL loci.

## 4. Extended figures

The Extended Figure scripts provide supporting analyses for the major conclusions of the manuscript.

They cover the complete sequence from variant detection and transcript discovery through ASE/ASTS, QTL characterization, fine-mapping, ancestry-dependent effects, and disease-associated regulatory variation.

## 5. Supplementary figures

The supplementary scripts are organized around technical validation and additional analyses.

Important groups include:

### Data quality and variant characterization

- WGS quality statistics.
- Short- and long-read RNA quality.
- ADMIXTURE ancestry components.
- Variant detection comparison.
- VEP annotation.

### Transcriptome characterization

- Consensus transcript discovery.
- Transcript filtering.
- Peptide support.
- Predicted protein structures.
- Cross-platform expression consistency.
- Iso-Seq sampling-depth analysis.
- Co-expression modules.

### Allele-specific analyses

- Genotype phasing.
- ASE and ASTS.
- Allele-specific splicing.
- Cross-tissue consistency.

### QTL analyses

- PCA-based covariate selection.
- sQTL consistency.
- SV/TR QTL tagging and representation.
- Genomic feature enrichment.
- Cross-tissue and cross-resource eQTL comparisons.
- eQTL portability.

### Fine-mapping

- SV/TR representation in credible sets.
- Coverage-filtered joint fine-mapping.
- Robustness of SV/TR prioritization.
- Representative SV fine-mapping loci.

### Disease association

- Heritability enrichment.
- GWAS colocalization.
- Characterization of colocalized/non-colocalized loci.
- Frequency-differentiated QTLs.

## 6. Input organization

Inputs are kept close to the plotting scripts:

```text
Figure_1/
├── Figure1.R
├── README-FIG.1.md
└── input/

Figure_2/
├── Figure2.R
├── README-FIG.2.md
└── input/

...
```

This organization is intentional: a figure can be reproduced without searching through the entire repository for its input tables.

## 7. Data flow

```text
analysis_pipeline/
        │
        │ processed analysis results
        ▼
figure-specific input/
        │
        ▼
Figure-specific R script
        │
        ▼
manuscript figure / panel
```

The plotting layer therefore depends on processed outputs rather than raw sequencing data.

## 8. Environment and paths

Several scripts contain study-specific working-directory settings or absolute paths. Before execution, replace these paths with locations appropriate to the local environment.

For example:

```r
setwd("/path/to/GTOP/figure_generation/Figure_4")
```

The scripts also rely on R packages that may need to be installed separately.

## 9. Recommended reproducibility workflow

For a single figure:

```text
Identify figure
    ↓
Open figure directory
    ↓
Check README-FIG.X.md
    ↓
Inspect input/
    ↓
Update paths
    ↓
Install required R packages
    ↓
Run FigureX.R
```

For the complete manuscript:

```text
Upstream analysis
       ↓
Processed result tables / R objects
       ↓
Figure-specific input directories
       ↓
Main figures
       +
Extended figures
       +
Supplementary figures
```

## 10. Important distinction from the analysis pipeline

The figure-generation scripts are **visualization workflows**, not the primary statistical analysis workflows.

The statistical analyses and generation of molecular QTLs, fine-mapping results, enrichment statistics, and colocalization results are documented under:

```text
analysis_pipeline/
```

The figure scripts consume those processed results and convert them into manuscript-ready visualizations.

## 11. Reuse

Individual figure scripts can be adapted for related datasets provided that the input objects retain the expected column names and data structures. For substantial reuse, inspect the relevant input-loading section of the script before substituting datasets.

## 12. Citation

Please cite the GTOP study when reusing these figure-generation workflows or results, and cite the original software packages used for the corresponding analyses.
