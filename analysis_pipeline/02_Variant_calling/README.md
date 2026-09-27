# Variant calling and downstream genomic analysis

To generate a high-confidence genome-wide variant dataset, we implemented a comprehensive variant calling and annotation workflow using both long-read sequencing (LRS) and short-read sequencing (SRS) whole-genome sequencing data. The workflow consists of four major components: (1) small variant calling from LRS WGS data, (2) small variant calling from SRS WGS data, (3) population structure analysis using PCA and ADMIXTURE, and (4) functional annotation of variants using Ensembl Variant Effect Predictor (VEP).

All analyses were performed using recommended parameters unless otherwise specified.

---

## 1. Small variant calling from long-read WGS data

Long-read WGS small variants, including single-nucleotide variants (SNVs) and small insertions/deletions (indels), were identified using a long-read-specific variant calling workflow.

Aligned PacBio HiFi reads were used as input for small variant calling. The workflow performs variant detection, genotyping, and downstream filtering to generate high-confidence LRS-derived SNV and indel call sets.

The output includes:

* Sample-level variant calls
* Joint genotyped VCF files
* Filtered high-confidence SNV/indel datasets


```bash
bash 01_LRS_WGS_Variant_Calling.sh
```

---

## 2.  Small variant calling from short-read WGS data


Short-read WGS data were processed independently to generate SRS-derived small variant call sets for comparison and integration with LRS-derived variants.

```bash
bash 02_SRS_WGS_Variant_Calling.sh
```

Illumina short-read aligned BAM files were used for SNV and indel calling. The workflow includes variant calling, genotype refinement, quality filtering, and generation of final SRS-derived variant datasets.

The output includes:

* Sample-level variant calls
* Joint genotyped VCF files
* Filtered high-confidence SNV/indel datasets

---

## 3. Population structure analysis using PCA and ADMIXTURE

Genome-wide genetic variation was used to evaluate population structure and ancestry composition among samples.

```bash
bash 03_PCA_ADMIXTURE.sh
```

---

## 4. Functional annotation of variants using VEP

Variant consequences and functional annotations were obtained using Ensembl Variant Effect Predictor (VEP).

```bash
bash 04_vep_annotation.sh
```
---

