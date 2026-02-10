# Transcript integration, annotation, filtering, and reference construction

To generate a unified and high-confidence transcript reference, transcript models identified by the three discovery pipelines were integrated based on intron-chain identity. Transcripts supported by multiple methods were prioritized, followed by annotation, quantification, and stringent filtering to produce the final long-read transcript reference used in downstream analyses.

---

```bash
cd script
```

## 1. Transcript integration and candidate transcript generation

```bash
merge_run.sh
```

Transcript models identified by the three methods were merged based on intron-chain identity. Transcripts sharing identical ordered splice junctions were grouped into a single intron-chain–defined model. For intron chains with heterogeneous 5′ or 3′ boundaries, the transcript spanning the most distal coordinates was selected as the representative model. Each intron chain was assigned a priority score corresponding to the number of independent methods supporting it. Transcripts supported by at least two methods were retained, and mono-exonic transcripts were excluded to generate a unified candidate transcript set.

---

## 2. Transcript-level read quantification using FLAIR

To improve computational efficiency, transcript quantification was performed in parallel at the sample level, followed by aggregation across samples.

### Sample-level quantification

```bash
python flair_quant_run.py flair_quant
```

### Merging quantification results

```bash
python flair_quant_run.py combine
```

---

## 3. Transcript annotation with SQANTI3

SQANTI3 was used to annotate candidate transcripts with respect to known gene models and splice junction features. To accelerate processing, annotation was performed in parallel at the chromosome level and subsequently merged.

### Chromosome-level annotation

```bash
python sqanti3_run.py run_sqanti3
```

### Merging SQANTI3 annotations

```bash
python sqanti3_run.py merge
```

---

## 4. Filtering of candidate transcripts

```bash
python sqanti3_run.py custom_filter
```

Candidate transcripts were filtered based on SQANTI3 annotation categories and transcript-level read support. Only transcripts meeting predefined structural and expression criteria were retained as high-confidence models.

---

## 5. Construction of enhanced transcript references

```bash
python enhanced_gtf_run.py enhanced_gtf
```

Two transcript reference annotations were generated: (i) a long-read–derived transcript reference (GTOP), and (ii) an enhanced reference combining GTOP novel transcripts with GENCODE v47 annotations. These references were used for downstream quantification and integrative analyses.
