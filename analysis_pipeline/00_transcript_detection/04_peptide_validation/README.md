# Proteomic validation of long-read RNA-seq transcripts

We used data-independent acquisition (DIA) proteomics data matched to long-read RNA sequencing (RNA-seq) samples to evaluate the translation potential of detected transcripts.  
For each tissue, we constructed protein sequence databases derived from predicted coding sequences (CDS) of highly expressed transcripts (median transcripts per million, TPM > 5) predicted to be protein-coding. Mass spectrometry data from each tissue were then searched against the corresponding protein database using DIA-NN (v1.8.1).

---

```bash
cd scripts
```

## 1. Construction of protein sequence databases from transcript coding sequences (CDS) for each tissue

We first selected transcripts with median TPM > 5 in each tissue and retained those predicted to be protein-coding. Protein-coding potential was predicted using TransDecoder2 (generated during the Structural and Quality Annotation of Novel Transcript Isoforms, SQANTI3, annotation step). Based on these transcripts, we constructed tissue-specific reference protein sequence databases.

In addition to the TPM > 5 threshold, we also constructed alternative protein databases using transcripts with TPM > 0.1 to include lower-abundance transcripts. While this increases sensitivity for detecting low-expression translation events, it also expands the search space.

```bash
python prepare_faa.py
```

## 2. Database search using DIA-NN

Mass spectrometry data from each tissue were searched separately against the corresponding protein database using DIA-NN (v1.8.1).

```bash
qsub run.diann.sh 
```

---
