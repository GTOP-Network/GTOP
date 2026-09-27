# Proteomic validation

Tissue-specific protein databases combine GENCODE translations with predicted coding sequences from the transcript annotation workflow. DIA-NN searches matched mass-spectrometry data against these databases.

## 1. Prepare protein databases

Complete transcript quantification and configure the tissue metadata and mass-spectrometry input paths.

```bash
cd scripts
python prepare_faa.py
```

## 2. Run DIA-NN

```bash
bash run_diann.sh
```

## 3. Export protein abundance

After all DIA-NN jobs finish:

```bash
python protein_abundance.py
```

Tissue protein databases are stored in `release/cds/tissue_based/`; abundance matrices are written to `release/molec_pheno/transcript_raw_protein/`. Peptide-level evidence is available in the DIA-NN reports.
