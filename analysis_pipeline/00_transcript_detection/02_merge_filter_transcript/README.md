# Transcript merging and annotation

Caller outputs are filtered by read support, merged within and across methods, and annotated with SQANTI3. Structural and expression evidence is used to construct the final transcript references.

```bash
cd scripts
```

## 1. Prepare transcript models

```bash
python 01_filter_tx_by_support_FLNC.py
python 02_raw_transcript_qc.py
python 03_split_gtf_by_chr_pipeline.py
python 04_trans_gtf2bed.py
```

Raw transcript QC summarizes per-sample transcript counts before and after filtering.

## 2. Merge within callers

```bash
python 05_tama_run_multi_submit.py submit-config
# After configuration jobs finish:
python 05_tama_run_multi_submit.py submit-merge
# After tissue-level jobs finish:
python 06_tama_run_tissue.py
# After cross-tissue jobs finish:
python 07_merge_bed2gtf.py
```

## 3. Integrate and annotate transcripts

Run these commands in order, waiting for submitted jobs to complete between steps:

```bash
python 08_merge_pipeline.py
python 09_sqanti3_run.py run_sqanti3
python 09_sqanti3_run.py merge
python 09_sqanti3_run.py filter_step1
python 10_gtf_process.py
```

## 4. Evaluate transcript support

Quantify candidate transcripts and combine the completed sample results:

```bash
python 11_flair_quant_run.py flair_quant
# After quantification jobs finish:
python 11_flair_quant_run.py combine
```

Prepare and submit first-exon and short-read junction checks, then merge each set of completed results:

```bash
bash 12_alt_first_exon_run.sh prepare
bash 12_alt_first_exon_run.sh submit
# After scan jobs finish:
bash 12_alt_first_exon_run.sh merge

bash 13_junction_srs_run.sh prepare
bash 13_junction_srs_run.sh submit
# After junction jobs finish:
bash 13_junction_srs_run.sh merge
```

## 5. Construct final references

```bash
python 09_sqanti3_run.py custom_filter_gt_1k
# After filtering finishes:
python 14_enhanced_gtf_run.py enhanced_gtf_gt_1k
```

Reference construction writes the GTOP and enhanced GENCODE transcript references, associated annotation tables, and predicted protein sequences to `release/` for subsequent analysis.
