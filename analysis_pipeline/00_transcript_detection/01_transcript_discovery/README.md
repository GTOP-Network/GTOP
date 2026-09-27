# Transcript discovery

Bambu, FLAIR, FLAMES, IsoQuant, Iso-Seq, IsoTools, and TALON independently identify transcript models from long-read RNA-seq data aligned to hg38.

## 1. Prepare reference objects

```bash
python bambu/bambu_run.py prepare
python isoquant/isoquant_run.py prepare
python talon/talon_run.py prepare
```

After preparation jobs finish:

```bash
python bambu/bambu_run.py prepare_check
python isoquant/isoquant_run.py prepare_check
python talon/talon_run.py prepare_check
```

FLAMES uses preprocessed reads from `FLAMES_FASTQ_DIR` and the configuration specified by `FLAMES_CONFIG`. The accompanying `exclude_IG_region_FASTQ.py` script prepares these reads from hg38 alignments.

## 2. Run transcript discovery

```bash
for tool in bambu flair flames isoquant isoseq isotools talon; do
  python "$tool/${tool}_run.py" run
done
```

## 3. Check and collect results

After discovery jobs finish, check the results before collecting outputs:

```bash
for tool in bambu flair flames isoquant isoseq isotools talon; do
  python "$tool/${tool}_run.py" check
done

for tool in bambu flair flames isoquant isoseq isotools talon; do
  python "$tool/${tool}_run.py" transfer
done
```

Collected annotations are stored in `output/assembly/LRS/sample_based/<sample>/hg38/<caller>/`. Complete collection before proceeding to merging.
