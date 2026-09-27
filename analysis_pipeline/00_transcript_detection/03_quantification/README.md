# Transcript and gene quantification

FLAIR quantifies FLNC reads against the final GTOP transcript reference. Complete reference construction before running this step.

```bash
cd scripts
python flair_quant_run.py flair_quant
```

After all sample jobs finish:

```bash
python flair_quant_run.py check
python flair_quant_run.py combine
```

The combined transcript and gene count/TPM matrices are written to `release/LRS_quant/GTOP/`.
