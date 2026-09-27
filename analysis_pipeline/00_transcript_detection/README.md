# Transcript discovery and validation

The workflow consists of four steps:

1. [Transcript discovery](01_transcript_discovery/README.md) using seven complementary tools.
2. [Merging and annotation](02_merge_filter_transcript/README.md) to construct the GTOP transcript references.
3. [Quantification](03_quantification/README.md) using FLAIR.
4. [Proteomic validation](04_peptide_validation/README.md) using tissue-specific protein databases and DIA-NN.

## Setup

Edit the server paths and software environments in `config.example.sh`, then load the configuration before running each step:

```bash
cp config.example.sh config.sh
source config.sh
```

The workflow uses Python, NumPy, pandas, Polars, and the bioinformatics tools invoked by each step. Jobs are submitted with SLURM. Complete each stage before starting the next.

Inputs include PacBio FLNC reads, indexed hg38 alignments, GENCODE v47 annotations, tissue metadata, short-read splice-junction evidence, and DIA mass-spectrometry data. Reference files and alignments must use consistent chromosome names. Input and output locations are configured in the scripts and shell configuration.
