#!/usr/bin/env bash
# Copy to config.sh, edit absolute paths, then source config.sh before every stage.
export GTOP_PROJECT_DIR=/path/to/GTOP-RNA/20260815
export GTOP_RUN_DIR="$GTOP_PROJECT_DIR/scratch"
export GTOP_REF_DIR="$GTOP_PROJECT_DIR/reference"
export HG38_FASTA="$GTOP_REF_DIR/hg38.fa"
export HG38_GTF="$GTOP_REF_DIR/gencode.v47.annotation.gtf"
export HG38_MMI="$GTOP_REF_DIR/hg38.isoseq.mmi"
export FLNC_DIR="$GTOP_PROJECT_DIR/input/flnc"
export FLAMES_FASTQ_DIR="$GTOP_PROJECT_DIR/input/flames_non_ig"
export TISSUE_META="$GTOP_PROJECT_DIR/input/tissue_code.csv"
export SRS_JUNCTION_DIR="$GTOP_PROJECT_DIR/input/srs_junctions"
export FLAMES_ROOT=/lustre/home/lhgong/2026-02-26-mengxin/2026-05-12-flames/FLAMES
# Server-side configuration from the supplied revision; override if relocated.
export FLAMES_CONFIG="$FLAMES_ROOT/gtop_LRS_RNA_config.json"
export SQANTI3_DIR=/path/to/SQANTI3
export TAMA_DIR=/path/to/tama
export BUILD_LOCI=/path/to/buildLoci/buildLoci.pl
export GENCODE_PROTEINS="$GTOP_REF_DIR/gencode.v47.pc_translations.fa.gz"
export MS_RAW_DIR="$GTOP_PROJECT_DIR/input/ms_raw"
export SLURM_PARTITION=cpuPartition
# Commands default to tools on PATH. Set environment activation commands as needed.
# export LOAD_BASE_ENV_CMD='source /path/to/miniforge/etc/profile.d/conda.sh'
# export LOAD_POLARS_ENVS_CMD='conda activate pipeline'
# export LOAD_FLAIR_ENVS_CMD='source /path/to/conda.sh && conda activate flair'
# export LOAD_ISOSEQ_ENVS_CMD='source /path/to/conda.sh && conda activate isoseq'
# export LOAD_TALON_ENV_CMD='source /path/to/conda.sh && conda activate talon'
# export LOAD_ISOTOOLS_ENV_CMD='source /path/to/conda.sh && conda activate isotools'
# export LOAD_SQANTI3_ENVS_CMD='source /path/to/conda.sh && conda activate sqanti3'
# export LOAD_TAMA_ENVS_CMD='source /path/to/conda.sh && conda activate tama_py2'
# export FLAMES_ENV='source /path/to/conda.sh && conda activate flames'
# export BAMBU_RSCRIPT=/path/to/bambu/bin/Rscript
# export ISOQUANT_PYTHON=/path/to/isoquant/bin/python
# export ISOQUANT_SCRIPT=/path/to/isoquant.py
# export ISOTOOLS_PYTHON=/path/to/isotools/bin/python
export DIANN_IMAGE=/path/to/diann.sif
export DIANN_BIN=/path/inside/container/diann-1.8.1
export DIANN_RUN_TAG=run
export DIANN_RESULTS="$GTOP_PROJECT_DIR/output/MS/diann/output/diann/tissues/$DIANN_RUN_TAG"
