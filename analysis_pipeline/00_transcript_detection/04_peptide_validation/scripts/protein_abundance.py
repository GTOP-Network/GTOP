# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/7/30 20:57
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_DIANN_RESULTS = os.environ.get('DIANN_RESULTS', str(Path(CONFIG_PROJECT_DIR) / 'output/MS/diann/output/diann/tissues/run'))

import logging
import pandas as pd


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

PROJ_DIR = os.environ.get("PROJECT_ROOT")

def generate_tissue_protein_abundance():
    """Generate per-tissue protein abundance matrix.
    """
    dia_nn_dir=CONFIG_DIANN_RESULTS
    main_dir=f'{CONFIG_PROJECT_DIR}'
    out_dir = f'{main_dir}/release/molec_pheno/transcript_raw_protein'
    os.makedirs(out_dir, exist_ok=True)
    for tissue in os.listdir(dia_nn_dir):
        tissue_path =f'{dia_nn_dir}/{tissue}/{tissue}.pg_matrix.tsv'
        if not os.path.exists(tissue_path):
            print(f'warning: tissue {tissue} does not exist')
            continue
        df=pd.read_csv(tissue_path, sep='\t')
        df.index=df['Protein.Group']
        intensity_cols = [c for c in df.columns if c.lower().endswith(('.mzml', '.raw', '.d'))]
        if not intensity_cols:
            raise ValueError(f'No recognizable MS sample columns in {tissue_path}')
        df=df[intensity_cols]
        df.columns=df.columns.map(lambda x: os.path.basename(x).split('.')[0])
        df = df[~df.index.astype(str).str.contains(";")]
        df.to_csv(f'{out_dir}/{tissue}.protein_abundance.tsv.gz', sep='\t')
        print(f'{tissue} done, {df.shape[0]} proteins x {df.shape[1]} samples')


if __name__ == '__main__':
    generate_tissue_protein_abundance()
