# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/10/6 01:42
@Desc    : Description
"""
import logging
import os
import re
import subprocess
import sys
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd

import pandas as pd

from util import PROJ_DIR, run_commands_threadpool

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

STAR_BIN_DIR='/lustre/home/cxue/software/STAR-2.7.3a/bin'
RSEM_BIN_DIR='/lustre/home/cxue/software/RSEM-1.3.3'

# ref_gtf_prefix = {'enhanced': f'{PROJ_ROOT_DIR}/project/GMTiP-RNA/20251031/long_read/HPC/output/enhanced_gtf/GTOP_novel-GENCODE_v47',
#             'gencode': f'{PROJ_ROOT_DIR}/raw_data/GMTiP/ref/LRS/gencode.v47.annotation'}

ref_gtf_prefix = {'enhanced': f'{PROJ_DIR}/project/GMTiP-RNA/20260131/output/LRS/isoform_discovery/merged/enhanced_gtf/GTOP_novel-GENCODE_v47',}

REF_GENOME_FA=f'{PROJ_DIR}/raw_data/GMTiP/ref/LRS/genome.fa'


ARGS = sys.argv[1:]
STEP = ARGS[0]
CONF_CSV = ARGS[1]
ROOT_OUTPUT_DIR=ARGS[2]
LOG_NAME = ARGS[3]
N_TASK = int(ARGS[4])
NT_PER_TASK = int(ARGS[5])

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

LOAD_BASE_ENV_CMD='module load anaconda && source ~/.bashrc'
LOAD_POLARS_ENV_CMD='source ~/.bashrc && mamba activate polars_env'

def prepare_rsem(n_task:int, log_name:str, nt_per_task):
    logging.info(f'prepare for rsem')
    parallel_cmds=[]
    for k,v in ref_gtf_prefix.items():
        cmd=f'''{RSEM_BIN_DIR}/rsem-prepare-reference 
            --gtf {v}.gtf
            --star
            --star-path {STAR_BIN_DIR}
            -p {nt_per_task}
        	{REF_GENOME_FA}
		    {v}.RSEM_ref/ref
        '''
        cmd=re.sub(r'\s+',' ',cmd)
        log_path=f'{v}.RSEM_ref/run_RSEM_ref.log'
        series_cmds=[]
        series_cmds.append({'cmd': cmd, 'log_path': log_path})
        parallel_cmds.append(series_cmds)
    logging.info(f'remain {len(parallel_cmds)} tasks')
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task)


def run_rsem_quant(task_conf_csv:str, n_task:int, log_name:str, nt_per_task:int):
    # run_label='enhanced'
    # 1. LIBTYPE is ISR; 2. added read in some samples.
    OUTPUT_DIR=f'{ROOT_OUTPUT_DIR}/sample_based'
    df=pd.read_csv(task_conf_csv)
    parallel_cmds=[]
    for run_label in ['enhanced']:
        RSEM_ref = f'{ref_gtf_prefix[run_label]}.RSEM_ref/ref'
        for i in df.index:
            sample_id = df.loc[i, 'sample_id']
            fqs1=df.loc[i,'fqs1']
            fqs2=df.loc[i,'fqs2']
            tissue_code=sample_id.split('-')[2]
            out_quant=f'{OUTPUT_DIR}/{sample_id}/RSEM_{run_label}'
            os.makedirs(os.path.dirname(out_quant), exist_ok=True)
            cmd=f'''
                {RSEM_BIN_DIR}/rsem-calculate-expression --paired-end -p {nt_per_task}
                --star
                --star-path {STAR_BIN_DIR}
                --star-gzipped-read-file
                --no-bam-output
                {fqs1} {fqs2} {RSEM_ref} {out_quant}
            '''
            log_path = f'{out_quant}/run_RSEM.log'
            series_cmds = []
            series_cmds.append({'cmd': cmd, 'log_path': log_path})
            parallel_cmds.append(series_cmds)
    logging.info(f'remain {len(parallel_cmds)} tasks')
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task)

def combine_rsem_quant(task_conf_csv:str, n_task:int, log_name:str, nt_per_task:int):
    res_dir = f'{ROOT_OUTPUT_DIR}/sample_based'
    out_dir = f'{ROOT_OUTPUT_DIR}/combine/rsem'
    parallel_cmds = []
    df = pd.read_csv(task_conf_csv)
    for i in df.index:
        ref_tag = df.loc[i, 'ref_tag']
        gene_type = df.loc[i, 'gene_type']
        quant_type = df.loc[i, 'quant_type']
        cmd = f'''{LOAD_POLARS_ENV_CMD} && python {CURRENT_DIR}/rsem_helper.py
           combine -s {res_dir} -t {ref_tag} -g {gene_type} -q {quant_type} -o {out_dir}/quant -nt {nt_per_task}
        '''
        log_path = f'{out_dir}/log/run_merge.{"_".join([ref_tag, gene_type, quant_type])}.log'
        series_cmds = []
        series_cmds.append({'cmd': cmd, 'log_path': log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task)
    logging.info(f'finish combined')


if __name__ == '__main__':
    if STEP == 'prepare':
        prepare_rsem(n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    if STEP == 'run':
        run_rsem_quant(task_conf_csv=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    if STEP == 'combine':
        combine_rsem_quant(task_conf_csv=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
