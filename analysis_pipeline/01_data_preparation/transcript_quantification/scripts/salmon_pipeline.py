# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/10/6 01:42
@Desc    : Description
"""
import logging
import os
import re
import sys

import pandas as pd

from util import run_commands_threadpool, PROJ_DIR

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LOAD_POLARS_ENV_CMD='source ~/.bashrc && mamba activate polars_env'

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


# ref_gtf_prefix = {'enhanced': f'{PROJ_ROOT_DIR}/project/GMTiP-RNA/20251031/long_read/HPC/output/enhanced_gtf/GTOP_novel-GENCODE_v47',
#             'gencode': f'{PROJ_ROOT_DIR}/raw_data/GMTiP/ref/LRS/gencode.v47.annotation'}

ref_gtf_prefix = {'enhanced': f'{PROJ_DIR}/project/GMTiP-RNA/20260131/output/LRS/isoform_discovery/merged/enhanced_gtf/GTOP_novel-GENCODE_v47',}
salmon_bin=f'{PROJ_DIR}/software/salmon-latest_linux_x86_64/bin/salmon'

ARGS = sys.argv[1:]
STEP = ARGS[0]
CONF_CSV = ARGS[1]
ROOT_OUTPUT_DIR=ARGS[2]
LOG_NAME = ARGS[3]
N_TASK = int(ARGS[4])
NT_PER_TASK = int(ARGS[5])

LOAD_BASE_ENV_CMD='module load anaconda && source ~/.bashrc'

def prepare_for_salmon(n_task:int, log_name:str, nt_per_task):
    logging.info(f'prepare for salmon')
    parallel_cmds=[]
    for k,v in ref_gtf_prefix.items():
        cmd = f'{salmon_bin} index -t {v}.fa -i {v}.salmon_index -k 31 -p {nt_per_task}'
        cmd=re.sub(r'\s+',' ',cmd)
        log_path=f'{v}.salmon_index/run_salmon_ref.log'
        series_cmds=[]
        series_cmds.append({'cmd': cmd, 'log_path': log_path})
        parallel_cmds.append(series_cmds)
    logging.info(f'remain {len(parallel_cmds)} tasks')
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task)


def run_salmon_quant(task_conf_csv:str, n_task:int, log_name:str, nt_per_task:int):
    # run_label='enhanced'
    # 1. LIBTYPE is ISR; 2. added read in some samples.
    OUTPUT_DIR=f'{ROOT_OUTPUT_DIR}/sample_based'
    df=pd.read_csv(task_conf_csv)
    parallel_cmds=[]
    for run_label in ['enhanced']:
        salmon_index = f'{ref_gtf_prefix[run_label]}.salmon_index'
        for i in df.index:
            sample_id = df.loc[i, 'sample_id']
            fqs1=df.loc[i,'fqs1']
            fqs2=df.loc[i,'fqs2']
            out_quant=f'{OUTPUT_DIR}/{sample_id}/salmon_{run_label}'
            os.makedirs(os.path.dirname(out_quant), exist_ok=True)
            cmd=f'{salmon_bin} quant -p {nt_per_task} -i {salmon_index} -l ISR -1 {fqs1} -2 {fqs2} --validateMappings --seqBias --gcBias --posBias -o {out_quant}'
            log_path = f'{out_quant}/run_salmon.log'
            series_cmds = []
            series_cmds.append({'cmd': cmd, 'log_path': log_path})
            parallel_cmds.append(series_cmds)
    logging.info(f'remain {len(parallel_cmds)} tasks')
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task)

def __isoform_to_gene_sum(df: pd.DataFrame, iso2gene: pd.Series) -> pd.DataFrame:
    df = df.copy()
    df["gene_id"] = iso2gene.loc[df.index].values
    sdf=df.groupby("gene_id").sum()
    return sdf

def merge_gene(task_conf_csv:str, n_task:int, log_name:str, nt_per_task:int):
    out_dir = f'{ROOT_OUTPUT_DIR}/combine/salmon'
    parallel_cmds=[]
    df = pd.read_csv(task_conf_csv)
    for i in df.index:
        ref_tag = df.loc[i, 'ref_tag']
        quant_type = df.loc[i, 'quant_type']
        cmd=(f'{LOAD_POLARS_ENV_CMD} && python {CURRENT_DIR}/salmon_helper.py merge_gene -o {out_dir}/quant -t {ref_tag} '
             f'-q {quant_type} -g {ref_gtf_prefix[ref_tag]}')
        log_path = f'{out_dir}/log/run_merge_gene.{"_".join([ref_tag, quant_type])}.log'
        series_cmds = []
        series_cmds.append({'cmd': cmd, 'log_path': log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task)
    logging.info(f'finish combined')


def combine_salmon_quant(task_conf_csv, n_task, log_name, nt_per_task):
    res_dir = f'{ROOT_OUTPUT_DIR}/sample_based'
    out_dir = f'{ROOT_OUTPUT_DIR}/combine/salmon'
    parallel_cmds = []
    df = pd.read_csv(task_conf_csv)
    for i in df.index:
        ref_tag = df.loc[i, 'ref_tag']
        quant_type = df.loc[i, 'quant_type']
        cmd = f'''{LOAD_POLARS_ENV_CMD} && python {CURRENT_DIR}/salmon_helper.py
           combine -s {res_dir} -t {ref_tag} -q {quant_type} -o {out_dir}/quant -nt {nt_per_task}
        '''
        log_path = f'{out_dir}/log/run_merge.{"_".join([ref_tag, quant_type])}.log'
        series_cmds = []
        series_cmds.append({'cmd': cmd, 'log_path': log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task)
    logging.info(f'finish combined')


if __name__ == '__main__':
    if STEP == 'prepare':
        prepare_for_salmon(n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    if STEP == 'run':
        run_salmon_quant(task_conf_csv=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    if STEP == 'combine':
        combine_salmon_quant(task_conf_csv=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    if STEP == 'merge_gene':
        merge_gene(task_conf_csv=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)