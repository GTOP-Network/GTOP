# -*- coding: utf-8 -*-

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_FLNC_DIR = os.environ.get('FLNC_DIR', str(Path(CONFIG_PROJECT_DIR) / 'input/flnc'))
CONFIG_EXCLUDED_SIDS = ['GTOP-CA241-5032-LN-YGE8', 'GTOP-BI281-0087-LN-L7VG', 'GTOP-CB271-4155-LN-A5YE']

import sys

import numpy as np
import pandas as pd


EXCLUDED_SIDS = CONFIG_EXCLUDED_SIDS

def get_LRS_correct_sample_id(sample_id):
    ind_codes = {'CC331': 'CC311'}
    tissue_merge = {'0260': '0262','FZDM':'4151','QZDM':'4150'}
    arr=sample_id.split('-')
    if len(arr) < 5:
        return sample_id
    ind_id=arr[1]
    if ind_id in ind_codes:
        ind_id=ind_codes[ind_id]
        arr[1]=ind_id
    tissue_id=arr[2]
    arr[0]='GTOP'
    if tissue_id in tissue_merge:
        tissue_id=tissue_merge[tissue_id]
        arr[2]=tissue_id
    return '-'.join(arr)

PROJ_DIR = os.environ.get("PROJECT_ROOT")
RUN_LOG_DIR=f'{CONFIG_PROJECT_DIR}/run_log'

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
py_name=os.path.basename(__file__)[:-3]

LOAD_BASE_ENV_CMD = os.environ.get('LOAD_BASE_ENV_CMD', 'true')
LOAD_PYSAM_ENV_CMD = os.environ.get('LOAD_PYSAM_ENV_CMD', 'true')
LOAD_POLARS_ENVS_CMD = os.environ.get('LOAD_POLARS_ENVS_CMD', 'true')

MAIN_DIR=f'{CONFIG_PROJECT_DIR}/output/assembly/LRS/tx_merge/merged_tools'
REF_FA=f'{MAIN_DIR}/sqanti3_filter_1/filtered.fa'
OUT_DIR=f'{MAIN_DIR}/sqanti3_filter_1/flair_quant'


def _load_flair_quant_conf_df():
    FLNC_dir = CONFIG_FLNC_DIR
    data=[]
    for f in os.listdir(FLNC_dir):
        sid = get_LRS_correct_sample_id(f)
        if sid.startswith("GTOP") and sid not in EXCLUDED_SIDS:
            data.append([f'{FLNC_dir}/{sid}/{sid}.flnc.fastq.gz', REF_FA, f'{OUT_DIR}/sample_based/{sid}'])
    df = pd.DataFrame(data, columns=['fastq_path', 'fa_path', 'out_dir'])
    return df


def submit_job(
        df,
        STEP,
        COMPUTER='HPC',
        PARTITION='cu-1',
        NODES = 1,
        CPU_PER_NODE = 30,
        NT_PER_TASK = 32,
        N_TASK_PER_NODE = 1,
        py_name='flair_quant_pipeline.py',
        conf_sep=',',has_header=True
    ):
    '''

    :param df: task dataframe. if none, use single node.
    :param STEP: keyword matching `iso-seq_pipline`
    :return:
    '''
    PARTITION = os.environ.get("SLURM_PARTITION", PARTITION)
    # assign task to nodes, build a sample list file as para for pipeline scripts per node.
    if df is not None:
        if df.empty:
            raise ValueError('No samples found for this step')
        if NODES>df.shape[0]:
            NODES=df.shape[0]
        sub_dfs = np.array_split(df, NODES)
    else:
        sub_dfs=[None,]
    cmds = []
    for i, sdf in enumerate(sub_dfs, 1):
        node_prefix = f'{RUN_LOG_DIR}/{COMPUTER}/conf/{py_name}_{STEP}/node{i}'
        os.makedirs(os.path.dirname(node_prefix), exist_ok=True)
        conf_path = f'{node_prefix}.conf.txt'
        hpc_job_shell = f'{node_prefix}.job_shell.sh'
        hpc_log_path = f'{node_prefix}.job_log.log'
        node_log_path = f'{node_prefix}.node_log.log'
        if sdf is not None:
            sdf.to_csv(conf_path, index=False, sep= conf_sep, header=has_header)
        # build hpc job submit list
        shell_lines = [
            f'#!/bin/bash',
            f'#SBATCH -J {STEP}.node{i}',
            f'#SBATCH -o {hpc_log_path}.out',
            f'#SBATCH -e {hpc_log_path}.err',
            f'#SBATCH -p {PARTITION} -N 1 -n {CPU_PER_NODE}',
            LOAD_BASE_ENV_CMD,
            LOAD_POLARS_ENVS_CMD,
            f'python {CURRENT_DIR}/{py_name} {STEP} {conf_path} {node_log_path} '
            f'{N_TASK_PER_NODE} {NT_PER_TASK}'
        ]
        with open(hpc_job_shell, "w") as f:
            for line in shell_lines:
                f.write(line + "\n")
        cmd = f'''
            sbatch {hpc_job_shell}
        '''
        cmd = ' '.join(cmd.split())
        cmds.append(cmd)
    for i, cmd in enumerate(cmds, 1):
        print(cmd)
        __import__('subprocess').run(cmd, shell=True, executable='/bin/bash', check=True)


if __name__ == '__main__':
    COMPUTER='HPC'
    ARGS=sys.argv[1:]
    STEP=ARGS[0]
    ## for step 1: sample-based multiple task.
    if STEP == 'flair_quant':
        df=_load_flair_quant_conf_df()
        PARTITION = os.environ.get("SLURM_PARTITION", 'cu-1,cpuPartition,fat-1')
        NODES = 400
        CPU_PER_NODE = 16
        NT_PER_TASK = 16
        N_TASK_PER_NODE = 1
        submit_job(df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE)

    ## for step 2: combine
    if STEP == 'combine':
        PARTITION = os.environ.get("SLURM_PARTITION", 'cu-1,cpuPartition,fat-1')
        NODES = 1
        CPU_PER_NODE = 10
        NT_PER_TASK = 1
        N_TASK_PER_NODE = 1
        submit_job(None,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE)
