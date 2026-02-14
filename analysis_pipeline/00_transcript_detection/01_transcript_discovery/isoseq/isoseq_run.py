# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/10/30 11:07
@Desc    : Run Iso-Seq pipeline in HPC or Single node.
"""
import os
import subprocess
import sys
from os import makedirs

import numpy as np
import pandas as pd

# config for computer evn
HIFI_BAM_DIR='/lustre/home/cxue/raw_data/GMTiP-LR-RNA-Seq/HiFi_read'
NEW_SAMPLE_ID='/lustre/home/cxue/project/GMTiP-RNA/20251031/long_read/HPC/conf/new_13_sample_ids.txt'
OUTPUT_DIR='/lustre/home/cxue/project/GMTiP-RNA/20251031/long_read'
RUN_LOG_DIR=f'{OUTPUT_DIR}/run_log'

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

LOAD_BASE_ENV_CMD='module load anaconda && source ~/.bashrc'
LOAD_PYSAM_ENV_CMD='module load anaconda && source ~/.bashrc && conda activate pysam'
LOAD_POLARS_ENV_CMD='source ~/.bashrc && mamba activate polars_env'


def _load_hifi_info():
    # HPC bam dir
    bam_dir = HIFI_BAM_DIR
    data=[]
    for f in os.listdir(bam_dir):
        if f.endswith(".bam"):
            data.append([f'{bam_dir}/{f}', f.split('.')[0]])
    df=pd.DataFrame(data=data,columns=['bam_path','sample_id'])
    return df

def _load_hifi_info_add():
    bam_dir = HIFI_BAM_DIR
    wpath = NEW_SAMPLE_ID
    conf_df=pd.read_csv(wpath,sep='\t',header=None)
    print(f'load {conf_df.shape[0]} samples')
    data=[]
    for sid in conf_df.iloc[:,0]:
        data.append([f'{bam_dir}/{sid}.HiFi_reads.bam',sid])
    df=pd.DataFrame(data=data,columns=['bam_path','sample_id'])
    print(f'remain {df.shape[0]} jobs.')
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
        pipeline_path = 'isoseq_pipeline.py',
    ):
    '''

    :param df: task dataframe. if none, use single node.
    :param STEP: keyword matching `iso-seq_pipline`
    :return:
    '''
    # assign task to nodes, build a sample list file as para for pipeline scripts per node.
    output_dir=f'{OUTPUT_DIR}/{COMPUTER}/output'
    if df is not None:
        sub_dfs = np.array_split(df, NODES)
    else:
        sub_dfs=[None,]
    cmds = []
    for i, sdf in enumerate(sub_dfs, 1):
        node_prefix = f'{RUN_LOG_DIR}/{COMPUTER}/conf/{STEP}/node{i}'
        makedirs(os.path.dirname(node_prefix), exist_ok=True)
        conf_path = f'{node_prefix}.conf.csv'
        hpc_job_shell = f'{node_prefix}.job_shell.sh'
        hpc_log_path = f'{node_prefix}.job_log.log'
        node_log_path = f'{node_prefix}.node_log.log'
        if sdf is not None:
            sdf.to_csv(conf_path, index=False)
        # build hpc job submit list
        shell_lines = [
            f'#!/bin/bash',
            f'#SBATCH -J {STEP}.node{i}',
            f'#SBATCH -o {hpc_log_path}.out',
            f'#SBATCH -e {hpc_log_path}.err',
            f'#SBATCH -p {PARTITION} -N 1 -n {CPU_PER_NODE}',
            LOAD_POLARS_ENV_CMD,
            f'python {CURRENT_DIR}/{pipeline_path} {STEP} {conf_path} {output_dir} {node_log_path} '
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
        os.system(cmd)



if __name__ == '__main__':
    ARGS=sys.argv[1:]
    COMPUTER='HPC'
    STEP=ARGS[0]
    df=_load_hifi_info()
    # df=_load_hifi_info_add()
    chr_df = pd.DataFrame([[chr] for chr in
                           "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM others".split()],
                          columns=['chr'])
    ## for step 1: sample-based multiple task.
    if STEP == 'hifi_to_mapped_bam':
        PARTITION='cu-1,cpuPartition,fat-1'
        NODES = 177
        CPU_PER_NODE = 1
        NT_PER_TASK = 1
        N_TASK_PER_NODE = 1
        submit_job(df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE)

    ## for step 2: chr-based multiple task.
    if STEP == 'bam_to_isoform':
        PARTITION='cu-1,cpuPartition,fat-1'
        NODES = 26
        CPU_PER_NODE = 10
        NT_PER_TASK = 32
        N_TASK_PER_NODE = 1
        submit_job(chr_df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE)

    ## for step 3: single task
    if STEP == 'merge':
        PARTITION='cu-1,cpuPartition'
        NODES = 1
        CPU_PER_NODE = 10
        NT_PER_TASK = 1
        N_TASK_PER_NODE = 1
        submit_job(None,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE)


