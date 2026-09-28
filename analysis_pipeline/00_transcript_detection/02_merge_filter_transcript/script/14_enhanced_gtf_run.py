# -*- coding: utf-8 -*-


import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))

import shlex
import subprocess
import sys
from os import makedirs

import numpy as np
import pandas as pd

# config for computer evn


RUN_LOG_DIR=f'{CONFIG_PROJECT_DIR}/run_log'
OUTPUT_DIR=f'{CONFIG_PROJECT_DIR}/output/assembly/LRS/tx_merge/merged_tools'

CLEAN_SQANTI3_TAG='sqanti3_clean'
FINAL_SQANTI3_TAG='sqanti3_filter_final'
ENHANCED_GTF_TAG='enhanced_gtf'

CLEAN_SQANTI3_GT_1K_TAG='sqanti3_clean_gt_1k'
FINAL_SQANTI3_GT_1K_TAG='sqanti3_filter_final_gt_1k'
ENHANCED_GTF_GT_1K_TAG='enhanced_gtf_gt_1k'

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LOAD_BASE_ENV_CMD = os.environ.get('LOAD_BASE_ENV_CMD', 'true')

def submit_job(
        df,
        STEP,
        COMPUTER='HPC',
        PARTITION='cu-1',
        NODES = 1,
        CPU_PER_NODE = 30,
        NT_PER_TASK = 32,
        N_TASK_PER_NODE = 1,
        pipeline_path = 'enhanced_gtf_pipeline.py',
        clean_sqanti3_tag = CLEAN_SQANTI3_TAG,
        final_sqanti3_tag = FINAL_SQANTI3_TAG,
        enhanced_gtf_tag = ENHANCED_GTF_TAG,
    ):
    '''

    :param df: task dataframe. if none, use single node.
    :param STEP: keyword matching `iso-seq_pipline`
    :return:
    '''
    PARTITION = os.environ.get("SLURM_PARTITION", PARTITION)
    # assign task to nodes, build a sample list file as para for pipeline scripts per node.
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
            f'#SBATCH -J {enhanced_gtf_tag}.node{i}',
            f'#SBATCH -o {hpc_log_path}.out',
            f'#SBATCH -e {hpc_log_path}.err',
            f'#SBATCH -p {PARTITION} -N 1 -n {CPU_PER_NODE}',
            LOAD_BASE_ENV_CMD,
            f'python {CURRENT_DIR}/{pipeline_path} {STEP} {conf_path} {OUTPUT_DIR} {node_log_path} '
            f'{N_TASK_PER_NODE} {NT_PER_TASK} '
            f'{shlex.quote(clean_sqanti3_tag)} {shlex.quote(final_sqanti3_tag)} {shlex.quote(enhanced_gtf_tag)}'
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
    ARGS=sys.argv[1:]
    COMPUTER='HPC'
    STEP=ARGS[0]
    ## single task
    if STEP in {'enhanced_gtf', 'enhanced_gtf_gt_1k'}:
        PARTITION = os.environ.get("SLURM_PARTITION", 'cu-1,cpuPartition,cu-debug')
        NODES = 1
        CPU_PER_NODE = 20
        NT_PER_TASK = 20
        N_TASK_PER_NODE = 1
        if STEP == 'enhanced_gtf':
            clean_sqanti3_tag = CLEAN_SQANTI3_TAG
            final_sqanti3_tag = FINAL_SQANTI3_TAG
            enhanced_gtf_tag = ENHANCED_GTF_TAG
        else:
            clean_sqanti3_tag = CLEAN_SQANTI3_GT_1K_TAG
            final_sqanti3_tag = FINAL_SQANTI3_GT_1K_TAG
            enhanced_gtf_tag = ENHANCED_GTF_GT_1K_TAG
        submit_job(
            None, 'enhanced_gtf', COMPUTER, PARTITION, NODES, CPU_PER_NODE,
            NT_PER_TASK, N_TASK_PER_NODE,
            clean_sqanti3_tag=clean_sqanti3_tag,
            final_sqanti3_tag=final_sqanti3_tag,
            enhanced_gtf_tag=enhanced_gtf_tag,
        )
    else:
        raise ValueError('STEP must be enhanced_gtf or enhanced_gtf_gt_1k')
