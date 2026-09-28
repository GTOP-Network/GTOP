# -*- coding: utf-8 -*-

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_FLNC_DIR = os.environ.get('FLNC_DIR', str(Path(CONFIG_PROJECT_DIR) / 'input/flnc'))
CONFIG_EXCLUDED_SIDS = ['GTOP-CA241-5032-LN-YGE8', 'GTOP-BI281-0087-LN-L7VG', 'GTOP-CB271-4155-LN-A5YE']

import shutil
import sys

import numpy as np
import pandas as pd

PROJ_DIR = os.environ.get("PROJECT_ROOT")

RUN_LOG_DIR=f'{CONFIG_PROJECT_DIR}/run_log'

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
py_name=os.path.basename(__file__)[:-3]

LOAD_BASE_ENV_CMD = os.environ.get('LOAD_BASE_ENV_CMD', 'true')
LOAD_PYSAM_ENV_CMD = os.environ.get('LOAD_PYSAM_ENV_CMD', 'true')
LOAD_POLARS_ENVS_CMD = os.environ.get('LOAD_POLARS_ENVS_CMD', 'true')

MAIN_RESULT_DIR=f'{CONFIG_PROJECT_DIR}/output/LRS/quantification'
REF_GTF_PREFIX={
    'GTOP':f'{CONFIG_PROJECT_DIR}/release/gtf/GTOP'
}
COMBINE_REF_NAME='GTOP'

EXCLUDED_SIDS = CONFIG_EXCLUDED_SIDS
fq_dir = CONFIG_FLNC_DIR


# ref_fas={'GTOP_enhanced':gtop_enhanced_fa}
# ref_fas={'gencode':gencode_fa}

def _load_flair_quant_conf_df():
    main_dir=MAIN_RESULT_DIR
    quant_dir=f'{main_dir}/flair_quant'
    ref_fas={k:f'{v}.fa' for k,v in REF_GTF_PREFIX.items()}
    conf_dir=f'{quant_dir}/conf'
    output_dir=f'{quant_dir}/output'
    if os.path.isdir(conf_dir):
        shutil.rmtree(conf_dir)
    os.makedirs(conf_dir, exist_ok=True)
    n=0
    for f in os.listdir(fq_dir):
        if f.startswith('GTOP') and os.path.isdir(os.path.join(fq_dir, f)) and f not in EXCLUDED_SIDS:
            fa=f'{fq_dir}/{f}/{f}.flnc.fastq.gz'
            if not os.path.exists(fa):
                raise Exception(f'{fa} does not exist')
            sample_id=f
            data=[[sample_id,'cond1','batch1',fa]]
            df=pd.DataFrame(data)
            df.to_csv(f'{conf_dir}/{sample_id}.flair_quant_fa.conf.txt', index=False, sep='\t', header=False)
            n+=1
    print(f'find {n} samples fq')
    # all task
    xdata = []
    for ref_name,ref_fa in ref_fas.items():
        for f in os.listdir(conf_dir):
            sid=f.split('.')[0]
            xdata.append([f'{conf_dir}/{f}',ref_fa,f'{output_dir}/{ref_name}/{sid}'])
    df=pd.DataFrame(xdata,columns=['fastq_path','fa_path','out_dir'])
    return df

def _load_combine_conf_df():
    main_dir = MAIN_RESULT_DIR
    ref_name = COMBINE_REF_NAME
    gtf_prefix = REF_GTF_PREFIX[ref_name]
    df = pd.DataFrame(
        [[main_dir, ref_name, gtf_prefix]],
        columns=['main_dir', 'ref_name', 'gtf_prefix']
    )
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
        if NODES > df.shape[0]:
            NODES = df.shape[0]
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
            f'#SBATCH -J flair_quant.{STEP}.node{i}',
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


def check():
    ref_key=COMBINE_REF_NAME
    res_dir=f'{MAIN_RESULT_DIR}/flair_quant/output/{ref_key}'
    fail_ids=[]
    succ_n=0
    for f in os.listdir(res_dir):
        quant_f=f'{res_dir}/{f}/flair_quant.counts.tsv'
        if os.path.exists(quant_f) and os.path.getsize(quant_f) > 1024:
            succ_n+=1
        else:
            fail_ids.append(f)
    print(f'succ_n={succ_n}; fail_n={len(fail_ids)}')
    print(f'failed id: {fail_ids}')

if __name__ == '__main__':
    COMPUTER='HPC'
    ARGS=sys.argv[1:]
    STEP=ARGS[0]
    ## for step 1: sample-based multiple task.
    if STEP == 'flair_quant':
        df=_load_flair_quant_conf_df()
        PARTITION = os.environ.get("SLURM_PARTITION", 'cu-1,cpuPartition,fat-1,cu-short')
        NODES = 400
        CPU_PER_NODE = 12
        NT_PER_TASK = 12
        N_TASK_PER_NODE = 1
        submit_job(df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE)

    ## for step 2: combine
    if STEP == 'combine':
        df = _load_combine_conf_df()
        PARTITION = os.environ.get("SLURM_PARTITION", 'cu-1,cpuPartition,fat-1')
        NODES = 1
        CPU_PER_NODE = 10
        NT_PER_TASK = 1
        N_TASK_PER_NODE = 1
        submit_job(df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE)

    if STEP == 'check':
        check()
