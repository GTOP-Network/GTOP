# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/11/23 21:43
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""
import os
import sys

import numpy as np
import pandas as pd

RUN_LOG_DIR=f'/lustre/home/cxue/project/GMTiP-RNA/20260131/run_log'
# OUTPUT_DIR='/media/dubai/home/xuechao/2025-10-29-GMTiP-RNA/data/output/long_read'
OUTPUT_DIR='/lustre/home/cxue/project/GMTiP-RNA/20260131/output/SRS/03_quantification'

LOAD_BASE_ENV_CMD='module load anaconda && source ~/.bashrc'
LOAD_PYSAM_ENV_CMD='module load anaconda && source ~/.bashrc && conda activate pysam'
LOAD_POLARS_ENV_CMD='source ~/.bashrc && mamba activate polars_env'

# batch 1
fastq_dirs=[
    '/flashfs1/scratch.global/rywangz/gtop_rna_fastq',
    '/lustre/home/xdzou/data/GMTiP_RNAseq',
    '/lustre/home/xdzou/2024-10-21-GTBMap/2025-02-11-RNA-mapping/input/fastq',
    '/flashfs1/scratch.global/xdzou/GMTiP_srRNA_fastq/GTOP_RNAseq_fq'
]
passed_srRNA_sample_path='/lustre/home/cxue/raw_data/GMTiP/meta/RNA/SRS_passed_sample_id.csv'
#batch 2
# fastq_dirs=[
# '/lustre/home/xdzou/2024-10-21-GTBMap/2025-02-11-RNA-mapping/input/fastq',
# '/flashfs1/scratch.global/xdzou/GTOP_RNAseq_fq',
# ]


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
CURRENT_PY = os.path.splitext(os.path.basename(__file__))[0]
def submit_job(
        df,
        STEP,
        COMPUTER='HPC',
        PARTITION='cu-1',
        NODES = 1,
        CPU_PER_NODE = 30,
        NT_PER_TASK = 32,
        N_TASK_PER_NODE = 1,
        pipeline_py='',
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
        node_prefix = f'{RUN_LOG_DIR}/{COMPUTER}/conf/{CURRENT_PY}_{STEP}/node{i}'
        os.makedirs(os.path.dirname(node_prefix), exist_ok=True)
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
            LOAD_BASE_ENV_CMD,
            f'python {CURRENT_DIR}/{pipeline_py} {STEP} {conf_path} {output_dir} {node_log_path} '
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


def _load_all_fastq():
    sample_files={}
    # load all sample fastq files.
    sid_prefix={}
    for fastq_dir in fastq_dirs:
        for root, dirs, files in os.walk(fastq_dir):
            for file in files:
                if file.endswith('.fq.gz'):
                    sample_id=file.split('_')[0]
                    sid_prefix[sample_id]=f'{root}/{sample_id}'
    k=0
    for sample_id,prefix in sid_prefix.items():
        read_1_fq=f'{prefix}_1.fq.gz'
        read_2_fq=f'{prefix}_2.fq.gz'
        if not (os.path.isfile(read_1_fq) and os.path.isfile(read_2_fq)):
            read_1_fq=f'{prefix}_1.merged.fq.gz'
            read_2_fq=f'{prefix}_2.merged.fq.gz'
        read_add_1_fq=f'{prefix}_1.add.fq.gz'
        read_add_2_fq=f'{prefix}_2.add.fq.gz'
        if not (os.path.isfile(read_1_fq) and os.path.isfile(read_2_fq)):
            print(f'ERROR: {sample_id} has no fastq files.')
        else:
            read_1s=[read_1_fq]
            read_2s=[read_2_fq]
            if os.path.isfile(read_add_1_fq) and os.path.isfile(read_add_2_fq):
                read_1s.append(read_add_1_fq)
                read_2s.append(read_add_2_fq)
                k+=1
            sample_files[sample_id]=[' '.join(read_1s),' '.join(read_2s)]
    print(f'{k} samples with add fq; all {len(sample_files)} samples.')
    data=[]
    passed_samples=pd.read_csv(passed_srRNA_sample_path)['sample_id'].tolist()
    for sid,fqs in sample_files.items():
        if sid in passed_samples:
            data.append([sid,fqs[0],fqs[1]])
    df=pd.DataFrame(data=data,columns=['sample_id','fqs1', 'fqs2'])
    print(f'submit {len(data)} samples.')
    return df

def _load_quant_conf():
    data=[]
    for ref_tag in ['enhanced']:
        for quant_type in ['count','tpm']:
            data.append([ref_tag,quant_type])
    df=pd.DataFrame(data=data,columns=['ref_tag','quant_type'])
    print(f'loaded {df.shape[0]} conf files')
    return df


if __name__ == '__main__':
    pipeline_py='salmon_pipeline.py'
    ARGS=sys.argv[1:]
    COMPUTER='HPC'
    STEP=ARGS[0]

    if STEP == 'prepare':
        PARTITION='cu-1,cpuPartition,fat-1,cu-short'
        NODES = 1
        CPU_PER_NODE = 30
        NT_PER_TASK = 30
        N_TASK_PER_NODE = 1
        submit_job(None,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE,pipeline_py)
    if STEP == 'run':
        PARTITION='cu-1,cpuPartition,fat-1'
        NODES = 300
        CPU_PER_NODE = 10
        NT_PER_TASK = 12
        N_TASK_PER_NODE = 1
        df = _load_all_fastq()
        submit_job(df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE,pipeline_py)
    if STEP == 'combine':
        PARTITION='cu-1,cpuPartition,fat-1'
        NODES = 2
        CPU_PER_NODE = 10
        NT_PER_TASK = 10
        N_TASK_PER_NODE = 1
        df=_load_quant_conf()
        submit_job(df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE,pipeline_py)
    if STEP == 'merge_gene':
        PARTITION='cu-1,cpuPartition,fat-1'
        NODES = 1
        CPU_PER_NODE = 40
        NT_PER_TASK = 10
        N_TASK_PER_NODE = 4
        df=_load_quant_conf()
        submit_job(df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE,pipeline_py)