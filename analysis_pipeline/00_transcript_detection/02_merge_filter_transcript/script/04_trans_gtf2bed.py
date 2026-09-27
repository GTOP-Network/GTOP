# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/5/18 21:26
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_TAMA_DIR = os.environ.get('TAMA_DIR', str(Path(CONFIG_PROJECT_DIR) / 'software/tama'))
CONFIG_TOOLS = ['bambu', 'flair', 'flames', 'isoquant', 'isoseq', 'isotools', 'talon']

import logging
import csv
import sys

import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_commands_threadpool(cmd_list, main_log, max_workers=4, main_log_prefix='batch'):
    if main_log:
        Path(main_log).parent.mkdir(parents=True, exist_ok=True)
    def worker(commands):
        for item in commands:
            log_path = Path(item['log_path']).resolve()
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with log_path.open('w') as handle:
                subprocess.run(item['cmd'], shell=True, executable='/bin/bash',
                               check=True, stdout=handle, stderr=subprocess.STDOUT,
                               cwd=log_path.parent)
    failures = []
    with ThreadPoolExecutor(max_workers=max(1, int(max_workers))) as pool:
        futures = [pool.submit(worker, commands) for commands in cmd_list]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as error:
                failures.append(str(error))
    if main_log:
        with open(main_log, 'a') as handle:
            handle.write(f'{main_log_prefix}: tasks={len(cmd_list)}, failures={len(failures)}\n')
            handle.writelines(error + '\n' for error in failures)
    if failures:
        raise RuntimeError(f'{len(failures)} task(s) failed; see {main_log}')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
MAIN_DIR = f'{CONFIG_PROJECT_DIR}/output/assembly/LRS'
SAMPLE_BASED_DIR = f'{MAIN_DIR}/sample_based'

TOOLS = CONFIG_TOOLS

RUN_LOG_DIR = f'{MAIN_DIR}/run_log'
LOG_DIR = f"{MAIN_DIR}/log/internal"
LOAD_POLARS_ENV_CMD = os.environ.get('LOAD_POLARS_ENV_CMD', 'true')

def trans_gtf2bed(gtf_tasks, n_task:int, main_log:str, nt_per_task:int):
    logger.info(f"load {len(gtf_tasks)} samples.")
    parallel_cmds=[]
    for gtf_path,bed_path in gtf_tasks:
        out_dir = os.path.dirname(bed_path)
        series_cmds=[]
        # run trans to bed
        log_path=f'{out_dir}/run_gtf_split.log'
        cmd=f'''
        python {CONFIG_TAMA_DIR}/tama_go/format_converter/tama_format_gtf_to_bed12_ncbi.py 
        {gtf_path} {bed_path}
        '''
        cmd = ' '.join(cmd.split())
        # tmp
        series_cmds.append({'cmd':cmd, 'log_path':log_path})

        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=main_log,max_workers=n_task)


def get_all_tasks():
    # tools = ['flair','bambu']
    # tools = ['isotools', 'isoseq','talon']
    tools = TOOLS
    gtf_tasks=[]
    for tool in tools:
        for sid in os.listdir(SAMPLE_BASED_DIR):
            if sid.startswith('GTOP'):
                cdir=f'{SAMPLE_BASED_DIR}/{sid}/hg38/{tool}/chr_gtf'
                if not os.path.isdir(cdir):
                    logger.warning(f"Skip missing chr_gtf dir: {cdir}")
                    continue
                bed_dir=f'{SAMPLE_BASED_DIR}/{sid}/hg38/{tool}/chr_bed'
                os.makedirs(bed_dir, exist_ok=True)
                for chr in os.listdir(cdir):
                    if chr.endswith('.gtf'):
                        chr_id=chr.split('.')[0]
                        # if chr_id not in ['chr14-1','chr14-2','chr1-2','chr2-2','chr22-1']:
                        #     continue
                        gtf_tasks.append([f"{cdir}/{chr}",f'{bed_dir}/{chr_id}.bed'])

    return gtf_tasks


def read_task_conf(conf_path):
    with open(conf_path, newline='') as f:
        reader = csv.DictReader(f)
        return [[row['gtf_path'], row['bed_path']] for row in reader]


def write_task_conf(gtf_tasks, conf_path):
    os.makedirs(os.path.dirname(conf_path), exist_ok=True)
    with open(conf_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['gtf_path', 'bed_path'])
        writer.writerows(gtf_tasks)


def split_tasks(tasks, n_chunks):
    n_chunks = max(1, min(n_chunks, len(tasks)))
    base, extra = divmod(len(tasks), n_chunks)
    chunks = []
    start = 0
    for i in range(n_chunks):
        size = base + (1 if i < extra else 0)
        chunks.append(tasks[start:start + size])
        start += size
    return chunks


def submit_job(
        gtf_tasks,
        step='trans_gtf2bed',
        computer='HPC',
        partition='cu-1,cpuPartition,fat-1',
        nodes=400,
        cpu_per_node=32,
        nt_per_task=1,
        n_task_per_node=32,
):
    partition = os.environ.get("SLURM_PARTITION", partition)
    if len(gtf_tasks) == 0:
        logger.warning("No trans_gtf2bed tasks found. Skip job submission.")
        return
    if nodes>len(gtf_tasks):
        nodes=len(gtf_tasks)
    sub_tasks = split_tasks(gtf_tasks, nodes)
    cmds = []
    for i, task_list in enumerate(sub_tasks, 1):
        node_prefix = f'{RUN_LOG_DIR}/{computer}/conf/{step}/node{i}'
        conf_path = f'{node_prefix}.conf.csv'
        hpc_job_shell = f'{node_prefix}.job_shell.sh'
        hpc_log_path = f'{node_prefix}.job_log.log'
        node_log_path = f'{node_prefix}.node_log.log'

        write_task_conf(task_list, conf_path)
        shell_lines = [
            '#!/bin/bash',
            f'#SBATCH -J {step}.node{i}',
            f'#SBATCH -o {hpc_log_path}.out',
            f'#SBATCH -e {hpc_log_path}.err',
            f'#SBATCH -p {partition} -N 1 -n {cpu_per_node}',
            LOAD_POLARS_ENV_CMD,
            f'python {CURRENT_DIR}/04_trans_gtf2bed.py run-conf {conf_path} {node_log_path} '
            f'{n_task_per_node} {nt_per_task}'
        ]
        os.makedirs(os.path.dirname(hpc_job_shell), exist_ok=True)
        with open(hpc_job_shell, "w") as f:
            for line in shell_lines:
                f.write(line + "\n")
        cmds.append(f'sbatch {hpc_job_shell}')

    for cmd in cmds:
        logger.info(cmd)
        __import__('subprocess').run(cmd, shell=True, executable='/bin/bash', check=True)


def run_local(n_task):
    trans_gtf2bed(get_all_tasks(), n_task, f"{LOG_DIR}/trans_gtf2bed.log", 1)


def run_conf(conf_path, main_log, n_task, nt_per_task):
    trans_gtf2bed(read_task_conf(conf_path), n_task, main_log, nt_per_task)


def submit_all(nodes=400, cpu_per_node=32, n_task_per_node=32, nt_per_task=1):
    submit_job(
        get_all_tasks(),
        nodes=nodes,
        cpu_per_node=cpu_per_node,
        n_task_per_node=n_task_per_node,
        nt_per_task=nt_per_task,
    )


def run():
    args = sys.argv[1:]
    if len(args) == 0:
        run_local(1)
        return

    if args[0] == 'submit':
        nodes = int(args[1]) if len(args) > 1 else 400
        cpu_per_node = int(args[2]) if len(args) > 2 else 8
        n_task_per_node = int(args[3]) if len(args) > 3 else 8
        nt_per_task = int(args[4]) if len(args) > 4 else 1
        submit_all(nodes, cpu_per_node, n_task_per_node, nt_per_task)
        return

    if args[0] == 'run-conf':
        run_conf(args[1], args[2], int(args[3]), int(args[4]))
        return

    run_local(int(args[0]))

if __name__ == '__main__':
    run()
