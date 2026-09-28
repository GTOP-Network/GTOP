# -*- coding: utf-8 -*-

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_ASSEMBLY_DIR = str(Path(CONFIG_PROJECT_DIR) / 'output/assembly/LRS')
CONFIG_TOOLS = ['bambu', 'flair', 'flames', 'isoquant', 'isoseq', 'isotools', 'talon']
CONFIG_EXCLUDED_SIDS = ['GTOP-CA241-5032-LN-YGE8', 'GTOP-BI281-0087-LN-L7VG', 'GTOP-CB271-4155-LN-A5YE']

import logging
import sys

import numpy as np
import pandas as pd
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


EXCLUDED_SIDS = CONFIG_EXCLUDED_SIDS

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

MAIN_DIR = f'{CONFIG_PROJECT_DIR}/output/assembly/LRS'
RUN_LOG_DIR = f'{MAIN_DIR}/run_log'
LOG_DIR = f"{MAIN_DIR}/log/internal"

TOOLS = CONFIG_TOOLS

LOAD_POLARS_ENV_CMD = os.environ.get('LOAD_POLARS_ENV_CMD', 'true')

def split_gtf(gtf_tasks, n_task:int, main_log:str, nt_per_task:int):
    logger.info(f"load {len(gtf_tasks)} samples.")
    parallel_cmds=[]
    for gtf_path,out_dir,tx_list_path in gtf_tasks:
        series_cmds=[]
        log_path=f'{out_dir}/run_gtf_split.log'
        tx_filter=''
        if tx_list_path != 'none':
            tx_filter=f'-t {tx_list_path}'
        cmd=f'''
        python {CURRENT_DIR}/split_gtf_by_chr_helper.py {gtf_path} {out_dir} {tx_filter}
        '''
        cmd = ' '.join(cmd.split())
        series_cmds.append({'cmd':cmd, 'log_path':log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=main_log,max_workers=n_task)

def read_task_conf(conf_path):
    df = pd.read_csv(conf_path)
    df = df.where(pd.notnull(df), None)
    return df[['gtf_path', 'out_dir', 'tx_list_path']].values.tolist()

def write_task_conf(gtf_tasks, conf_path):
    os.makedirs(os.path.dirname(conf_path), exist_ok=True)
    df = pd.DataFrame(gtf_tasks, columns=['gtf_path', 'out_dir', 'tx_list_path'])
    df.to_csv(conf_path, index=False)

def submit_job(
        gtf_tasks,
        step='split_gtf',
        computer='HPC',
        partition='cu-1,cpuPartition,fat-1',
        nodes=400,
        cpu_per_node=32,
        nt_per_task=1,
        n_task_per_node=32,
):
    partition = os.environ.get("SLURM_PARTITION", partition)
    if len(gtf_tasks) == 0:
        logger.warning("No split tasks found. Skip job submission.")
        return

    if nodes > len(gtf_tasks):
        nodes = len(gtf_tasks)
    sub_tasks = np.array_split(np.array(gtf_tasks, dtype=object), nodes)
    cmds = []

    for i, task_arr in enumerate(sub_tasks, 1):
        node_prefix = f'{RUN_LOG_DIR}/{computer}/conf/{step}/node{i}'
        conf_path = f'{node_prefix}.conf.csv'
        hpc_job_shell = f'{node_prefix}.job_shell.sh'
        hpc_log_path = f'{node_prefix}.job_log.log'
        node_log_path = f'{node_prefix}.node_log.log'

        task_list = task_arr.tolist()
        write_task_conf(task_list, conf_path)

        shell_lines = [
            '#!/bin/bash',
            f'#SBATCH -J {step}.node{i}',
            f'#SBATCH -o {hpc_log_path}.out',
            f'#SBATCH -e {hpc_log_path}.err',
            f'#SBATCH -p {partition} -N 1 -n {cpu_per_node}',
            LOAD_POLARS_ENV_CMD,
            f'python {CURRENT_DIR}/03_split_gtf_by_chr_pipeline.py run-conf {conf_path} {node_log_path} '
            f'{n_task_per_node} {nt_per_task}'
        ]
        os.makedirs(os.path.dirname(hpc_job_shell), exist_ok=True)
        with open(hpc_job_shell, "w") as f:
            for line in shell_lines:
                f.write(line + "\n")

        cmd = f'sbatch {hpc_job_shell}'
        cmds.append(cmd)

    for cmd in cmds:
        logger.info(cmd)
        __import__('subprocess').run(cmd, shell=True, executable='/bin/bash', check=True)

def get_tasks(tool_name):
    from pathlib import Path
    sample_based_dir = Path(CONFIG_ASSEMBLY_DIR) / 'sample_based'
    names = {
        'bambu': 'discoveryOnly.ndr.default.gtf',
        'flair': '{sid}.isoforms.gtf',
        'flames': 'isoform_annotated.gff3',
        'isoquant': 'OUT/OUT.transcript_models.gtf',
        'isoseq': '{sid}.collapsed.gff',
        'isotools': '{sid}.isotools.gtf',
        'talon': 'talon_observed*.gtf',
    }
    tasks = []
    for sample_dir in sorted(sample_based_dir.iterdir()):
        sid = sample_dir.name
        if not sample_dir.is_dir() or not sid.startswith('GTOP') or sid in EXCLUDED_SIDS:
            continue
        tool_dir = sample_dir / 'hg38' / tool_name
        matches = sorted(tool_dir.glob(names[tool_name].format(sid=sid)))
        if len(matches) != 1:
            raise FileNotFoundError(f'Expected one {tool_name} annotation for {sid}: {matches}')
        tx_list = 'none'
        if tool_name in {'isoseq', 'isotools', 'talon'}:
            tx_list = str(tool_dir / 'filter_tx' / f'{sid}.txt')
            if not os.path.isfile(tx_list):
                raise FileNotFoundError(tx_list)
        tasks.append([str(matches[0]), str(tool_dir / 'chr_gtf'), tx_list])
    if not tasks:
        raise ValueError(f'No samples for {tool_name} in {sample_based_dir}')
    return tasks

def get_all_tasks():
    tools = TOOLS
    gtf_tasks=[]
    for tool_name in tools:
        gtf_tasks+=get_tasks(tool_name)
    return gtf_tasks

def run_local(n_task):
    split_gtf(get_all_tasks(), n_task, f"{LOG_DIR}/split_gtf_by_chr.log", 1)

def run_conf(conf_path, main_log, n_task, nt_per_task):
    split_gtf(read_task_conf(conf_path), n_task, main_log, nt_per_task)

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
        cpu_per_node = int(args[2]) if len(args) > 2 else 4
        n_task_per_node = int(args[3]) if len(args) > 3 else 4
        nt_per_task = int(args[4]) if len(args) > 4 else 1
        submit_all(nodes, cpu_per_node, n_task_per_node, nt_per_task)
        return

    if args[0] == 'run-conf':
        run_conf(args[1], args[2], int(args[3]), int(args[4]))
        return

    run_local(int(args[0]))

if __name__ == '__main__':
    run()
