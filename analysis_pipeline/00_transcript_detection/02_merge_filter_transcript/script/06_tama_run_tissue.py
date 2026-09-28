# -*- coding: utf-8 -*-

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_TAMA_DIR = os.environ.get('TAMA_DIR', str(Path(CONFIG_PROJECT_DIR) / 'software/tama'))
CONFIG_TISSUE_META = os.environ.get('TISSUE_META', str(Path(CONFIG_PROJECT_DIR) / 'input/tissue_code.csv'))
CONFIG_TOOLS = ['bambu', 'flair', 'flames', 'isoquant', 'isoseq', 'isotools', 'talon']

import logging
import subprocess
import sys
from typing import List

import pandas as pd

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

def is_file_valid(file_path: str) -> bool:
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            first_line = f.readline()
        return bool(first_line.strip())
    except (FileNotFoundError, PermissionError):
        return False

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LOAD_TAMA_ENVS_CMD = os.environ.get('LOAD_TAMA_ENVS_CMD', 'true')
PARTITION = os.environ.get("SLURM_PARTITION", 'cpuPartition,cu-1,fat-1,cu-short,cu-debug')
JOB_TIME = '72:00:00'

TISSUE_META_PATH = CONFIG_TISSUE_META

TOOLS = CONFIG_TOOLS

result_ver=''

def tama_pipeline(conf_files, n_task:int, main_log:str, nt_per_task:int):
    logger.info(f"load {len(conf_files)} tasks.")
    parallel_cmds=[]
    for conf_file in conf_files:
        out_dir = os.path.dirname(conf_file)
        series_cmds=[]
        log_path=f'{out_dir}/run_gtf_split.log'
        # py_cmd=f'''
        # '''
        py_cmd=f'''
        {LOAD_TAMA_ENVS_CMD} &&
        python {CONFIG_TAMA_DIR}/tama_merge.py
        '''
        cmd=f'''
        {py_cmd}
        -f {conf_file}
          -p {out_dir}/merged{result_ver} 
          -e common_ends 
          -a 50 
          -m 10 
          -z 100 
          -d merge_dup
        '''

        cmd = ' '.join(cmd.split())
        series_cmds.append({'cmd':cmd, 'log_path':log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=main_log,max_workers=n_task)


def submit_tissue_job(job_name: str, conf_files: List[str], cpu_per_task: int, main_dir: str):
    submit_dir = f'{main_dir}/log/tama_submit{result_ver}/{job_name}'
    os.makedirs(submit_dir, exist_ok=True)
    conf_list_path = f'{submit_dir}/conf_files.txt'
    with open(conf_list_path, 'w') as f:
        for conf_file in conf_files:
            f.write(conf_file + '\n')

    shell_path = f'{submit_dir}/submit_job.sh'
    hpc_out = f'{submit_dir}/job.out'
    hpc_err = f'{submit_dir}/job.err'
    node_log = f'{submit_dir}/node.log'
    shell_lines = [
        '#!/bin/bash',
        f'#SBATCH --job-name=tama_{job_name}',
        '#SBATCH --nodes=1',
        '#SBATCH --ntasks=1',
        f'#SBATCH --cpus-per-task={cpu_per_task}',
        f'#SBATCH --time={JOB_TIME}',
        f'#SBATCH --output={hpc_out}',
        f'#SBATCH --error={hpc_err}',
        f'#SBATCH --partition={PARTITION}',
        '',
        f'cd {CURRENT_DIR}',
        (
            f'python {CURRENT_DIR}/06_tama_run_tissue.py run_tissue '
            f'{conf_list_path} {node_log} {cpu_per_task}'
        ),
    ]
    with open(shell_path, 'w') as f:
        for line in shell_lines:
            f.write(line + '\n')

    logger.info(f'submit tissue job: {job_name}')
    subprocess.run(['sbatch', shell_path], check=True)


def run_tissue(conf_list_path: str, main_log: str, n_task: int):
    with open(conf_list_path) as f:
        conf_files = [line.strip() for line in f if line.strip()]
    tama_pipeline(conf_files, n_task, main_log, 1)


def make_config():
    main_dir = f'{CONFIG_PROJECT_DIR}/output/assembly/LRS/tx_merge'
    tools = TOOLS
    cpu_per_task=26
    for tool in tools:
        tool_main_dir = f'{main_dir}/tool_based/{tool}'
        chr_dir = f'{tool_main_dir}/merged/chr_bed'
        conf_files=[]
        for chr in os.listdir(chr_dir):
            wdir=f'{chr_dir}/{chr}/merged'
            os.makedirs(wdir, exist_ok=True)
            conf_file = f'{wdir}/filelist.txt'
            tissue_dir=f'{chr_dir}/{chr}/tissue_based'
            with open(conf_file, 'w') as f:
                for tissue in os.listdir(tissue_dir):
                    bed_path = f'{tissue_dir}/{tissue}/merged.bed'
                    if is_file_valid(bed_path):
                        f.write(
                            '\t'.join([
                                bed_path,
                                'no_cap',
                                '1,1,1',
                                tissue
                            ]) + '\n'
                        )
                    elif is_file_valid(f'{tissue_dir}/{tissue}/filelist.txt'):
                        raise FileNotFoundError(f'Incomplete tissue merge: {bed_path}')
            logger.info(f'save conf to {conf_file}')
            if is_file_valid(conf_file):
                conf_files.append(conf_file)
        if not conf_files:
            raise ValueError(f'No completed TAMA tissue outputs for {tool}')
        job_name = f'{tool}'
        submit_tissue_job(
            job_name,
            conf_files,
            cpu_per_task,
            main_dir
        )


if __name__ == '__main__':
    args = sys.argv[1:]
    if args and args[0] == 'run_tissue':
        _, conf_list_path, main_log, n_task = args
        run_tissue(conf_list_path, main_log, int(n_task))
    else:
        make_config()
