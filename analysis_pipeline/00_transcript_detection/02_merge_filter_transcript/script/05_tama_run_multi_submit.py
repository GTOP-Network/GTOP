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
CONFIG_TISSUE_META = os.environ.get('TISSUE_META', str(Path(CONFIG_PROJECT_DIR) / 'input/tissue_code.csv'))
CONFIG_TOOLS = ['bambu', 'flair', 'flames', 'isoquant', 'isoseq', 'isotools', 'talon']
CONFIG_EXCLUDED_SIDS = ['GTOP-CA241-5032-LN-YGE8', 'GTOP-BI281-0087-LN-L7VG', 'GTOP-CB271-4155-LN-A5YE']

import logging
import csv
import shlex
import subprocess
import sys
from typing import Dict, List, Tuple

import pandas as pd


EXCLUDED_SIDS = CONFIG_EXCLUDED_SIDS

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
LOAD_POLARS_ENV_CMD = os.environ.get('LOAD_POLARS_ENV_CMD', 'true')
PARTITION = os.environ.get("SLURM_PARTITION", 'cpuPartition,cu-1,fat-1,cu-short,cu-debug')
JOB_TIME = '72:00:00'


TISSUE_META_PATH = CONFIG_TISSUE_META

TOOLS = CONFIG_TOOLS

result_ver=''
MAIN_DIR = f'{CONFIG_PROJECT_DIR}/output/assembly/LRS'

SAMPLE_BASED_DIR = f'{MAIN_DIR}/sample_based'
MERGE_MAIN_DIR = f'{MAIN_DIR}/tx_merge'
RUN_LOG_DIR = f'{MAIN_DIR}/run_log'
LOG_DIR = f"{MAIN_DIR}/log/internal"

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
    if len(conf_files) == 0:
        logger.warning(f'No conf files found for {job_name}. Skip job submission.')
        return

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
            f'python {CURRENT_DIR}/05_tama_run_multi_submit.py run_tissue '
            f'{conf_list_path} {node_log} {cpu_per_task}'
        ),
    ]
    with open(shell_path, 'w') as f:
        for line in shell_lines:
            f.write(line + '\n')

    logger.info(f'submit merge job: {job_name}')
    subprocess.run(['sbatch', shell_path], check=True)


def run_tissue(conf_list_path: str, main_log: str, n_task: int):
    with open(conf_list_path) as f:
        conf_files = [line.strip() for line in f if line.strip()]
    tama_pipeline(conf_files, n_task, main_log, 1)


def load_tissue_sids() -> Tuple[Dict[str, List[str]], List[str]]:
    tissue_meta=pd.read_csv(TISSUE_META_PATH,dtype=str)
    code_tissue=dict(zip(tissue_meta['Tissue_Code'], tissue_meta['Tissue']))
    sample_ids=[]
    chrs=set()
    for sid in sorted(os.listdir(SAMPLE_BASED_DIR)):
        if sid.startswith('GTOP') and sid not in EXCLUDED_SIDS:
            sample_ids.append(sid)
            for tool in TOOLS:
                chr_dir = os.path.join(SAMPLE_BASED_DIR, sid, 'hg38', tool, 'chr_gtf')
                if not os.path.isdir(chr_dir):
                    raise FileNotFoundError(chr_dir)
                chrs.update(f[:-4] for f in os.listdir(chr_dir) if f.endswith('.gtf'))
    chrs = sorted(chrs)
    logger.info(f"load {len(sample_ids)} samples.")
    logger.info(f'load {len(chrs)} chrs.')
    tissue_sids={}
    for sid in sample_ids:
        tissue_code=sid.split('-')[2]
        tissue=code_tissue[tissue_code]
        if tissue not in tissue_sids:
            tissue_sids[tissue]=[]
        tissue_sids[tissue].append(sid)

    tissue_counts = sorted(
        ((tissue, len(sids)) for tissue, sids in tissue_sids.items()),
        key=lambda item: (-item[1], item[0])
    )
    for tissue, count in tissue_counts:
        print(f'{tissue}\t{count}')

    return tissue_sids, chrs


def get_config_tasks() -> List[List[str]]:
    tissue_sids, chrs = load_tissue_sids()
    # tissues=['Adrenal_Gland','Aortic_Arch','Uterus','Common_Iliac_Artery']
    # chr_ids=['chr1','chr17','chr13']
    tissues = sorted(tissue_sids.keys())
    tasks=[]
    for tool in TOOLS:
        tool_main_dir = f'{MERGE_MAIN_DIR}/tool_based/{tool}'
        out_dir = f'{tool_main_dir}/merged/chr_bed'
        for chrom in chrs:
            for tissue in tissues:
                wdir = f'{out_dir}/{chrom}/tissue_based/{tissue}'
                conf_file = f'{wdir}/filelist.txt'
                tasks.append([tool, tissue, chrom, conf_file, ';'.join(tissue_sids[tissue])])
    logger.info(f'load {len(tasks)} config tasks.')
    return tasks


def write_task_conf(tasks: List[List[str]], conf_path: str):
    os.makedirs(os.path.dirname(conf_path), exist_ok=True)
    with open(conf_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['tool', 'tissue', 'chrom', 'conf_file', 'sample_ids'])
        writer.writerows(tasks)


def read_task_conf(conf_path: str) -> List[List[str]]:
    with open(conf_path, newline='') as f:
        reader = csv.DictReader(f)
        return [
            [row['tool'], row['tissue'], row['chrom'], row['conf_file'], row['sample_ids']]
            for row in reader
        ]


def split_tasks(tasks: List[List[str]], n_chunks: int) -> List[List[List[str]]]:
    if len(tasks) == 0:
        return []
    n_chunks = max(1, min(n_chunks, len(tasks)))
    base, extra = divmod(len(tasks), n_chunks)
    chunks = []
    start = 0
    for i in range(n_chunks):
        size = base + (1 if i < extra else 0)
        chunks.append(tasks[start:start + size])
        start += size
    return chunks


def make_config_from_tasks(tasks: List[List[str]]):
    for tool, tissue, chrom, conf_file, sample_ids in tasks:
        os.makedirs(os.path.dirname(conf_file), exist_ok=True)
        n = 0
        with open(conf_file, 'w') as f:
            for sid in sample_ids.split(';'):
                bed_path = f'{SAMPLE_BASED_DIR}/{sid}/hg38/{tool}/chr_bed/{chrom}.bed'
                if is_file_valid(bed_path):
                    f.write(
                        '\t'.join([
                            bed_path,
                            'no_cap',
                            '1,1,1',
                            sid
                        ]) + '\n'
                    )
                    n += 1
                else:
                    logger.warning(f'no such file: {bed_path}')
        logger.info(f'save conf to {conf_file}, records={n}')


def make_config():
    make_config_from_tasks(get_config_tasks())


def submit_config_jobs(
        tasks: List[List[str]],
        step='tama_make_config',
        computer='HPC',
        partition=PARTITION,
        nodes=400,
        cpu_per_node=8,
        n_task_per_node=8,
):
    partition = os.environ.get("SLURM_PARTITION", partition)
    if len(tasks) == 0:
        logger.warning('No config tasks found. Skip job submission.')
        return

    cmds = []
    for i, task_list in enumerate(split_tasks(tasks, nodes), 1):
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
            f'python {CURRENT_DIR}/05_tama_run_multi_submit.py run-config-conf {conf_path} {node_log_path} '
            f'{n_task_per_node}'
        ]
        os.makedirs(os.path.dirname(hpc_job_shell), exist_ok=True)
        with open(hpc_job_shell, "w") as f:
            for line in shell_lines:
                f.write(line + "\n")
        cmds.append(f'sbatch {hpc_job_shell}')

    for cmd in cmds:
        logger.info(cmd)
        subprocess.run(cmd.split(), check=True)


def run_config_conf(conf_path: str, main_log: str, n_task: int):
    tasks = read_task_conf(conf_path)
    parallel_cmds = []
    for tool, tissue, chrom, conf_file, sample_ids in tasks:
        out_dir = os.path.dirname(conf_file)
        cmd = (
            f'python {CURRENT_DIR}/05_tama_run_multi_submit.py make-config-one '
            f'{shlex.quote(tool)} {shlex.quote(tissue)} {shlex.quote(chrom)} '
            f'{shlex.quote(conf_file)} {shlex.quote(sample_ids)}'
        )
        parallel_cmds.append([{'cmd': cmd, 'log_path': f'{out_dir}/make_config.log'}])
    run_commands_threadpool(parallel_cmds, main_log=main_log, max_workers=n_task)


def get_merge_conf_files() -> Dict[str, List[str]]:
    tissue_sids, chrs = load_tissue_sids()
    tissues = sorted(tissue_sids.keys())
    merge_conf_files = {}
    for tool in TOOLS:
        tool_main_dir = f'{MERGE_MAIN_DIR}/tool_based/{tool}'
        out_dir = f'{tool_main_dir}/merged/chr_bed'
        for tissue in tissues:
            job_name = f'{tool}_{tissue}'
            conf_files = []
            for chrom in chrs:
                conf_file = f'{out_dir}/{chrom}/tissue_based/{tissue}/filelist.txt'
                if is_file_valid(conf_file):
                    conf_files.append(conf_file)
                else:
                    logger.warning(f'no valid conf file: {conf_file}')
            if conf_files:
                merge_conf_files[job_name] = conf_files
    return merge_conf_files


def submit_merge_jobs(cpu_per_task: int = 26):
    merge_conf_files = get_merge_conf_files()
    logger.info(
        f'submit {len(merge_conf_files)} tool-tissue merge jobs '
        f'with {cpu_per_task} cpus per job.'
    )
    for job_name, conf_files in merge_conf_files.items():
        submit_tissue_job(
            job_name,
            conf_files,
            cpu_per_task,
            MERGE_MAIN_DIR
        )


def submit_all_config(nodes=400, cpu_per_node=8, n_task_per_node=8):
    submit_config_jobs(
        get_config_tasks(),
        nodes=nodes,
        cpu_per_node=cpu_per_node,
        n_task_per_node=n_task_per_node,
    )


def run():
    args = sys.argv[1:]
    if len(args) == 0:
        make_config()
        return

    if args[0] == 'submit-config':
        nodes = int(args[1]) if len(args) > 1 else 400
        cpu_per_node = int(args[2]) if len(args) > 2 else 8
        n_task_per_node = int(args[3]) if len(args) > 3 else 8
        submit_all_config(nodes, cpu_per_node, n_task_per_node)
        return

    if args[0] == 'submit-merge':
        cpu_per_task = int(args[1]) if len(args) > 1 else 26
        submit_merge_jobs(cpu_per_task)
        return

    if args[0] == 'run_tissue':
        _, conf_list_path, main_log, n_task = args
        run_tissue(conf_list_path, main_log, int(n_task))
        return

    if args[0] == 'run-config-conf':
        run_config_conf(args[1], args[2], int(args[3]))
        return

    if args[0] == 'make-config-one':
        make_config_from_tasks([args[1:6]])
        return

    make_config()


if __name__ == '__main__':
    run()
