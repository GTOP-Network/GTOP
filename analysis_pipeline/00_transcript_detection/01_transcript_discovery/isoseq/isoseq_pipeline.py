# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/10/30 11:07
@Desc    : Run Iso-Seq pipeline in HPC or Single node.
"""

import shutil
import shlex
import os
import sys
from pathlib import Path

import pandas as pd

from datetime import datetime

def log(msg, log_file=None):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line)
    if log_file:
        with open(log_file, 'a') as f:
            f.write(line + '\n')

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

ARGS = sys.argv[1:]
STEP = ARGS[0]
CONF_CSV = ARGS[1]
LOG_NAME = ARGS[2]
N_TASK = int(ARGS[3])
NT_PER_TASK = int(ARGS[4])

LOAD_ISOSEQ_ENVS_CMD = os.environ.get('LOAD_ISOSEQ_ENVS_CMD', 'true')

REQUIRED_OUTPUT_NAMES = [
    '{sample_id}.collapsed.gff',
    '{sample_id}.collapsed.abundance.txt'
]
OPTIONAL_TRANSFER_NAMES = [
]

def make_dir(*dirname):
    for sdir in dirname:
        os.makedirs(sdir,exist_ok=True)

def _quote_cmd(parts):
    return ' '.join(shlex.quote(str(x)) for x in parts)


def _sample_label(row):
    sample_id = row.get('sample_id', 'NA')
    ref_genome_tag = row.get('ref_genome_tag', '')
    return f'{sample_id}:{ref_genome_tag}' if ref_genome_tag else sample_id


def _is_nonempty_file(file_path):
    return os.path.exists(file_path) and os.path.getsize(file_path) > 0


def _format_names(names, sample_id):
    return [name.format(sample_id=sample_id) for name in names]


def cluster_cmd(sample_id, flnc_dir, out_dir, nt):
    make_dir(out_dir)
    cmd = (
        f'{LOAD_ISOSEQ_ENVS_CMD} && ' +
        _quote_cmd([
            'isoseq', 'cluster2',
            f'{flnc_dir}/{sample_id}.flnc.bam',
            f'{out_dir}/{sample_id}.clustered.bam',
            '--sort-threads', nt,
            '--num-threads', nt,
        ])
    )
    return cmd

def pbmm2_cmd(sample_id, out_dir, ref_genome_index, nt):
    make_dir(out_dir)
    cmd = (
        f'{LOAD_ISOSEQ_ENVS_CMD} && ' +
        _quote_cmd([
            'pbmm2', 'align',
            ref_genome_index,
            f'{out_dir}/{sample_id}.clustered.bam',
            f'{out_dir}/{sample_id}.mapped.bam',
            '--preset', 'ISOSEQ',
            '--unmapped',
            '--sort',
            '-j', nt,
            '--sample', sample_id,
            '--log-level', 'INFO',
        ])
    )
    return cmd

def collapse_cmd(sample_id, out_dir, nt):
    make_dir(out_dir)
    inner_log_path = f'{out_dir}/{sample_id}.inner_collapse.log'
    out_gff = f'{out_dir}/{sample_id}.collapsed.gff'
    cmd = (
        f'{LOAD_ISOSEQ_ENVS_CMD} && ' +
        _quote_cmd([
            'isoseq', 'collapse',
            f'{out_dir}/{sample_id}.mapped.bam',
            out_gff,
            '--num-threads', nt,
            '--log-file', inner_log_path,
        ])
    )
    return cmd


def _build_required_files(out_dir, sample_id):
    return [f'{out_dir}/{name}' for name in _format_names(REQUIRED_OUTPUT_NAMES, sample_id)]


def _build_transfer_files(out_dir, sample_id):
    names = _format_names(REQUIRED_OUTPUT_NAMES + OPTIONAL_TRANSFER_NAMES, sample_id)
    return [f'{out_dir}/{name}' for name in names]


def _build_final_required_files(final_out_dir, sample_id):
    names = _format_names(REQUIRED_OUTPUT_NAMES, sample_id)
    return [f'{final_out_dir}/{name}' for name in names]


def _check_sample(out_dir, sample_id):
    missing_files = []
    for file_path in _build_required_files(out_dir, sample_id):
        if not _is_nonempty_file(file_path):
            missing_files.append(file_path)
    return len(missing_files) == 0, missing_files


def _check_final_transfer(final_out_dir, sample_id):
    missing_files = []
    for file_path in _build_final_required_files(final_out_dir, sample_id):
        if not _is_nonempty_file(file_path):
            missing_files.append(file_path)
    return len(missing_files) == 0, missing_files


def run_isoform_discovery(task_conf_csv:str, n_task:int, log_name:str, nt_per_task:int):
    '''
    Sample based multiple task.
    From HiFi read (.bam file) to chr-seperated mapped genome (.bam file).
    :return:
    '''
    df=pd.read_csv(task_conf_csv)
    parallel_cmds=[]
    for i in df.index:
        series_cmds=[]
        sample_id=df.loc[i,'sample_id']
        flnc_dir = df.loc[i, 'flnc_dir']
        out_dir = df.loc[i, 'out_dir']
        ref_genome_index = df.loc[i, 'ref_genome_index']

        # isoseq cluster
        cmd=cluster_cmd(sample_id, flnc_dir, out_dir, nt=nt_per_task)
        log_path=f'{out_dir}/{sample_id}.run_cluster.log'
        series_cmds.append({'cmd':cmd, 'log_path':log_path})

        # pbmm2 mapping
        cmd=pbmm2_cmd(sample_id, out_dir, ref_genome_index, nt=nt_per_task)
        log_path=f'{out_dir}/{sample_id}.run_pbmm2.log'
        series_cmds.append({'cmd':cmd, 'log_path':log_path})

        # isoseq collapse
        cmd=collapse_cmd(sample_id, out_dir, nt_per_task)
        log_path=f'{out_dir}/{sample_id}.run_collapse.log'
        series_cmds.append({'cmd':cmd, 'log_path':log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=log_name,max_workers=n_task)


def _is_relative_to(child_path, parent_path):
    try:
        child_path.relative_to(parent_path)
        return True
    except ValueError:
        return False


def _validate_transfer_paths(out_dir, final_out_dir, run_out_main_dir=None):
    out_path = Path(out_dir).resolve()
    final_path = Path(final_out_dir).resolve()
    if out_path == final_path:
        raise Exception(f'out_dir and final_out_dir are the same: {out_path}')
    if _is_relative_to(final_path, out_path):
        raise Exception(f'final_out_dir is inside out_dir; refuse to delete: {final_path}')
    if run_out_main_dir is not None:
        run_out_path = Path(run_out_main_dir).resolve()
        if not _is_relative_to(out_path, run_out_path):
            raise Exception(f'out_dir is not inside run_out_main_dir; refuse to delete: {out_path}')


def _transfer_one_row(row, main_log=None):
    sample_id = row['sample_id']
    sample_label = _sample_label(row)
    out_dir = row['out_dir']
    final_out_dir = row['final_out_dir']
    run_out_main_dir = row.get('run_out_main_dir')
    is_success, missing_files = _check_sample(out_dir, sample_id)
    if not is_success:
        log(f'[WARN] Skip transfer sample={sample_label}; missing files: {missing_files}', main_log)
        return sample_label, False

    _validate_transfer_paths(out_dir, final_out_dir, run_out_main_dir)
    os.makedirs(final_out_dir, exist_ok=True)
    copied_files = []
    for source_file in _build_transfer_files(out_dir, sample_id):
        if not os.path.exists(source_file):
            continue
        target_file = f'{final_out_dir}/{os.path.basename(source_file)}'
        shutil.copy2(source_file, target_file)
        copied_files.append(target_file)

    is_transferred, final_missing_files = _check_final_transfer(final_out_dir, sample_id)
    if not is_transferred:
        log(
            f'[WARN] Skip remove sample={sample_label}; '
            f'final missing files after copy: {final_missing_files}',
            main_log,
        )
        return sample_label, False

    # Keep working outputs and logs after copying.
    log(
        f'[INFO] Transferred sample={sample_label}; copied_files={len(copied_files)}; '
        f'retained_work_dir={out_dir}',
        main_log,
    )
    return sample_label, True


def transfer_one(task_conf_csv, row_idx, log_name):
    df = pd.read_csv(task_conf_csv)
    row = df.iloc[int(row_idx)]
    _, is_transferred = _transfer_one_row(row, log_name)
    if not is_transferred:
        sys.exit(1)


def transfer(task_conf_csv, n_task=1, log_name=None):
    df = pd.read_csv(task_conf_csv)
    if 'final_out_dir' not in df.columns:
        raise Exception('final_out_dir column is required for transfer step')

    n_task = max(1, int(n_task))
    parallel_cmds = []
    for i in df.index:
        sample_id = df.loc[i, 'sample_id']
        out_dir = df.loc[i, 'out_dir']
        cmd = _quote_cmd([
            sys.executable,
            os.path.abspath(__file__),
            'transfer_one',
            task_conf_csv,
            log_name,
            1,
            1,
            i,
        ])
        log_path = f'{out_dir}/{sample_id}.run_transfer.log'
        parallel_cmds.append([{'cmd': cmd, 'log_path': log_path}])
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task, main_log_prefix='transfer')


if __name__ == '__main__':
    # sample based task
    if STEP == 'run':
        run_isoform_discovery(task_conf_csv=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    if STEP == 'transfer':
        transfer(CONF_CSV, n_task=N_TASK, log_name=LOG_NAME)
    if STEP == 'transfer_one':
        transfer_one(CONF_CSV, row_idx=ARGS[5], log_name=LOG_NAME)
