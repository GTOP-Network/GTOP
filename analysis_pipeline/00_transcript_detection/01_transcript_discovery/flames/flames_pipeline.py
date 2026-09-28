# -*- coding: utf-8 -*-

import os
import re
import shlex
import shutil
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

REQUIRED_OUTPUT_NAMES = [
    'isoform_annotated.gff3',
]
OPTIONAL_TRANSFER_NAMES = [
    'run_flames.log',
]

FLAMES_PYTHON = os.environ.get('FLAMES_PYTHON', 'python')
FLAMES_ENV = os.environ.get('FLAMES_ENV', 'true')


def _clean_cmd(cmd):
    return ' '.join(re.split(r'\s+', cmd.strip()))


def _quote_cmd(parts):
    return ' '.join(shlex.quote(str(x)) for x in parts)


def _sample_label(row):
    sample_id = row.get('sample_id', 'NA')
    ref_genome_tag = row.get('ref_genome_tag', '')
    return f'{sample_id}:{ref_genome_tag}' if ref_genome_tag else sample_id


def _is_nonempty_file(file_path):
    return os.path.exists(file_path) and os.path.getsize(file_path) > 0


def _get_ref_gtf(row):
    return row.get('ref_gtf_path', row.get('ref_gtf'))


def _get_ref_genome(row):
    return row.get('ref_genome_path', row.get('ref_genome'))


def _build_sample_commands(row):
    fastq_dir = row['fastq_dir']
    ref_gtf = _get_ref_gtf(row)
    ref_genome = _get_ref_genome(row)
    out_dir = row['out_dir']
    config_file = row['config_file']
    flames_python_dir = row['flames_python_dir']

    os.makedirs(out_dir, exist_ok=True)
    if not os.path.isdir(fastq_dir):
        raise FileNotFoundError(fastq_dir)
    for path in [ref_gtf, ref_genome, config_file]:
        if not _is_nonempty_file(path):
            raise Exception(f'{path} does not exist or is empty')

    cmd = _clean_cmd(f'''
        {FLAMES_ENV} &&  
        export PATH={flames_python_dir}:$PATH &&  
        bulk_long_pipeline.py 
        --gff3 {ref_gtf}
        --genomefa {ref_genome}
        -i {fastq_dir}
        --outdir {out_dir}
        --config_file {config_file}
    ''')
    return [{'cmd': cmd, 'log_path': f'{out_dir}/run_flames.log'}]


def flames(conf_path, n_task, log_name):
    parallel_cmds = []
    df = pd.read_csv(conf_path)
    for i, row in df.iterrows():
        parallel_cmds.append(_build_sample_commands(row))
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task)


def _build_required_files(out_dir):
    return [f'{out_dir}/{name}' for name in REQUIRED_OUTPUT_NAMES]


def _build_transfer_files(out_dir):
    names = REQUIRED_OUTPUT_NAMES + OPTIONAL_TRANSFER_NAMES
    return [f'{out_dir}/{name}' for name in names]


def _build_final_required_files(final_out_dir):
    return [f'{final_out_dir}/{name}' for name in REQUIRED_OUTPUT_NAMES]


def _check_sample(out_dir):
    missing_files = []
    for file_path in _build_required_files(out_dir):
        if not _is_nonempty_file(file_path):
            missing_files.append(file_path)
    return len(missing_files) == 0, missing_files


def _check_final_transfer(final_out_dir):
    missing_files = []
    for file_path in _build_final_required_files(final_out_dir):
        if not _is_nonempty_file(file_path):
            missing_files.append(file_path)
    return len(missing_files) == 0, missing_files


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
    sample_label = _sample_label(row)
    out_dir = row['out_dir']
    final_out_dir = row['final_out_dir']
    run_out_main_dir = row.get('run_out_main_dir')
    is_success, missing_files = _check_sample(out_dir)
    if not is_success:
        log(f'[WARN] Skip transfer sample={sample_label}; missing files: {missing_files}', main_log)
        return sample_label, False

    _validate_transfer_paths(out_dir, final_out_dir, run_out_main_dir)
    os.makedirs(final_out_dir, exist_ok=True)
    copied_files = []
    for source_file in _build_transfer_files(out_dir):
        if not os.path.exists(source_file):
            continue
        target_file = f'{final_out_dir}/{os.path.basename(source_file)}'
        shutil.copy2(source_file, target_file)
        copied_files.append(target_file)

    is_transferred, final_missing_files = _check_final_transfer(final_out_dir)
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


def transfer_one(conf_path, row_idx, log_name):
    df = pd.read_csv(conf_path)
    row = df.iloc[int(row_idx)]
    _, is_transferred = _transfer_one_row(row, log_name)
    if not is_transferred:
        sys.exit(1)


def transfer(conf_path, n_task=1, log_name=None):
    df = pd.read_csv(conf_path)
    if 'final_out_dir' not in df.columns:
        raise Exception('final_out_dir column is required for transfer step')

    n_task = max(1, int(n_task))
    parallel_cmds = []
    for i, row in df.iterrows():
        out_dir = row['out_dir']
        sample_id = row['sample_id']
        cmd = _quote_cmd([
            sys.executable,
            os.path.abspath(__file__),
            'transfer_one',
            conf_path,
            log_name,
            1,
            1,
            i,
        ])
        log_path = f'{out_dir}/{sample_id}.run_transfer.log'
        parallel_cmds.append([{'cmd': cmd, 'log_path': log_path}])
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task, main_log_prefix='transfer')


if __name__ == '__main__':
    if STEP == 'run':
        flames(conf_path=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME)
    if STEP == 'transfer':
        transfer(CONF_CSV, n_task=N_TASK, log_name=LOG_NAME)
    if STEP == 'transfer_one':
        transfer_one(CONF_CSV, row_idx=ARGS[5], log_name=LOG_NAME)
