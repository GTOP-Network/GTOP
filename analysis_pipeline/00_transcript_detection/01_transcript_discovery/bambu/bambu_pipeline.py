# -*- coding: utf-8 -*-


import os.path
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

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

ARGS = sys.argv[1:]
STEP = ARGS[0]
CONF_CSV = ARGS[1]
LOG_NAME = ARGS[2]
N_TASK = int(ARGS[3])
NT_PER_TASK = int(ARGS[4])

R_env = os.environ.get('BAMBU_RSCRIPT', 'Rscript')

REQUIRED_OUTPUT_NAMES = [
    'discoveryOnly.ndr.default.gtf',
]
OPTIONAL_TRANSFER_NAMES = [
]
def _sample_label(row):
    sample_id = row.get('sample_id', row.get('mapped_bam_path', 'NA'))
    ref_genome_tag = row.get('ref_genome_tag', '')
    return f'{sample_id}:{ref_genome_tag}' if ref_genome_tag else sample_id


def _is_nonempty_file(file_path):
    return os.path.exists(file_path) and os.path.getsize(file_path) > 0


def _rscript_cmd(script_path, args):
    return ' '.join(shlex.quote(str(x)) for x in [R_env, script_path] + list(args))


def _quote_cmd(parts):
    return ' '.join(shlex.quote(str(x)) for x in parts)


def _deduplicate_prepare_jobs(df):
    if 'ref_gtf_path' not in df.columns or 'refanno_path' not in df.columns:
        raise Exception('ref_gtf_path and refanno_path columns are required for prepare step')
    return df.drop_duplicates(subset=['refanno_path']).copy()


def prepare(conf_path, n_task, log_name):
    parallel_cmds=[]
    df = _deduplicate_prepare_jobs(pd.read_csv(conf_path))
    skip_n = 0
    for _, row in df.iterrows():
        ref_gtf_path = row['ref_gtf_path']
        refanno_path = row['refanno_path']
        if _is_nonempty_file(refanno_path):
            skip_n += 1
            print(f'[INFO] Skip existing bambu annotation: {refanno_path}')
            continue

        out_dir = os.path.dirname(refanno_path)
        os.makedirs(out_dir, exist_ok=True)
        log_path = f'{out_dir}/prepare_bambu_annotations.log'
        r_cmd = _rscript_cmd(
            f'{CURRENT_DIR}/preparing.R',
            [ref_gtf_path, refanno_path]
        )
        cmd = r_cmd
        parallel_cmds.append([{'cmd': cmd, 'log_path': log_path}])

    print(f'prepare_n={len(parallel_cmds)}; skip_existing_n={skip_n}')
    if not parallel_cmds:
        return
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task)


def bambu(conf_path, n_task, log_name, nt_per_task):
    '''
    Flair quant.
    :return:
    '''
    parallel_cmds=[]
    df=pd.read_csv(conf_path)
    for i,row in df.iterrows():
        mapped_bam_path=row['mapped_bam_path']
        out_dir=row['out_dir']
        ref_genome_path=row['ref_genome_path']
        ref_gtf_path=row['ref_gtf_path']
        refanno_path=row.get('refanno_path', f'{out_dir}/bambuAnnotations.rds')
        os.makedirs(out_dir, exist_ok=True)
        log_path=f'{out_dir}/run_bambu.log'
        series_cmds=[]
        r_cmd = _rscript_cmd(
            f'{CURRENT_DIR}/bambu.R',
            [
                str(nt_per_task),
                out_dir,
                mapped_bam_path,
                ref_genome_path,
                ref_gtf_path,
                refanno_path,
            ]
        )
        cmd=r_cmd
        series_cmds.append({'cmd':cmd, 'log_path':log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=log_name,max_workers=n_task)


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


def _check_prepare_sample(refanno_path):
    missing_files = []
    if not _is_nonempty_file(refanno_path):
        missing_files.append(refanno_path)
    return len(missing_files) == 0, missing_files


def prepare_check(conf_path):
    df = _deduplicate_prepare_jobs(pd.read_csv(conf_path))
    fail_ids = []
    succ_n = 0
    for _, row in df.iterrows():
        sample_label = _sample_label(row)
        refanno_path = row['refanno_path']
        is_success, missing_files = _check_prepare_sample(refanno_path)
        if is_success:
            succ_n += 1
        else:
            fail_ids.append(sample_label)
            print(f'[WARN] prepare sample={sample_label} missing files: {missing_files}')
    print(f'prepare_succ_n={succ_n}; prepare_fail_n={len(fail_ids)}')
    print(f'prepare_failed sample_hap: {fail_ids}')


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
        sample_id = row.get('sample_id', os.path.basename(out_dir))
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
    if STEP == 'prepare':
        prepare(conf_path=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME)
    if STEP == 'prepare_check':
        prepare_check(CONF_CSV)
    if STEP == 'run':
        bambu(conf_path=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    if STEP == 'transfer':
        transfer(CONF_CSV, n_task=N_TASK, log_name=LOG_NAME)
    if STEP == 'transfer_one':
        transfer_one(CONF_CSV, row_idx=ARGS[5], log_name=LOG_NAME)
