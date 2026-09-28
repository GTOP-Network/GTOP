# -*- coding: utf-8 -*-

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_RUN_DIR = os.environ.get('GTOP_RUN_DIR', str(Path(CONFIG_PROJECT_DIR) / 'scratch'))
CONFIG_REF_DIR = os.environ.get('GTOP_REF_DIR', str(Path(CONFIG_PROJECT_DIR) / 'reference'))
CONFIG_HG38_FASTA = os.environ.get('HG38_FASTA', str(Path(CONFIG_REF_DIR) / 'hg38.fa'))
CONFIG_HG38_GTF = os.environ.get('HG38_GTF', str(Path(CONFIG_REF_DIR) / 'gencode.v47.annotation.gtf'))
CONFIG_TALON_DB = os.environ.get('TALON_DB', str(Path(CONFIG_REF_DIR) / 'gencode.v47.talon.db'))
CONFIG_EXCLUDED_SIDS = ['GTOP-CA241-5032-LN-YGE8', 'GTOP-BI281-0087-LN-L7VG', 'GTOP-CB271-4155-LN-A5YE']

import shlex
import sys

import numpy as np
import pandas as pd

def get_LRS_correct_sample_id(sample_id):
    ind_codes = {'CC331': 'CC311'}
    tissue_merge = {'0260': '0262','FZDM':'4151','QZDM':'4150'}
    arr=sample_id.split('-')
    ind_id=arr[1]
    if ind_id in ind_codes:
        ind_id=ind_codes[ind_id]
        arr[1]=ind_id
    tissue_id=arr[2]
    if tissue_id in tissue_merge:
        tissue_id=tissue_merge[tissue_id]
        arr[2]=tissue_id
    return '-'.join(arr)

PROJ_DIR = os.environ.get("PROJECT_ROOT")
EXE_PROJ_DIR = os.environ.get("EXE_PROJECT_ROOT")
RUN_OUT_MAIN_DIR = CONFIG_RUN_DIR
FINAL_OUT_MAIN_DIR = CONFIG_PROJECT_DIR
RUN_LOG_DIR = f'{FINAL_OUT_MAIN_DIR}/run_log'
METHOD_NAME = 'talon'
REF_GENOME_TAGS = ['hg38']

TALON_TMP_ROOT = os.environ.get('TALON_TMP_ROOT', '/tmp')

FAILED_SAMPLES_PATH = os.environ.get(
    'FAILED_SAMPLES_PATH',
    f'{RUN_OUT_MAIN_DIR}/task_conf/talon/failed_samples.tsv'
)

REQUIRED_OUTPUT_NAMES = [
    'talon.db',
    'TALON_AND_GTF_SUCCESS.done',
    'tx_abundance/read_cov.tsv.gz',
]
REQUIRED_OUTPUT_GLOBS = [
    'talon_observed*.gtf',
]

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
py_name = os.path.basename(__file__)[:-3]

LOAD_BASE_ENV_CMD = os.environ.get('LOAD_BASE_ENV_CMD', 'true')
LOAD_TALON_ENV_CMD = os.environ.get('LOAD_TALON_ENV_CMD', 'true')
LOAD_POLARS_ENVS_CMD = os.environ.get('LOAD_POLARS_ENVS_CMD', 'true')
JOB_SAMPLE_ID_PATH = None


def _get_ref_paths(sample_id, ref_genome_tag):
    if ref_genome_tag != "hg38":
        raise ValueError("This release supports hg38 only")
    return CONFIG_HG38_FASTA, CONFIG_HG38_GTF, CONFIG_TALON_DB, "hg38", "gencode_v47"


def _load_sample_id_set(sample_path):
    sam_df = pd.read_csv(sample_path, sep='\t')
    if 'sample_id' not in sam_df.columns:
        raise Exception(f'sample_id column is required in sample list: {sample_path}')
    return set(sam_df['sample_id'].dropna().astype(str).unique().tolist())


def _load_sample_hap_set(sample_path):
    sam_df = pd.read_csv(sample_path, sep='\t')
    required_cols = ['sample_id', 'ref_genome_tag']
    missing_cols = [col for col in required_cols if col not in sam_df.columns]
    if missing_cols:
        raise Exception(
            f'{required_cols} columns are required in sample-hap list: {sample_path}; '
            f'missing={missing_cols}'
        )
    pair_df = sam_df[required_cols].dropna().astype(str)
    return set(zip(pair_df['sample_id'], pair_df['ref_genome_tag']))


def _load_target_filters(failed_only=False):
    sample_ids = None
    sample_haps = None
    if JOB_SAMPLE_ID_PATH:
        sample_ids = _load_sample_id_set(JOB_SAMPLE_ID_PATH)

    if failed_only:
        if not FAILED_SAMPLES_PATH:
            raise Exception('FAILED_SAMPLES_PATH is required when failed_only=True')
        if not os.path.exists(FAILED_SAMPLES_PATH):
            raise Exception(f'FAILED_SAMPLES_PATH does not exist: {FAILED_SAMPLES_PATH}')
        sample_haps = _load_sample_hap_set(FAILED_SAMPLES_PATH)

    return sample_ids, sample_haps


def _is_failed_only_arg(args):
    if not args:
        return False
    failed_only_values = {'--failed-only', 'failed_only', 'failed', 'true', '1'}
    value = str(args[0]).strip().lower()
    if value in failed_only_values:
        return True
    if value in {'--all', 'all', 'false', '0'}:
        return False
    raise ValueError(
        f'Unsupported run sample filter argument: {args[0]}; '
        'use --failed-only to run FAILED_SAMPLES_PATH only'
    )


def _is_nonempty_file(file_path):
    return os.path.exists(file_path) and os.path.getsize(file_path) > 0


def _has_required_glob(out_dir, pattern):
    import glob
    return any(_is_nonempty_file(path) for path in glob.glob(f'{out_dir}/{pattern}'))


def _check_talon_row_success(row):
    missing_files = []
    for name in REQUIRED_OUTPUT_NAMES:
        file_path = f"{row['out_dir']}/{name}"
        if not _is_nonempty_file(file_path):
            missing_files.append(file_path)
    for pattern in REQUIRED_OUTPUT_GLOBS:
        if not _has_required_glob(row['out_dir'], pattern):
            missing_files.append(f"{row['out_dir']}/{pattern}")
    return len(missing_files) == 0, missing_files


def _load_conf_df(apply_sample_filter=True, failed_only=False):
    input_sample_based_dir = f'{FINAL_OUT_MAIN_DIR}/output/assembly/LRS/sample_based'
    run_sample_based_dir = f'{RUN_OUT_MAIN_DIR}/output/assembly/LRS/sample_based'
    final_sample_based_dir = f'{FINAL_OUT_MAIN_DIR}/output/assembly/LRS/sample_based'
    excluded_sids = CONFIG_EXCLUDED_SIDS
    data = []
    job_samples, job_sample_haps = (
        _load_target_filters(failed_only=failed_only)
        if apply_sample_filter else (None, None)
    )
    for f in os.listdir(input_sample_based_dir):
        sid = get_LRS_correct_sample_id(f)
        if sid.startswith("GTOP") and sid not in excluded_sids:
            if job_samples is not None and sid not in job_samples:
                continue
            for ref_genome_tag in REF_GENOME_TAGS:
                if (
                    job_sample_haps is not None
                    and (sid, ref_genome_tag) not in job_sample_haps
                ):
                    continue
                mapped_bam_path = (
                    f'{input_sample_based_dir}/{sid}/{ref_genome_tag}/'
                    f'{sid}.flnc_minimap2.bam'
                )
                out_dir = f'{run_sample_based_dir}/{sid}/{ref_genome_tag}/{METHOD_NAME}'
                final_out_dir = f'{final_sample_based_dir}/{sid}/{ref_genome_tag}/{METHOD_NAME}'
                ref_genome_path, ref_gtf_path, clean_db_path, genome_build, annot_name = (
                    _get_ref_paths(sid, ref_genome_tag)
                )
                data.append([
                    sid, ref_genome_tag, mapped_bam_path, out_dir, final_out_dir,
                    RUN_OUT_MAIN_DIR, ref_genome_path, ref_gtf_path, clean_db_path,
                    genome_build, annot_name
                ])
    return pd.DataFrame(
        data,
        columns=[
            'sample_id', 'ref_genome_tag', 'mapped_bam_path', 'out_dir',
            'final_out_dir', 'run_out_main_dir', 'ref_genome_path',
            'ref_gtf_path', 'clean_db_path', 'genome_build', 'annot_name'
        ]
    )


def _load_prepare_conf_df():
    df = _load_conf_df()
    if df.empty:
        return df
    return df.drop_duplicates(subset=['clean_db_path']).copy()


def check_failed_samples():
    if not FAILED_SAMPLES_PATH:
        raise Exception('FAILED_SAMPLES_PATH is required for check_failed_samples')

    df = _load_conf_df(apply_sample_filter=False)
    failed_rows = []
    succ_n = 0
    for _, row in df.iterrows():
        sample_id = row['sample_id']
        ref_genome_tag = row['ref_genome_tag']
        is_success, missing_files = _check_talon_row_success(row)
        if not is_success:
            failed_rows.append({
                'sample_id': sample_id,
                'ref_genome_tag': ref_genome_tag,
                'sample_hap': f'{sample_id}:{ref_genome_tag}',
            })
            print(f'[WARN] sample={sample_id}; ref={ref_genome_tag}; missing_files={missing_files}')
        else:
            succ_n += 1

    failed_df = pd.DataFrame(failed_rows, columns=['sample_id', 'ref_genome_tag', 'sample_hap'])
    failed_samples_dir = os.path.dirname(FAILED_SAMPLES_PATH)
    if failed_samples_dir:
        os.makedirs(failed_samples_dir, exist_ok=True)
    failed_df.to_csv(FAILED_SAMPLES_PATH, sep='\t', index=False)
    print(f'succ_n={succ_n}; fail_n={failed_df.shape[0]}')
    print(f'failed_samples_path={FAILED_SAMPLES_PATH}')


def submit_job(
    df,
    STEP,
    COMPUTER='HPC',
    PARTITION='cu-1',
    NODES=1,
    CPU_PER_NODE=30,
    NT_PER_TASK=32,
    N_TASK_PER_NODE=1,
    py_name='talon_pipeline.py',
    conf_sep=',',
    has_header=True,
):
    PARTITION = os.environ.get("SLURM_PARTITION", PARTITION)
    if df is not None:
        if df.shape[0] == 0:
            print(f'[WARN] No tasks for step={STEP}')
            return
        if NODES > df.shape[0]:
            NODES = df.shape[0]
        sub_dfs = np.array_split(df, NODES)
    else:
        sub_dfs = [None]

    cmds = []
    for i, sdf in enumerate(sub_dfs, 1):
        node_prefix = f'{RUN_LOG_DIR}/{COMPUTER}/conf/{py_name}_{STEP}/node{i}'
        os.makedirs(os.path.dirname(node_prefix), exist_ok=True)
        conf_path = f'{node_prefix}.conf.txt'
        hpc_job_shell = f'{node_prefix}.job_shell.sh'
        hpc_log_path = f'{node_prefix}.job_log.log'
        node_log_path = f'{node_prefix}.node_log.log'
        if sdf is not None:
            sdf.to_csv(conf_path, index=False, sep=conf_sep, header=has_header)
        shell_lines = [
            '#!/bin/bash',
            f'#SBATCH -J {METHOD_NAME}.{STEP}.node{i}',
            f'#SBATCH -o {hpc_log_path}.out',
            f'#SBATCH -e {hpc_log_path}.err',
            f'#SBATCH -p {PARTITION} -N 1 -n {CPU_PER_NODE}',
            LOAD_BASE_ENV_CMD,
            LOAD_POLARS_ENVS_CMD,
            f'python {CURRENT_DIR}/{py_name} {STEP} {conf_path} {node_log_path} '
            f'{N_TASK_PER_NODE} {NT_PER_TASK}'
        ]
        if STEP in {'prepare', 'run'}:
            shell_lines = [
                '#!/bin/bash',
                f'#SBATCH -J {METHOD_NAME}.{STEP}.node{i}',
                f'#SBATCH -o {hpc_log_path}.out',
                f'#SBATCH -e {hpc_log_path}.err',
                f'#SBATCH -p {PARTITION} -N 1 -n {CPU_PER_NODE}',
                LOAD_TALON_ENV_CMD,
                f'python {CURRENT_DIR}/{py_name} {STEP} {conf_path} {node_log_path} '
                f'{N_TASK_PER_NODE} {NT_PER_TASK}'
            ]
            if STEP == 'run':
                shell_lines.insert(6, f'export TALON_TMP_ROOT={shlex.quote(TALON_TMP_ROOT)}')
        with open(hpc_job_shell, "w") as f:
            for line in shell_lines:
                f.write(line + "\n")
        cmd = f'''
            sbatch {hpc_job_shell}
        '''
        cmd = ' '.join(cmd.split())
        cmds.append(cmd)
    for cmd in cmds:
        print(cmd)
        __import__('subprocess').run(cmd, shell=True, executable='/bin/bash', check=True)


if __name__ == '__main__':
    COMPUTER = 'HPC'
    ARGS = sys.argv[1:]
    STEP = ARGS[0]
    STEP_ARGS = ARGS[1:]

    if STEP == 'prepare':
        df = _load_prepare_conf_df()
        # df = df.iloc[:2,:]
        # print(df)
        submit_job(
            df, STEP, COMPUTER,
            PARTITION='cu-1,cpuPartition,fat-1',
            NODES=400,
            CPU_PER_NODE=4,
            NT_PER_TASK=1,
            N_TASK_PER_NODE=1,
        )

    if STEP == 'prepare_check':
        df = _load_prepare_conf_df()
        submit_job(
            df, STEP, COMPUTER,
            PARTITION='cu-1,cpuPartition,fat-1',
            NODES=1,
            CPU_PER_NODE=2,
            NT_PER_TASK=1,
            N_TASK_PER_NODE=1,
        )

    if STEP == 'run':
        failed_only = _is_failed_only_arg(STEP_ARGS)
        df = _load_conf_df(failed_only=failed_only)
        # df = df.iloc[:2,:]
        # print(df)
        if failed_only:
            print(f'[INFO] Run FAILED_SAMPLES_PATH only: {FAILED_SAMPLES_PATH}')
        submit_job(
            df, STEP, COMPUTER,
            PARTITION='cu-1,cpuPartition,fat-1',
            NODES=1000,
            CPU_PER_NODE=16,
            NT_PER_TASK=16,
            N_TASK_PER_NODE=1,
        )

    if STEP == 'check':
        check_failed_samples()

    if STEP == 'transfer':
        df = _load_conf_df()
        submit_job(
            df, STEP, COMPUTER,
            PARTITION='cu-1,cpuPartition,fat-1',
            NODES=1,
            CPU_PER_NODE=32,
            NT_PER_TASK=1,
            N_TASK_PER_NODE=32,
        )
