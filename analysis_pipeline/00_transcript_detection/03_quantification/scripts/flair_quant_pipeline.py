# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/12/29 21:19
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_FLNC_DIR = os.environ.get('FLNC_DIR', str(Path(CONFIG_PROJECT_DIR) / 'input/flnc'))
CONFIG_EXCLUDED_SIDS = ['GTOP-CA241-5032-LN-YGE8', 'GTOP-BI281-0087-LN-L7VG', 'GTOP-CB271-4155-LN-A5YE']

import os.path
import sys
import polars as pl
import subprocess
from datetime import datetime
import numpy as np
import pandas as pd
import re
from concurrent.futures import ThreadPoolExecutor, as_completed


ARGS = sys.argv[1:]
STEP = ARGS[0]
CONF_CSV = ARGS[1]
LOG_NAME = ARGS[2]
N_TASK = int(ARGS[3])
NT_PER_TASK = int(ARGS[4])

# LOAD_FLAIR_ENVS_CMD='module load anaconda && source ~/.bashrc && mamba activate flair'
LOAD_FLAIR_ENVS_CMD = os.environ.get('LOAD_FLAIR_ENVS_CMD', 'true')

def log(msg, log_file=None):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line)
    if log_file:
        with open(log_file, 'a') as f:
            f.write(line + '\n')
def format_to_decimal(df, n=6):
    def format_func(x):
        if isinstance(x, (int, float, np.number)):
            return str(round(float(x), n))
        else:
            return str(x)

    new_df=df.map(format_func)
    return new_df


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


def flair_quant(conf_path, n_task, log_name, nt_per_task):
    '''
    Flair quant.
    :return:
    '''
    parallel_cmds=[]
    df=pd.read_csv(conf_path)
    for i,row in df.iterrows():
        fa_path=row['fa_path']
        fastq_path=row['fastq_path']
        out_dir=row['out_dir']
        log_path=f'{out_dir}/run_flair_quant.log'
        series_cmds=[]
        os.makedirs(out_dir, exist_ok=True)
        cmd=f'''
            {LOAD_FLAIR_ENVS_CMD} && 
            flair quantify -r {fastq_path} -i {fa_path} --threads {nt_per_task} --sample_id_only 
            --output {out_dir}/flair_quant 
        '''
        cmd=' '.join(re.split(r'\s+', cmd))
        series_cmds.append({'cmd':cmd, 'log_path':log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=log_name,max_workers=n_task)


def combine_flair(main_dir,ref_name='GTOP'):
    quant_dir=f'{main_dir}/flair_quant'
    flair_task_out_dir=f'{quant_dir}/output/{ref_name}'
    flair_combined_out_dir=f'{quant_dir}/combined'
    os.makedirs(flair_combined_out_dir, exist_ok=True)

    merged_df = None
    expected_samples = sorted(s for s in os.listdir(CONFIG_FLNC_DIR)
                              if s.startswith('GTOP') and s not in CONFIG_EXCLUDED_SIDS
                              and os.path.isdir(os.path.join(CONFIG_FLNC_DIR, s)))
    for task in expected_samples:
        print(f'start {task}')
        count_path=f'{flair_task_out_dir}/{task}/flair_quant.counts.tsv'
        if not os.path.exists(count_path):
            raise FileNotFoundError(count_path)
        df = pl.read_csv(count_path,separator='\t',has_header=True)
        gene_col = df.columns[0]
        df = df.with_columns(pl.col(gene_col).cast(pl.String))
        df = df.unique(subset=[gene_col], keep="first")
        if merged_df is None:
            merged_df = df
        else:
            merged_gene_col = merged_df.columns[0]
            df_renamed = df.rename({gene_col: merged_gene_col})
            merged_df = merged_df.join(
                df_renamed,
                on=merged_gene_col,
                how="full",
                coalesce=True
            )
    if merged_df is None:
        raise ValueError('No quantification results found')
    if merged_df is not None:
        merged_df = merged_df.fill_null(0)
        merged_df = merged_df.sort(merged_df.columns[0])
    merged_df.write_csv(f'{flair_combined_out_dir}/{ref_name}.transcript.count.flair.tsv',separator='\t')
    print(f'save {ref_name}: {merged_df.shape}')


def __count_to_tpm(counts: pd.DataFrame, lengths: pd.Series) -> pd.DataFrame:
    lengths = lengths.loc[counts.index]
    rpk = counts.div(lengths / 1000, axis=0)
    tpm = rpk.div(rpk.sum(axis=0).replace(0, np.nan), axis=1).fillna(0) * 1e6
    return tpm

def __isoform_to_gene_sum(df: pd.DataFrame, iso2gene: pd.Series) -> pd.DataFrame:
    df = df.copy()
    df["gene_id"] = iso2gene.loc[df.index].values
    sdf=df.groupby("gene_id").sum()
    return sdf

def gene_expr(main_dir, gtf_prefix):
    quant_dir=f'{main_dir}/flair_quant'
    flair_combined_out_dir=f'{quant_dir}/combined'
    for ref,gtf_pre in gtf_prefix.items():
        isoform_len_df = pd.read_csv(f'{gtf_pre}.transcript_length', sep='\t')
        isoform_lengths = isoform_len_df.set_index("isoform")["length"]
        gene_annot_df = pd.read_csv(f'{gtf_pre}.gene_transcript_map.txt', sep='\t')
        isoform_gene = gene_annot_df.set_index("transcript_id")["gene_id"]
        suffix='flair.tsv'
        df=pd.read_csv(f'{flair_combined_out_dir}/{ref}.transcript.count.{suffix}',sep='\t',index_col=0)
        # calculate tpm
        tpm_df=__count_to_tpm(df, isoform_lengths)
        format_to_decimal(tpm_df).to_csv(f'{flair_combined_out_dir}/{ref}.transcript.tpm.{suffix}',sep='\t',float_format='%.6f')
        # gene level count
        gene_count_df=__isoform_to_gene_sum(df, isoform_gene)
        gene_count_df.to_csv(f'{flair_combined_out_dir}/{ref}.gene.count.{suffix}',sep='\t')
        # gene level tpm
        gene_tpm_df=__isoform_to_gene_sum(tpm_df, isoform_gene)
        format_to_decimal(gene_tpm_df).to_csv(f'{flair_combined_out_dir}/{ref}.gene.tpm.{suffix}',sep='\t',float_format='%.6f')

def combine(conf_path):
    df = pd.read_csv(conf_path)
    for i, row in df.iterrows():
        main_dir = row['main_dir']
        ref_name = row['ref_name']
        gtf_prefix = {ref_name: row['gtf_prefix']}
        combine_flair(main_dir, ref_name)
        gene_expr(main_dir, gtf_prefix)
        import shutil
        destination = Path(CONFIG_PROJECT_DIR) / 'release/LRS_quant' / ref_name
        destination.mkdir(parents=True, exist_ok=True)
        for level in ['transcript', 'gene']:
            for unit in ['count', 'tpm']:
                name = f'{ref_name}.{level}.{unit}.flair.tsv'
                source = Path(main_dir) / 'flair_quant/combined' / name
                shutil.copy2(source, destination / name)

if __name__ == '__main__':
    if STEP == 'flair_quant':
        flair_quant(conf_path=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    if STEP == 'combine':
        combine(conf_path=CONF_CSV)
