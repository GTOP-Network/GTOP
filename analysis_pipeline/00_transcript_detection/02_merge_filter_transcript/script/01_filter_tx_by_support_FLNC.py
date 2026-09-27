# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/5/19 16:37
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_ASSEMBLY_DIR = str(Path(CONFIG_PROJECT_DIR) / 'output/assembly/LRS')
CONFIG_EXCLUDED_SIDS = ['GTOP-CA241-5032-LN-YGE8', 'GTOP-BI281-0087-LN-L7VG', 'GTOP-CB271-4155-LN-A5YE']

import logging
import sys

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

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)
EXCLUDED_SIDS = CONFIG_EXCLUDED_SIDS

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

def filter_tx(input_path, tool_name, output_path, min_read_count):
    min_read_count=int(min_read_count)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    transcript_list=[]

    if tool_name == 'isotools':
        df=pd.read_csv(input_path,sep='\t')
        transcript_list=df.loc[(df['read_count'] >= min_read_count) &
                               (df['in_IG_region'] == 0),'transcript_id'].tolist()
        # add an extra filter for The IGH locus (Immunoglobulin Heavy Locus)
        # because isotools assembly large number of transcripts in the region, especially in
        # immune-relevant tissues, such as spleen, colon, et, al.

    if tool_name == 'talon':
        df=pd.read_csv(input_path,sep='\t')
        transcript_list=df.loc[(df['read_count'] >= min_read_count) &
                               (df['in_IG_region'] == 0),'transcript_id'].tolist()

    if tool_name == 'isoseq':
        df = pd.read_csv(input_path, sep=None, engine='python', comment='#')
        id_col = next((c for c in ('pbid', 'id', 'isoform') if c in df), None)
        count_col = next((c for c in ('count_fl', 'fl_count', 'sample') if c in df), None)
        if id_col is None or count_col is None:
            raise ValueError(f'Unrecognized Iso-Seq abundance columns: {list(df.columns)}')
        transcript_list = df.loc[pd.to_numeric(df[count_col]) >= min_read_count, id_col].tolist()


    with open(output_path, "w", encoding="utf-8") as f:
        for tid in transcript_list:
            f.write(f"{tid}\n")


def run():
    from pathlib import Path
    sample_root = Path(CONFIG_ASSEMBLY_DIR) / 'sample_based'
    for sample_dir in sorted(sample_root.iterdir()):
        sid = sample_dir.name
        if not sample_dir.is_dir() or not sid.startswith('GTOP') or sid in EXCLUDED_SIDS:
            continue
        inputs = {
            'isotools': f'{sid}.read_cov.tsv.gz',
            'talon': 'tx_abundance/read_cov.tsv.gz',
            'isoseq': f'{sid}.collapsed.abundance.txt',
        }
        for tool, name in inputs.items():
            base = sample_dir / 'hg38' / tool
            filter_tx(str(base / name), tool, str(base / 'filter_tx' / f'{sid}.txt'), 3)


if __name__ == '__main__':
    args = sys.argv[1:]
    if args and args[0] == 'filter_tx':
        _, input_path, tool_name, output_path, min_read_count = args
        filter_tx(input_path, tool_name, output_path, min_read_count)
    else:
        run()
