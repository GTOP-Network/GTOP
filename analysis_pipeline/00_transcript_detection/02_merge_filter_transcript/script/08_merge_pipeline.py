# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/1/29 21:19
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_REF_DIR = os.environ.get('GTOP_REF_DIR', str(Path(CONFIG_PROJECT_DIR) / 'reference'))
CONFIG_HG38_FASTA = os.environ.get('HG38_FASTA', str(Path(CONFIG_REF_DIR) / 'hg38.fa'))
CONFIG_TOOLS = ['bambu', 'flair', 'flames', 'isoquant', 'isoseq', 'isotools', 'talon']

import logging

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


REF_DIR = CONFIG_REF_DIR
REF_GENOME_FASTA = CONFIG_HG38_FASTA

## common para
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LOAD_cerberus_ENVS_CMD = os.environ.get('LOAD_cerberus_ENVS_CMD', 'true')
gffread_PATH = os.environ.get('GFFREAD_BIN', 'gffread')

PROJ_DIR = os.environ.get("PROJECT_ROOT")
main_dir = f'{CONFIG_PROJECT_DIR}/output/assembly/LRS/tx_merge'
merge_main_dir=f'{main_dir}/merged_tools'
gtf_dir=f'{merge_main_dir}/raw_isoform'
gtf_conf_path = f'{gtf_dir}/gtfs_path.csv'

raw_merged_gtf = f'{gtf_dir}/raw_merged.gtf'
raw_merged_meta_tsv=f'{gtf_dir}/raw_merged.meta.tsv'
tool_col='sources'
n_tool_filter_gtf = f'{gtf_dir}/tool_gt3_filter.gtf'
# tool2_filter_gtf_chr_dir = f'{merge_main_dir}/chr_based'
gtf_chr_dir = f'{merge_main_dir}/chr_based'

TOOLS = CONFIG_TOOLS


def make_gtf_conf():
    gtfs={}
    new_tools=TOOLS
    for t in new_tools:
        gtfs[t]=f'{main_dir}/tool_based/{t}/merged/gtf/merged.gtf'
    data=[]
    for k,v in gtfs.items():
        data.append([k, v])
    df=pd.DataFrame(data,columns=['tool_name','gtf_path'])
    os.makedirs(os.path.dirname(gtf_conf_path), exist_ok=True)
    df.to_csv(gtf_conf_path, index=False)
    pass

def merge_gtf_by_intron_chain():
    cmd=(f'{LOAD_cerberus_ENVS_CMD} '
         f'&& python {CURRENT_DIR}/merge_gtf_by_ic.py {gtf_conf_path} {raw_merged_gtf}')
    print(cmd)
    __import__('subprocess').run(cmd, shell=True, executable='/bin/bash', check=True)
    pass


def split_gtf_chr():
    input_gtf = raw_merged_gtf
    outdir = gtf_chr_dir
    os.makedirs(outdir, exist_ok=True)
    valid = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
    handles = {}
    def get_handle(name):
        if name not in handles:
            odir=f'{outdir}/{name}/gtf'
            os.makedirs(odir, exist_ok=True)
            handles[name] = open(f'{odir}/merged.gtf', "w")
        return handles[name]

    with open(input_gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chr_ = line.split("\t", 1)[0]
            if chr_ in valid:
                get_handle(chr_).write(line)
            else:
                get_handle("chr_other").write(line)
    for h in handles.values():
        h.close()


import re
from collections import defaultdict


if __name__ == '__main__':
    make_gtf_conf()
    merge_gtf_by_intron_chain()
    split_gtf_chr()
    # prefilter_by_n_tools(n_tools=3)
    # make_fa()
    # make_bed()


