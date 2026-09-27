# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/12/24 09:16
@Email   : xuechao@szbl.ac.cn
@Desc    :
"""

# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/10/30 11:07
@Desc    : Run Iso-Seq pipeline in HPC or Single node.
"""

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_REF_DIR = os.environ.get('GTOP_REF_DIR', str(Path(CONFIG_PROJECT_DIR) / 'reference'))
CONFIG_HG38_FASTA = os.environ.get('HG38_FASTA', str(Path(CONFIG_REF_DIR) / 'hg38.fa'))
CONFIG_HG38_GTF = os.environ.get('HG38_GTF', str(Path(CONFIG_REF_DIR) / 'gencode.v47.annotation.gtf'))
CONFIG_SQANTI3_DIR = os.environ.get('SQANTI3_DIR', str(Path(CONFIG_PROJECT_DIR) / 'software/SQANTI3'))

import re
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import subprocess
import sys
from statistics import median

import numpy as np
import pandas as pd

# config for computer evn

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))


REF_DIR = CONFIG_REF_DIR
PRIMER_PATH=f'{REF_DIR}/IsoSeq_v2_primers_12.fasta'
REF_GENOME_INDEX=f'{REF_DIR}/genome_ISOSEQ.mmi'
CAGE_PEAK = f"{REF_DIR}/human.refTSS_v3.1.hg38.bed"
POLYA = f"{REF_DIR}/mouse_and_human.polyA_motif.txt"
REF_GENOME_FASTA = CONFIG_HG38_FASTA
REF_GTF = CONFIG_HG38_GTF
REF_SORTED_GTF = f"{REF_DIR}/gencode.v47.annotation.sorted.gtf.gz"
REF_GTF_FASTA = f"{REF_DIR}/gencode.v47.annotation.fa"

ARGS = sys.argv[1:]
STEP = ARGS[0]
CONF_CSV = ARGS[1]
ROOT_OUTPUT_DIR = ARGS[2]
LOG_NAME = ARGS[3]
N_TASK = int(ARGS[4])
NT_PER_TASK = int(ARGS[5])

LOAD_ISOSEQ_ENVS_CMD = os.environ.get('LOAD_ISOSEQ_ENVS_CMD', 'true')
LOAD_PBINDEX_ENVS_CMD = os.environ.get('LOAD_PBINDEX_ENVS_CMD', 'true')
LOAD_SAMTOOLS_ENVS_CMD = os.environ.get('LOAD_SAMTOOLS_ENVS_CMD', 'true')

LOAD_SQANTI3_ENVS_CMD = os.environ.get('LOAD_SQANTI3_ENVS_CMD', 'true')
SQANTI3_DIR = CONFIG_SQANTI3_DIR

# LOAD_SQANTI3_ENVS_CMD='source ~/.bashrc && mamba activate SQANTI3.env'

LOAD_isoLASER_ENVS_CMD = os.environ.get('LOAD_isoLASER_ENVS_CMD', 'true')
LOAD_PICARD_ENVS_CMD = os.environ.get('LOAD_PICARD_ENVS_CMD', 'true')
LOAD_SALMON_ENVS_CMD = os.environ.get('LOAD_SALMON_ENVS_CMD', 'true')
LOAD_POLARS_ENVS_CMD = os.environ.get('LOAD_POLARS_ENVS_CMD', 'true')


def make_dir(*dirname):
    for sdir in dirname:
        os.makedirs(sdir,exist_ok=True)

def log(msg, log_file=None):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line)
    if log_file:
        with open(log_file, 'a') as f:
            f.write(line + '\n')


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


def sqanti3_annot_filter(task_conf_csv:str, n_task:int, log_name:str, nt_per_task:int):
    '''
    :return:
    '''

    OUTPUT_DIR=f'{ROOT_OUTPUT_DIR}/chr_based'
    df=pd.read_csv(task_conf_csv)
    parallel_cmds=[]
    SQANTI3_QC = f"{SQANTI3_DIR}/sqanti3_qc.py"
    SQANTI3_FILTER = f"{SQANTI3_DIR}/sqanti3_filter.py"
    for i in df.index:
        series_cmds=[]
        chr_id=df.loc[i,'chr']
        root_out_dir=f'{OUTPUT_DIR}/{chr_id}'
        gtf_path=f'{root_out_dir}/gtf/merged.gtf'
        # sqanti3 qc and filter.
        log(f'start sqanti3')
        out_dir = f'{root_out_dir}/sqanti3'
        out_prefix = 'sqanti3'
        # rerun
        # if os.path.isfile(f'{out_dir}/{out_prefix}_corrected.gtf'):
        #     log(f'{chr_id} exist, skip')
        #     continue
        if os.path.isdir(out_dir):
            shutil.rmtree(out_dir)
        filter_out_dir = f'{out_dir}/filter'
        make_dir(out_dir,filter_out_dir)
        log_path = f'{out_dir}/run_sqanti3.log'

        chunk_size=1
        if chr_id in ['chrM','chr_other','other']:
            chunk_size=1
        cmd = f'''
                {LOAD_SQANTI3_ENVS_CMD} && 
                python {SQANTI3_QC}
                --isoforms {gtf_path}
                --refGTF {REF_GTF}
                --refFasta {REF_GENOME_FASTA}
                --CAGE_peak {CAGE_PEAK}
                --polyA_motif_list {POLYA}
                --saturation
                --report skip
                --output {out_prefix}
                --dir {out_dir}
                --cpus {nt_per_task}
                --chunks {chunk_size}
                &&
                python {SQANTI3_FILTER} rules 
                --sqanti_class {out_dir}/{out_prefix}_classification.txt 
                --filter_gtf {out_dir}/{out_prefix}_corrected.gtf
                --filter_faa {out_dir}/{out_prefix}_corrected.faa
                --skip_report
                --filter_mono_exonic
                --output {out_prefix} 
                --dir {filter_out_dir}
                --cpus {nt_per_task}
                '''
        cmd=' '.join(cmd.split())
        series_cmds.append({'cmd':cmd, 'log_path':log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=log_name,max_workers=n_task)


def merge_sqanti3_annotation():
    chr_dir=f'{ROOT_OUTPUT_DIR}/chr_based'
    out_prefix=f'{ROOT_OUTPUT_DIR}/sqanti3_merged/filtered'
    make_dir(os.path.dirname(out_prefix))
    # Process per-chromosome directories
    dfs=[]
    sorted_chrs=sorted([x for x in os.listdir(chr_dir) if x.startswith('chr')])

    # merge annotation
    for d in sorted_chrs:
        print(f"Processing class: {d}")
        # locate files
        class_path = f'{chr_dir}/{d}/sqanti3/filter/sqanti3_RulesFilter_result_classification.txt'
        if not os.path.exists(class_path):
            print(f'bad {d}')
            continue
        df = pd.read_csv(class_path, sep="\t", dtype=str, low_memory=False)
        dfs.append(df)
    mdf=pd.concat(dfs, ignore_index=True)
    mdf.to_csv(f'{out_prefix}.RulesFilter_result_classification.txt', sep="\t", index=False)
    log(f'save merged annotation.')

    # merge gtf
    gtf_output_lines=[]
    for d in sorted_chrs:
        print(f"Processing gtf: {d}")
        gtf_path = f'{chr_dir}/{d}/sqanti3/filter/sqanti3.filtered.gtf'
        with open(gtf_path, "r") as gfh:
            for ln in gfh:
                if ln.startswith("#") or not ln.strip():
                    continue
                # split to attrs
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                gtf_output_lines.append("\t".join(parts) + "\n")
    with open(f'{out_prefix}.gtf', 'w') as fout:
        fout.writelines(gtf_output_lines)
    log(f'save merged gtf.')

    # merge cds
    gtf_output_lines = []
    for d in sorted_chrs:
        print(f"Processing gtf: {d}")
        gtf_path = f'{chr_dir}/{d}/sqanti3/sqanti3_corrected.cds.gff3'
        with open(gtf_path, "r") as gfh:
            for ln in gfh:
                if ln.startswith("#") or not ln.strip():
                    continue
                # split to attrs
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                gtf_output_lines.append("\t".join(parts) + "\n")
    with open(f'{out_prefix}.full.cds.gff3', 'w') as fout:
        fout.writelines(gtf_output_lines)
    log(f'save merged CDS.')

    # merge predicted protein sequencing
    gtf_output_lines = []
    for d in sorted_chrs:
        print(f"Processing gtf: {d}")
        gtf_path = f'{chr_dir}/{d}/sqanti3/sqanti3_corrected.faa'
        with open(gtf_path, "r") as gfh:
            for ln in gfh:
                gtf_output_lines.append(ln.strip()+'\n')
    with open(f'{out_prefix}.faa', 'w') as fout:
        fout.writelines(gtf_output_lines)
    log(f'save merged protein sequence.')


def sqanti3_filter_step1():
    sqanti3_prefix=f'{ROOT_OUTPUT_DIR}/sqanti3_merged/filtered'
    filter_prefix=f'{ROOT_OUTPUT_DIR}/sqanti3_filter_1/filtered'
    tools_support_path=f'{ROOT_OUTPUT_DIR}/raw_isoform/raw_merged.meta.tsv'

    def get_transcript_id(attr):
        tid_match = re.search(r'transcript_id\s+"([^"]+)"', attr)
        if tid_match:
            return tid_match.group(1)
        return None

    def get_exon_bounds(gtf_path):
        exon_bounds = {}
        with open(gtf_path, 'r') as in_f:
            for line in in_f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9 or parts[2] != 'exon':
                    continue
                tid = get_transcript_id(parts[8])
                if not tid:
                    continue
                start, end = int(parts[3]), int(parts[4])
                if tid not in exon_bounds:
                    exon_bounds[tid] = [start, end]
                else:
                    exon_bounds[tid][0] = min(exon_bounds[tid][0], start)
                    exon_bounds[tid][1] = max(exon_bounds[tid][1], end)
        return exon_bounds

    def write_filtered_gtf(gtf_path, out_gtf, transcripts, exon_bounds):
        with open(gtf_path, 'r') as in_f, open(out_gtf, 'w') as out_f:
            for line in in_f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                tid = get_transcript_id(parts[8])
                if tid not in transcripts:
                    continue
                if parts[2] == 'transcript' and tid in exon_bounds:
                    parts[3] = str(exon_bounds[tid][0])
                    parts[4] = str(exon_bounds[tid][1])
                    line = '\t'.join(parts) + '\n'
                out_f.write(line)

    # 1. add supported tools
    df=pd.read_csv(f'{sqanti3_prefix}.RulesFilter_result_classification.txt',sep='\t')
    tool_df = pd.read_csv(tools_support_path, sep='\t')
    tool_count_df = tool_df.set_index('transcript_id').apply(pd.to_numeric, errors='coerce').fillna(0)
    tx_n_tools = tool_count_df.sum(axis=1).to_dict()
    df['n_support_tools']=pd.to_numeric(df['isoform'].map(lambda x: tx_n_tools.get(x)), errors='coerce').fillna(0).astype(int)
    n_tools_gt0_transcripts=set(df.loc[(df['n_support_tools']>0) & (df['filter_result']=='Isoform'),'isoform'].unique().tolist())
    # update gtf and annotation file
    make_dir(os.path.dirname(filter_prefix))
    out_gtf=f'{filter_prefix}.gtf'
    out_n_tools_gt0_gtf=f'{filter_prefix}.n_tools_gt0.gtf'
    out_cds_gtf = f'{filter_prefix}.cds.gff3'
    out_annot = f'{filter_prefix}.RulesFilter_result_classification.txt'
    #  with >= 3 tools
    df.loc[(df['n_support_tools']<3),'filter_result']='Artifact'
    # update annotation file
    df.to_csv(out_annot,sep='\t',index=False)
    log(f'save updated annotation file')
    transcripts=set(df.loc[df['filter_result']=='Isoform','isoform'].unique().tolist())
    log(f'isoform: {len(transcripts)}')
    # update gtf
    gtf_path=f'{sqanti3_prefix}.gtf'
    exon_bounds = get_exon_bounds(gtf_path)
    write_filtered_gtf(gtf_path, out_gtf, transcripts, exon_bounds)
    log(f'save updated gtf')
    write_filtered_gtf(gtf_path, out_n_tools_gt0_gtf, n_tools_gt0_transcripts, exon_bounds)
    log(f'save n_tools > 0 gtf')
    # update cds
    cds_gtf_path=f'{sqanti3_prefix}.full.cds.gff3'
    with open(cds_gtf_path, 'r') as in_f, open(out_cds_gtf, 'w') as out_f:
        for line in in_f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            tid_match = re.search(r'transcript_id\s+"([^"]+)"', parts[8])
            if tid_match and tid_match.group(1) in transcripts:
                out_f.write(line)
    log(f'save updated cds gff')


def custom_sqanti3_filter(exclude_sample_ids=[]):
    min_read_a_sample = 5
    min_support_samples = 2
    filter_dir=f'{ROOT_OUTPUT_DIR}/sqanti3_filter_1'
    raw_sqanti3_prefix=f'{filter_dir}/filtered'
    out_filter_dir = f'{ROOT_OUTPUT_DIR}/sqanti3_filter_final'
    if len(exclude_sample_ids)>0:
        out_filter_dir=f'{ROOT_OUTPUT_DIR}/sqanti3_filter_final_gt_1k'
    os.makedirs(out_filter_dir, exist_ok=True)
    sqanti3_prefix=f'{out_filter_dir}/filtered'

    flnc_count_path=f'{filter_dir}/flair_quant/combined/transcript.count.flair.tsv'
    first_exon_path=f'{filter_dir}/AF_filter/final_novel_first_exon.transcript_qc.tsv'
    srs_junction_path=f'{filter_dir}/SRS_junction_filter/srs_jobs/merged/srs.transcript_qc.tsv'

    # 1. add supported all reads
    sq_df=pd.read_csv(f'{raw_sqanti3_prefix}.RulesFilter_result_classification.txt',sep='\t',index_col=0)
    expr_df = pd.read_csv(flnc_count_path, index_col=0, sep='\t')
    # Add SQANTI3 Isoform transcripts missing from the expression matrix and assign 0 counts.
    raw_quant_tx_ids = sq_df.index[sq_df['filter_result'] == 'Isoform']
    expr_df = expr_df.reindex(expr_df.index.union(raw_quant_tx_ids), fill_value=0)
    # Add the remaining SQANTI3 transcripts missing from the expression matrix and assign NA.
    expr_df = expr_df.reindex(expr_df.index.union(sq_df.index), fill_value=np.nan)
    # if excluding low length samples
    if len(exclude_sample_ids) > 0:
        expr_df=expr_df.loc[:,~expr_df.columns.isin(exclude_sample_ids)]
        log(f'exclude {len(exclude_sample_ids)} samples; remain {expr_df.shape[1]} samples')
        pass
    sum_read = expr_df.sum(axis=1, min_count=1)
    support_samples = (expr_df >= min_read_a_sample).sum(axis=1).astype(float)
    support_samples[expr_df.notna().sum(axis=1) == 0] = np.nan
    sq_df['support_read'] = sum_read.reindex(sq_df.index)
    # 2. add supported samples (with read >= min_read_a_sample)
    sq_df['support_samples'] = support_samples.reindex(sq_df.index)
    # 3. add first-exon QC
    fe_df = None
    fe_df=pd.read_csv(first_exon_path,sep='\t')
    fail_af_txs=set(fe_df.loc[fe_df['qc_status']!='PASS','transcript_id'].tolist())
    sq_df['AF_filter']=sq_df.index.isin(fail_af_txs)
    # 4. add short-read RNA-seq junction QC
    srj_df=pd.read_csv(srs_junction_path,sep='\t')
    fail_sj_txs = set(srj_df.loc[srj_df['qc_pass'] != 'PASS', 'transcript_id'].tolist())
    sq_df['SRS_junction_filter']=sq_df.index.isin(fail_sj_txs) | ~sq_df.index.isin(srj_df['transcript_id'])

    # save
    advanced_filter_path=f'{sqanti3_prefix}.advanced_filter.txt'
    sq_df.to_csv(advanced_filter_path,sep='\t')

    # update gtf and annotation file
    out_gtf=f'{sqanti3_prefix}.gtf'
    out_cds_gtf = f'{sqanti3_prefix}.cds.gff3'
    out_annot = f'{sqanti3_prefix}.RulesFilter_result_classification.txt'
    df=pd.read_csv(advanced_filter_path,sep='\t')
    # all isoform with >= 10 supported reads
    df.loc[(df['support_read']<10),'filter_result']='Artifact'
    # non-FSM
    # df.loc[((df['associated_transcript']=='novel') & (df["polyA_motif_found"] == False)),'filter_result']='Artifact'
    df.loc[((df['associated_transcript']=='novel') & (df['support_samples']<min_support_samples)),'filter_result']='Artifact'
    # novel: first exon filter
    df.loc[((df['associated_transcript']=='novel') & (df['AF_filter'])),'filter_result']='Artifact'
    # SRS junction QC
    df.loc[((df['associated_transcript']=='novel') & (df['SRS_junction_filter'])),'filter_result']='Artifact'

    isoform=df.loc[df['filter_result']=='Isoform',:].shape[0]
    log(f'isoform: {isoform}')
    # update annotation file
    df.to_csv(out_annot,sep='\t',index=False)
    log(f'save updated annotation file')
    transcripts=set(df.loc[df['filter_result']=='Isoform','isoform'].unique().tolist())
    log(f'isoform: {len(transcripts)}')
    # update gtf
    gtf_path=f'{raw_sqanti3_prefix}.gtf'
    with open(gtf_path, 'r') as in_f, open(out_gtf, 'w') as out_f:
        for line in in_f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            tid_match = re.search(r'transcript_id\s+"([^"]+)"', parts[8])
            if tid_match and tid_match.group(1) in transcripts:
                out_f.write(line)
    log(f'save updated gtf')
    # update cds
    cds_gtf_path=f'{raw_sqanti3_prefix}.cds.gff3'
    with open(cds_gtf_path, 'r') as in_f, open(out_cds_gtf, 'w') as out_f:
        for line in in_f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            tid_match = re.search(r'transcript_id\s+"([^"]+)"', parts[8])
            if tid_match and tid_match.group(1) in transcripts:
                out_f.write(line)
    log(f'save updated cds gff')


if __name__ == '__main__':
    if STEP == 'run_sqanti3':
        sqanti3_annot_filter(task_conf_csv=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    if STEP == 'merge':
        merge_sqanti3_annotation()
    if STEP == 'filter_step1':
        sqanti3_filter_step1()
    if STEP == 'custom_filter':
        custom_sqanti3_filter()
    if STEP == 'custom_filter_gt_1k':
        Low_than_1k_sample_IDs = ['GTOP-BF221-0378-LN-2ZKY', 'GTOP-BG171-0378-LN-TDV5', 'GTOP-CF242-0378-LN-F6UE', 'GTOP-BL021-0378-LN-RZ42', 'GTOP-AJ221-1053-LN-DP1V', 'GTOP-AJ221-4019-LN-D5AL', 'GTOP-CF241-1236-LN-NJ1W', 'GTOP-BA131-2099-LN-H19X', 'GTOP-CF242-2099-LN-6Y37', 'GTOP-CA081-1392-LN-20Q9', 'GTOP-CB161-1392-LN-I067', 'GTOP-CB271-1584-LN-N1C6', 'GTOP-CG041-5156-LN-75CF']
        custom_sqanti3_filter(Low_than_1k_sample_IDs)
