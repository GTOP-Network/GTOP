# -*- coding: utf-8 -*-

import os
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from util import PROJ_DIR


LOAD_BASE_ENV_CMD='module load anaconda && source ~/.bashrc'
LOAD_PYSAM_ENV_CMD='module load anaconda && source ~/.bashrc && conda activate pysam'


RUN_LOG_DIR=f'{PROJ_DIR}/path/to/run_log'

OUTPUT_DIR='/path/to/out'

passed_srRNA_sample_path=f'{PROJ_DIR}/path/to/SRS_passed_sample_id.csv'
# REF_TAG='DSA_to_hg38'
# REF_TAG='hg38'
REF_TAG='enhanced'

## real data
fastq_dirs=[
    '/path/to/raw_data/GTOP/RNA/SRS/fastq',
    '/path/to/raw_data/GTOP/RNA/SRS/fastq_merged',
    '/path/to/input/fastq',
]

# downsampling data
# fastq_dirs=[
#     '/path/to/output_80m'
# ]
# RUN_LOG_DIR=f'{PROJ_DIR}/path/to/run_log/SRS_downsampling'
# OUTPUT_DIR='/path/to/SRS_downsampling_quant'

check_dir = f'{RUN_LOG_DIR}/HPC/check'
FAILED_SAMPLES_LIST_PATH = f'{check_dir}/failed_sample_ids.txt'
FAILED_SAMPLES_DETAIL_PATH = f'{check_dir}/failed_sample_details.tsv'
BAD_TRANSCRIPT_ID_REPORT_PATH = f'{check_dir}/rsem_bad_transcript_ids.tsv'
DUP_TRANSCRIPT_ID_REPORT_PATH = f'{check_dir}/rsem_duplicate_transcript_ids.tsv'
TRANSCRIPT_ID_SUMMARY_REPORT_PATH = f'{check_dir}/rsem_transcript_id_check_summary.tsv'


NASFS1_FASTQ_COPY_DIR='/path/to/raw_data/GTOP/RNA/SRS/fastq'
COPY_CHUNK_SIZE=64*1024*1024
COPY_WORKERS=4
COPY_PROGRESS_INTERVAL=5
UNPASSED_SRRNA_SAMPLE_TABLE=f'{RUN_LOG_DIR}/HPC/check/unpassed_srRNA_sample_ids.csv'
MIN_PASS_FASTQ_SIZE=1024*1024*1024

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
CURRENT_PY = os.path.splitext(os.path.basename(__file__))[0]


def _ensure_check_dir():
    os.makedirs(check_dir, exist_ok=True)


def submit_job(
        df,
        STEP,
        COMPUTER='HPC',
        PARTITION='cu-1',
        NODES = 1,
        CPU_PER_NODE = 30,
        NT_PER_TASK = 32,
        N_TASK_PER_NODE = 1,
        pipeline_py='',
    ):
    '''

    :param df: task dataframe. if none, use single node.
    :param STEP: keyword matching `iso-seq_pipline`
    :return:
    '''
    # assign task to nodes, build a sample list file as para for pipeline scripts per node.
    output_dir=f'{OUTPUT_DIR}/{COMPUTER}/output'
    if df is not None:
        if NODES>df.shape[0]:
            NODES=df.shape[0]
        sub_dfs = np.array_split(df, NODES)
    else:
        sub_dfs=[None,]
    cmds = []
    for i, sdf in enumerate(sub_dfs, 1):
        node_prefix = f'{RUN_LOG_DIR}/{COMPUTER}/conf/{CURRENT_PY}_{STEP}/node{i}'
        os.makedirs(os.path.dirname(node_prefix), exist_ok=True)
        conf_path = f'{node_prefix}.conf.csv'
        hpc_job_shell = f'{node_prefix}.job_shell.sh'
        hpc_log_path = f'{node_prefix}.job_log.log'
        node_log_path = f'{node_prefix}.node_log.log'
        if sdf is not None:
            sdf.to_csv(conf_path, index=False)
        # build hpc job submit list
        shell_lines = [
            f'#!/bin/bash',
            f'#SBATCH -J RSEM.{STEP}.node{i}',
            f'#SBATCH -o {hpc_log_path}.out',
            f'#SBATCH -e {hpc_log_path}.err',
            f'#SBATCH -p {PARTITION} -N 1 -n {CPU_PER_NODE}',
            LOAD_BASE_ENV_CMD,
            f'python {CURRENT_DIR}/{pipeline_py} {STEP} {conf_path} {output_dir} {node_log_path} '
            f'{N_TASK_PER_NODE} {NT_PER_TASK}'
        ]
        with open(hpc_job_shell, "w") as f:
            for line in shell_lines:
                f.write(line + "\n")
        cmd = f'''
            sbatch {hpc_job_shell}
        '''
        cmd = ' '.join(cmd.split())
        cmds.append(cmd)
    for i, cmd in enumerate(cmds, 1):
        print(cmd)
        os.system(cmd)

def _format_size(size):
    units=['B','KB','MB','GB','TB']
    value=float(size)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f'{value:.2f}{unit}'
        value/=1024


def _load_all_unique_fastq_gt_1g(must_in_passed_id=True):
    passed_samples=pd.read_csv(passed_srRNA_sample_path)['sample_id'].tolist()
    sample_files={}
    size_too_small_count=0
    missing_pair_count=0
    sid_prefix={}
    sample_prefixes={}
    dir_sample_counts={}
    dir_sample_ids={}
    for fastq_dir in fastq_dirs:
        dir_sample_counts[fastq_dir]=0
        dir_sample_ids[fastq_dir]=set()
        dir_seen_prefixes=set()
        for root, dirs, files in os.walk(fastq_dir):
            for file in files:
                if file.endswith('.fq.gz'):
                    sample_id='_'.join(file.split('_')[:-1])
                    prefix=f'{root}/{sample_id}'
                    if prefix in dir_seen_prefixes:
                        continue
                    dir_seen_prefixes.add(prefix)
                    dir_sample_counts[fastq_dir]+=1
                    dir_sample_ids[fastq_dir].add(sample_id)
                    sample_prefixes.setdefault(sample_id,[]).append((fastq_dir,prefix))
                    if sample_id not in sid_prefix:
                        sid_prefix[sample_id]=prefix
    duplicate_sample_ids={
        sample_id:prefixes
        for sample_id,prefixes in sample_prefixes.items()
        if len(prefixes) > 1
    }
    duplicate_prefix_count=sum(len(prefixes)-1 for prefixes in duplicate_sample_ids.values())
    print(f'duplicate sample IDs: {len(duplicate_sample_ids)}; duplicate files skipped: {duplicate_prefix_count}.')
    for sample_id,prefixes in sorted(duplicate_sample_ids.items()):
        kept_prefix=prefixes[0][1]
        skipped_prefixes=[prefix for fastq_dir,prefix in prefixes[1:]]
        print(f'DUPLICATE: {sample_id}; keep {kept_prefix}; skip {len(skipped_prefixes)} files: {"; ".join(skipped_prefixes)}')
    for fastq_dir in fastq_dirs:
        sample_ids=dir_sample_ids[fastq_dir]
        multi_dir_sample_count=0
        for sample_id in sample_ids:
            sample_dirs={sample_fastq_dir for sample_fastq_dir,prefix in sample_prefixes[sample_id]}
            if len(sample_dirs) > 1:
                multi_dir_sample_count+=1
        print(
            f'fastq_dir: {fastq_dir}; '
            f'samples: {dir_sample_counts[fastq_dir]}; '
            f'unique_samples: {len(sample_ids)}; '
            f'multi_dir_samples: {multi_dir_sample_count}'
        )
    k=0
    for sample_id, prefix in sid_prefix.items():
        read_1_fq=f'{prefix}_1.fq.gz'
        read_2_fq=f'{prefix}_2.fq.gz'
        if not (os.path.isfile(read_1_fq) and os.path.isfile(read_2_fq)):
            read_1_fq=f'{prefix}_1.merged.fq.gz'
            read_2_fq=f'{prefix}_2.merged.fq.gz'
        read_add_1_fq=f'{prefix}_1.add.fq.gz'
        read_add_2_fq=f'{prefix}_2.add.fq.gz'
        if not (os.path.isfile(read_1_fq) and os.path.isfile(read_2_fq)):
            print(f'ERROR: {sample_id} ({prefix}) has no fastq files.')
            missing_pair_count+=1
            continue
        if os.path.getsize(read_1_fq) <= MIN_PASS_FASTQ_SIZE or os.path.getsize(read_2_fq) <= MIN_PASS_FASTQ_SIZE:
            print(
                f'SKIP: {sample_id} fq size <= {_format_size(MIN_PASS_FASTQ_SIZE)}; '
                f'R1={_format_size(os.path.getsize(read_1_fq))}; '
                f'R2={_format_size(os.path.getsize(read_2_fq))}'
            )
            size_too_small_count+=1
            continue
        read_1s=[read_1_fq]
        read_2s=[read_2_fq]
        if os.path.isfile(read_add_1_fq) and os.path.isfile(read_add_2_fq):
            read_1s.append(read_add_1_fq)
            read_2s.append(read_add_2_fq)
            k+=1
        sample_files[sample_id]=[' '.join(read_1s),' '.join(read_2s)]
    data=[]

    failed_samples=None
    if os.path.isfile(FAILED_SAMPLES_LIST_PATH):
        failed_samples=[]
        with open(FAILED_SAMPLES_LIST_PATH) as br:
            for line in br:
                failed_samples.append(line.strip())
        print(f'load failed samples: {len(failed_samples)}')

    for sid,fqs in sample_files.items():
        sid='-'.join(['GTOP']+sid.split('-')[1:])
        if must_in_passed_id and sid not in passed_samples:
            continue
        # only failed samples
        if failed_samples is not None and sid not in failed_samples:
            continue

        data.append([sid,fqs[0],fqs[1]])
    df=pd.DataFrame(data=data,columns=['sample_id','fqs1', 'fqs2'])
    print(f'{k} samples with add fq; all {len(sample_files)} samples.')
    print(f'missing paired fq samples: {missing_pair_count}')
    print(f'fq size <= {_format_size(MIN_PASS_FASTQ_SIZE)} samples: {size_too_small_count}')
    print(f'submit {len(data)} samples.')
    return df

def remain_sr_task():
    sid_df=_load_all_unique_fastq_gt_1g()
    fq_ids=set(sid_df['sample_id'].tolist())
    quant_methods=['salmon','rsem','stringtie']
    sr_dir=f'{PROJ_DIR}/path/to/short_read/HPC/output'
    for quant_method in ['rsem']:
        df=pd.read_csv(f'{sr_dir}/{quant_method}/quant/gencode.gene.count.{quant_method}.tsv',sep='\t',
                       index_col=0,nrows=0)
        sid=set(df.columns.tolist())
        # remain=sample_ids-sid
        remain=fq_ids-sid
        remain_df=sid_df.loc[sid_df['sample_id'].isin(remain),:]
        print(f'[{quant_method}]: quant {len(sid)} samples; remain {remain_df.shape[0]} samples.')
        return remain_df

def _load_quant_conf():
    data=[]
    for ref_tag in [REF_TAG]:
        for gene_type in ['transcript','gene']:
            for quant_type in ['count','tpm']:
                data.append([ref_tag,gene_type,quant_type])
    df=pd.DataFrame(data=data,columns=['ref_tag','gene_type','quant_type'])
    print(f'loaded {df.shape[0]} conf files')
    return df

def _check_success_samples(must_in_db=True, ref_key=REF_TAG):
    passed_samples=pd.read_csv(passed_srRNA_sample_path)['sample_id'].tolist()
    base_dir = f'{OUTPUT_DIR}/HPC/output/sample_based'
    success_count = 0
    total_count = 0
    failed_samples = {}
    if not os.path.isdir(base_dir):
        print(f'ERROR: sample base dir not found: {base_dir}')
        return
    if not must_in_db:
        passed_samples = os.listdir(base_dir)
    _ensure_check_dir()

    def _read_header(path):
        with open(path, 'r') as handle:
            line = handle.readline().strip()
        if not line:
            return []
        return line.split('\t')

    def _check_rsem_table(path, required_cols):
        reasons = []
        if not os.path.isfile(path):
            return [f'missing file: {path}']
        if os.path.getsize(path) <= 1024:
            reasons.append(f'file size <= 1KB: {path}')
        try:
            header = _read_header(path)
        except Exception as exc:
            reasons.append(f'read header failed: {path}; {type(exc).__name__}: {exc}')
            return reasons
        missing_cols = [col for col in required_cols if col not in header]
        if missing_cols:
            reasons.append(
                f'missing columns in {path}: {",".join(missing_cols)}; '
                f'valid columns: {",".join(header)}'
            )
        return reasons

    for sample_id in passed_samples:
        sample_dir = f'{base_dir}/{sample_id}'
        # if not os.path.isdir(sample_dir):
        #     continue
        total_count += 1
        sample_reasons = []
        table_checks = [
            (
                f'{sample_dir}/RSEM_{ref_key}.isoforms.results',
                ['transcript_id', 'expected_count', 'TPM']
            ),
            (
                f'{sample_dir}/RSEM_{ref_key}.genes.results',
                ['gene_id', 'expected_count', 'TPM']
            ),
        ]
        for quant_sf, required_cols in table_checks:
            sample_reasons.extend(_check_rsem_table(quant_sf, required_cols))
        if not sample_reasons:
            success_count += 1
        else:
            failed_samples[sample_id] = sample_reasons

    failed_out = FAILED_SAMPLES_LIST_PATH
    with open(failed_out, 'w') as f:
        for sample_id in sorted(failed_samples):
            f.write(sample_id + '\n')
            # out_res_dir=f'{OUTPUT_DIR}/HPC/output/sample_based/{sample_id}'
            # if os.path.isdir(out_res_dir):
            #     shutil.rmtree(out_res_dir)
            #     print(f'rm {sample_id}')

    with open(FAILED_SAMPLES_DETAIL_PATH, 'w') as f:
        f.write('sample_id\treason\n')
        for sample_id in sorted(failed_samples):
            for reason in failed_samples[sample_id]:
                f.write(f'{sample_id}\t{reason}\n')
                print(f'ERROR: {sample_id}: {reason}')

    print(f'total_samples: {total_count}')
    print(f'success_samples: {success_count}')
    print(f'failed_samples: {len(failed_samples)}')
    print(f'failed_sample_id_file: {failed_out}')
    print(f'failed_sample_detail_file: {FAILED_SAMPLES_DETAIL_PATH}')


def _check_rsem_transcript_ids(ref_key=REF_TAG):
    base_dir = f'{OUTPUT_DIR}/HPC/output/sample_based'
    if not os.path.isdir(base_dir):
        print(f'ERROR: sample base dir not found: {base_dir}')
        return

    import polars as pl

    _ensure_check_dir()
    sample_ids = sorted([
        sample_id for sample_id in os.listdir(base_dir)
        if os.path.isdir(f'{base_dir}/{sample_id}')
    ])
    n_threads = min(32, max(1, os.cpu_count() or 1))

    def _scan_sample(sample_id):
        path = f'{base_dir}/{sample_id}/RSEM_{ref_key}.isoforms.results'
        if not os.path.isfile(path):
            return {
                'sample_id': sample_id,
                'path': path,
                'status': 'missing',
                'reason': 'missing file',
                'total_rows': 0,
                'empty_id_rows': 0,
                'bad_counts': {},
                'duplicate_counts': {},
            }
        try:
            df = pl.read_csv(
                path,
                separator='\t',
                has_header=True,
                columns=['transcript_id'],
                schema_overrides={'transcript_id': pl.Utf8},
            ).with_columns(
                pl.col('transcript_id').fill_null('').str.strip_chars()
            )
        except Exception as exc:
            return {
                'sample_id': sample_id,
                'path': path,
                'status': 'error',
                'reason': f'{type(exc).__name__}: {exc}',
                'total_rows': 0,
                'empty_id_rows': 0,
                'bad_counts': {},
                'duplicate_counts': {},
            }

        total_rows = df.height
        empty_id_rows = df.filter(pl.col('transcript_id') == '').height
        bad_df = (
            df
            .filter(
                ~pl.col('transcript_id').str.starts_with('ENST')
                & ~pl.col('transcript_id').str.starts_with('GTOPT')
            )
            .group_by('transcript_id')
            .len()
            .sort('transcript_id')
        )
        duplicate_df = (
            df
            .filter(pl.col('transcript_id') != '')
            .group_by('transcript_id')
            .len()
            .filter(pl.col('len') > 1)
            .sort('transcript_id')
        )
        bad_counts = dict(zip(bad_df['transcript_id'].to_list(), bad_df['len'].to_list()))
        duplicate_counts = dict(zip(duplicate_df['transcript_id'].to_list(), duplicate_df['len'].to_list()))
        return {
            'sample_id': sample_id,
            'path': path,
            'status': 'problem' if bad_counts or duplicate_counts else 'ok',
            'reason': '',
            'total_rows': total_rows,
            'empty_id_rows': empty_id_rows,
            'bad_counts': bad_counts,
            'duplicate_counts': duplicate_counts,
        }

    print(f'checking {len(sample_ids)} samples with {n_threads} threads')
    with ThreadPoolExecutor(max_workers=n_threads) as pool:
        scan_results = list(pool.map(_scan_sample, sample_ids))

    missing_files = sum(1 for result in scan_results if result['status'] == 'missing')
    samples_with_bad_ids = sum(1 for result in scan_results if result['bad_counts'])
    samples_with_duplicate_ids = sum(1 for result in scan_results if result['duplicate_counts'])
    problem_sample_ids = sorted({
        result['sample_id']
        for result in scan_results
        if result['status'] != 'ok'
    })

    with open(BAD_TRANSCRIPT_ID_REPORT_PATH, 'w') as bad_out, \
            open(DUP_TRANSCRIPT_ID_REPORT_PATH, 'w') as dup_out, \
            open(TRANSCRIPT_ID_SUMMARY_REPORT_PATH, 'w') as summary_out, \
            open(FAILED_SAMPLES_LIST_PATH, 'w') as failed_out:
        bad_out.write('sample_id\ttranscript_id\trow_count\tpath\n')
        dup_out.write('sample_id\ttranscript_id\trow_count\tpath\n')
        summary_out.write(
            'sample_id\tstatus\ttotal_rows\tempty_id_rows\tbad_id_count\t'
            'bad_id_row_count\tduplicate_id_count\tduplicate_id_row_count\tpath\treason\n'
        )
        for sample_id in problem_sample_ids:
            failed_out.write(sample_id + '\n')

        for result in scan_results:
            sample_id = result['sample_id']
            path = result['path']
            total_rows = result['total_rows']
            empty_id_rows = result['empty_id_rows']
            bad_counts = result['bad_counts']
            duplicate_counts = result['duplicate_counts']
            bad_id_count = len(bad_counts)
            bad_id_row_count = sum(bad_counts.values())
            duplicate_id_count = len(duplicate_counts)
            duplicate_id_row_count = sum(duplicate_counts.values())
            if result['status'] in {'missing', 'error'}:
                summary_out.write(
                    f'{sample_id}\t{result["status"]}\t0\t0\t0\t0\t0\t0\t{path}\t{result["reason"]}\n'
                )
                continue
            if bad_id_count:
                for transcript_id, count in sorted(bad_counts.items()):
                    bad_out.write(f'{sample_id}\t{transcript_id}\t{count}\t{path}\n')
            if duplicate_id_count:
                for transcript_id, count in sorted(duplicate_counts.items()):
                    dup_out.write(f'{sample_id}\t{transcript_id}\t{count}\t{path}\n')

            summary_out.write(
                f'{sample_id}\t{result["status"]}\t{total_rows}\t{empty_id_rows}\t'
                f'{bad_id_count}\t{bad_id_row_count}\t'
                f'{duplicate_id_count}\t{duplicate_id_row_count}\t{path}\t\n'
            )

    print(f'total_samples: {len(sample_ids)}')
    print(f'missing_files: {missing_files}')
    print(f'failed_samples: {len(problem_sample_ids)}')
    print(f'samples_with_bad_transcript_ids: {samples_with_bad_ids}')
    print(f'samples_with_duplicate_transcript_ids: {samples_with_duplicate_ids}')
    print(f'failed_sample_id_file: {FAILED_SAMPLES_LIST_PATH}')
    print(f'bad_transcript_id_report: {BAD_TRANSCRIPT_ID_REPORT_PATH}')
    print(f'duplicate_transcript_id_report: {DUP_TRANSCRIPT_ID_REPORT_PATH}')
    print(f'summary_report: {TRANSCRIPT_ID_SUMMARY_REPORT_PATH}')


if __name__ == '__main__':
    pipeline_py='rsem_pipeline.py'
    ARGS=sys.argv[1:]
    COMPUTER='HPC'
    STEP=ARGS[0]
    # df=_load_all_fastq()
    ## for step 1: sample-based multiple task.
    if STEP == 'prepare':
        PARTITION='cu-1,cpuPartition,fat-1,cu-short'
        NODES = 1
        CPU_PER_NODE = 30
        NT_PER_TASK = 30
        N_TASK_PER_NODE = 1
        submit_job(None,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE,pipeline_py)
    ## for step 1: sample-based multiple task.
    if STEP == 'run':
        PARTITION='cu-1,fat-1,cu-short'
        NODES = 400
        CPU_PER_NODE = 30
        NT_PER_TASK = 30
        N_TASK_PER_NODE = 1
        df=_load_all_unique_fastq_gt_1g(True)
        submit_job(df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE,pipeline_py)
    if STEP == 'combine':
        PARTITION='cu-1,cpuPartition,fat-1,cu-short'
        NODES = 4
        CPU_PER_NODE = 24
        NT_PER_TASK = 24
        N_TASK_PER_NODE = 1
        df=_load_quant_conf()
        submit_job(df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE,pipeline_py)

    if STEP == 'check':
        _check_success_samples(True)
    if STEP == 'check_tx_ids':
        _check_rsem_transcript_ids()
