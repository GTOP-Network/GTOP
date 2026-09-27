# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/11/23 21:43
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""
import os
import shlex
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

PROJ_DIR = os.environ.get("PROJECT_ROOT")
LOAD_BASE_ENV_CMD='source ~/.bashrc'

RUN_LOG_DIR=f'{PROJ_DIR}/project/GTOP-RNA/20260815/run_log'
OUTPUT_DIR=f'{PROJ_DIR}/project/GTOP-RNA/20260815/output/SRS/quantification'
# all samples
passed_srRNA_sample_path=f'{PROJ_DIR}/raw_data/GMTiP/meta/RNA/SRS_passed_sample_id.csv'
# add skin samples
# passed_srRNA_sample_path=f'{PROJ_DIR}/raw_data/GMTiP/meta/RNA/SRS_add_skin_sample_ids.csv'

# REF_TAG='DSA_to_hg38'
# REF_TAG='hg38'
REF_TAG='enhanced'

## real data
fastq_dirs=[
    '/flashfs1/scratch.global/cxue/raw_data/GTOP/RNA/SRS/fastq',
    '/flashfs1/scratch.global/cxue/raw_data/GTOP/RNA/SRS/fastq_merged',
    '/lustre/home/xdzou/2024-10-21-GTBMap/2025-02-11-RNA-mapping/input/fastq',
]

## downsampling data
fastq_dirs=[
    '/flashfs1/scratch.global/lhgong/mengxin/output_80m'
]

RUN_LOG_DIR=f'{PROJ_DIR}/project/GTOP-RNA/20260815/run_log/SRS_downsampling'
OUTPUT_DIR=f'{PROJ_DIR}/project/GTOP-RNA/20260815/output/SRS_downsampling_quant'


NASFS1_FASTQ_COPY_DIR='/flashfs1/scratch.global/cxue/raw_data/GTOP/RNA/SRS/fastq'
COPY_CHUNK_SIZE=64*1024*1024
COPY_WORKERS=4
COPY_PROGRESS_INTERVAL=5
UNPASSED_SRRNA_SAMPLE_TABLE=f'{RUN_LOG_DIR}/HPC/check/unpassed_srRNA_sample_ids.csv'
MIN_PASS_FASTQ_SIZE=1024*1024*1024

LOCAL_WDIR='/media/london_A/xuechao/common_code/sync_data/srsRNA_fq'
LOCAL_UNPASSED_SRRNA_SAMPLE_TABLE=f'{LOCAL_WDIR}/{os.path.basename(UNPASSED_SRRNA_SAMPLE_TABLE)}'
LOCAL_TEMP_LINK_DIR=f'{LOCAL_WDIR}/temp_fq'

SRC_FASTQ_UPLOAD_DIR='/media/london_B/mengxin/GTOP_RNA_FQ'

LFTP_REMOTE='sftp://szbl_hpc'
LFTP_UPLOAD_PARALLEL=4
SRC_FASTQ_SIZE_TABLE=f'{RUN_LOG_DIR}/HPC/check/src_fastq_file_sizes.csv'


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
CURRENT_PY = os.path.splitext(os.path.basename(__file__))[0]
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
        MEM_GB=4,
    ):
    '''

    :param df: task dataframe. if none, use single node.
    :param STEP: keyword matching `iso-seq_pipline`
    :return:
    '''
    # assign task to nodes, build a sample list file as para for pipeline scripts per node.
    output_dir=f'{OUTPUT_DIR}/{COMPUTER}/output'
    if df is not None:
        if NODES > df.shape[0]:
            NODES = df.shape[0]
        sub_dfs = [df.iloc[part] for part in np.array_split(range(len(df)), NODES)]
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
            f'#SBATCH -J Salmon.{STEP}.node{i}',
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


def _load_all_fastq():
    sample_files={}
    size_table,conflict_size_files=_load_src_fastq_size_table()
    size_mismatch_count=0
    size_missing_count=0
    # load all sample fastq files.
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
    for sample_id,prefix in sid_prefix.items():
        read_1_fq=f'{prefix}_1.fq.gz'
        read_2_fq=f'{prefix}_2.fq.gz'
        if not (os.path.isfile(read_1_fq) and os.path.isfile(read_2_fq)):
            read_1_fq=f'{prefix}_1.merged.fq.gz'
            read_2_fq=f'{prefix}_2.merged.fq.gz'
        read_add_1_fq=f'{prefix}_1.add.fq.gz'
        read_add_2_fq=f'{prefix}_2.add.fq.gz'
        if not (os.path.isfile(read_1_fq) and os.path.isfile(read_2_fq)):
            print(f'ERROR: {sample_id} ({prefix}) has no fastq files.')
        else:
            read_1s=[read_1_fq]
            read_2s=[read_2_fq]
            if os.path.isfile(read_add_1_fq) and os.path.isfile(read_add_2_fq):
                read_1s.append(read_add_1_fq)
                read_2s.append(read_add_2_fq)
                k+=1
            size_check_ok=True
            for fq_path in read_1s+read_2s:
                check_status=_check_fastq_size_with_src_table(fq_path,size_table,conflict_size_files)
                if check_status == 'mismatch':
                    size_mismatch_count+=1
                    size_check_ok=False
                elif check_status == 'missing':
                    size_missing_count+=1
            if not size_check_ok:
                print(f'ERROR: {sample_id} skipped because fq size check failed.')
                continue
            sample_files[sample_id]=[' '.join(read_1s),' '.join(read_2s)]
    print(f'{k} samples with add fq; all {len(sample_files)} samples.')
    if size_table is not None:
        print(f'fq size mismatches: {size_mismatch_count}')
        print(f'fq files not found in size table: {size_missing_count}')
    data=[]
    passed_samples=pd.read_csv(passed_srRNA_sample_path)['sample_id'].tolist()
    print(f'load {len(passed_samples)} passed samples. ')
    # passed_samples=[f'GTOP{x[5:]}' for x in passed_samples]
    # passed_samples=[]
    check_dir = f'{RUN_LOG_DIR}/HPC/check'
    failed_out = f'{check_dir}/failed_sample_ids.txt'
    if os.path.isfile(failed_out):
        with open(failed_out, 'r') as f:
            for line in f:
                passed_samples.append(line.strip())
        for sid,fqs in sample_files.items():
            if sid in passed_samples:
                data.append([sid,fqs[0],fqs[1]])
    else:
        for sid,fqs in sample_files.items():
            # temp code start
            # if os.path.isdir(f'{result_dir}/{sid}'):
            #     continue
            # temp code end
            data.append([sid,fqs[0],fqs[1]])
    df=pd.DataFrame(data=data,columns=['sample_id','fqs1', 'fqs2'])
    print(f'submit {len(data)} samples.')
    return df

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
    for sid,fqs in sample_files.items():
        sid='-'.join(['GTOP']+sid.split('-')[1:])
        if must_in_passed_id and sid not in passed_samples:
            continue
        data.append([sid,fqs[0],fqs[1]])
    df=pd.DataFrame(data=data,columns=['sample_id','fqs1', 'fqs2'])
    print(f'{k} samples with add fq; all {len(sample_files)} samples.')
    print(f'missing paired fq samples: {missing_pair_count}')
    print(f'fq size <= {_format_size(MIN_PASS_FASTQ_SIZE)} samples: {size_too_small_count}')
    print(f'submit {len(data)} samples.')
    return df

def _write_src_fastq_size_table(
        out_path=SRC_FASTQ_SIZE_TABLE,
        src_dir=SRC_FASTQ_UPLOAD_DIR,
    ):
    data=[]
    for root, dirs, files in os.walk(src_dir):
        for file in files:
            if not file.endswith('.fq.gz'):
                continue
            fq_path=f'{root}/{file}'
            real_path=os.path.realpath(fq_path)
            raw_sample_id='_'.join(file.split('_')[:-1])
            sample_id=_normalize_srRNA_sample_id(raw_sample_id)
            read_type=_get_fastq_read_type(file)
            data.append([
                sample_id,
                raw_sample_id,
                file,
                read_type,
                fq_path,
                real_path,
                os.path.getsize(real_path),
                int(os.path.islink(fq_path)),
            ])
    df=pd.DataFrame(
        data=data,
        columns=['sample_id','raw_sample_id','fq_file','read_type','src_path','real_path','real_size','is_symlink']
    )
    os.makedirs(os.path.dirname(out_path),exist_ok=True)
    df.to_csv(out_path,index=False)
    print(f'src fq files: {df.shape[0]}')
    print(f'src fq samples: {df["sample_id"].nunique()}')
    print(f'src symlink fq files: {df["is_symlink"].sum()}')
    print(f'src fq total real size: {_format_size(df["real_size"].sum())}')
    print(f'src fq size table: {out_path}')
    return df

def _load_src_fastq_size_table(size_table_path=SRC_FASTQ_SIZE_TABLE):
    if not os.path.isfile(size_table_path):
        print(f'WARNING: src fastq size table not found, skip size check: {size_table_path}')
        return None,set()
    df=pd.read_csv(size_table_path)
    size_table={}
    conflict_files=set()
    for fq_file,sdf in df.groupby('fq_file'):
        sizes=set(sdf['real_size'].astype(int).tolist())
        if len(sizes) == 1:
            size_table[fq_file]=next(iter(sizes))
        else:
            conflict_files.add(fq_file)
    print(f'loaded src fastq size table: {size_table_path}')
    print(f'src fq size records: {df.shape[0]}')
    print(f'src fq unique basenames: {len(size_table)}')
    print(f'src fq basename size conflicts: {len(conflict_files)}')
    return size_table,conflict_files

def _check_fastq_size_with_src_table(fq_path,size_table,conflict_size_files):
    if size_table is None:
        return 'skipped'
    fq_file=os.path.basename(fq_path)
    if fq_file in conflict_size_files:
        print(f'ERROR: {fq_file} has conflicting sizes in src size table.')
        return 'mismatch'
    if fq_file not in size_table:
        print(f'WARNING: {fq_file} not found in src size table; size check skipped.')
        return 'missing'
    local_size=os.path.getsize(fq_path)
    expected_size=size_table[fq_file]
    if local_size != expected_size:
        print(f'ERROR: fq size mismatch: {fq_path}; local={local_size}; src={expected_size}')
        return 'mismatch'
    return 'match'

def _normalize_srRNA_sample_id(sample_id):
    if '-' in sample_id:
        return sample_id.split('-',1)[1]
    return sample_id

def _get_fastq_read_type(file):
    if file.endswith('_1.fq.gz'):
        return 'read1'
    if file.endswith('_2.fq.gz'):
        return 'read2'
    if file.endswith('_1.merged.fq.gz'):
        return 'read1_merged'
    if file.endswith('_2.merged.fq.gz'):
        return 'read2_merged'
    if file.endswith('_1.add.fq.gz'):
        return 'read1_add'
    if file.endswith('_2.add.fq.gz'):
        return 'read2_add'
    return ''

def _load_passed_srRNA_sample_ids():
    passed_df=pd.read_csv(passed_srRNA_sample_path)
    passed_samples=[]
    for sample_id in passed_df['sample_id'].astype(str).tolist():
        if sample_id.startswith('GTOP'):
            passed_samples.append(_normalize_srRNA_sample_id(sample_id))
        elif len(sample_id) > 5:
            passed_samples.append(_normalize_srRNA_sample_id(f'GTOP{sample_id[5:]}'))
    return set(passed_samples)

def _write_unpassed_srRNA_sample_table(out_path=UNPASSED_SRRNA_SAMPLE_TABLE):
    existing_files={}
    sample_dirs={}
    for fastq_dir in fastq_dirs:
        for root, dirs, files in os.walk(fastq_dir):
            for file in files:
                if not file.endswith('.fq.gz'):
                    continue
                raw_sample_id='_'.join(file.split('_')[:-1])
                sample_id=_normalize_srRNA_sample_id(raw_sample_id)
                read_type=_get_fastq_read_type(file)
                if read_type == '':
                    continue
                existing_files.setdefault(sample_id,{})[read_type]={
                    'raw_sample_id':raw_sample_id,
                    'fq_file':file,
                    'size':os.path.getsize(f'{root}/{file}'),
                }
                sample_dirs.setdefault(sample_id,set()).add(fastq_dir)

    sample_files={}
    for sample_id,files in existing_files.items():
        if (
                'read1' in files and
                'read2' in files and
                files['read1']['size'] > MIN_PASS_FASTQ_SIZE and
                files['read2']['size'] > MIN_PASS_FASTQ_SIZE
        ):
            sample_files[sample_id]=[files['read1']['fq_file'],files['read2']['fq_file']]
            continue
        if (
                'read1_merged' in files and
                'read2_merged' in files and
                files['read1_merged']['size'] > MIN_PASS_FASTQ_SIZE and
                files['read2_merged']['size'] > MIN_PASS_FASTQ_SIZE
        ):
            sample_files[sample_id]=[files['read1_merged']['fq_file'],files['read2_merged']['fq_file']]

    passed_samples=_load_passed_srRNA_sample_ids()
    data=[]
    missing_sample_ids=set()
    for sample_id in sorted(passed_samples):
        if sample_id in sample_files:
            continue
        files=existing_files.get(sample_id,{})
        missing_sample_ids.add(sample_id)
        expected_files=[]
        if 'read1' in files or 'read2' in files:
            expected_files.extend([
                (f'{sample_id}_1.fq.gz','read1'),
                (f'{sample_id}_2.fq.gz','read2'),
            ])
        elif 'read1_merged' in files or 'read2_merged' in files:
            expected_files.extend([
                (f'{sample_id}_1.merged.fq.gz','read1_merged'),
                (f'{sample_id}_2.merged.fq.gz','read2_merged'),
            ])
        else:
            expected_files.extend([
                (f'{sample_id}_1.fq.gz','read1'),
                (f'{sample_id}_2.fq.gz','read2'),
            ])
        for fq_file,read_type in expected_files:
            existing_size=files[read_type]['size'] if read_type in files else 0
            fail_reason='missing'
            if read_type in files and existing_size <= MIN_PASS_FASTQ_SIZE:
                fail_reason='size_le_1G'
            if read_type in files and existing_size > MIN_PASS_FASTQ_SIZE:
                continue
            data.append([
                sample_id,
                fq_file,
                read_type,
                existing_size,
                fail_reason,
                ';'.join(sorted(sample_dirs.get(sample_id,set()))),
            ])
    df=pd.DataFrame(
        data=data,
        columns=['sample_id','fq_file','read_type','existing_size','fail_reason','existing_fastq_dirs']
    )
    os.makedirs(os.path.dirname(out_path),exist_ok=True)
    df.to_csv(out_path,index=False)
    print(f'passed srRNA sample IDs: {len(passed_samples)}')
    print(f'sample_files sample IDs: {len(sample_files)}')
    print(f'passed sample IDs missing in sample_files: {len(missing_sample_ids)}')
    print(f'missing fq files: {df.shape[0]}')
    print(f'min pass fq size: {_format_size(MIN_PASS_FASTQ_SIZE)}')
    print(f'missing fq table: {out_path}')
    return df


def _make_upload_temp_link():
    os.makedirs(LOCAL_TEMP_LINK_DIR,exist_ok=True)
    df=pd.read_csv(LOCAL_UNPASSED_SRRNA_SAMPLE_TABLE)
    # don't consider add fa
    for fq_file in df['fq_file']:
        has_fq=False
        for prefix in ['AGTEX','GTOP']:
            src_fq=f'{SRC_FASTQ_UPLOAD_DIR}/{prefix}-{fq_file}'
            if os.path.exists(src_fq):
                has_fq=True
                # add soft link to temp dir
                target = os.readlink(src_fq)
                dst=f'{LOCAL_TEMP_LINK_DIR}/{prefix}-{fq_file}'
                os.symlink(target, dst)
                break

        if not has_fq:
            print(f'WARNING: no {fq_file} in {SRC_FASTQ_UPLOAD_DIR}')

        pass


def _format_size(size):
    units=['B','KB','MB','GB','TB']
    value=float(size)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f'{value:.2f}{unit}'
        value/=1024

def _update_copy_progress(progress_state,progress_lock,src_path,copied_size,total_size,speed,avg_speed,eta):
    if progress_state is None:
        return
    with progress_lock:
        progress_state[src_path]={
            'copied_size':copied_size,
            'total_size':total_size,
            'speed':speed,
            'avg_speed':avg_speed,
            'eta':eta,
            'basename':os.path.basename(src_path),
        }

def _print_copy_progress(progress_state,copy_totals,progress_lock,stop_event,total_size,total_files,interval=COPY_PROGRESS_INTERVAL):
    while not stop_event.is_set():
        with progress_lock:
            active_items=list(progress_state.values())
            completed_size=copy_totals['completed_size']
            completed_files=copy_totals['completed_files']
        copied_size=completed_size+sum(item['copied_size'] for item in active_items)
        copied_files=completed_files+len(active_items)
        speed=sum(item['speed'] for item in active_items)
        percent=copied_size/total_size*100 if total_size > 0 else 100
        active_files=', '.join(item['basename'] for item in active_items[:3])
        if len(active_items) > 3:
            active_files+=f' +{len(active_items)-3}'
        print(
            f'\rCOPY_PROGRESS files={copied_files}/{total_files} '
            f'{percent:6.2f}% copied={_format_size(copied_size)}/{_format_size(total_size)} '
            f'speed={_format_size(speed)}/s active=[{active_files}]',
            end='',
            flush=True
        )
        time.sleep(interval)
    print()

def _copy_file_with_progress(src_path,target_path,chunk_size=COPY_CHUNK_SIZE,progress_state=None,progress_lock=None):
    src_size=os.path.getsize(src_path)
    part_path=f'{target_path}.part'
    copied_size=0
    if os.path.exists(target_path):
        target_size=os.path.getsize(target_path)
        if target_size == src_size:
            return 'skipped'
        return 'failed'
    if os.path.exists(part_path):
        part_size=os.path.getsize(part_path)
        if part_size < src_size:
            copied_size=part_size
        else:
            os.remove(part_path)

    start_time=time.time()
    last_update_time=start_time
    last_update_size=copied_size
    _update_copy_progress(progress_state,progress_lock,src_path,copied_size,src_size,0,0,0)
    mode='ab' if copied_size > 0 else 'wb'
    with open(src_path,'rb') as src_f, open(part_path,mode) as dst_f:
        if copied_size > 0:
            src_f.seek(copied_size)
        while copied_size < src_size:
            chunk=src_f.read(chunk_size)
            if not chunk:
                break
            dst_f.write(chunk)
            copied_size+=len(chunk)
            now=time.time()
            if now-last_update_time >= COPY_PROGRESS_INTERVAL or copied_size == src_size:
                interval=max(now-last_update_time,1e-6)
                speed=(copied_size-last_update_size)/interval
                total_interval=max(now-start_time,1e-6)
                avg_speed=copied_size/total_interval
                remain_size=max(src_size-copied_size,0)
                eta=remain_size/avg_speed if avg_speed > 0 else 0
                _update_copy_progress(progress_state,progress_lock,src_path,copied_size,src_size,speed,avg_speed,eta)
                last_update_time=now
                last_update_size=copied_size
        dst_f.flush()
        os.fsync(dst_f.fileno())
    if progress_state is not None:
        with progress_lock:
            progress_state.pop(src_path,None)
    part_size=os.path.getsize(part_path)
    if part_size != src_size:
        print(f'ERROR: incomplete copy: {part_path}; source={src_size}; copied={part_size}')
        return 'failed'
    os.replace(part_path,target_path)
    shutil.copystat(src_path,target_path)
    return 'copied'

def _copy_nasfs1_only_fastq(output_dir=NASFS1_FASTQ_COPY_DIR,workers=COPY_WORKERS):
    sample_dirs={}
    nasfs1_sample_files={}
    for fastq_dir in fastq_dirs:
        is_nasfs1_dir=fastq_dir.startswith('/nasfs1/')
        for root, dirs, files in os.walk(fastq_dir):
            for file in files:
                if not file.endswith('.fq.gz'):
                    continue
                sample_id='_'.join(file.split('_')[:-1])
                sample_dirs.setdefault(sample_id,set()).add(fastq_dir)
                if is_nasfs1_dir:
                    fq_path=f'{root}/{file}'
                    nasfs1_sample_files.setdefault(sample_id,[]).append(fq_path)

    os.makedirs(output_dir,exist_ok=True)
    copy_items=[]
    copy_sample_ids=set()
    seen_target_paths=set()
    duplicate_target_count=0
    for sample_id,fq_paths in sorted(nasfs1_sample_files.items()):
        non_nasfs1_dirs=[
            fastq_dir
            for fastq_dir in sample_dirs[sample_id]
            if not fastq_dir.startswith('/nasfs1/')
        ]
        if non_nasfs1_dirs:
            continue
        for fq_path in sorted(set(fq_paths)):
            target_path=f'{output_dir}/{os.path.basename(fq_path)}'
            if target_path in seen_target_paths:
                duplicate_target_count+=1
                continue
            seen_target_paths.add(target_path)
            copy_items.append((sample_id,fq_path,target_path))
            copy_sample_ids.add(sample_id)

    total_copy_size=sum(os.path.getsize(fq_path) for sample_id,fq_path,target_path in copy_items)
    print(f'nasfs1-only samples to copy: {len(copy_sample_ids)}')
    print(f'fq files to copy/check: {len(copy_items)}')
    print(f'total fq size: {_format_size(total_copy_size)}')
    print(f'copy output dir: {output_dir}')
    print(f'copy workers: {workers}')
    print(f'duplicate target files skipped: {duplicate_target_count}')

    copy_count=0
    skip_exists_count=0
    failed_count=0
    progress_state={}
    copy_totals={'completed_size':0,'completed_files':0}
    progress_lock=threading.Lock()
    stop_event=threading.Event()
    progress_thread=threading.Thread(
        target=_print_copy_progress,
        args=(progress_state,copy_totals,progress_lock,stop_event,total_copy_size,len(copy_items)),
        daemon=True
    )
    progress_thread.start()
    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_item={
            executor.submit(_copy_file_with_progress,fq_path,target_path,COPY_CHUNK_SIZE,progress_state,progress_lock):(sample_id,fq_path,target_path)
            for sample_id,fq_path,target_path in copy_items
        }
        for future in as_completed(future_to_item):
            sample_id,fq_path,target_path=future_to_item[future]
            try:
                status=future.result()
            except Exception as e:
                print(f'ERROR: copy failed: {fq_path} -> {target_path}; {e}')
                status='failed'
            if status == 'copied':
                copy_count+=1
            elif status == 'skipped':
                skip_exists_count+=1
            else:
                failed_count+=1
            with progress_lock:
                copy_totals['completed_size']+=os.path.getsize(fq_path)
                copy_totals['completed_files']+=1
    stop_event.set()
    progress_thread.join()
    print(
        f'\rCOPY_PROGRESS files={len(copy_items)}/{len(copy_items)} '
        f'100.00% copied={_format_size(total_copy_size)}/{_format_size(total_copy_size)} '
        f'speed=0.00B/s active=[]',
        flush=True
    )
    print(f'nasfs1-only samples copied: {len(copy_sample_ids)}')
    print(f'fq files copied: {copy_count}')
    print(f'fq files skipped because target exists: {skip_exists_count}')
    print(f'fq files failed: {failed_count}')
    print(f'copy output dir: {output_dir}')

def _load_quant_conf():
    data=[]
    for ref_tag in [REF_TAG]:
        for quant_type in ['count','tpm']:
            data.append([ref_tag,quant_type])
    df=pd.DataFrame(data=data,columns=['ref_tag','quant_type'])
    print(f'loaded {df.shape[0]} conf files')
    return df

def _check_success_samples(ref_key=REF_TAG):
    passed_samples=pd.read_csv(passed_srRNA_sample_path)['sample_id'].tolist()
    base_dir = f'{OUTPUT_DIR}/HPC/output/sample_based'
    success_count = 0
    total_count = 0
    failed_sample_ids = []
    if not os.path.isdir(base_dir):
        print(f'ERROR: sample base dir not found: {base_dir}')
        return
    for sample_id in passed_samples:
        sample_dir = f'{base_dir}/{sample_id}'
        # if not os.path.isdir(sample_dir):
        #     continue
        total_count += 1
        quant_sf = f'{sample_dir}/salmon_{ref_key}/quant.sf'
        # quant_sf = f'{sample_dir}/RSEM_{ref_key}.genes.results'
        if os.path.isfile(quant_sf) and os.path.getsize(quant_sf) > 1024:
            success_count += 1
        else:
            failed_sample_ids.append(sample_id)

    check_dir = f'{RUN_LOG_DIR}/HPC/check'
    os.makedirs(check_dir, exist_ok=True)
    failed_out = f'{check_dir}/failed_sample_ids.txt'
    with open(failed_out, 'w') as f:
        for sample_id in sorted(failed_sample_ids):
            f.write(sample_id + '\n')

    print(f'total_samples: {total_count}')
    print(f'success_samples(quant.sf > 1KB): {success_count}')
    print(f'failed_samples: {len(failed_sample_ids)}')
    print(f'failed_sample_id_file: {failed_out}')


def _common_check_success_samples(ref_key=REF_TAG):
    base_dir = f'{OUTPUT_DIR}/HPC/output/sample_based'
    success_count = 0
    if not os.path.isdir(base_dir):
        print(f'ERROR: sample base dir not found: {base_dir}')
        return

    for sample_id in sorted(os.listdir(base_dir)):
        sample_dir = f'{base_dir}/{sample_id}'
        if not os.path.isdir(sample_dir):
            continue
        quant_sf = f'{sample_dir}/salmon_{ref_key}/quant.sf'
        if os.path.isfile(quant_sf) and os.path.getsize(quant_sf) > 1024:
            success_count += 1
            # print(sample_id)
        else:
            pass
            print(sample_id)

    print(f'success_samples(quant.sf > 1KB): {success_count}')


if __name__ == '__main__':
    pipeline_py='salmon_pipeline.py'
    ARGS=sys.argv[1:]
    COMPUTER='HPC'
    STEP=ARGS[0]
    partition = 'cu-1,cpuPartition,fat-1,cu-short'
    if STEP == 'prepare':
        PARTITION=partition
        NODES = 1
        CPU_PER_NODE = 30
        NT_PER_TASK = 30
        N_TASK_PER_NODE = 1
        submit_job(None,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE,pipeline_py)
    if STEP == 'run':
        PARTITION=partition
        NODES = 500
        CPU_PER_NODE = 10
        NT_PER_TASK = 10
        N_TASK_PER_NODE = 1
        df = _load_all_unique_fastq_gt_1g(False)
        # print(df)
        submit_job(df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE,pipeline_py)
    if STEP == 'combine':
        PARTITION=partition
        NODES = 2
        CPU_PER_NODE = 30
        NT_PER_TASK = 30
        N_TASK_PER_NODE = 1
        df=_load_quant_conf()
        submit_job(df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE,pipeline_py)
    if STEP == 'merge_gene':
        PARTITION=partition
        NODES = 2
        CPU_PER_NODE = 30
        NT_PER_TASK = 30
        N_TASK_PER_NODE = 2
        df=_load_quant_conf()
        submit_job(df,STEP,COMPUTER,PARTITION,NODES,CPU_PER_NODE,NT_PER_TASK,N_TASK_PER_NODE,pipeline_py)
    if STEP == 'check':
        _check_success_samples()
    if STEP == 'common_check':
        _common_check_success_samples()
    if STEP == 'copy_nasfs1_fastq':
        workers=COPY_WORKERS
        if len(ARGS) > 1:
            workers=max(1,int(ARGS[1]))
        _copy_nasfs1_only_fastq(workers=workers)
    if STEP == 'write_unpassed_srRNA_sample_table':
        _write_unpassed_srRNA_sample_table()

    # run in local machine.
    if STEP == 'upload_miss_srRNA_fq':
        _make_upload_temp_link()

