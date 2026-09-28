# -*- coding: utf-8 -*-
"""
@Desc    : Filter rules for alternative first exon:
    1. Novel transcripts with new first exons should only be accepted if the exon is at least 30 bases long
    2. and has no more than one SNV/mismatch per 30 nucleotides.
"""

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_EXCLUDED_SIDS = ['GTOP-CA241-5032-LN-YGE8', 'GTOP-BI281-0087-LN-L7VG', 'GTOP-CB271-4155-LN-A5YE']

import re

EXCLUDED_SIDS = CONFIG_EXCLUDED_SIDS

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR=f'{CONFIG_PROJECT_DIR}/output/assembly/LRS/tx_merge/merged_tools/sqanti3_filter_1'
sample_based_dir = f'{CONFIG_PROJECT_DIR}/output/assembly/LRS/sample_based'


PARTITION = os.environ.get("SLURM_PARTITION", 'fat-1,cpuPartition,cu-1')

CPU_PER_NODE = 2
NT_PER_TASK = 2
N_TASK_PER_NODE = 1

def submit_scan_task():
    STEP='scan'
    py_name='alt_novel_first_exon_qc.py'
    root_dir=f'{OUTPUT_DIR}/AF_filter/scan_results'
    wdir = sample_based_dir
    bams=[]
    n=0
    for f in os.listdir(wdir):
        sid=f
        if sid.startswith('GTOP') and sid not in EXCLUDED_SIDS:
            n+=1
            node_dir=f'{root_dir}/{sid}'
            os.makedirs(node_dir,exist_ok=True)
            bam_list_path=f'{node_dir}/bam_list.txt'
            with open(bam_list_path,'w') as bw:
                bw.write(f'{wdir}/{f}/hg38/{f}.flnc_minimap2.bam\n')
            print(f'save {len(bams)} bam to {bam_list_path}')
            # build hpc job submit list
            hpc_log_path = f'{node_dir}/job_log'
            hpc_job_shell = f'{node_dir}/job_shell.sh'
            shell_lines = [
                f'#!/bin/bash',
                f'#SBATCH -J {STEP}.{sid}',
                f'#SBATCH -o {hpc_log_path}.out',
                f'#SBATCH -e {hpc_log_path}.err',
                f'#SBATCH -p {PARTITION} -N 1 -n {CPU_PER_NODE}',
                f'cd {OUTPUT_DIR}/AF_filter && python {CURRENT_DIR}/{py_name} '
                f'  --log-file {node_dir}/program.log '
                  f'scan '
                  f'--regions novel_first_exon.regions.tsv  '
                  f'--bam-list {bam_list_path}  '
                 f' --out-prefix {node_dir}/out '
                  f'--sample-processes {NT_PER_TASK} '
                f'--merge-gap 10000 '
                 f' --min-mapq 0 '
            ]
            with open(hpc_job_shell, "w") as f:
                for line in shell_lines:
                    f.write(line + "\n")
            cmd = f'''
                sbatch {hpc_job_shell}
            '''
            cmd = ' '.join(cmd.split())
            print(cmd)
            __import__('subprocess').run(cmd, shell=True, executable='/bin/bash', check=True)
            bams.append(cmd)
    print(f'submit {len(bams)} tasks')


if __name__ == '__main__':
    submit_scan_task()
