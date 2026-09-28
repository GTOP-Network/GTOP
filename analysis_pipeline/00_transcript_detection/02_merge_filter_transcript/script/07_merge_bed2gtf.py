# -*- coding: utf-8 -*-

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_TOOLS = ['bambu', 'flair', 'flames', 'isoquant', 'isoseq', 'isotools', 'talon']

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

def is_file_valid(file_path: str) -> bool:
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            first_line = f.readline()
        return bool(first_line.strip())
    except (FileNotFoundError, PermissionError):
        return False

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

TOOLS = CONFIG_TOOLS

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

def bed12_list_to_gtf(bed_files_path, output_gtf):
    """
    Convert a list of BED12 files to one GTF file.

    Parameters
    ----------
    bed_files : list[str]
        BED12 file path list. BED name field should be like: gene_id;transcript_id
        Example: G1;G1.1

    output_gtf : str
        Output GTF path.

    Notes
    -----
    BED is 0-based half-open.
    GTF is 1-based closed.

    Because the same gene_id/transcript_id in different chromosome/strand BED files
    may not represent the same transcript, this function remaps IDs globally using:
        original_gene_id + original_transcript_id + chrom + strand

    Output IDs will be:
        G1, G2, ...
        G1.1, G1.2, ...
    """

    bed_files=[]
    with open(bed_files_path,'r') as br:
        for line in br:
            bed_files.append(line.strip())

    gene_id_map = {}
    transcript_id_map = {}
    gene_transcript_count = {}

    global_gene_index = 0

    def parse_bed12_line(line, bed_path, line_no):
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 12:
            raise ValueError(
                f"Invalid BED12 line: {bed_path}:{line_no}, "
                f"expected at least 12 columns, got {len(fields)}"
            )

        chrom = fields[0]
        chrom_start = int(fields[1])
        chrom_end = int(fields[2])
        name = fields[3]
        score = fields[4]
        strand = fields[5]
        thick_start = int(fields[6])
        thick_end = int(fields[7])
        item_rgb = fields[8]
        block_count = int(fields[9])
        block_sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
        block_starts = [int(x) for x in fields[11].rstrip(",").split(",") if x]

        if len(block_sizes) != block_count:
            raise ValueError(
                f"Invalid BED12 line: {bed_path}:{line_no}, "
                f"blockCount={block_count}, but blockSizes has {len(block_sizes)} values"
            )

        if len(block_starts) != block_count:
            raise ValueError(
                f"Invalid BED12 line: {bed_path}:{line_no}, "
                f"blockCount={block_count}, but blockStarts has {len(block_starts)} values"
            )

        if ";" not in name:
            raise ValueError(
                f"Invalid BED name field: {bed_path}:{line_no}, "
                f"expected gene_id;transcript_id, got {name!r}"
            )

        old_gene_id, old_transcript_id = name.split(";", 1)

        return {
            "chrom": chrom,
            "chrom_start": chrom_start,
            "chrom_end": chrom_end,
            "name": name,
            "score": score,
            "strand": strand,
            "thick_start": thick_start,
            "thick_end": thick_end,
            "item_rgb": item_rgb,
            "block_count": block_count,
            "block_sizes": block_sizes,
            "block_starts": block_starts,
            "old_gene_id": old_gene_id,
            "old_transcript_id": old_transcript_id,
        }

    def get_new_ids(record):
        nonlocal global_gene_index

        old_gene_id = record["old_gene_id"]
        old_transcript_id = record["old_transcript_id"]
        chrom = record["chrom"]
        strand = record["strand"]

        # 关键点：
        # 同一个 G1;G1.1 如果来自不同染色体/不同链，要视为不同转录本/基因
        gene_key = (chrom, strand, old_gene_id)
        transcript_key = (chrom, strand, old_gene_id, old_transcript_id)

        if gene_key not in gene_id_map:
            global_gene_index += 1
            new_gene_id = f"G{global_gene_index}"
            gene_id_map[gene_key] = new_gene_id
            gene_transcript_count[new_gene_id] = 0

        new_gene_id = gene_id_map[gene_key]

        if transcript_key not in transcript_id_map:
            gene_transcript_count[new_gene_id] += 1
            transcript_index = gene_transcript_count[new_gene_id]
            new_transcript_id = f"{new_gene_id}.{transcript_index}"
            transcript_id_map[transcript_key] = new_transcript_id

        return new_gene_id, transcript_id_map[transcript_key]

    def gtf_attributes(gene_id, transcript_id):
        return f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'

    with open(output_gtf, "w") as out:
        for bed_path in bed_files:
            with open(bed_path) as bed:
                for line_no, line in enumerate(bed, start=1):
                    if not line.strip() or line.startswith("#"):
                        continue

                    record = parse_bed12_line(line, bed_path, line_no)
                    gene_id, transcript_id = get_new_ids(record)

                    chrom = record["chrom"]
                    strand = record["strand"]
                    score = record["score"]
                    chrom_start = record["chrom_start"]
                    chrom_end = record["chrom_end"]
                    block_sizes = record["block_sizes"]
                    block_starts = record["block_starts"]

                    source = "bed12"
                    attributes = gtf_attributes(gene_id, transcript_id)

                    # transcript feature
                    transcript_start = chrom_start + 1
                    transcript_end = chrom_end

                    out.write(
                        "\t".join(
                            [
                                chrom,
                                source,
                                "transcript",
                                str(transcript_start),
                                str(transcript_end),
                                score,
                                strand,
                                ".",
                                attributes,
                            ]
                        )
                        + "\n"
                    )

                    # exon features
                    for exon_number, (block_size, block_start) in enumerate(
                        zip(block_sizes, block_starts),
                        start=1,
                    ):
                        exon_start = chrom_start + block_start + 1
                        exon_end = chrom_start + block_start + block_size

                        exon_attributes = (
                            f'gene_id "{gene_id}"; '
                            f'transcript_id "{transcript_id}"; '
                            f'exon_number "{exon_number}";'
                        )

                        out.write(
                            "\t".join(
                                [
                                    chrom,
                                    source,
                                    "exon",
                                    str(exon_start),
                                    str(exon_end),
                                    score,
                                    strand,
                                    ".",
                                    exon_attributes,
                                ]
                            )
                            + "\n"
                        )


def run():
    tools=TOOLS
    n_task=len(tools)
    main_dir = f'{CONFIG_PROJECT_DIR}/output/assembly/LRS/tx_merge/tool_based'
    main_log=f'{main_dir}/log/internal/merge_bed2gtf.log'
    parallel_cmds = []
    for tool in tools:
        wdir=f'{main_dir}/{tool}/merged/gtf'
        os.makedirs(wdir,exist_ok=True)
        beds_path=f'{wdir}/beds_path.txt'
        with open(beds_path, 'w') as bw:
            chr_dir=f'{main_dir}/{tool}/merged/chr_bed'
            for chr in os.listdir(chr_dir):
                bed_path=f'{chr_dir}/{chr}/merged/merged.bed'
                if is_file_valid(bed_path):
                    bw.write(f'{bed_path}\n')
                elif is_file_valid(f'{chr_dir}/{chr}/merged/filelist.txt'):
                    raise FileNotFoundError(f'Incomplete cross-tissue merge: {bed_path}')

        output_gtf=f'{wdir}/merged.gtf'
        series_cmds = []
        log_path = f'{wdir}/run_merge_bed2gtf.log'
        cmd = f'''
                python {CURRENT_DIR}/07_merge_bed2gtf.py merge_bed
                {beds_path} {output_gtf} 
                '''
        cmd = ' '.join(cmd.split())
        series_cmds.append({'cmd': cmd, 'log_path': log_path})
        parallel_cmds.append(series_cmds)
    run_commands_threadpool(parallel_cmds, main_log=main_log, max_workers=n_task)


if __name__ == '__main__':
    args = sys.argv[1:]
    if args and args[0] == 'merge_bed':
        _, beds_path, output_gtf = args
        bed12_list_to_gtf(beds_path, output_gtf)
    else:
        run()
