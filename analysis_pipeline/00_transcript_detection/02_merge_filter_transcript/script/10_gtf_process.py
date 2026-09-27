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

import logging
import re
from collections import defaultdict


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
n_tool_filter_gtf=f'{main_dir}/merged_tools/sqanti3_filter_1/filtered.gtf'

def make_fa():
    cmd=f'{gffread_PATH} {n_tool_filter_gtf} -g {REF_GENOME_FASTA} -w {n_tool_filter_gtf[:-4]}.fa'
    print(cmd)
    __import__('subprocess').run(cmd, shell=True, executable='/bin/bash', check=True)


def parse_gtf_attributes(attr_string):
    """
    解析 GTF 第 9 列 attributes
    例如: gene_id "GENE1"; transcript_id "TX1";
    """
    attrs = {}
    for key, value in re.findall(r'(\S+)\s+"([^"]+)"', attr_string):
        attrs[key] = value
    return attrs


def gtf_to_bed12(
    gtf_file,
    bed_file,
    transcript_attr="transcript_id",
    gene_attr="gene_id",
    name_format="transcript_id",
):
    """
    从 GTF 生成 BED12。

    参数
    ----
    gtf_file : str
        输入 GTF 文件路径

    bed_file : str
        输出 BED 文件路径

    transcript_attr : str
        用哪个 attribute 作为 transcript ID，默认 transcript_id

    gene_attr : str
        用哪个 attribute 作为 gene ID，默认 gene_id

    name_format : str
        BED 第 4 列命名方式：
        - "transcript_id"：只用 transcript_id
        - "gene_id|transcript_id"：gene_id|transcript_id
    """

    transcripts = defaultdict(list)
    tx_info = {}

    with open(gtf_file) as fin:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs_raw = fields

            if feature != "exon":
                continue

            attrs = parse_gtf_attributes(attrs_raw)

            if transcript_attr not in attrs:
                continue

            tx_id = attrs[transcript_attr]
            gene_id = attrs.get(gene_attr, "NA")

            start = int(start)
            end = int(end)

            # GTF: 1-based closed
            # BED: 0-based half-open
            bed_start = start - 1
            bed_end = end

            transcripts[tx_id].append((bed_start, bed_end))

            if tx_id not in tx_info:
                tx_info[tx_id] = {
                    "chrom": chrom,
                    "strand": strand,
                    "gene_id": gene_id,
                }

    with open(bed_file, "w") as fout:
        for tx_id, exons in transcripts.items():
            if not exons:
                continue

            info = tx_info[tx_id]
            chrom = info["chrom"]
            strand = info["strand"]
            gene_id = info["gene_id"]

            # BED12 blocks 必须按基因组坐标升序
            exons = sorted(exons, key=lambda x: x[0])

            chrom_start = min(x[0] for x in exons)
            chrom_end = max(x[1] for x in exons)

            block_count = len(exons)
            block_sizes = [end - start for start, end in exons]
            block_starts = [start - chrom_start for start, end in exons]

            if name_format == "gene_id|transcript_id":
                name = f"{gene_id}|{tx_id}"
            else:
                name = tx_id

            score = "0"
            thick_start = chrom_start
            thick_end = chrom_end
            item_rgb = "0"

            bed_fields = [
                chrom,
                str(chrom_start),
                str(chrom_end),
                name,
                score,
                strand,
                str(thick_start),
                str(thick_end),
                item_rgb,
                str(block_count),
                ",".join(map(str, block_sizes)) + ",",
                ",".join(map(str, block_starts)) + ",",
            ]

            fout.write("\t".join(bed_fields) + "\n")

def make_bed():
    gtf_to_bed12(f'{n_tool_filter_gtf}',f'{n_tool_filter_gtf[:-4]}.bed')


if __name__ == '__main__':
    make_fa()
    make_bed()

    pass
    # plot_stat()


