# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/1/5 09:13
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""
import argparse
import csv
import logging
import os
import re
import subprocess
from collections import defaultdict

import pandas as pd
import polars as pl

from util import PROJ_DIR

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

import os
import polars as pl
import logging
from concurrent.futures import ThreadPoolExecutor

logger = logging.getLogger(__name__)


def combine_tables(
    tab_paths,
    gene_col,
    value_col,
    new_col_names,
    sep="\t",
    out_path=None,
    n_threads=8
):
    """
    Merge one value column from each table by gene/transcript ID

    Parameters
    ----------
    tab_paths : list[str]
        Input tables
    gene_col : str
        Join key column
    value_col : str
        Column to merge (single column per file)
    new_col_names : list[str]
        New column name for each table (must match tab_paths)
    sep : str
        Separator
    out_path : str or None
        Output path
    n_threads : int
        Thread number
    """
    # check
    if len(tab_paths) != len(new_col_names):
        raise ValueError("new_col_names length must match tab_paths")
    fine_table_paths = []
    fine_new_col_names=[]
    for i in range(len(tab_paths)):
        if os.path.exists(tab_paths[i]):
            fine_table_paths.append(tab_paths[i])
            fine_new_col_names.append(new_col_names[i])
    logger.info(f'remain {len(fine_table_paths)} tables; raw {len(tab_paths)}')
    tab_paths=fine_table_paths
    new_col_names=fine_new_col_names
    # ---------- parallel read ----------
    def _read_one(path, new_col):
        if not os.path.exists(path):
            logger.info(f"skip {path}")
            return None

        df = pl.read_csv(
            path,
            separator=sep,
            has_header=True,
            columns=[gene_col, value_col]
        )

        return df.rename({value_col: new_col})

    with ThreadPoolExecutor(max_workers=n_threads) as pool:
        dfs = list(pool.map(_read_one, tab_paths, new_col_names))

    dfs = [df for df in dfs if df is not None]
    if not dfs:
        return None

    # ---------- parallel tree-reduce join ----------
    def _merge_two(df1, df2):
        return df1.join(df2, on=gene_col, how="full", coalesce=True)

    while len(dfs) > 1:
        merged = []
        with ThreadPoolExecutor(max_workers=n_threads) as pool:
            it = iter(dfs)
            futures = []
            for df1 in it:
                df2 = next(it, None)
                if df2 is None:
                    merged.append(df1)
                else:
                    futures.append(pool.submit(_merge_two, df1, df2))

            for f in futures:
                merged.append(f.result())

        dfs = merged

    merged_df = (
        dfs[0]
        .fill_null(0)
        .sort(gene_col)
    )

    if out_path:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        merged_df.write_csv(out_path, separator='\t')

    return merged_df



def combine_samples(sample_dir, ref_tag, gene_type, quant_type,out_dir,nt):
    logger.info(f"combining {ref_tag}, {gene_type} samples into {out_dir}")
    paths=[]
    gene_type_alias={'gene':'genes','transcript':'isoforms'}
    gene_col_name={'gene':'gene_id','transcript':'transcript_id'}
    quant_type_alias={'count':'expected_count','tpm':'TPM'}
    sample_ids=[]
    for sample_id in os.listdir(sample_dir):
        paths.append(f'{sample_dir}/{sample_id}/RSEM_{ref_tag}.{gene_type_alias[gene_type]}.results')
        sample_ids.append(sample_id)
    combine_tables(paths,gene_col_name[gene_type],quant_type_alias[quant_type],sample_ids,
                   sep='\t',out_path=f'{out_dir}/{ref_tag}.{gene_type}.{quant_type}.rsem.tsv',
                   n_threads=nt)
    logger.info('done')


def main():
    parser = argparse.ArgumentParser(
        description='RSEM helper.',
    )
    subparsers = parser.add_subparsers(
        dest='command',
        required=True,
        help='Available commands'
    )

    # combine samples
    parser_com = subparsers.add_parser(
        'combine',
        help='RSEM combine'
    )
    parser_com.add_argument(
        '-s', '--sample-dir',
        required=True,
        type=str,
        help='Sample based dir'
    )
    parser_com.add_argument(
        '-t', '--tag',
        required=True,
        type=str,
        help='Ref tag'
    )
    parser_com.add_argument(
        '-g', '--gene-type',
        required=True,
        type=str,
        help='Ref tag'
    )
    parser_com.add_argument(
        '-q', '--quant-type',
        required=True,
        type=str,
        help='Ref tag'
    )
    parser_com.add_argument(
        '-o', '--out-dir',
        required=True,
        type=str,
        help='Output dir'
    )
    parser_com.add_argument(
        '-nt', '--nt',
        required=True,
        type=int,
        help='Thread'
    )
    args = parser.parse_args()
    if args.command == 'combine':
        combine_samples(args.sample_dir, args.tag, args.gene_type,args.quant_type, args.out_dir, args.nt)
    else:
        logger.error(f"Unsupported command: {args.command}")


if __name__ == '__main__':
    main()
