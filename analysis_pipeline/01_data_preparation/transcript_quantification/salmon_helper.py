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
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import polars as pl

from util import PROJ_DIR

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
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
    quant_type_alias={'count':'NumReads','tpm':'TPM'}
    sample_ids=[]
    for sample_id in os.listdir(sample_dir):
        paths.append(f'{sample_dir}/{sample_id}/salmon_{ref_tag}/quant.sf')
        sample_ids.append(sample_id)
    combine_tables(paths,'Name',quant_type_alias[quant_type],sample_ids,
                   sep='\t',out_path=f'{out_dir}/{ref_tag}.{gene_type}.{quant_type}.salmon.tsv',
                   n_threads=nt)
    logger.info('done')


def __isoform_to_gene_sum(df: pd.DataFrame, iso2gene: pd.Series) -> pd.DataFrame:
    df = df.copy()
    df["gene_id"] = iso2gene.loc[df.index].values
    sdf=df.groupby("gene_id").sum()
    return sdf

def merge_gene(out_dir,ref_tag,quant_type,gtf_prefix):
    gene_annot_df=pd.read_csv(f'{gtf_prefix}.gene_transcript_map.txt',sep='\t')
    isoform_gene = gene_annot_df.set_index("transcript_id")["gene_id"]
    logger.info(f'start {ref_tag}, {quant_type}')
    df=pd.read_csv(f'{out_dir}/{ref_tag}.transcript.{quant_type}.salmon.tsv',sep='\t',index_col=0)
    # gene level
    gene_tpm_df=__isoform_to_gene_sum(df, isoform_gene)
    gene_tpm_df.to_csv(f'{out_dir}/{ref_tag}.gene.{quant_type}.salmon.tsv',sep='\t')
    logger.info(f'done {ref_tag}, {quant_type}')


def main():
    parser = argparse.ArgumentParser(
        description='Stringtie2 helper.',
    )
    subparsers = parser.add_subparsers(
        dest='command',
        required=True,
        help='Available commands'
    )

    # combine samples
    parser_com = subparsers.add_parser(
        'combine',
        help='Salmon combine'
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
        help='Output dir'
    )

    # merge gene
    parser_gene = subparsers.add_parser(
        'merge_gene',
        help='Merge'
    )
    parser_gene.add_argument(
        '-o', '--out-dir',
        required=True,
        type=str,
        help=''
    )
    parser_gene.add_argument(
        '-t', '--ref-tag',
        required=True,
        type=str,
        help='Ref tag'
    )
    parser_gene.add_argument(
        '-q', '--quant-type',
        required=True,
        type=str,
        help='Ref tag'
    )

    parser_gene.add_argument(
        '-g', '--gtf-prefix',
        required=True,
        type=str,
        help='Ref tag'
    )

    args = parser.parse_args()
    if args.command == 'combine':
        combine_samples(args.sample_dir, args.tag, 'transcript',args.quant_type, args.out_dir, args.nt)
    elif args.command == 'merge_gene':
        merge_gene(args.out_dir, args.ref_tag, args.quant_type, args.gtf_prefix)
    else:
        logger.error(f"Unsupported command: {args.command}")

if __name__ == '__main__':
    main()
