# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/1/5 09:13
@Email   : xuechao@szbl.ac.cn
@Desc    :  
"""
import argparse
import logging
import os

import polars as pl

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

from concurrent.futures import ThreadPoolExecutor


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

    n_threads = max(1, n_threads)
    problem_tables = []
    logger.info(f'checking {len(tab_paths)} tables')

    def _read_header(path):
        with open(path, "r") as handle:
            line = handle.readline().strip()
        if not line:
            return []
        return line.split(sep)

    def _write_problem_report():
        if not out_path or not problem_tables:
            return
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        problem_path = f"{out_path}.problems.tsv"
        with open(problem_path, "w") as handle:
            handle.write("sample_id\tpath\treason\n")
            for sample_id, path, reason in problem_tables:
                handle.write(f"{sample_id}\t{path}\t{reason}\n")
        logger.error(f"wrote problem table report: {problem_path}")

    def _raise_if_problems():
        if not problem_tables:
            return
        logger.error(f"found {len(problem_tables)} problem tables; stop combining")
        for sample_id, path, reason in problem_tables:
            logger.error(f"problem sample={sample_id}, path={path}, reason={reason}")
        _write_problem_report()
        raise RuntimeError(
            f"Found {len(problem_tables)} problem tables. "
            "Fix them before combining; see log for all sample IDs and reasons."
        )

    for path, sample_id in zip(tab_paths, new_col_names):
        if not os.path.exists(path):
            problem_tables.append((sample_id, path, "missing file"))
            continue
        try:
            header = _read_header(path)
        except Exception as exc:
            problem_tables.append(
                (sample_id, path, f"read header failed: {type(exc).__name__}: {exc}")
            )
            continue
        missing_cols = [col for col in (gene_col, value_col) if col not in header]
        if missing_cols:
            problem_tables.append((
                sample_id,
                path,
                f"missing columns: {','.join(missing_cols)}; valid columns: {','.join(header)}"
            ))

    _raise_if_problems()

    # ---------- parallel read ----------
    def _read_one(path, new_col):
        try:
            df = pl.read_csv(
                path,
                separator=sep,
                has_header=True,
                columns=[gene_col, value_col]
            )
        except Exception as exc:
            return None, (new_col, path, f"read table failed: {type(exc).__name__}: {exc}")

        valid_gene_expr = (
            pl.col(gene_col).is_not_null()
            & (pl.col(gene_col).cast(pl.Utf8).str.strip_chars() != "")
        )
        invalid_row_count = df.filter(~valid_gene_expr).height
        if invalid_row_count:
            logger.warning(
                f"sample={new_col}, path={path}, "
                f"dropped {invalid_row_count} rows with empty {gene_col}"
            )
            df = df.filter(valid_gene_expr)

        return df.rename({value_col: new_col}), None

    with ThreadPoolExecutor(max_workers=n_threads) as pool:
        read_results = list(pool.map(_read_one, tab_paths, new_col_names))

    dfs = []
    for df, problem in read_results:
        if df is not None:
            dfs.append(df)
        if problem is not None:
            problem_tables.append(problem)

    _raise_if_problems()

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



def combine_samples(sample_dir, ref_tag, gene_type, quant_type, out_dir, nt):
    logger.info(f"combining {ref_tag}, {gene_type} samples into {out_dir}")
    paths = []
    gene_type_alias = {'gene': 'genes', 'transcript': 'isoforms'}
    gene_col_name = {'gene': 'gene_id', 'transcript': 'transcript_id'}
    quant_type_alias = {'count': 'expected_count', 'tpm': 'TPM'}
    sample_ids = []
    for sample_id in sorted(os.listdir(sample_dir)):
        sample_path = os.path.join(sample_dir, sample_id)
        if not os.path.isdir(sample_path):
            continue
        paths.append(f'{sample_dir}/{sample_id}/RSEM_{ref_tag}.{gene_type_alias[gene_type]}.results')
        sample_ids.append(sample_id)
    combine_tables(paths, gene_col_name[gene_type], quant_type_alias[quant_type], sample_ids,
                   sep='\t', out_path=f'{out_dir}/{ref_tag}.{gene_type}.{quant_type}.rsem.tsv',
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
        choices=['gene', 'transcript'],
        type=str,
        help='Ref tag'
    )
    parser_com.add_argument(
        '-q', '--quant-type',
        required=True,
        choices=['count', 'tpm'],
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
