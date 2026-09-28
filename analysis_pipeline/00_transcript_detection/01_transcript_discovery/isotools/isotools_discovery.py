# -*- coding: utf-8 -*-

import argparse
import os

import isotools
import pandas as pd


IG_REGIONS = {
    'chr14': [(103000000, 107500000)],
    'chr2': [(88000000, 90000000)],
    'chr22': [(22000000, 23500000)],
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref-gtf', required=True)
    parser.add_argument('--bam', required=True)
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--out-gtf', required=True)
    parser.add_argument('--out-table', required=True)
    parser.add_argument('--out-read-cov', required=True)
    parser.add_argument('--min-coverage', type=int, default=3)
    return parser.parse_args()


def get_ig_transcript_ids(tx_gtf):
    transcript_ids = set()
    ig_chroms = set(IG_REGIONS)
    usecols = [0, 2, 3, 4, 8]
    names = ['chrom', 'feature', 'start', 'end', 'attrs']
    for chunk in pd.read_csv(
            tx_gtf,
            sep='\t',
            comment='#',
            header=None,
            usecols=usecols,
            names=names,
            compression='infer',
            chunksize=200000,
    ):
        chunk = chunk[
            (chunk['feature'] == 'transcript') &
            (chunk['chrom'].isin(ig_chroms))
        ]
        if chunk.empty:
            continue
        region_mask = pd.Series(False, index=chunk.index)
        for chrom, regions in IG_REGIONS.items():
            chrom_mask = chunk['chrom'] == chrom
            for start, end in regions:
                region_mask |= chrom_mask & (chunk['start'] <= end) & (chunk['end'] >= start)
        chunk = chunk[region_mask]
        if chunk.empty:
            continue
        ids = chunk['attrs'].str.extract(r'transcript_id "([^"]+)"', expand=False)
        transcript_ids.update(ids.dropna())
    return transcript_ids


def write_read_coverage(tx_gtf, transcript_info_table, output_table):
    df = pd.read_csv(transcript_info_table, sep='\t')
    count_col = ''
    for col in df.columns:
        if col.endswith('_coverage'):
            count_col = col
            break
    if not count_col:
        raise Exception(f'No *_coverage column found in {transcript_info_table}')

    df['transcript_id'] = df['gene_id'].astype(str) + '_' + df['transcript_nr'].astype(str)
    df['read_count'] = df[count_col]
    ig_transcript_ids = get_ig_transcript_ids(tx_gtf)
    df['in_IG_region'] = df['transcript_id'].isin(ig_transcript_ids).astype('int8')
    os.makedirs(os.path.dirname(output_table), exist_ok=True)
    df[['transcript_id', 'read_count', 'in_IG_region']].to_csv(
        output_table,
        sep='\t',
        index=False,
        compression='infer',
    )


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.out_gtf), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_table), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_read_cov), exist_ok=True)

    tr = isotools.Transcriptome.from_reference(args.ref_gtf)
    tr.add_sample_from_bam(args.bam, sample_name=args.sample_id)
    tr.make_index()
    tr.write_gtf(args.out_gtf)

    transcript_tab = tr.transcript_table(
        coverage=True,
        tpm=True,
        min_coverage=args.min_coverage,
    )
    transcript_tab.to_csv(args.out_table, sep='\t', index=False)
    write_read_coverage(args.out_gtf, args.out_table, args.out_read_cov)
    print(f'Finished {args.sample_id}')


if __name__ == '__main__':
    main()
