# -*- coding: utf-8 -*-
"""
@Desc    :  
Extract novel isoforms from SQANTI3 filtered GTF based on classification report.

Enhanced version:
 - Extracts novel isoforms (Isoform + not known categories)
 - Replaces gene_id and transcript_id in GTF using associated_gene / associated_transcript mappings
"""

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_REF_DIR = os.environ.get('GTOP_REF_DIR', str(Path(CONFIG_PROJECT_DIR) / 'reference'))
CONFIG_HG38_FASTA = os.environ.get('HG38_FASTA', str(Path(CONFIG_REF_DIR) / 'hg38.fa'))
CONFIG_HG38_GTF = os.environ.get('HG38_GTF', str(Path(CONFIG_REF_DIR) / 'gencode.v47.annotation.gtf'))
CONFIG_BUILD_LOCI = os.environ.get('BUILD_LOCI', str(Path(CONFIG_PROJECT_DIR) / 'software/buildLoci/buildLoci.pl'))

import csv
from contextlib import contextmanager
import logging
import re
import shlex
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# all samples
CLEAN_SQANTI3_TAG='sqanti3_clean'
FINAL_SQANTI3_TAG='sqanti3_filter_final'
ENHANCED_GTF_TAG='enhanced_gtf'

# samples without averaged aligned read length < 1000bp.
# CLEAN_SQANTI3_TAG='sqanti3_clean_gt_1k'
# FINAL_SQANTI3_TAG='sqanti3_filter_final_gt_1k'
# ENHANCED_GTF_TAG='enhanced_gtf_gt_1k'


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_DIR = os.environ.get("PROJECT_ROOT")

REF_DIR = CONFIG_REF_DIR
REF_GENOME_FASTA = CONFIG_HG38_FASTA
REF_GTF = CONFIG_HG38_GTF
GFFREAD_BIN = os.environ.get('GFFREAD_BIN', 'gffread')
LOAD_HTSLIB_ENV_CMD = os.environ.get('LOAD_HTSLIB_ENV_CMD', 'true')
GTF_NESTED_SORT_PROGRESS_INTERVAL = int(os.environ.get("GTOP_GTF_SORT_PROGRESS_INTERVAL", "500000"))
GTF_TRANSCRIPT_FEATURES = {'transcript', 'mRNA'}
GTF_GENE_ID_RE = re.compile(r'(?:^|;\s*)gene_id "([^"]*)"')
GTF_TRANSCRIPT_ID_RE = re.compile(r'(?:^|;\s*)transcript_id "([^"]*)"')


def configure_output_tags(clean_sqanti3_tag=None, final_sqanti3_tag=None, enhanced_gtf_tag=None):
    global CLEAN_SQANTI3_TAG, FINAL_SQANTI3_TAG, ENHANCED_GTF_TAG

    if clean_sqanti3_tag:
        CLEAN_SQANTI3_TAG = clean_sqanti3_tag
    if final_sqanti3_tag:
        FINAL_SQANTI3_TAG = final_sqanti3_tag
    if enhanced_gtf_tag:
        ENHANCED_GTF_TAG = enhanced_gtf_tag
    logger.info(
        "Using output tags: "
        f"CLEAN_SQANTI3_TAG={CLEAN_SQANTI3_TAG}, "
        f"FINAL_SQANTI3_TAG={FINAL_SQANTI3_TAG}, "
        f"ENHANCED_GTF_TAG={ENHANCED_GTF_TAG}"
    )


class NestedGTFRecord:
    __slots__ = ('line', 'sort_key', 'child_sort_key')

    def __init__(self, line, sort_key, child_sort_key):
        self.line = line
        self.sort_key = sort_key
        self.child_sort_key = child_sort_key

    def to_line(self):
        return self.line + "\n"


def _extract_gtf_gene_transcript_id(attr_text):
    gene_match = GTF_GENE_ID_RE.search(attr_text)
    transcript_match = GTF_TRANSCRIPT_ID_RE.search(attr_text)
    gene_id = gene_match.group(1) if gene_match else None
    transcript_id = transcript_match.group(1) if transcript_match else None
    return gene_id, transcript_id


def _gtf_nested_chrom_sort_key(chrom):
    c = chrom.strip()
    c_no_chr = c[3:] if c.lower().startswith("chr") else c
    c_lower = c_no_chr.lower()
    if c_lower.isdigit():
        return 0, int(c_lower)
    if c_lower == "x":
        return 0, 23
    if c_lower == "y":
        return 0, 24
    if c_lower in {"m", "mt"}:
        return 0, 25
    return 1, c


def _gtf_nested_feature_rank(feature):
    ranks = {
        'exon': 0,
        'CDS': 1,
        'start_codon': 2,
        'stop_codon': 3,
        'UTR': 4,
        'five_prime_utr': 5,
        'three_prime_utr': 6,
        'Selenocysteine': 7,
    }
    return ranks.get(feature, 50)


def _gtf_nested_record_sort_key(record):
    return record.sort_key


def _gtf_nested_child_sort_key(record):
    return record.child_sort_key


def _update_gtf_nested_anchor(anchor_map, key, sort_key):
    current = anchor_map.get(key)
    if current is None or sort_key < current:
        anchor_map[key] = sort_key


def _log_gtf_nested_sort_progress(stage, count, start_time, extra=''):
    elapsed = time.perf_counter() - start_time
    rate = count / elapsed if elapsed > 0 else 0.0
    suffix = f"; {extra}" if extra else ""
    logger.info(f"{stage} {count:,} records in {elapsed:.1f}s ({rate:,.0f} records/s){suffix}")


def _load_and_group_gtf_for_nested_sort(gtf_path, progress_interval):
    comments = []
    gene_records = defaultdict(list)
    transcript_records = defaultdict(list)
    transcript_children = defaultdict(list)
    gene_to_transcripts = defaultdict(set)
    orphan_records = []
    gene_anchor_keys = {}
    transcript_anchor_keys = {}
    record_count = 0
    start_time = time.perf_counter()

    logger.info(f"Reading and grouping GTF for nested sort: {gtf_path}")
    with open(gtf_path) as fin:
        for line_no, line in enumerate(fin, start=1):
            if line.startswith("#"):
                comments.append(line)
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Line {line_no} does not have 9 GTF columns: {gtf_path}")
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError as exc:
                raise ValueError(f"Line {line_no} has invalid start/end coordinates: {gtf_path}") from exc

            feature = fields[2]
            gene_id, transcript_id = _extract_gtf_gene_transcript_id(fields[8])
            chrom_key = _gtf_nested_chrom_sort_key(fields[0])
            sort_key = (chrom_key, start, end, line_no)
            record = NestedGTFRecord(
                line,
                sort_key,
                (chrom_key, start, end, _gtf_nested_feature_rank(feature), line_no),
            )

            if feature == 'gene' and gene_id:
                gene_records[gene_id].append(record)
                _update_gtf_nested_anchor(gene_anchor_keys, gene_id, sort_key)
            elif feature in GTF_TRANSCRIPT_FEATURES and transcript_id:
                transcript_records[transcript_id].append(record)
                _update_gtf_nested_anchor(transcript_anchor_keys, transcript_id, sort_key)
                if gene_id:
                    gene_to_transcripts[gene_id].add(transcript_id)
                    _update_gtf_nested_anchor(gene_anchor_keys, gene_id, sort_key)
            elif transcript_id:
                transcript_children[transcript_id].append(record)
                _update_gtf_nested_anchor(transcript_anchor_keys, transcript_id, sort_key)
                if gene_id:
                    gene_to_transcripts[gene_id].add(transcript_id)
                    _update_gtf_nested_anchor(gene_anchor_keys, gene_id, sort_key)
            elif gene_id:
                gene_records[gene_id].append(record)
                _update_gtf_nested_anchor(gene_anchor_keys, gene_id, sort_key)
            else:
                orphan_records.append(record)

            record_count += 1
            if progress_interval > 0 and record_count % progress_interval == 0:
                _log_gtf_nested_sort_progress(
                    'Read/grouped',
                    record_count,
                    start_time,
                    f"genes={len(gene_anchor_keys):,}, transcripts={len(transcript_anchor_keys):,}",
                )

    _log_gtf_nested_sort_progress(
        'Finished reading/grouping',
        record_count,
        start_time,
        f"comments={len(comments):,}, genes={len(gene_anchor_keys):,}, transcripts={len(transcript_anchor_keys):,}",
    )
    return (
        comments,
        gene_records,
        transcript_records,
        transcript_children,
        gene_to_transcripts,
        orphan_records,
        gene_anchor_keys,
        transcript_anchor_keys,
        record_count,
    )


def _write_nested_sorted_gtf(input_gtf, output_gtf, progress_interval=GTF_NESTED_SORT_PROGRESS_INTERVAL):
    total_start = time.perf_counter()
    (
        comments,
        gene_records,
        transcript_records,
        transcript_children,
        gene_to_transcripts,
        orphan_records,
        gene_anchor_keys,
        transcript_anchor_keys,
        record_count,
    ) = _load_and_group_gtf_for_nested_sort(input_gtf, progress_interval)

    sort_start = time.perf_counter()
    logger.info("Sorting nested GTF groups")
    sorted_gene_ids = sorted(gene_anchor_keys, key=gene_anchor_keys.__getitem__)
    all_transcript_ids = set(transcript_records) | set(transcript_children)
    logger.info(
        f"Finished sorting top-level groups in {time.perf_counter() - sort_start:.1f}s; "
        f"genes={len(sorted_gene_ids):,}, transcripts={len(all_transcript_ids):,}"
    )

    emitted_transcripts = set()
    emitted_lines = 0
    os.makedirs(os.path.dirname(os.path.abspath(output_gtf)), exist_ok=True)
    write_start = time.perf_counter()
    logger.info(f"Writing nested sorted GTF: {output_gtf}")
    with open(output_gtf, "w") as fout:
        for comment in comments:
            fout.write(comment)

        next_write_log = progress_interval if progress_interval > 0 else None
        for gene_id in sorted_gene_ids:
            for gene_record in sorted(gene_records.get(gene_id, []), key=_gtf_nested_record_sort_key):
                fout.write(gene_record.to_line())
                emitted_lines += 1

            transcript_ids = sorted(
                gene_to_transcripts.get(gene_id, set()),
                key=transcript_anchor_keys.__getitem__,
            )
            for transcript_id in transcript_ids:
                if transcript_id in emitted_transcripts:
                    continue
                for tx_record in sorted(transcript_records.get(transcript_id, []), key=_gtf_nested_record_sort_key):
                    fout.write(tx_record.to_line())
                    emitted_lines += 1
                for child_record in sorted(
                    transcript_children.get(transcript_id, []),
                    key=_gtf_nested_child_sort_key,
                ):
                    fout.write(child_record.to_line())
                    emitted_lines += 1
                emitted_transcripts.add(transcript_id)
                if next_write_log is not None and emitted_lines >= next_write_log:
                    _log_gtf_nested_sort_progress('Written', emitted_lines, write_start)
                    while emitted_lines >= next_write_log:
                        next_write_log += progress_interval

        remaining_transcripts = sorted(
            all_transcript_ids - emitted_transcripts,
            key=transcript_anchor_keys.__getitem__,
        )
        for transcript_id in remaining_transcripts:
            for tx_record in sorted(transcript_records.get(transcript_id, []), key=_gtf_nested_record_sort_key):
                fout.write(tx_record.to_line())
                emitted_lines += 1
            for child_record in sorted(
                transcript_children.get(transcript_id, []),
                key=_gtf_nested_child_sort_key,
            ):
                fout.write(child_record.to_line())
                emitted_lines += 1
            if next_write_log is not None and emitted_lines >= next_write_log:
                _log_gtf_nested_sort_progress('Written', emitted_lines, write_start)
                while emitted_lines >= next_write_log:
                    next_write_log += progress_interval

        for orphan_record in sorted(orphan_records, key=_gtf_nested_record_sort_key):
            fout.write(orphan_record.to_line())
            emitted_lines += 1
            if next_write_log is not None and emitted_lines >= next_write_log:
                _log_gtf_nested_sort_progress('Written', emitted_lines, write_start)
                while emitted_lines >= next_write_log:
                    next_write_log += progress_interval

    _log_gtf_nested_sort_progress('Finished writing', emitted_lines, write_start)
    logger.info(
        f"Nested sort complete for {input_gtf}: input={record_count:,}; output={emitted_lines:,}; "
        f"genes={len(sorted_gene_ids):,}; transcripts={len(all_transcript_ids):,}; "
        f"elapsed={time.perf_counter() - total_start:.1f}s"
    )
    return record_count, emitted_lines


def sort_gtf_nested_inplace(gtf_path, progress_interval=GTF_NESTED_SORT_PROGRESS_INTERVAL):
    tmp_path = f"{gtf_path}.nested_sort.{os.getpid()}.tmp"
    try:
        logger.info(f"Nested-sort and replace GTF: {gtf_path}")
        input_records, output_records = _write_nested_sorted_gtf(gtf_path, tmp_path, progress_interval)
        if input_records != output_records:
            raise RuntimeError(
                f"Nested sort record count mismatch for {gtf_path}: "
                f"input={input_records:,}; output={output_records:,}"
            )
        os.replace(tmp_path, gtf_path)
        logger.info(f"Replaced original GTF with nested sorted file: {gtf_path}")
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
            logger.info(f"Deleted temporary nested-sort file: {tmp_path}")


@contextmanager
def timed_step(label):
    start = time.perf_counter()
    logger.info(f"START {label}")
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        logger.info(f"DONE {label} in {elapsed:.2f}s")


class GTFParser:
    def __init__(self, gtf_file=None):
        self.gtf_file = gtf_file
        self.gtf_df:pd.DataFrame
        self.comment = None
        self.cols = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
        self._attribute_patterns = {}
        self._cache = {}

    @staticmethod
    def parse_attributes(attr_str):
        attr_dict = {}
        for attr in attr_str.strip().split(";"):
            attr = attr.strip()
            if not attr:
                continue
            if " " in attr:
                key, val = attr.split(" ", 1)
                attr_dict[key] = val.strip('"')
        return attr_dict

    @staticmethod
    def build_attributes(attr_dict):
        return "; ".join(f'{k} "{v}"' for k, v in attr_dict.items()) + ";"

    @classmethod
    def extract_attribute(cls, attributes, key):
        pattern = getattr(cls, '_shared_attribute_patterns', None)
        if pattern is None:
            cls._shared_attribute_patterns = {}
            pattern = cls._shared_attribute_patterns
        if key not in pattern:
            pattern[key] = re.compile(r'(?:^|;\s*)' + re.escape(key) + r' "([^"]*)"')
        return attributes.astype(str).str.extract(pattern[key], expand=False)

    @classmethod
    def get_attr_value(cls, attr_str, key):
        pattern = getattr(cls, '_shared_attribute_patterns', None)
        if pattern is None:
            cls._shared_attribute_patterns = {}
            pattern = cls._shared_attribute_patterns
        if key not in pattern:
            pattern[key] = re.compile(r'(?:^|;\s*)' + re.escape(key) + r' "([^"]*)"')
        match = pattern[key].search(attr_str)
        return match.group(1) if match else None

    @classmethod
    def add_id_columns(cls, df, include_gene_name=False):
        if 'transcript_id' not in df.columns:
            df.loc[:, 'transcript_id'] = cls.extract_attribute(df['attributes'], 'transcript_id')
        if 'gene_id' not in df.columns:
            df.loc[:, 'gene_id'] = cls.extract_attribute(df['attributes'], 'gene_id')
        if include_gene_name and 'gene_name' not in df.columns:
            df.loc[:, 'gene_name'] = cls.extract_attribute(df['attributes'], 'gene_name')
        return df

    def load_comment(self):
        comment_arr=[]
        with open(self.gtf_file) as f:
            for line in f:
                if line.startswith("#"):
                    comment_arr.append(line.strip())
                else:
                    break
        return comment_arr

    def create_from_arr(self,arr):
        self._cache.clear()
        if isinstance(arr, (list)):
            cols = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
            df = pd.DataFrame(arr, columns=cols)
        else:
            df=arr.copy()
        df = self.add_id_columns(df)
        self.gtf_df = df
        self.gtf_summary()
        return self

    def write_file(self, out_file, comment_str):
        df = self.gtf_df
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        logger.info(f'Writing {len(df):,} GTF records to {out_file}')
        with open(out_file, "w") as fout:
            fout.write(comment_str or '')
            df.to_csv(
                fout,
                sep='\t',
                columns=self.cols,
                header=False,
                index=False,
                quoting=csv.QUOTE_NONE,
                escapechar='\\',
            )
        logger.info(f"Wrote corrected GTF to {out_file} ({len(df):,} records)")

    def load_df(self):
        self._cache.clear()
        cols = self.cols
        df = pd.read_csv(self.gtf_file, sep='\t', comment='#', names=cols, dtype=str, memory_map=True)
        df['start'] = df['start'].astype('int64')
        df['end'] = df['end'].astype('int64')
        df = self.add_id_columns(df, include_gene_name=True)
        self.gtf_df = df
        self.gtf_summary()
        return self

    def get_unique_gene_id(self, only_gene_future=False):
        cache_key = ('unique_gene_id', only_gene_future)
        if cache_key in self._cache:
            return self._cache[cache_key]
        anno = self.add_id_columns(self.gtf_df.copy())
        if only_gene_future:
            anno_gene_mask = anno['feature'] == 'gene'
        else:
            anno_gene_mask = pd.Series(True, index=anno.index)
        result = set(anno.loc[anno_gene_mask, 'gene_id'].dropna())
        self._cache[cache_key] = result
        return result

    def get_transcript_gene_id_map(self,gene_item='gene_id'):
        cache_key = ('transcript_gene_id_map', gene_item)
        if cache_key in self._cache:
            return self._cache[cache_key]
        anno = self.gtf_df
        anno_gene_mask = anno['feature'] == 'transcript'
        sub = anno.loc[anno_gene_mask, ['attributes', 'transcript_id']].copy()
        if gene_item in anno.columns:
            sub.loc[:, gene_item] = anno.loc[anno_gene_mask, gene_item].values
        else:
            sub.loc[:, gene_item] = self.extract_attribute(sub['attributes'], gene_item)
        sub = sub.dropna(subset=['transcript_id', gene_item])
        result = dict(zip(sub['transcript_id'].values, sub[gene_item].values))
        self._cache[cache_key] = result
        return result

    def get_transcript_id_attr_map(self):
        cache_key = ('transcript_id_attr_map',)
        if cache_key in self._cache:
            return self._cache[cache_key]
        anno = self.gtf_df
        anno_gene_mask = anno['feature'] == 'transcript'
        sub = anno.loc[anno_gene_mask, ['transcript_id', 'attributes']].dropna(subset=['transcript_id'])
        result = dict(zip(sub['transcript_id'].values, sub['attributes'].values))
        self._cache[cache_key] = result
        return result

    def get_transcript_length_map(self):
        cache_key = ('transcript_length_map',)
        if cache_key in self._cache:
            return self._cache[cache_key]
        anno = self.gtf_df
        anno_gene_mask = anno['feature'] == 'transcript'
        sub = anno.loc[anno_gene_mask, ['transcript_id', 'start', 'end']].dropna(subset=['transcript_id'])
        lengths = (sub['start'] - sub['end']).abs()
        result = dict(zip(sub['transcript_id'].values, lengths.values))
        self._cache[cache_key] = result
        return result

    def get_attribute_map(self,feature='gene',id='gene_id'):
        cache_key = ('attribute_map', feature, id)
        if cache_key in self._cache:
            return self._cache[cache_key]
        anno = self.gtf_df
        attr_map = {}
        anno_gene_mask = anno['feature'] == feature
        for attr in anno.loc[anno_gene_mask, 'attributes'].values:
            d = self.parse_attributes(attr)
            attr_map[d[id]] = d
        self._cache[cache_key] = attr_map
        return attr_map

    def get_position_map(self,feature='gene',id='gene_id'):
        cache_key = ('position_map', feature, id)
        if cache_key in self._cache:
            return self._cache[cache_key]
        anno = self.gtf_df
        anno_gene_mask = anno['feature'] == feature
        sub = anno.loc[anno_gene_mask, ['attributes', 'start', 'end']].copy()
        if id in anno.columns:
            sub.loc[:, id] = anno.loc[anno_gene_mask, id].values
        else:
            sub.loc[:, id] = self.extract_attribute(sub['attributes'], id)
        sub = sub.dropna(subset=[id])
        result = dict(zip(sub[id].values, zip(sub['start'].values, sub['end'].values)))
        self._cache[cache_key] = result
        return result

    def gtf_summary(self):
        df  = self.gtf_df
        name = ''
        if self.gtf_file is not None:
            name = os.path.basename(self.gtf_file)
        gene_ids = df['gene_id'] if 'gene_id' in df.columns else self.extract_attribute(df['attributes'], 'gene_id')
        transcript_ids = df['transcript_id'] if 'transcript_id' in df.columns else self.extract_attribute(df['attributes'], 'transcript_id')
        genes = gene_ids.nunique(dropna=True)
        transcripts = transcript_ids.nunique(dropna=True)
        exons = int((df['feature'] == 'exon').sum())
        logger.info(f"[{name} GTF summary] genes={genes:,}; transcripts={transcripts:,}; exons={exons:,}")

    def correct_transcript_boundaries(self):
        df = self.gtf_df.copy()
        exons = df[df["feature"] == "exon"]
        exon_bounds = exons.groupby("transcript_id", sort=False).agg(
            new_start=("start", "min"),
            new_end=("end", "max"),
        )
        transcript_mask = df["feature"] == "transcript"
        mapped_start = df.loc[transcript_mask, "transcript_id"].map(exon_bounds["new_start"])
        mapped_end = df.loc[transcript_mask, "transcript_id"].map(exon_bounds["new_end"])
        has_bounds = mapped_start.notna() & mapped_end.notna()
        tx_index = df.loc[transcript_mask].index
        df.loc[tx_index[has_bounds.values], "start"] = mapped_start.loc[has_bounds].astype('int64').values
        df.loc[tx_index[has_bounds.values], "end"] = mapped_end.loc[has_bounds].astype('int64').values
        ## correct start/end
        mask = df['start'] > df['end']
        if mask.any():
            logger.info(f"Warning: found {mask.sum()} lines with start > end; fixing...")
            df.loc[mask, ['start', 'end']] = df.loc[mask, ['end', 'start']].values
        logger.info(f" * Corrected transcript start/end using exon boundaries")
        self.gtf_df = df
        self._cache.clear()
        return

    def update_gene_coor(self):
        df_updated = self.add_id_columns(self.gtf_df.copy())
        gene_mask = df_updated['feature'] == 'gene'
        transcript_rows = df_updated[df_updated['feature'] == 'transcript']
        if not transcript_rows.empty:
            transcript_stats = transcript_rows.groupby('gene_id').agg(
                min_transcript_start=('start', 'min'),
                max_transcript_end=('end', 'max')
            )
            new_start = df_updated.loc[gene_mask, 'gene_id'].map(transcript_stats['min_transcript_start'])
            new_end = df_updated.loc[gene_mask, 'gene_id'].map(transcript_stats['max_transcript_end'])
            has_stats = new_start.notna() & new_end.notna()
            gene_index = df_updated.loc[gene_mask].index
            df_updated.loc[gene_index[has_stats.values], 'start'] = new_start.loc[has_stats].astype('int64').values
            df_updated.loc[gene_index[has_stats.values], 'end'] = new_end.loc[has_stats].astype('int64').values
        self.gtf_df = df_updated
        self._cache.clear()
        return self

    def filter_transcripts_with_genes(self, transcript_ids):
        gtf_df = self.add_id_columns(self.gtf_df.copy())
        transcript_ids = set(transcript_ids)
        target_genes = gtf_df.loc[gtf_df['transcript_id'].isin(transcript_ids), 'gene_id'].dropna().unique()
        mask = gtf_df['transcript_id'].isin(transcript_ids) | (
                    (gtf_df['feature'] == 'gene') & gtf_df['gene_id'].isin(target_genes))
        return gtf_df.loc[mask, self.cols]

    def filter_genes(self, gene_list):
        gtf_df = self.add_id_columns(self.gtf_df.copy())
        gene_list = set(gene_list)
        mask = (gtf_df['feature'] == 'gene') & gtf_df['gene_id'].isin(gene_list)
        return gtf_df.loc[mask, self.cols]


def _fields_to_record(fields):
    return {
        'chrom': fields[0],
        'source': fields[1],
        'feature': fields[2],
        'start': int(fields[3]),
        'end': int(fields[4]),
        'score': fields[5],
        'strand': fields[6],
        'frame': fields[7],
        'attributes': fields[8],
    }


def _records_to_gtf_df(records):
    cols = GTFParser().cols
    if not records:
        return pd.DataFrame(columns=cols)
    df = pd.DataFrame.from_records(records, columns=cols)
    df = GTFParser.add_id_columns(df)
    return df


def _append_gtf_dataframe(fout, df):
    if df is None or df.empty:
        return
    df.to_csv(
        fout,
        sep='\t',
        columns=GTFParser().cols,
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE,
        escapechar='\\',
    )


def _concat_gtf_frames(frames):
    cols = GTFParser().cols
    valid_frames = [df.loc[:, cols] for df in frames if df is not None and not df.empty]
    if not valid_frames:
        return pd.DataFrame(columns=cols)
    return pd.concat(valid_frames, ignore_index=True, sort=False)


class ReferenceGTFIndex:
    """Small in-memory index for REF_GTF; full reference records are streamed from disk."""

    def __init__(self, gtf_file):
        self.gtf_file = gtf_file
        self.cols = GTFParser().cols
        self.gene_attr = {}
        self.gene_records = {}
        self.gene_order = []
        self.tx_gene = {}
        self.transcript_bounds_by_gene = {}
        self._cache = {}

    def load_index(self):
        self.gene_attr.clear()
        self.gene_records.clear()
        self.gene_order.clear()
        self.tx_gene.clear()
        self.transcript_bounds_by_gene.clear()
        self._cache.clear()

        with open(self.gtf_file) as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                feature = fields[2]
                attr = fields[8]
                if feature == 'gene':
                    gid = GTFParser.get_attr_value(attr, 'gene_id')
                    if gid is None:
                        continue
                    self.gene_records[gid] = _fields_to_record(fields)
                    self.gene_attr[gid] = GTFParser.parse_attributes(attr)
                    self.gene_order.append(gid)
                elif feature == 'transcript':
                    tid = GTFParser.get_attr_value(attr, 'transcript_id')
                    gid = GTFParser.get_attr_value(attr, 'gene_id')
                    if tid is None or gid is None:
                        continue
                    self.tx_gene[tid] = gid
                    start = int(fields[3])
                    end = int(fields[4])
                    if gid in self.transcript_bounds_by_gene:
                        old_start, old_end = self.transcript_bounds_by_gene[gid]
                        self.transcript_bounds_by_gene[gid] = (min(old_start, start), max(old_end, end))
                    else:
                        self.transcript_bounds_by_gene[gid] = (start, end)

        logger.info(
            f"[{os.path.basename(self.gtf_file)} reference index] "
            f"genes={len(self.gene_attr):,}; transcripts={len(self.tx_gene):,}"
        )
        return self

    def get_transcript_gene_id_map(self, gene_item='gene_id'):
        if gene_item != 'gene_id':
            raise ValueError("ReferenceGTFIndex only supports transcript -> gene_id")
        return self.tx_gene

    def get_attribute_map(self, feature='gene', id='gene_id'):
        if feature != 'gene' or id != 'gene_id':
            raise ValueError("ReferenceGTFIndex only supports gene attribute map by gene_id")
        return self.gene_attr

    def filter_genes(self, gene_list):
        gene_set = set(gene_list)
        records = [self.gene_records[gid] for gid in self.gene_order if gid in gene_set and gid in self.gene_records]
        return _records_to_gtf_df(records)

    def filter_transcripts_with_genes(self, transcript_ids):
        transcript_ids = set(transcript_ids)
        target_genes = {self.tx_gene[tid] for tid in transcript_ids if tid in self.tx_gene}
        records = []
        with open(self.gtf_file) as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                attr = fields[8]
                keep = False
                if fields[2] == 'gene':
                    gid = GTFParser.get_attr_value(attr, 'gene_id')
                    keep = gid in target_genes
                else:
                    tid = GTFParser.get_attr_value(attr, 'transcript_id')
                    keep = tid in transcript_ids
                if keep:
                    records.append(_fields_to_record(fields))
        return _records_to_gtf_df(records)

    def write_enhanced_gtf(self, out_file, comment, gene_bounds, append_dfs):
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        logger.info(f"Streaming reference GTF into {out_file}")
        ref_records = 0
        updated_gene_records = 0
        with open(self.gtf_file) as fin, open(out_file, "w") as fout:
            fout.write(comment or '')
            for line in fin:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                if fields[2] == 'gene':
                    gid = GTFParser.get_attr_value(fields[8], 'gene_id')
                    if gid in gene_bounds:
                        fields[3], fields[4] = map(str, gene_bounds[gid])
                        updated_gene_records += 1
                fout.write("\t".join(fields) + "\n")
                ref_records += 1
            for df in append_dfs:
                _append_gtf_dataframe(fout, df)
        appended = sum(0 if df is None else len(df) for df in append_dfs)
        logger.info(
            f"Wrote enhanced GTF by streaming: ref_records={ref_records:,}; "
            f"updated_genes={updated_gene_records:,}; appended={appended:,}"
        )


def _predicted_nmd_type(value):
    if pd.isna(value):
        return 'noncoding_RNA'
    if isinstance(value, (bool, np.bool_)):
        return 'nonsense_mediated_decay' if value else 'protein_coding'
    value = str(value).strip().lower()
    if value in {'true', 't', '1', 'yes'}:
        return 'nonsense_mediated_decay'
    if value in {'false', 'f', '0', 'no'}:
        return 'protein_coding'
    return 'noncoding_RNA'


def _load_isoform_annotation(class_file, ref_gtf):
    annot_df = pd.read_csv(class_file, sep='\t')
    for col in ["isoform", "structural_category", "filter_result", "associated_transcript", "associated_gene"]:
        if col not in annot_df.columns:
            raise ValueError(f"Missing required column '{col}' in class file")

    annot_df = annot_df.loc[annot_df["filter_result"] == "Isoform", :].copy()
    ref_isoform_gene_map = ref_gtf.get_transcript_gene_id_map()
    known_mask = annot_df['associated_transcript'] != 'novel'
    if known_mask.any():
        mapped_gene = annot_df.loc[known_mask, 'associated_transcript'].map(ref_isoform_gene_map)
        if mapped_gene.isna().any():
            missing_tx = annot_df.loc[known_mask, 'associated_transcript'].loc[mapped_gene.isna()].iloc[0]
            raise KeyError(f"associated_transcript not found in reference GTF: {missing_tx}")
        annot_df.loc[known_mask, 'associated_gene'] = mapped_gene.values

    if 'predicted_NMD' not in annot_df.columns:
        annot_df['predicted_NMD'] = np.nan
    annot_df['predicted_NMD.full'] = annot_df['predicted_NMD'].map(_predicted_nmd_type)
    return annot_df


def _isoform_type_maps(annot_df, gtf):
    isoform_ids = set(annot_df['isoform'])
    transcript_rows = gtf.gtf_df.loc[
        (gtf.gtf_df['feature'] == 'transcript') & (gtf.gtf_df['transcript_id'].isin(isoform_ids)),
        ['transcript_id', 'gene_id', 'start', 'end']
    ].copy()
    missing_isoforms = isoform_ids.difference(set(transcript_rows['transcript_id']))
    if missing_isoforms:
        raise KeyError(f"isoform not found in SQANTI3 GTF: {next(iter(missing_isoforms))}")

    length_map = dict(zip(
        transcript_rows['transcript_id'].values,
        (transcript_rows['start'] - transcript_rows['end']).abs().values,
    ))
    isoform_gene = dict(zip(transcript_rows['transcript_id'].values, transcript_rows['gene_id'].values))

    formal_isoform_type = {}
    gene_best_priority = {}
    gene_best_type = {}
    priority = {'protein_coding': 0, 'lncRNA': 1, 'nonsense_mediated_decay': 2, 'other_ncRNA': 3}
    for iso, raw_type in zip(annot_df['isoform'].values, annot_df['predicted_NMD.full'].values):
        formal_type = raw_type
        if formal_type == 'noncoding_RNA':
            formal_type = 'lncRNA' if length_map[iso] >= 200 else 'other_ncRNA'
        formal_isoform_type[iso] = formal_type
        gene = isoform_gene[iso]
        rank = priority[formal_type]
        if gene not in gene_best_priority or rank < gene_best_priority[gene]:
            gene_best_priority[gene] = rank
            gene_best_type[gene] = formal_type
    return formal_isoform_type, gene_best_type


def _rewrite_isoform_attributes(sdf, isoform_assoc_gene, isoform_type, known_gene_attr, gene_type, add_transcript_type):
    known_gene_ids = set(known_gene_attr)
    attrs = []
    for attr_str in sdf['attributes'].values:
        attr = GTFParser.parse_attributes(attr_str)
        tid = attr["transcript_id"]
        gid = attr["gene_id"]
        if add_transcript_type:
            attr['transcript_type'] = isoform_type[tid]
        annot_gid = isoform_assoc_gene[tid]
        if annot_gid in known_gene_ids:
            attr.update(known_gene_attr[annot_gid])
        else:
            attr['gene_id'] = gid
            attr['gene_name'] = gid
            attr['gene_type'] = gene_type[gid]
        attrs.append(GTFParser.build_attributes(attr))
    out = sdf.copy()
    out.loc[:, 'attributes'] = attrs
    for col in ('transcript_id', 'gene_id', 'gene_name'):
        if col in out.columns:
            del out[col]
    return GTFParser.add_id_columns(out)


def _build_novel_gene_rows(isoform_df, known_gene_attr, gene_type):
    if isoform_df.empty:
        return pd.DataFrame(columns=GTFParser().cols)
    known_gene_ids = set(known_gene_attr)
    novel_gene_df = isoform_df.loc[~isoform_df['gene_id'].isin(known_gene_ids), :]
    if novel_gene_df.empty:
        return pd.DataFrame(columns=GTFParser().cols)

    grouped = novel_gene_df.groupby('gene_id', sort=False).agg(
        chrom=('chrom', 'first'),
        start=('start', 'min'),
        end=('end', 'max'),
        strand=('strand', 'first'),
    ).reset_index()
    grouped.loc[:, 'source'] = 'GTOP'
    grouped.loc[:, 'feature'] = 'gene'
    grouped.loc[:, 'score'] = '.'
    grouped.loc[:, 'frame'] = '.'
    grouped.loc[:, 'attributes'] = [
        GTFParser.build_attributes({'gene_id': gid, 'gene_name': gid, 'gene_type': gene_type[gid]})
        for gid in grouped['gene_id'].values
    ]
    return grouped.loc[:, GTFParser().cols]


def _deduplicate_gene_records(merge_df):
    gene_dup = (merge_df['feature'] == 'gene') & merge_df.loc[merge_df['feature'] == 'gene', 'attributes'].duplicated(keep='first').reindex(merge_df.index, fill_value=False)
    return merge_df.loc[~gene_dup, :].copy()


def _gene_bounds_for_enhanced(ref_gtf, novel_isoform_df):
    if hasattr(ref_gtf, 'transcript_bounds_by_gene'):
        gene_bounds = dict(ref_gtf.transcript_bounds_by_gene)
    else:
        transcript_rows = ref_gtf.gtf_df.loc[ref_gtf.gtf_df['feature'] == 'transcript', ['gene_id', 'start', 'end']]
        gene_bounds = {
            gid: (int(sub['start'].min()), int(sub['end'].max()))
            for gid, sub in transcript_rows.groupby('gene_id', sort=False)
        }

    novel_transcripts = novel_isoform_df.loc[
        novel_isoform_df['feature'] == 'transcript',
        ['gene_id', 'start', 'end']
    ]
    for gid, sub in novel_transcripts.groupby('gene_id', sort=False):
        start = int(sub['start'].min())
        end = int(sub['end'].max())
        if gid in gene_bounds:
            old_start, old_end = gene_bounds[gid]
            gene_bounds[gid] = (min(old_start, start), max(old_end, end))
        else:
            gene_bounds[gid] = (start, end)
    return gene_bounds


def __make_enhanced_gtf(class_file, gtf:GTFParser, ref_gtf:GTFParser, enhanced_gtf_path, enhanced_comment,
                      novel_gtf_path=None, novel_gtf_comment='', gtop_gtf_path=None, gtop_gtf_comment=''):
    annot_df = _load_isoform_annotation(class_file, ref_gtf)
    # known = ["full-splice_match", "incomplete-splice_match"]
    annot_df_novel = annot_df.loc[annot_df['associated_transcript'] == 'novel', :]
    novel_isoform_ids = annot_df_novel['isoform'].unique()
    known_isoform_ids_in_all_isoforms=annot_df.loc[annot_df['associated_transcript'] != 'novel', 'associated_transcript'].unique()

    # get known gene id and other info, e.g., gene name, gene type.
    known_gene_attr=ref_gtf.get_attribute_map(feature='gene',id='gene_id')
    isoform_assoc_gene = dict(zip(annot_df_novel['isoform'], annot_df_novel['associated_gene']))
    logger.info(f'start make gene/transcript type')
    formal_novel_isoform_type, gene_type = _isoform_type_maps(annot_df_novel, gtf)
    logger.info(f'start make novel isoform gtf')

    ldf=gtf.gtf_df
    sdf=ldf.loc[ldf['transcript_id'].isin(novel_isoform_ids), :].copy()
    known_gene_ids = set(known_gene_attr)
    annot_geneid_in_novel_isoforms=set(sdf['transcript_id'].map(isoform_assoc_gene).dropna()).intersection(known_gene_ids)
    sdf = _rewrite_isoform_attributes(
        sdf, isoform_assoc_gene, formal_novel_isoform_type, known_gene_attr, gene_type, add_transcript_type=True
    )

    novel_isoform_gtf=GTFParser()
    novel_isoform_gtf.gtf_df=sdf
    # add novel gene feature line
    ndf = novel_isoform_gtf.gtf_df.copy()
    logger.info(f'start make novel gene gtf record in novel isoform ')
    all_gene_count = ndf['gene_id'].nunique(dropna=True)
    novel_genes_gtf_record_in_novel_isoform_df = _build_novel_gene_rows(ndf, known_gene_attr, gene_type)
    logger.info(f'Novel isoform: all gene: {all_gene_count}; novel gene: {len(novel_genes_gtf_record_in_novel_isoform_df)}')
    # save novel gtf
    logger.info(f'Annotated genes: {len(annot_geneid_in_novel_isoforms)} in GTOP novel isoforms.')
    annot_genes_gtf_record_in_novel_isoform_df = ref_gtf.filter_genes(annot_geneid_in_novel_isoforms)
    logger.info(f'Get {annot_genes_gtf_record_in_novel_isoform_df.shape[0]} annotated genes record for ref gtf.')
    if novel_gtf_path is not None:
        merge_df = _concat_gtf_frames([
            annot_genes_gtf_record_in_novel_isoform_df,
            novel_genes_gtf_record_in_novel_isoform_df,
            novel_isoform_gtf.gtf_df,
        ])
        out_gtf=GTFParser()
        out_gtf.gtf_df=merge_df
        logger.info(f'Novel: start writing gtf')
        out_gtf.write_file(novel_gtf_path,novel_gtf_comment)
        logger.info(f"Novel: Writing GTF: {novel_gtf_path}")
    # save gtop gtf: novel + known (replace with GENCODE v47 transcript annotation).
    if gtop_gtf_path is not None:
        annot_isoform_df=ref_gtf.filter_transcripts_with_genes(known_isoform_ids_in_all_isoforms)
        logger.info(f'Known isoforms ID: {len(known_isoform_ids_in_all_isoforms)} in all GTOP isoforms.')
        merge_df = _concat_gtf_frames([
            annot_isoform_df,
            annot_genes_gtf_record_in_novel_isoform_df,
            novel_genes_gtf_record_in_novel_isoform_df,
            novel_isoform_gtf.gtf_df,
        ])
        ## remove duplicate gene record
        merge_df = _deduplicate_gene_records(merge_df)
        out_gtf=GTFParser()
        out_gtf.gtf_df=merge_df
        logger.info(f'GTOP: start writing gtf')
        out_gtf.write_file(gtop_gtf_path,gtop_gtf_comment)
        logger.info(f"GTOP: Writing GTF: {gtop_gtf_path}")
    # save enhanced gtf: novel + GENCODE v47 transcript annotation.
    if enhanced_gtf_path is not None:
        logger.info(f'Enhanced: start writing streamed gtf')
        if hasattr(ref_gtf, 'write_enhanced_gtf'):
            gene_bounds = _gene_bounds_for_enhanced(ref_gtf, novel_isoform_gtf.gtf_df)
            ref_gtf.write_enhanced_gtf(
                enhanced_gtf_path,
                enhanced_comment,
                gene_bounds,
                [novel_genes_gtf_record_in_novel_isoform_df, novel_isoform_gtf.gtf_df],
            )
        else:
            merge_df = _concat_gtf_frames([
                ref_gtf.gtf_df,
                novel_genes_gtf_record_in_novel_isoform_df,
                novel_isoform_gtf.gtf_df,
            ])
            out_gtf=GTFParser()
            out_gtf.gtf_df=merge_df
            out_gtf.update_gene_coor()
            out_gtf.write_file(enhanced_gtf_path,enhanced_comment)
        logger.info(f"Enhanced: Writing GTF: {enhanced_gtf_path}")

def __make_GTOP_gtf(class_file, gtf:GTFParser, ref_gtf:GTFParser, gtop_gtf_path=None, gtop_gtf_comment=''):
    annot_df = _load_isoform_annotation(class_file, ref_gtf)
    # known = ["full-splice_match", "incomplete-splice_match"]
    isoform_ids = annot_df['isoform'].unique()
    logger.info(f'load {len(isoform_ids)} isoforms.')
    # get known gene id and other info, e.g., gene name, gene type.
    known_gene_attr=ref_gtf.get_attribute_map(feature='gene',id='gene_id')
    isoform_assoc_gene = dict(zip(annot_df['isoform'], annot_df['associated_gene']))
    logger.info(f'start make gene/transcript type')
    formal_isoform_type, gene_type = _isoform_type_maps(annot_df, gtf)
    logger.info(f'start make novel isoform gtf')
    ldf=gtf.gtf_df
    sdf=ldf.loc[ldf['transcript_id'].isin(isoform_ids), :].copy()
    known_gene_ids = set(known_gene_attr)
    annot_geneid_in_isoforms=set(sdf['transcript_id'].map(isoform_assoc_gene).dropna()).intersection(known_gene_ids)
    sdf = _rewrite_isoform_attributes(
        sdf, isoform_assoc_gene, formal_isoform_type, known_gene_attr, gene_type, add_transcript_type=False
    )

    isoform_gtf=GTFParser()
    isoform_gtf.gtf_df=sdf
    # add novel gene feature line
    ndf = isoform_gtf.gtf_df.copy()
    logger.info(f'start make novel gene gtf record in novel isoform ')
    all_gene_count = ndf['gene_id'].nunique(dropna=True)
    novel_genes_gtf_record_in_novel_isoform_df = _build_novel_gene_rows(ndf, known_gene_attr, gene_type)
    logger.info(f'GTOP isoform: all gene: {all_gene_count}; novel gene: {len(novel_genes_gtf_record_in_novel_isoform_df)}')
    # save novel gtf
    logger.info(f'Annotated genes: {len(annot_geneid_in_isoforms)} in GTOP isoforms.')
    annot_genes_gtf_record_in_isoform_df = ref_gtf.filter_genes(annot_geneid_in_isoforms)
    logger.info(f'Get {annot_genes_gtf_record_in_isoform_df.shape[0]} annotated genes record for ref gtf.')
    if gtop_gtf_path is not None:
        merge_df = _concat_gtf_frames([
            annot_genes_gtf_record_in_isoform_df,
            novel_genes_gtf_record_in_novel_isoform_df,
            isoform_gtf.gtf_df,
        ])
        out_gtf=GTFParser()
        out_gtf.gtf_df=merge_df
        logger.info(f'GTOP: start writing gtf')
        out_gtf.write_file(gtop_gtf_path,gtop_gtf_comment)
        logger.info(f"GTOP: Writing GTF: {gtop_gtf_path}")

def __merge_cds_to_enhanced_gtf(cds_path,enhanced_gtf_path,out_path,enhanced_with_CDS_comment):
    trans_attr = {}
    with open(enhanced_gtf_path) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != 'transcript':
                continue
            tid = GTFParser.get_attr_value(fields[8], 'transcript_id')
            if tid and not tid.startswith('ENS'):
                trans_attr[tid] = fields[8]

    logger.info(f'load {len(trans_attr)} novel isoforms for CDS merge')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    copied = 0
    cds_written = 0
    with open(out_path, "w") as fout:
        fout.write(enhanced_with_CDS_comment or '')
        with open(enhanced_gtf_path) as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                fout.write(line)
                copied += 1
        with open(cds_path) as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != 'CDS':
                    continue
                tid = GTFParser.get_attr_value(fields[8], 'transcript_id')
                if tid not in trans_attr:
                    continue
                fields[8] = trans_attr[tid]
                fout.write("\t".join(fields) + "\n")
                cds_written += 1
    logger.info(f'Wrote {out_path}: copied={copied:,}; CDS_added={cds_written:,}')


class GTFAnnotation:
    def __init__(self, gtf_file):
        self.gtf_file = gtf_file

    @staticmethod
    def normalize_chr(c):
        if pd.isna(c):
            return c
        c = str(c).strip()
        if c.startswith("chr"):
            return c
        if re.match(r"^\d+$", c):
            return "chr" + c
        if c.lower() in ["m", "mt", "chrM", "chrMT", "mitochondria", "mitochondrial"]:
            return "chrM"
        if c.lower() == "x":
            return "chrX"
        if c.lower() == "y":
            return "chrY"
        return c

    def make_annotation_tables(self):
        gtf_file = self.gtf_file
        wdir = os.path.dirname(gtf_file)
        f = os.path.basename(gtf_file)
        gtf_file = f'{wdir}/{f}'
        gene_annot_file = f'{wdir}/{f[:-4]}.gene_annot.txt'
        gene_transcript_file = f'{wdir}/{f[:-4]}.gene_transcript_map.txt'
        autosome_file = f'{wdir}/{f[:-4]}.autosome_isoform.csv'
        autosomes = {f"chr{i}" for i in range(1, 23)}
        seen_tx = set()
        seen_auto = set()

        with open(gtf_file) as gtf, \
                open(gene_annot_file, "w") as gene_out, \
                open(gene_transcript_file, "w") as tx_out, \
                open(autosome_file, "w") as auto_out:
            tx_out.write("transcript_id\tchr\tstart\tend\tgene_id\tgene_name\ttranscript_type\tgene_type\n")
            auto_out.write("transcript_id,gene_id,chromosome\n")
            for line in gtf:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                feature = fields[2].strip()
                if feature not in {'gene', 'transcript', 'mRNA'}:
                    continue
                attr_dict = GTFParser.parse_attributes(fields[8])
                gid = attr_dict.get("gene_id")

                if feature == 'gene':
                    gene_out.write(
                        f"{attr_dict.get('gene_name')}\t{attr_dict.get('gene_type')}\t{gid}\t"
                        f"{fields[0]}\t{fields[3]}\t{fields[4]}\t{fields[6]}\n"
                    )
                    continue

                tid = attr_dict.get("transcript_id")
                if not tid or not gid:
                    continue
                if tid not in seen_tx:
                    tx_out.write(
                        f"{tid}\t{fields[0]}\t{fields[3]}\t{fields[4]}\t{gid}\t"
                        f"{attr_dict.get('gene_name')}\t{attr_dict.get('transcript_type')}\t"
                        f"{attr_dict.get('gene_type')}\n"
                    )
                    seen_tx.add(tid)
                chromosome = self.normalize_chr(fields[0])
                auto_key = (tid, gid, chromosome)
                if chromosome in autosomes and auto_key not in seen_auto:
                    auto_out.write(f"{tid},{gid},{chromosome}\n")
                    seen_auto.add(auto_key)
        logger.info(f'saved annotation tables for {gtf_file}')

    def conduct_gene_length(self):
        self.make_length_tables()

    def make_length_tables(self):
        gtf_file = self.gtf_file
        wdir = os.path.dirname(gtf_file)
        f = os.path.basename(gtf_file)
        gene_length_file = f'{wdir}/{f[:-4]}.gene_length'
        transcript_length_file = f'{wdir}/{f[:-4]}.transcript_length'
        gene_exons = {}
        gene2iso = {}
        iso2gene = {}
        iso_length = {}

        with open(gtf_file) as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != 'exon':
                    continue
                attr = fields[8]
                gid = GTFParser.get_attr_value(attr, 'gene_id')
                tid = GTFParser.get_attr_value(attr, 'transcript_id')
                if not gid or not tid:
                    continue
                start0 = int(fields[3]) - 1
                end = int(fields[4])
                length = end - start0
                iso2gene[tid] = gid
                gene2iso.setdefault(gid, set()).add(tid)
                iso_length[tid] = iso_length.get(tid, 0) + length
                gene_exons.setdefault(gid, []).append((fields[0], start0, end, fields[6]))

        def merged_exon_length(exons):
            if not exons:
                return 0
            exons = sorted(exons, key=lambda x: (x[0], x[3], x[1], x[2]))
            total = 0
            cur_chr, cur_start, cur_end, cur_strand = exons[0]
            for chrom, start, end, strand in exons[1:]:
                if chrom == cur_chr and strand == cur_strand and start <= cur_end:
                    cur_end = max(cur_end, end)
                else:
                    total += cur_end - cur_start
                    cur_chr, cur_start, cur_end, cur_strand = chrom, start, end, strand
            total += cur_end - cur_start
            return total

        with open(transcript_length_file, "w") as out:
            out.write('isoform\tgene\tlength\n')
            for tid, gid in iso2gene.items():
                out.write(f'{tid}\t{gid}\t{iso_length[tid]}\n')

        with open(gene_length_file, "w") as out:
            out.write('gene\tmean\tmedian\tlongest_isoform\tmerged\n')
            for gid, isoforms in gene2iso.items():
                lengths = [iso_length[tid] for tid in isoforms]
                out.write(
                    f'{gid}\t{int(np.mean(lengths))}\t{int(np.median(lengths))}\t'
                    f'{max(lengths)}\t{merged_exon_length(gene_exons.get(gid, []))}\n'
                )
        logger.info(f'saved length tables for {gtf_file}')

    def make_gene_transcript_tab(self):
        gtf_file = self.gtf_file
        wdir = os.path.dirname(gtf_file)
        f = os.path.basename(gtf_file)
        gtf_file = f'{wdir}/{f}'
        out_file = f'{wdir}/{f[:-4]}.gene_transcript_map.txt'
        seen = set()
        with open(gtf_file) as gtf, open(out_file, "w") as out:
            out.write("transcript_id\tchr\tstart\tend\tgene_id\tgene_name\ttranscript_type\tgene_type\n")
            for line in gtf:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue
                attr_field = fields[8]
                attr_dict = {}
                for attr in attr_field.strip().split(";"):
                    attr = attr.strip()
                    if not attr:
                        continue
                    try:
                        key, value = attr.split(" ", 1)
                        attr_dict[key] = value.strip('"')
                    except ValueError:
                        continue
                tid = attr_dict.get("transcript_id")
                gid = attr_dict.get("gene_id")
                gname = attr_dict.get("gene_name")
                transcript_type = attr_dict.get("transcript_type")
                gene_type = attr_dict.get("gene_type")
                chr = fields[0]
                start = fields[3]
                end = fields[4]
                if tid and gid and tid not in seen:
                    out.write(f"{tid}\t{chr}\t{start}\t{end}\t{gid}\t{gname}\t{transcript_type}\t{gene_type}\n")
                    seen.add(tid)

    def make_gene_annot_tab(self):
        gtf_file = self.gtf_file
        wdir = os.path.dirname(gtf_file)
        f = os.path.basename(gtf_file)
        gtf_file = f'{wdir}/{f}'
        out_file = f'{wdir}/{f[:-4]}.gene_annot.txt'
        seen = set()
        with open(gtf_file) as gtf, open(out_file, "w") as out:
            for line in gtf:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue
                if fields[2].strip() != 'gene':
                    continue
                attr_field = fields[8]
                attr_dict = {}
                for attr in attr_field.strip().split(";"):
                    attr = attr.strip()
                    if not attr:
                        continue
                    try:
                        key, value = attr.split(" ", 1)
                        attr_dict[key] = value.strip('"')
                    except ValueError:
                        continue
                gid = attr_dict.get("gene_id")
                gname = attr_dict.get("gene_name")
                gtype = attr_dict.get("gene_type")
                chr = fields[0]
                st = fields[3]
                ed = fields[4]
                chain = fields[6]
                out.write(f"{gname}\t{gtype}\t{gid}\t{chr}\t{st}\t{ed}\t{chain}\n")
        logger.info(f'complete {gtf_file}')

    def make_autosome_transcript_tab(self):
        keep_autosomes_only = True
        protein_coding_only = False
        # wdir = f'/media/dubai/home/xuechao/project/GMTiP-RNA/20251031/long_read/HPC/output/enhanced_gtf'
        # files=['GTOP_novel-GENCODE_v47.gtf','GTOP_all.gtf']
        gtf_file = self.gtf_file
        wdir = os.path.dirname(gtf_file)
        f = os.path.basename(gtf_file)
        gtf_file = f'{wdir}/{f}'
        out_tab = f'{wdir}/{f[:-4]}.autosome_isoform.csv'
        gtf = pd.read_csv(
            gtf_file,
            sep="\t",
            comment="#",
            header=None,
            names=["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
        )
        transcripts = gtf[gtf["feature"].isin(["transcript", "mRNA"])]  #
        if protein_coding_only:
            transcripts = transcripts[transcripts["attribute"].str.contains('gene_biotype "protein_coding"')]
        transcripts["gene_id"] = transcripts["attribute"].str.extract('gene_id "([^"]+)"')
        transcripts["transcript_id"] = transcripts["attribute"].str.extract('transcript_id "([^"]+)"')

        def normalize_chr(c):
            if pd.isna(c):
                return c
            c = str(c).strip()
            if c.startswith("chr"):
                return c
            if re.match(r"^\d+$", c):
                return "chr" + c
            if c.lower() in ["m", "mt", "chrM", "chrMT", "mitochondria", "mitochondrial"]:
                return "chrM"
            # sex chromosomes
            if c.lower() in ["x"]:
                return "chrX"
            if c.lower() in ["y"]:
                return "chrY"
            return c

        transcripts["chromosome"] = transcripts["chr"].apply(normalize_chr)
        if keep_autosomes_only:
            autosomes = [f"chr{i}" for i in range(1, 23)]
            transcripts = transcripts[transcripts["chromosome"].isin(autosomes)]
        output_df = transcripts[["transcript_id", "gene_id", "chromosome"]].drop_duplicates()
        output_df.to_csv(out_tab, index=False)
        logger.info(f'save to {out_tab}')

    def index_gtf(self):
        gtf_file = self.gtf_file
        wdir = os.path.dirname(gtf_file)
        f = os.path.basename(gtf_file)
        gtf_file = f'{wdir}/{f}'
        gtf_kw = f[:-4]
        gz_gtf = f'{wdir}/{gtf_kw}.sorted.gtf.gz'
        q_gtf = shlex.quote(gtf_file)
        q_gz = shlex.quote(gz_gtf)
        logger.info(f"Sorting and compressing {gtf_file}...")
        subprocess.run(
            f"{LOAD_HTSLIB_ENV_CMD} && sort -k1,1 -k4,4n {q_gtf} | bgzip -c > {q_gz}",
            shell=True,
            executable="/bin/bash",
            check=True,
        )
        logger.info(f"Indexing {gz_gtf}...")
        subprocess.run(
            f"{LOAD_HTSLIB_ENV_CMD} && tabix -p gff {q_gz}",
            shell=True,
            executable="/bin/bash",
            check=True,
        )
        logger.info("Done!")

    def make_gtf_fa(self):
        gtf_file = self.gtf_file
        wdir = os.path.dirname(gtf_file)
        f = os.path.basename(gtf_file)
        f_prefix = f[:-4]
        cmd = f'{GFFREAD_BIN} {gtf_file} -g {REF_GENOME_FASTA} -w {wdir}/{f_prefix}.fa'
        logger.info(f"cmd: {cmd}")
        subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)

    def annotate_gtf_task(self):
        self.make_gtf_fa()
        self.make_annotation_tables()
        self.conduct_gene_length()
        self.index_gtf()


def make_enhanced_gtf_and_annotation_task(LRS_out_dir):
    enhanced_gtf_dir=f'{LRS_out_dir}/{ENHANCED_GTF_TAG}'
    sqanti_dir=f'{LRS_out_dir}/{CLEAN_SQANTI3_TAG}'

    sqanti3_gtf=f'{sqanti_dir}/filter.clean.gtf'
    class_file=f'{sqanti_dir}/filter.clean.RulesFilter_result_classification.txt'
    cds_gtf=f'{sqanti_dir}/filter.clean.cds.gff3'
    ref_gtf=REF_GTF
    novel_gtf=f'{enhanced_gtf_dir}/GTOP_novel.gtf'
    novel_gtf_comment = f'# GTOP novel isoform annotation \n'
    gtop_gtf=f'{enhanced_gtf_dir}/GTOP.gtf'
    gtop_comment = f'# GTOP isoform annotation: GTOP novel + annotated \n'
    gtop_all_gtf=f'{enhanced_gtf_dir}/GTOP_all.gtf'
    gtop_all_comment = f'# GTOP isoform annotation: GTOP all \n'
    enhanced_gtf = f'{enhanced_gtf_dir}/GTOP_novel-GENCODE_v47.gtf'
    enhanced_comment = f'# Enhanced isoform annotation: GTOP novel + GENCODE v47 \n'
    enhanced_with_CDS_gtf = f'{enhanced_gtf_dir}/GTOP_novel-GENCODE_v47.with_CDS.gtf'
    enhanced_with_CDS_comment = f'# Enhanced isoform annotation: GTOP novel + GENCODE v47 with CDS\n'
    os.makedirs(os.path.dirname(enhanced_gtf), exist_ok=True)
    #
    with timed_step(f"load SQANTI3 GTF {sqanti3_gtf}"):
        gtf_par=GTFParser(sqanti3_gtf).load_df()
    # correct gtf generated by SQANTI3, because the start/end position of transcript in "-" chain is error.
    with timed_step("correct SQANTI3 transcript boundaries"):
        gtf_par.correct_transcript_boundaries()
    #
    with timed_step(f"index reference GTF {ref_gtf}"):
        ref_gtf_par=ReferenceGTFIndex(ref_gtf).load_index()
    with timed_step("make enhanced/novel/GTOP GTF"):
        __make_enhanced_gtf(class_file,gtf_par,ref_gtf_par,enhanced_gtf,enhanced_comment,
                          novel_gtf ,novel_gtf_comment,
                          gtop_gtf,gtop_comment)
    # generate GTOP raw isoform with known genes.
    with timed_step("make GTOP_all GTF"):
        __make_GTOP_gtf(class_file,gtf_par,ref_gtf_par, gtop_all_gtf,gtop_all_comment)
    logger.info(f"Done. Novel isoforms GTF written to: {gtop_all_gtf}")
    # merge CDS to enhanced gtf.
    with timed_step("merge CDS into enhanced GTF"):
        __merge_cds_to_enhanced_gtf(cds_gtf,enhanced_gtf,enhanced_with_CDS_gtf,enhanced_with_CDS_comment)
    logger.info(f"Done. CDS written to: {enhanced_with_CDS_gtf}")
    gtf_files=[novel_gtf,gtop_gtf, gtop_all_gtf, enhanced_gtf, enhanced_with_CDS_gtf]
    with timed_step(f"nested-sort and replace {len(gtf_files)} generated GTF files"):
        for gtf_file in gtf_files:
            sort_gtf_nested_inplace(gtf_file)
    # annotate GTF
    # gtf_files.append(REF_GTF)
    workers = int(os.environ.get("GTOP_ANNOTATE_WORKERS", str(min(3, len(gtf_files)))))
    workers = max(1, min(workers, len(gtf_files)))
    with timed_step(f"annotate {len(gtf_files)} GTF files with {workers} worker(s)"):
        if workers == 1:
            for gtf in gtf_files:
                logger.info(f"start to annotate {gtf}")
                GTFAnnotation(gtf).annotate_gtf_task()
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                future_map = {executor.submit(GTFAnnotation(gtf).annotate_gtf_task): gtf for gtf in gtf_files}
                for future in as_completed(future_map):
                    gtf = future_map[future]
                    future.result()
                    logger.info(f"completed annotation {gtf}")


def extract_transcripts_gtf(gtf_path, transcript_ids, out_path):
    ids = set(transcript_ids)
    all_tids = set()
    kept_tids = set()
    with open(gtf_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            tid = None
            for item in parts[8].split(";"):
                item = item.strip()
                if item.startswith("transcript_id"):
                    tid = item.split(" ", 1)[1].strip('"')
                    break
            if tid is None:
                continue
            all_tids.add(tid)
            if tid in ids:
                kept_tids.add(tid)
                fout.write(line)
    logger.info(f"Total transcripts in GTF: {len(all_tids)}")
    logger.info(f"Requested transcripts: {len(ids)}")
    logger.info(f"Successfully extracted transcripts: {len(kept_tids)}")


def fix_transcript_coords(gtf_in, gtf_out):
    coords = {}
    records = []

    with open(gtf_in) as f:
        for line in f:
            if line.startswith("#"):
                records.append((None, line))
                continue
            p = line.rstrip("\n").split("\t")
            tid = None
            for x in p[8].split(";"):
                x = x.strip()
                if x.startswith("transcript_id"):
                    tid = x.split(" ", 1)[1].strip('"')
                    break
            records.append((tid, p))
            if p[2] == "exon" and tid is not None:
                coords.setdefault(tid, []).extend((int(p[3]), int(p[4])))

    with open(gtf_out, "w") as f:
        for tid, rec in records:
            if tid is None:
                f.write(rec)
            else:
                if rec[2] == "transcript" and tid in coords:
                    rec[3], rec[4] = map(str, (min(coords[tid]), max(coords[tid])))
                f.write("\t".join(rec) + "\n")


def transfer_gene_id(ref_gtf, in_gtf, out_gtf):
    logger = logging.getLogger(__name__)
    tx2gene = {}
    genes = set()

    with open(ref_gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            tid = gid = None
            for x in p[8].split(";"):
                x = x.strip()
                if x.startswith("transcript_id"):
                    tid = x.split(" ", 1)[1].strip('"')
                elif x.startswith("gene_id"):
                    gid = x.split(" ", 1)[1].strip('"')
            if tid and gid and tid not in tx2gene:
                tx2gene[tid] = gid
                genes.add(gid)

    logger.info(f"Transcripts found in reference GTF: {len(tx2gene)}")
    logger.info(f"Genes found in reference GTF: {len(genes)}")

    updated = set()
    with open(in_gtf) as fin, open(out_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                fout.write(line)
                continue

            tid = None
            attrs = []
            old_gene = None

            for x in p[8].split(";"):
                x = x.strip()
                if not x:
                    continue
                if x.startswith("transcript_id"):
                    tid = x.split(" ", 1)[1].strip('"')
                    attrs.append(x)
                elif x.startswith("gene_id"):
                    old_gene = x
                else:
                    attrs.append(x)

            if tid in tx2gene:
                attrs.insert(0, f'gene_id "{tx2gene[tid]}"')
                updated.add(tid)
            elif old_gene is not None:
                attrs.insert(0, old_gene)

            p[8] = "; ".join(attrs) + ";"
            p[1] = 'GTOP'
            fout.write("\t".join(p) + "\n")

    logger.info(f"Transcripts updated or filled: {len(updated)}")


def load_gene_ids_from_gtf(gtf_path):
    gene_ids = set()
    with open(gtf_path) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            gid = GTFParser.get_attr_value(fields[8], 'gene_id')
            if gid:
                gene_ids.add(gid)
    return gene_ids


def replace_gtf_ids(in_gtf, out_gtf):
    """
    Replace transcript_id and gene_id in a GTF file.

    Rules:
    - transcript_id: transcript_number -> GTOPTXXXXXXXX
    - gene_id:
        1) transcript_number -> same as transcript_id
        2) LOC_number -> GTOPGXXXXXXXX
    """

    def format_id(prefix, number):
        return f"{prefix}{int(number):09d}"

    tx_re = re.compile(r'transcript_id "([^"]+)"')
    gene_re = re.compile(r'gene_id "([^"]+)"')
    tx_sub_re = re.compile(r'transcript_id "[^"]+"')
    gene_sub_re = re.compile(r'gene_id "[^"]+"')
    with open(in_gtf) as fin, open(out_gtf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                fout.write(line)
                continue

            attr = fields[8]

            # ---- transcript_id ----
            m_tx = tx_re.search(attr)
            tx_new = None
            if m_tx:
                tx_id = m_tx.group(1)
                if tx_id.startswith("transcript_"):
                    num = tx_id.split("_", 1)[1]
                    tx_new = format_id("GTOPT", num)
                    attr = tx_sub_re.sub(
                        f'transcript_id "{tx_new}"',
                        attr
                    )

            # ---- gene_id ----
            m_gene = gene_re.search(attr)
            if m_gene:
                gene_id = m_gene.group(1)

                if gene_id.startswith("transcript_") and tx_new is not None:
                    # rule 1: same as transcript_id
                    attr = gene_sub_re.sub(
                        f'gene_id "{tx_new}"',
                        attr
                    )

                elif gene_id.startswith("LOC_"):
                    # rule 2: LOC_number -> GTOPGXXXXXXXX
                    num = gene_id.split("_", 1)[1]
                    gene_new = format_id("GTOPG", num)
                    attr = gene_sub_re.sub(
                        f'gene_id "{gene_new}"',
                        attr
                    )

            fields[8] = attr
            fout.write("\t".join(fields) + "\n")

def replace_sqanti3_annot_txid(in_annot, out_annot, tool_support_path):
    os.makedirs(os.path.dirname(out_annot),exist_ok=True)
    def format_id(prefix, number):
        return f"{prefix}{int(number):09d}"
    df=pd.read_csv(in_annot, sep='\t')
    # add tool support info
    tool_df=pd.read_csv(tool_support_path, sep='\t', index_col=0)
    tool_cols = tool_df.columns.tolist()
    support_tools = {
        idx: ",".join([tool for tool, val in zip(tool_cols, row) if val == 1])
        for idx, row in zip(tool_df.index, tool_df[tool_cols].to_numpy())
    }
    df['support_tools']=df['isoform'].map(support_tools)
    if df['support_tools'].isna().any():
        missing_isoform = df.loc[df['support_tools'].isna(), 'isoform'].iloc[0]
        raise KeyError(f"isoform not found in support table: {missing_isoform}")
    df['isoform']=df['isoform'].apply(lambda x:format_id('GTOPT',x.split('_')[1]))
    df.to_csv(out_annot, sep='\t', index=False)

def replace_faa_ids(in_faa, out_faa):
    os.makedirs(os.path.dirname(out_faa),exist_ok=True)
    tid_re = re.compile(r'transcript_(\d+)')
    def replace_tid(line):
        return tid_re.sub(
            lambda m: f'GTOPT{int(m.group(1)):09d}',
            line
        )
    with open(in_faa, 'r') as fin, open(out_faa, "w") as fout:
        for line in fin:
            val=line.rstrip()
            if val.startswith(">"):
                val=replace_tid(val)
            fout.write(val+'\n')
    logger.info(f'save faa to {out_faa}')


def replace_flair_quant_ids(in_quant, clean_quant):
    df=pd.read_csv(in_quant, sep='\t',index_col=0)
    def format_id(prefix, number):
        return f"{prefix}{int(number):09d}"
    df.index=df.index.map(lambda x: format_id('GTOPT',x.split('_')[1]))
    os.makedirs(os.path.dirname(clean_quant),exist_ok=True)
    df.to_csv(clean_quant, sep='\t')
    pass


def assign_new_gene_loci(LRS_out_dir):
    clean_sqanti3_dir=f'{LRS_out_dir}/{CLEAN_SQANTI3_TAG}'
    sqant_name= FINAL_SQANTI3_TAG

    tag=''
    in_faa=f'{LRS_out_dir}/sqanti3_merged/filtered.faa'
    sqanti3_annot_path=f'{LRS_out_dir}/{sqant_name}/filtered{tag}.RulesFilter_result_classification.txt'
    in_gtf=f'{LRS_out_dir}/{sqant_name}/filtered{tag}.gtf'
    in_cds_gtf=f'{LRS_out_dir}/{sqant_name}/filtered{tag}.cds.gff3'
    fixed_in_gtf=f'{LRS_out_dir}/{sqant_name}/filtered{tag}.fixed.gtf'
    isoform_in_novel_gene_gtf=f'{LRS_out_dir}/{sqant_name}/filtered{tag}.fixed.isoform_in_novel_gene.gtf'
    novel_gene_gtf=f'{LRS_out_dir}/{sqant_name}/filtered{tag}.fixed.buildLoci_novel_gene.gtf'
    out_gtf=f'{LRS_out_dir}/{sqant_name}/filtered{tag}.final.gtf'
    in_quant=f'{LRS_out_dir}/sqanti3_filter_1/flair_quant/combined/transcript.count.flair.tsv'
    #
    # fix gtf, because transcript coor from SQANTI3 is wrong.
    logger.info(f'start fix gtf')
    fix_transcript_coords(in_gtf, fixed_in_gtf)
    # get known gene id
    gids=load_gene_ids_from_gtf(REF_GTF)
    logger.info(f'load {len(gids)} genes from GENCODE')
    # extract isoforms in novel gene
    sq_df = pd.read_csv(sqanti3_annot_path,sep='\t',index_col=0)
    sq_df = sq_df.loc[sq_df['filter_result']=='Isoform',:]
    print(sq_df.shape)
    sq_df=sq_df.loc[sq_df['associated_transcript']=='novel', :]
    isoform_ING_ids=sq_df.loc[~sq_df['associated_gene'].isin(gids), :].index.tolist()
    logger.info(f'load {len(isoform_ING_ids)} isoforms in novel gene')
    # save to tmp gtf
    extract_transcripts_gtf(fixed_in_gtf,isoform_ING_ids,isoform_in_novel_gene_gtf)
    # build gene loci
    buildLoci_cmd=f'perl {CONFIG_BUILD_LOCI}'
    bedtools_cmd='bedtools'
    cmd=f'{bedtools_cmd} intersect -s -wao -a {isoform_in_novel_gene_gtf} -b {isoform_in_novel_gene_gtf} | {buildLoci_cmd} - > {novel_gene_gtf}'
    print(cmd)
    __import__('subprocess').run(cmd, shell=True, executable='/bin/bash', check=True)
    # assign novel gene id to gtf
    transfer_gene_id(novel_gene_gtf,fixed_in_gtf,out_gtf)
    #
    # # rename gene id/transcript id
    tool_support_path=f'{LRS_out_dir}/raw_isoform/raw_merged.meta.tsv'
    os.makedirs(clean_sqanti3_dir,exist_ok=True)
    clean_gtf=f'{clean_sqanti3_dir}/filter.clean.gtf'
    clean_cds_gtf=f'{clean_sqanti3_dir}/filter.clean.cds.gff3'
    clean_annot=f'{clean_sqanti3_dir}/filter.clean.RulesFilter_result_classification.txt'
    clean_faa=f'{clean_sqanti3_dir}/filter.clean.faa'
    clean_quant=f'{clean_sqanti3_dir}/filter.clean.flair_quant.tsv'
    replace_gtf_ids(out_gtf,clean_gtf)
    replace_gtf_ids(in_cds_gtf,clean_cds_gtf)
    replace_sqanti3_annot_txid(sqanti3_annot_path,clean_annot,tool_support_path)
    replace_faa_ids(in_faa,clean_faa)
    replace_flair_quant_ids(in_quant, clean_quant)
    pass


def export_reference_files(root_output_dir):
    project = Path(CONFIG_PROJECT_DIR)
    merged = Path(root_output_dir)
    pairs = []
    for prefix in ['GTOP', 'GTOP_novel-GENCODE_v47']:
        for suffix in ['.gtf', '.fa', '.gene_transcript_map.txt', '.transcript_length']:
            pairs.append((merged / ENHANCED_GTF_TAG / (prefix + suffix),
                          project / 'release/gtf' / (prefix + suffix)))
    pairs.append((merged / CLEAN_SQANTI3_TAG / 'filter.clean.faa',
                  project / 'release/LRS_assembly/sqanti3/filter.clean.faa'))
    missing = [str(src) for src, _ in pairs if not src.is_file() or not src.stat().st_size]
    if missing:
        raise FileNotFoundError('Missing reference outputs: ' + ', '.join(missing))
    for src, dst in pairs:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


if __name__ == '__main__':
    ARGS = sys.argv[1:]
    STEP = ARGS[0]
    CONF_CSV = ARGS[1]
    ROOT_OUTPUT_DIR = ARGS[2]
    LOG_NAME = ARGS[3]
    N_TASK = int(ARGS[4])
    NT_PER_TASK = int(ARGS[5])
    if len(ARGS) >= 9:
        configure_output_tags(
            clean_sqanti3_tag=ARGS[6],
            final_sqanti3_tag=ARGS[7],
            enhanced_gtf_tag=ARGS[8],
        )
    else:
        configure_output_tags()
    # single task
    if STEP == 'enhanced_gtf':
        assign_new_gene_loci(ROOT_OUTPUT_DIR)
        make_enhanced_gtf_and_annotation_task(ROOT_OUTPUT_DIR)
        if ENHANCED_GTF_TAG == "enhanced_gtf":
            export_reference_files(ROOT_OUTPUT_DIR)

    # assign_new_gene_loci(LRS_out_dir)
    # make_enhanced_gtf_and_annotation_task(LRS_out_dir)
