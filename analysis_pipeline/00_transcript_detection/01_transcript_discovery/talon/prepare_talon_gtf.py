# -*- coding: utf-8 -*-
"""
Create a TALON initialization GTF by adding only missing gene records.

The output follows nested GTF order: gene -> transcript -> exon/children.
Existing records are written unchanged. Synthetic gene records are added only
for gene_ids that do not already have a gene feature.
"""

import logging
import os
import re
import sys
import time
from collections import defaultdict


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

GTF_NESTED_SORT_PROGRESS_INTERVAL = int(os.environ.get("GTOP_GTF_SORT_PROGRESS_INTERVAL", "500000"))
GTF_TRANSCRIPT_FEATURES = {'transcript', 'mRNA'}
GTF_GENE_ID_RE = re.compile(r'(?:^|;\s*)gene_id "([^"]*)"')
GTF_TRANSCRIPT_ID_RE = re.compile(r'(?:^|;\s*)transcript_id "([^"]*)"')
GTF_ATTR_RE = re.compile(r'([A-Za-z_][A-Za-z0-9_]*) "([^"]*)"')
GTF_GENE_ATTR_KEYS = {
    'gene_id',
    'gene_name',
    'gene_type',
    'gene_biotype',
    'gene_version',
}


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


def _parse_gene_attrs(attr_text, gene_id):
    seen = set()
    pairs = []
    for key, value in GTF_ATTR_RE.findall(attr_text):
        if key not in GTF_GENE_ATTR_KEYS or key in seen:
            continue
        pairs.append((key, value))
        seen.add(key)
    if 'gene_id' not in seen:
        pairs.insert(0, ('gene_id', gene_id))
    return ' '.join(f'{key} "{value}";' for key, value in pairs)


def _update_missing_gene_info(missing_gene_info, gene_id, fields):
    start = int(fields[3])
    end = int(fields[4])
    info = missing_gene_info.setdefault(
        gene_id,
        {
            'chrom': fields[0],
            'source': fields[1],
            'start': start,
            'end': end,
            'strand': fields[6],
            'attr_text': fields[8],
        },
    )
    if start < info['start']:
        info['start'] = start
    if end > info['end']:
        info['end'] = end


def _make_synthetic_gene_record(gene_id, gene_info, line_no):
    attrs = _parse_gene_attrs(gene_info['attr_text'], gene_id)
    line = '\t'.join([
        gene_info['chrom'],
        gene_info['source'],
        'gene',
        str(gene_info['start']),
        str(gene_info['end']),
        '.',
        gene_info['strand'],
        '.',
        attrs,
    ])
    chrom_key = _gtf_nested_chrom_sort_key(gene_info['chrom'])
    sort_key = (chrom_key, gene_info['start'], gene_info['end'], line_no)
    return NestedGTFRecord(line, sort_key, sort_key)


def _load_group_and_patch_gtf(gtf_path, progress_interval):
    comments = []
    gene_records = defaultdict(list)
    transcript_records = defaultdict(list)
    transcript_children = defaultdict(list)
    gene_to_transcripts = defaultdict(set)
    orphan_records = []
    gene_anchor_keys = {}
    transcript_anchor_keys = {}
    missing_gene_info = {}
    observed_gene_ids = set()
    record_count = 0
    start_time = time.perf_counter()

    logger.info(f"Reading and grouping GTF for TALON nested sort: {gtf_path}")
    with open(gtf_path, encoding='utf-8') as fin:
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
                observed_gene_ids.add(gene_id)
                gene_records[gene_id].append(record)
                _update_gtf_nested_anchor(gene_anchor_keys, gene_id, sort_key)
            elif feature in GTF_TRANSCRIPT_FEATURES and transcript_id:
                transcript_records[transcript_id].append(record)
                _update_gtf_nested_anchor(transcript_anchor_keys, transcript_id, sort_key)
                if gene_id:
                    gene_to_transcripts[gene_id].add(transcript_id)
                    _update_gtf_nested_anchor(gene_anchor_keys, gene_id, sort_key)
                    _update_missing_gene_info(missing_gene_info, gene_id, fields)
            elif transcript_id:
                transcript_children[transcript_id].append(record)
                _update_gtf_nested_anchor(transcript_anchor_keys, transcript_id, sort_key)
                if gene_id:
                    gene_to_transcripts[gene_id].add(transcript_id)
                    _update_gtf_nested_anchor(gene_anchor_keys, gene_id, sort_key)
                    _update_missing_gene_info(missing_gene_info, gene_id, fields)
            elif gene_id:
                gene_records[gene_id].append(record)
                _update_gtf_nested_anchor(gene_anchor_keys, gene_id, sort_key)
                _update_missing_gene_info(missing_gene_info, gene_id, fields)
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

    synthetic_line_no = -1
    added_gene_n = 0
    for gene_id, gene_info in missing_gene_info.items():
        if gene_id in observed_gene_ids:
            continue
        synthetic_record = _make_synthetic_gene_record(gene_id, gene_info, synthetic_line_no)
        synthetic_line_no -= 1
        gene_records[gene_id].append(synthetic_record)
        _update_gtf_nested_anchor(gene_anchor_keys, gene_id, synthetic_record.sort_key)
        added_gene_n += 1

    _log_gtf_nested_sort_progress(
        'Finished reading/grouping',
        record_count,
        start_time,
        (
            f"comments={len(comments):,}, genes={len(gene_anchor_keys):,}, "
            f"transcripts={len(transcript_anchor_keys):,}, added_genes={added_gene_n:,}"
        ),
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
        added_gene_n,
    )


def write_talon_init_gtf(input_gtf, output_gtf, progress_interval=GTF_NESTED_SORT_PROGRESS_INTERVAL):
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
        added_gene_n,
    ) = _load_group_and_patch_gtf(input_gtf, progress_interval)

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
    logger.info(f"Writing TALON init nested GTF: {output_gtf}")
    with open(output_gtf, "w", encoding='utf-8') as fout:
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
        f"TALON init nested GTF complete for {input_gtf}: "
        f"input={record_count:,}; output={emitted_lines:,}; "
        f"added_genes={added_gene_n:,}; elapsed={time.perf_counter() - total_start:.1f}s"
    )
    return record_count, emitted_lines, added_gene_n


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit(f'Usage: python {sys.argv[0]} input.gtf output.gtf')
    write_talon_init_gtf(sys.argv[1], sys.argv[2])
