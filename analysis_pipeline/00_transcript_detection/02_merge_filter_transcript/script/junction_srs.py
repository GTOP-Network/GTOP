# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/5/21 10:46
@Email   : xuechao@szbl.ac.cn
@Desc    : 从SRS样本中提取转录本的junction支持情况。
"""

#!/usr/bin/env python3

import argparse
import gzip
import multiprocessing as mp
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed


_WORKER_JUNCTION_MAP = None
_WORKER_OUT_DIR = None
_WORKER_COUNT_MODE = None


def open_text(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def infer_sample_id(path):
    base = os.path.basename(path)
    if base.endswith(".gz"):
        base = base[:-3]
    for suffix in (".SJ.out.tab", ".sj.out.tab", ".tsv", ".tab", ".txt"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def parse_gtf_attr(attr_text, key):
    for item in attr_text.rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(None, 1)
        if not parts or parts[0] != key:
            continue
        if len(parts) == 1:
            return ""
        value = parts[1].strip()
        if value.startswith('"'):
            end = value.find('"', 1)
            if end != -1:
                return value[1:end]
        return value.split()[0]
    return None


def junction_id(chrom, intron_start, intron_end, strand):
    return f"{chrom}:{intron_start}-{intron_end}:{strand}"


def build_index(args):
    os.makedirs(args.out_dir, exist_ok=True)

    exons_by_tx = defaultdict(list)

    with open_text(args.gtf) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue

            tid = parse_gtf_attr(fields[8], args.transcript_attr)
            if not tid:
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            if start > end:
                start, end = end, start

            exons_by_tx[tid].append((chrom, start, end, strand))

    unique_junctions = {}
    junction_to_transcripts = defaultdict(list)
    transcript_rows = []
    seen_tx_junction = set()

    for tid, exons in exons_by_tx.items():
        exons = sorted(set(exons), key=lambda x: (x[1], x[2], x[0], x[3]))
        if len(exons) < 2:
            continue

        chroms = {x[0] for x in exons}
        strands = {x[3] for x in exons}
        if len(chroms) != 1 or len(strands) != 1:
            print(f"[WARN] skip transcript with inconsistent chrom/strand: {tid}", file=sys.stderr)
            continue

        strand = exons[0][3]

        if strand == "-":
            tx_order_exons = sorted(exons, key=lambda x: (x[1], x[2]), reverse=True)
        else:
            tx_order_exons = sorted(exons, key=lambda x: (x[1], x[2]))

        for idx in range(len(tx_order_exons) - 1):
            exon_a = tx_order_exons[idx]
            exon_b = tx_order_exons[idx + 1]

            left_exon, right_exon = sorted((exon_a, exon_b), key=lambda x: (x[1], x[2]))
            chrom = left_exon[0]
            intron_start = left_exon[2] + 1
            intron_end = right_exon[1] - 1

            if intron_start > intron_end:
                continue

            jid = junction_id(chrom, intron_start, intron_end, strand)
            tx_junction_id = f"{tid}.J{idx + 1:04d}"
            key = (tid, jid, idx + 1)

            if key in seen_tx_junction:
                continue
            seen_tx_junction.add(key)

            unique_junctions[jid] = (chrom, intron_start, intron_end, strand)
            junction_to_transcripts[jid].append((tid, tx_junction_id, idx + 1))
            transcript_rows.append((tid, tx_junction_id, jid, chrom, intron_start, intron_end, strand, idx + 1))

    unique_path = os.path.join(args.out_dir, "unique_junctions.tsv")
    map_path = os.path.join(args.out_dir, "junction_to_transcript.tsv")
    tx_path = os.path.join(args.out_dir, "transcript_junctions.tsv")

    with open(unique_path, "w") as out:
        out.write("junction_id\tchrom\tintron_start\tintron_end\tstrand\ttranscript_count\n")
        for jid, meta in sorted(unique_junctions.items(), key=lambda x: (x[1][0], x[1][1], x[1][2], x[1][3])):
            chrom, intron_start, intron_end, strand = meta
            out.write(
                f"{jid}\t{chrom}\t{intron_start}\t{intron_end}\t{strand}\t"
                f"{len(junction_to_transcripts[jid])}\n"
            )

    with open(map_path, "w") as out:
        out.write("junction_id\ttranscript_id\ttranscript_junction_id\tjunction_order\n")
        for jid in sorted(junction_to_transcripts):
            for tid, tx_junction_id, order in sorted(junction_to_transcripts[jid], key=lambda x: (x[0], x[2])):
                out.write(f"{jid}\t{tid}\t{tx_junction_id}\t{order}\n")

    with open(tx_path, "w") as out:
        out.write("transcript_id\ttranscript_junction_id\tjunction_id\tchrom\tintron_start\tintron_end\tstrand\tjunction_order\n")
        for row in sorted(transcript_rows, key=lambda x: (x[0], x[7], x[2])):
            out.write("\t".join(map(str, row)) + "\n")

    print(f"[DONE] transcripts: {len(exons_by_tx)}", file=sys.stderr)
    print(f"[DONE] unique junctions: {len(unique_junctions)}", file=sys.stderr)
    print(f"[DONE] write: {unique_path}", file=sys.stderr)
    print(f"[DONE] write: {map_path}", file=sys.stderr)
    print(f"[DONE] write: {tx_path}", file=sys.stderr)


def load_junction_map(unique_junctions_path):
    junction_map = {}

    with open_text(unique_junctions_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: idx for idx, name in enumerate(header)}

        required = ["junction_id", "chrom", "intron_start", "intron_end", "strand"]
        missing = [x for x in required if x not in col]
        if missing:
            raise ValueError(f"missing columns in {unique_junctions_path}: {missing}")

        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            jid = fields[col["junction_id"]]
            chrom = fields[col["chrom"]]
            intron_start = int(fields[col["intron_start"]])
            intron_end = int(fields[col["intron_end"]])
            strand = fields[col["strand"]]
            junction_map[(chrom, intron_start, intron_end, strand)] = jid

    return junction_map


def parse_sj_list(path):
    samples = []

    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) == 1:
                sj_path = fields[0]
                sample_id = infer_sample_id(sj_path)
            else:
                sample_id = fields[0]
                sj_path = fields[1]

            samples.append((sample_id, sj_path))

    return samples


def count_one_sj(sample_id, sj_path, junction_map, out_dir, count_mode):
    counts = defaultdict(int)
    matched_sj_rows = 0
    scanned_sj_rows = 0

    with open_text(sj_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                fields = line.split()
            if len(fields) < 7:
                continue

            scanned_sj_rows += 1

            chrom = fields[0]
            intron_start = int(fields[1])
            intron_end = int(fields[2])
            strand_code = fields[3]

            if strand_code == "1":
                strand = "+"
            elif strand_code == "2":
                strand = "-"
            else:
                strand = "."

            unique_reads = int(fields[6])

            if count_mode == "total":
                multi_reads = int(fields[7]) if len(fields) > 7 else 0
                read_count = unique_reads + multi_reads
            else:
                read_count = unique_reads

            if read_count <= 0:
                continue

            jid = junction_map.get((chrom, intron_start, intron_end, strand))
            if jid is None:
                continue

            counts[jid] += read_count
            matched_sj_rows += 1

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{safe_name(sample_id)}.junction_counts.tsv")

    with open(out_path, "w") as out:
        out.write("sample_id\tjunction_id\tunique_reads\n")
        for jid in sorted(counts):
            out.write(f"{sample_id}\t{jid}\t{counts[jid]}\n")

    return {
        "sample_id": sample_id,
        "sj_path": sj_path,
        "out_path": out_path,
        "scanned_sj_rows": scanned_sj_rows,
        "matched_sj_rows": matched_sj_rows,
        "matched_junctions": len(counts),
        "total_reads": sum(counts.values()),
    }


def init_count_worker(junction_map, out_dir, count_mode):
    global _WORKER_JUNCTION_MAP
    global _WORKER_OUT_DIR
    global _WORKER_COUNT_MODE

    _WORKER_JUNCTION_MAP = junction_map
    _WORKER_OUT_DIR = out_dir
    _WORKER_COUNT_MODE = count_mode


def count_worker(sample):
    sample_id, sj_path = sample
    return count_one_sj(sample_id, sj_path, _WORKER_JUNCTION_MAP, _WORKER_OUT_DIR, _WORKER_COUNT_MODE)


def count_command(args):
    junction_map = load_junction_map(args.junctions)

    if args.sj_file:
        sample_id = args.sample_id or infer_sample_id(args.sj_file)
        samples = [(sample_id, args.sj_file)]
    else:
        samples = parse_sj_list(args.sj_list)

    if not samples:
        raise ValueError("no SJ.out.tab input found")

    os.makedirs(args.out_dir, exist_ok=True)

    print(f"[INFO] target junctions: {len(junction_map)}", file=sys.stderr)
    print(f"[INFO] samples: {len(samples)}", file=sys.stderr)
    print(f"[INFO] count mode: {args.count_mode}", file=sys.stderr)

    if args.sample_processes <= 1 or len(samples) == 1:
        for sample_id, sj_path in samples:
            result = count_one_sj(sample_id, sj_path, junction_map, args.out_dir, args.count_mode)
            print(
                f"[DONE] {result['sample_id']} matched_junctions={result['matched_junctions']} "
                f"total_reads={result['total_reads']} out={result['out_path']}",
                file=sys.stderr,
            )
    else:
        processes = min(args.sample_processes, len(samples))
        try:
            ctx = mp.get_context("fork")
        except ValueError:
            ctx = None

        with ProcessPoolExecutor(
            max_workers=processes,
            mp_context=ctx,
            initializer=init_count_worker,
            initargs=(junction_map, args.out_dir, args.count_mode),
        ) as pool:
            futures = [pool.submit(count_worker, sample) for sample in samples]
            for future in as_completed(futures):
                result = future.result()
                print(
                    f"[DONE] {result['sample_id']} matched_junctions={result['matched_junctions']} "
                    f"total_reads={result['total_reads']} out={result['out_path']}",
                    file=sys.stderr,
                )


def find_count_files(counts_dir):
    count_files = []
    for root, _, files in os.walk(counts_dir):
        for name in files:
            if name.endswith(".junction_counts.tsv") or name.endswith(".junction_counts.tsv.gz"):
                count_files.append(os.path.join(root, name))
    return sorted(count_files)


def parse_count_list(path):
    count_files = []

    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) == 1:
                count_files.append(fields[0])
            else:
                count_files.append(fields[1])

    return count_files


def load_transcript_junctions(path):
    rows = []

    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: idx for idx, name in enumerate(header)}

        required = [
            "transcript_id",
            "transcript_junction_id",
            "junction_id",
            "chrom",
            "intron_start",
            "intron_end",
            "strand",
            "junction_order",
        ]
        missing = [x for x in required if x not in col]
        if missing:
            raise ValueError(f"missing columns in {path}: {missing}")

        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            rows.append(
                (
                    fields[col["transcript_id"]],
                    fields[col["transcript_junction_id"]],
                    fields[col["junction_id"]],
                    fields[col["chrom"]],
                    fields[col["intron_start"]],
                    fields[col["intron_end"]],
                    fields[col["strand"]],
                    fields[col["junction_order"]],
                )
            )

    return rows


def merge_command(args):
    if args.count_list:
        count_files = parse_count_list(args.count_list)
    else:
        count_files = find_count_files(args.counts_dir)

    if not count_files:
        raise ValueError("no count files found")

    transcript_junctions = load_transcript_junctions(args.transcript_junctions)
    stats = defaultdict(lambda: [0, 0, 0, 0])
    scanned_files = 0

    for count_path in count_files:
        if not os.path.exists(count_path):
            if args.skip_missing:
                print(f"[WARN] missing count file, skip: {count_path}", file=sys.stderr)
                continue
            raise FileNotFoundError(count_path)

        scanned_files += 1
        print(count_path)
        with open_text(count_path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col = {name: idx for idx, name in enumerate(header)}

            if "junction_id" not in col:
                raise ValueError(f"missing junction_id column: {count_path}")

            if "unique_reads" in col:
                read_col = col["unique_reads"]
            elif "read_count" in col:
                read_col = col["read_count"]
            else:
                raise ValueError(f"missing unique_reads/read_count column: {count_path}")

            for line in fh:
                if not line.strip():
                    continue

                fields = line.rstrip("\n").split("\t")
                jid = fields[col["junction_id"]]
                reads = int(fields[read_col])

                item = stats[jid]
                item[0] += reads
                if reads >= 5:
                    item[1] += 1
                if reads >= 10:
                    item[2] += 1
                if reads >= 20:
                    item[3] += 1

    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    junction_summary_path = f"{args.out_prefix}.junction_summary.tsv"
    transcript_summary_path = f"{args.out_prefix}.transcript_junction_summary.tsv"

    with open(junction_summary_path, "w") as out:
        out.write("junction_id\ttotal_unique_reads\tsamples_ge_5\tsamples_ge_10\tsamples_ge_20\n")
        for jid in sorted(stats):
            total_reads, ge5, ge10, ge20 = stats[jid]
            out.write(f"{jid}\t{total_reads}\t{ge5}\t{ge10}\t{ge20}\n")

    with open(transcript_summary_path, "w") as out:
        out.write(
            "transcript_id\ttranscript_junction_id\tjunction_id\tchrom\tintron_start\tintron_end\tstrand\t"
            "junction_order\ttotal_unique_reads\tsamples_ge_5\tsamples_ge_10\tsamples_ge_20\n"
        )

        for row in transcript_junctions:
            tid, tx_jid, jid, chrom, intron_start, intron_end, strand, order = row
            total_reads, ge5, ge10, ge20 = stats.get(jid, [0, 0, 0, 0])
            out.write(
                f"{tid}\t{tx_jid}\t{jid}\t{chrom}\t{intron_start}\t{intron_end}\t{strand}\t"
                f"{order}\t{total_reads}\t{ge5}\t{ge10}\t{ge20}\n"
            )

    print(f"[DONE] scanned count files: {scanned_files}", file=sys.stderr)
    print(f"[DONE] write: {junction_summary_path}", file=sys.stderr)
    print(f"[DONE] write: {transcript_summary_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Build transcript junction index, count STAR SJ.out.tab unique reads, and merge SRS metrics."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    p = subparsers.add_parser("build-index", help="build unique junction index from transcript GTF")
    p.add_argument("--gtf", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--transcript-attr", default="transcript_id")
    p.set_defaults(func=build_index)

    p = subparsers.add_parser("count", help="count STAR SJ.out.tab support reads for indexed junctions")
    p.add_argument("--junctions", required=True, help="unique_junctions.tsv from build-index")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--sj-list", help="one line: SJ_path OR sample_id<TAB>SJ_path")
    group.add_argument("--sj-file", help="single STAR SJ.out.tab")
    p.add_argument("--sample-id")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sample-processes", type=int, default=1)
    p.add_argument("--count-mode", choices=["unique", "total"], default="unique")
    p.set_defaults(func=count_command)

    p = subparsers.add_parser("merge", help="merge sample count files and expand to transcript junctions")
    p.add_argument("--transcript-junctions", required=True, help="transcript_junctions.tsv from build-index")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--count-list", help="one line: count_path OR sample_id<TAB>count_path")
    group.add_argument("--counts-dir", help="directory containing *.junction_counts.tsv")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--skip-missing", action="store_true")
    p.set_defaults(func=merge_command)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
