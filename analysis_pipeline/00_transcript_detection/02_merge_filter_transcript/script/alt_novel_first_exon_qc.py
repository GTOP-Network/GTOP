#!/usr/bin/env python3
from __future__ import annotations


import argparse
import csv
import gzip
import logging
import multiprocessing as mp
import os
import re
import sys
import time
from bisect import bisect_left
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


LOGGER = logging.getLogger("novel_first_exon_qc")

DEFAULT_BLOCKS_PER_TASK = 20
DEFAULT_PROGRESS_INTERVAL = 30


@dataclass(frozen=True)
class ExonRecord:
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    start_1based: int
    end_1based: int
    start_0: int
    end_0: int
    length: int
    region_id: str


@dataclass(frozen=True)
class RegionRecord:
    region_id: str
    chrom: str
    strand: str
    start_1based: int
    end_1based: int
    start_0: int
    end_0: int
    length: int


@dataclass(frozen=True)
class ScanBlock:
    block_id: str
    chrom: str
    start_0: int
    end_0: int
    regions: tuple[RegionRecord, ...]


def setup_logging(log_file: str | None, verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] [pid=%(process)d] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    LOGGER.setLevel(level)
    LOGGER.handlers.clear()

    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(level)
    stderr_handler.setFormatter(formatter)
    LOGGER.addHandler(stderr_handler)

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        LOGGER.addHandler(file_handler)


def open_text(path: str, mode: str = "rt"):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")


def parse_gtf_attributes(attribute_field: str) -> dict[str, str]:
    attributes = {}
    for item in attribute_field.rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if " " in item:
            key, value = item.split(" ", 1)
            attributes[key] = value.strip().strip('"')
        elif "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value.strip().strip('"')
    return attributes


def normalize_transcript_id(transcript_id: str, strip_version: bool) -> str:
    transcript_id = transcript_id.strip()
    if strip_version and "." in transcript_id:
        head, tail = transcript_id.rsplit(".", 1)
        if tail.isdigit():
            return head
    return transcript_id


def make_region_id(chrom: str, strand: str, start_1based: int, end_1based: int) -> str:
    safe_chrom = re.sub(r"[^A-Za-z0-9_.-]+", "_", chrom)
    strand_name = "plus" if strand == "+" else "minus" if strand == "-" else "unknown"
    return f"{safe_chrom}:{start_1based}-{end_1based}:{strand_name}"


def sample_name_from_bam(bam_path: str) -> str:
    name = Path(bam_path).name
    if name.endswith(".bam"):
        name = name[:-4]
    return name


def is_better_first_exon(candidate, current) -> bool:
    candidate_chrom, candidate_start, candidate_end, candidate_strand = candidate
    current_chrom, current_start, current_end, current_strand = current

    if candidate_chrom != current_chrom or candidate_strand != current_strand:
        return False

    if candidate_strand == "-":
        if candidate_end != current_end:
            return candidate_end > current_end
        return candidate_start > current_start

    if candidate_start != current_start:
        return candidate_start < current_start
    return candidate_end < current_end


def load_sqanti_novel_ids(sqanti_path: str, strip_version: bool) -> set[str]:
    LOGGER.info("Loading SQANTI3 annotation: %s", sqanti_path)

    novel_ids = set()
    total_rows = 0

    with open_text(sqanti_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty SQANTI3 file: {sqanti_path}")

        required = {"isoform", "associated_transcript"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"SQANTI3 file missing columns: {', '.join(sorted(missing))}")

        for row in reader:
            total_rows += 1
            associated_transcript = (row.get("associated_transcript") or "").strip().lower()
            if associated_transcript != "novel":
                continue

            isoform = (row.get("isoform") or "").strip()
            if isoform:
                novel_ids.add(normalize_transcript_id(isoform, strip_version))

    LOGGER.info("SQANTI3 rows=%d, novel transcript IDs=%d", total_rows, len(novel_ids))
    return novel_ids


def load_first_exons(
    gtf_path: str,
    transcript_filter: set[str] | None,
    strip_version: bool,
) -> list[ExonRecord]:
    LOGGER.info("Loading first exons from GTF: %s", gtf_path)

    first_exons: dict[str, tuple[str, int, int, str, str]] = {}
    total_exon_lines = 0
    kept_exon_lines = 0

    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue

            total_exon_lines += 1

            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue

            transcript_id = normalize_transcript_id(transcript_id, strip_version)
            if transcript_filter is not None and transcript_id not in transcript_filter:
                continue

            kept_exon_lines += 1

            gene_id = attrs.get("gene_id", "NA")
            chrom = fields[0]
            start_1based = int(fields[3])
            end_1based = int(fields[4])
            strand = fields[6]

            candidate = (chrom, start_1based, end_1based, strand)
            current = first_exons.get(transcript_id)

            if current is None or is_better_first_exon(candidate, current[:4]):
                first_exons[transcript_id] = (chrom, start_1based, end_1based, strand, gene_id)

    records = []
    for transcript_id, (chrom, start_1based, end_1based, strand, gene_id) in first_exons.items():
        region_id = make_region_id(chrom, strand, start_1based, end_1based)
        records.append(
            ExonRecord(
                transcript_id=transcript_id,
                gene_id=gene_id,
                chrom=chrom,
                strand=strand,
                start_1based=start_1based,
                end_1based=end_1based,
                start_0=start_1based - 1,
                end_0=end_1based,
                length=end_1based - start_1based + 1,
                region_id=region_id,
            )
        )

    records.sort(key=lambda x: (x.chrom, x.start_1based, x.end_1based, x.strand, x.transcript_id))

    LOGGER.info(
        "GTF=%s total_exon_lines=%d kept_exon_lines=%d transcripts_with_first_exon=%d",
        gtf_path,
        total_exon_lines,
        kept_exon_lines,
        len(records),
    )
    return records


def load_all_exon_coordinate_keys(gtf_path: str) -> set[tuple[str, str, int, int]]:
    LOGGER.info("Loading all reference exon coordinates for exact comparison: %s", gtf_path)

    keys = set()
    total_exon_lines = 0

    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue

            total_exon_lines += 1
            chrom = fields[0]
            start_0 = int(fields[3]) - 1
            end_0 = int(fields[4])
            strand = fields[6]
            keys.add((chrom, strand, start_0, end_0))

    LOGGER.info("Reference all-exon coordinate keys=%d from exon lines=%d", len(keys), total_exon_lines)
    return keys


def write_tsv(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    LOGGER.info("Writing table: %s rows=%d", path, len(rows))
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: str) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def iter_tsv(path: str):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield row


def prepare_regions(args) -> None:
    start_time = time.time()

    novel_ids = load_sqanti_novel_ids(args.sqanti, args.strip_version)
    lrs_first_exons = load_first_exons(args.lrs_gtf, transcript_filter=novel_ids, strip_version=args.strip_version)

    if args.compare_to == "ref-first-exon":
        ref_first_exons = load_first_exons(args.ref_gtf, transcript_filter=None, strip_version=args.strip_version)
        reference_coordinate_keys = {
            (record.chrom, record.strand, record.start_0, record.end_0)
            for record in ref_first_exons
        }
        LOGGER.info("Comparison mode=ref-first-exon, reference first-exon keys=%d", len(reference_coordinate_keys))
    elif args.compare_to == "ref-any-exon-exact":
        reference_coordinate_keys = load_all_exon_coordinate_keys(args.ref_gtf)
        LOGGER.info("Comparison mode=ref-any-exon-exact")
    else:
        raise ValueError(f"Unsupported --compare-to: {args.compare_to}")

    novel_first_exons = [
        record
        for record in lrs_first_exons
        if (record.chrom, record.strand, record.start_0, record.end_0) not in reference_coordinate_keys
    ]

    unique_regions: dict[str, RegionRecord] = {}
    region_to_transcripts = defaultdict(list)
    region_to_genes = defaultdict(set)

    for record in novel_first_exons:
        unique_regions[record.region_id] = RegionRecord(
            region_id=record.region_id,
            chrom=record.chrom,
            strand=record.strand,
            start_1based=record.start_1based,
            end_1based=record.end_1based,
            start_0=record.start_0,
            end_0=record.end_0,
            length=record.length,
        )
        region_to_transcripts[record.region_id].append(record.transcript_id)
        region_to_genes[record.region_id].add(record.gene_id)

    LOGGER.info("Novel first-exon transcript records=%d", len(novel_first_exons))
    LOGGER.info("Unique novel first-exon regions=%d", len(unique_regions))

    if novel_first_exons:
        dedup_ratio = 100.0 * (1.0 - len(unique_regions) / len(novel_first_exons))
    else:
        dedup_ratio = 0.0
    LOGGER.info("Deduplicated BAM-search targets by %.2f%%", dedup_ratio)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    transcript_rows = [
        {
            "transcript_id": record.transcript_id,
            "gene_id": record.gene_id,
            "region_id": record.region_id,
            "chrom": record.chrom,
            "strand": record.strand,
            "first_exon_start": record.start_1based,
            "first_exon_end": record.end_1based,
            "first_exon_length": record.length,
        }
        for record in novel_first_exons
    ]

    region_rows = []
    for region in sorted(unique_regions.values(), key=lambda x: (x.chrom, x.start_1based, x.end_1based, x.strand)):
        transcripts = sorted(region_to_transcripts[region.region_id])
        genes = sorted(region_to_genes[region.region_id])
        region_rows.append(
            {
                "region_id": region.region_id,
                "chrom": region.chrom,
                "strand": region.strand,
                "start": region.start_1based,
                "end": region.end_1based,
                "length": region.length,
                "transcript_count": len(transcripts),
                "transcript_ids": ",".join(transcripts),
                "gene_ids": ",".join(genes),
            }
        )

    write_tsv(
        f"{args.out_prefix}.transcripts.tsv",
        transcript_rows,
        [
            "transcript_id",
            "gene_id",
            "region_id",
            "chrom",
            "strand",
            "first_exon_start",
            "first_exon_end",
            "first_exon_length",
        ],
    )

    write_tsv(
        f"{args.out_prefix}.regions.tsv",
        region_rows,
        [
            "region_id",
            "chrom",
            "strand",
            "start",
            "end",
            "length",
            "transcript_count",
            "transcript_ids",
            "gene_ids",
        ],
    )

    LOGGER.info("Prepare finished in %.2f sec", time.time() - start_time)


def load_regions(region_tsv: str) -> list[RegionRecord]:
    LOGGER.info("Loading unique regions: %s", region_tsv)

    regions = []
    for row in read_tsv(region_tsv):
        start_1based = int(row["start"])
        end_1based = int(row["end"])
        regions.append(
            RegionRecord(
                region_id=row["region_id"],
                chrom=row["chrom"],
                strand=row["strand"],
                start_1based=start_1based,
                end_1based=end_1based,
                start_0=start_1based - 1,
                end_0=end_1based,
                length=end_1based - start_1based + 1,
            )
        )

    regions.sort(key=lambda x: (x.chrom, x.start_0, x.end_0, x.strand, x.region_id))
    LOGGER.info("Loaded unique regions=%d", len(regions))
    return regions


def read_bam_list(args) -> list[str]:
    bam_paths = []

    if args.bam_list:
        with open(args.bam_list, encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if line and not line.startswith("#"):
                    bam_paths.append(line)

    bam_paths.extend(args.bams or [])
    bam_paths = list(dict.fromkeys(bam_paths))
    return bam_paths


def get_bam_contigs_and_check_index(bam_path: str) -> set[str]:
    import pysam
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if not bam.has_index():
            raise ValueError(f"BAM has no index. Please create .bai/.csi first: {bam_path}")
        return set(bam.references)


def build_scan_blocks(regions: list[RegionRecord], merge_gap: int) -> list[ScanBlock]:
    grouped = defaultdict(list)
    for region in regions:
        grouped[region.chrom].append(region)

    blocks = []
    for chrom, chrom_regions in grouped.items():
        chrom_regions.sort(key=lambda x: (x.start_0, x.end_0, x.strand, x.region_id))

        current_start = None
        current_end = None
        current_regions = []

        for region in chrom_regions:
            if current_start is None:
                current_start = region.start_0
                current_end = region.end_0
                current_regions = [region]
                continue

            if region.start_0 <= current_end + merge_gap:
                current_end = max(current_end, region.end_0)
                current_regions.append(region)
            else:
                block_id = f"{chrom}:{current_start + 1}-{current_end}"
                blocks.append(
                    ScanBlock(
                        block_id=block_id,
                        chrom=chrom,
                        start_0=current_start,
                        end_0=current_end,
                        regions=tuple(current_regions),
                    )
                )
                current_start = region.start_0
                current_end = region.end_0
                current_regions = [region]

        if current_start is not None:
            block_id = f"{chrom}:{current_start + 1}-{current_end}"
            blocks.append(
                ScanBlock(
                    block_id=block_id,
                    chrom=chrom,
                    start_0=current_start,
                    end_0=current_end,
                    regions=tuple(current_regions),
                )
            )

    blocks.sort(key=lambda x: (x.chrom, x.start_0, x.end_0))
    return blocks


def chunk_blocks(blocks: list[ScanBlock], blocks_per_task: int) -> list[list[ScanBlock]]:
    chunks = []
    current_chunk = []

    sorted_blocks = sorted(
        blocks,
        key=lambda x: max(1, x.end_0 - x.start_0) + len(x.regions) * 1000,
        reverse=True,
    )

    for block in sorted_blocks:
        current_chunk.append(block)
        if len(current_chunk) >= blocks_per_task:
            chunks.append(current_chunk)
            current_chunk = []

    if current_chunk:
        chunks.append(current_chunk)

    return chunks


def md_has_snv(md: str) -> bool:
    in_deletion = False

    for char in md:
        if char == "^":
            in_deletion = True
            continue
        if char.isdigit():
            in_deletion = False
            continue
        if char.isalpha() and not in_deletion:
            return True

    return False


def extract_md_snv_positions_0based(read, left_0: int, right_0: int) -> list[int] | None:
    if not read.has_tag("MD"):
        return None

    try:
        md = read.get_tag("MD")
    except Exception:
        return None

    if not md_has_snv(md):
        return []

    try:
        reference_positions = read.get_reference_positions(full_length=False)
    except Exception:
        return None

    mismatch_positions = []
    aligned_index = 0
    index = 0

    while index < len(md):
        char = md[index]

        if char.isdigit():
            start = index
            while index < len(md) and md[index].isdigit():
                index += 1
            aligned_index += int(md[start:index])
            continue

        if char == "^":
            index += 1
            while index < len(md) and md[index].isalpha():
                index += 1
            continue

        if char.isalpha():
            if aligned_index < len(reference_positions):
                position_0based = reference_positions[aligned_index]
                if position_0based is not None and left_0 <= position_0based < right_0:
                    mismatch_positions.append(position_0based)
            aligned_index += 1
            index += 1
            continue

        index += 1

    return mismatch_positions


def worker_scan_blocks(task):
    import pysam
    bam_path, blocks, min_mapq, ignore_strand = task

    start_time = time.time()

    region_reads = defaultdict(int)
    region_md_reads = defaultdict(int)
    region_md_aligned_bases = defaultdict(int)
    region_mismatch_observations = defaultdict(int)
    region_mismatch_positions = defaultdict(set)

    total_fetched_reads = 0
    total_read_region_overlaps = 0
    failed_fetch_blocks = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for block in blocks:
            regions = list(block.regions)
            region_starts = [region.start_0 for region in regions]

            try:
                iterator = bam.fetch(block.chrom, block.start_0, block.end_0)
            except ValueError:
                failed_fetch_blocks += 1
                continue

            for read in iterator:
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                    continue
                if read.mapping_quality < min_mapq:
                    continue
                if read.reference_start is None or read.reference_end is None:
                    continue

                read_start = read.reference_start
                read_end = read.reference_end

                if read_start >= block.end_0 or read_end <= block.start_0:
                    continue

                total_fetched_reads += 1
                read_strand = "-" if read.is_reverse else "+"

                candidate_regions = []
                limit = bisect_left(region_starts, read_end)

                for index in range(limit):
                    region = regions[index]
                    if region.end_0 <= read_start:
                        continue
                    if not ignore_strand and region.strand != read_strand:
                        continue
                    candidate_regions.append(region)

                if not candidate_regions:
                    continue

                for region in candidate_regions:
                    region_reads[region.region_id] += 1
                    total_read_region_overlaps += 1

                mismatch_positions = extract_md_snv_positions_0based(
                    read,
                    block.start_0,
                    block.end_0,
                )
                if mismatch_positions is None:
                    continue

                for region in candidate_regions:
                    aligned_bases = read.get_overlap(region.start_0, region.end_0)

                    if aligned_bases <= 0:
                        continue

                    region_md_reads[region.region_id] += 1
                    region_md_aligned_bases[region.region_id] += aligned_bases

                    left_mismatch = bisect_left(mismatch_positions, region.start_0)
                    right_mismatch = bisect_left(mismatch_positions, region.end_0)
                    mismatch_count = right_mismatch - left_mismatch

                    if mismatch_count <= 0:
                        continue

                    region_mismatch_observations[region.region_id] += mismatch_count
                    for position_0based in mismatch_positions[left_mismatch:right_mismatch]:
                        region_mismatch_positions[region.region_id].add(position_0based + 1)

    return {
        "region_reads": dict(region_reads),
        "region_md_reads": dict(region_md_reads),
        "region_md_aligned_bases": dict(region_md_aligned_bases),
        "region_mismatch_observations": dict(region_mismatch_observations),
        "region_mismatch_positions": {key: set(value) for key, value in region_mismatch_positions.items()},
        "blocks": len(blocks),
        "failed_fetch_blocks": failed_fetch_blocks,
        "fetched_reads": total_fetched_reads,
        "read_region_overlaps": total_read_region_overlaps,
        "seconds": time.time() - start_time,
    }


def merge_worker_result(target, source) -> None:
    for key, value in source["region_reads"].items():
        target["region_reads"][key] += int(value)

    for key, value in source["region_md_reads"].items():
        target["region_md_reads"][key] += int(value)

    for key, value in source["region_md_aligned_bases"].items():
        target["region_md_aligned_bases"][key] += int(value)

    for key, value in source["region_mismatch_observations"].items():
        target["region_mismatch_observations"][key] += int(value)

    for key, values in source["region_mismatch_positions"].items():
        target["region_mismatch_positions"][key].update(values)

    target["blocks"] += int(source["blocks"])
    target["failed_fetch_blocks"] += int(source["failed_fetch_blocks"])
    target["fetched_reads"] += int(source["fetched_reads"])
    target["read_region_overlaps"] += int(source["read_region_overlaps"])
    target["worker_seconds"] += float(source["seconds"])


def log_sample_progress(
    sample: str,
    completed_blocks: int,
    total_blocks: int,
    completed_tasks: int,
    total_tasks: int,
    sample_start_time: float,
) -> None:
    now = time.time()
    elapsed = now - sample_start_time
    rate = completed_blocks / elapsed if elapsed > 0 else 0.0
    remaining_blocks = max(0, total_blocks - completed_blocks)
    eta = remaining_blocks / rate if rate > 0 else -1.0

    LOGGER.info(
        "Sample=%s progress blocks=%d/%d tasks=%d/%d elapsed=%.1f sec rate=%.2f blocks/sec eta=%.1f sec",
        sample,
        completed_blocks,
        total_blocks,
        completed_tasks,
        total_tasks,
        elapsed,
        rate,
        eta,
    )


def scan_one_bam(
    bam_path: str,
    regions: list[RegionRecord],
    sample_processes: int,
    merge_gap: int,
    min_mapq: int,
    ignore_strand: bool,
) -> tuple[list[dict], dict]:
    sample = sample_name_from_bam(bam_path)
    sample_start_time = time.time()

    bam_contigs = get_bam_contigs_and_check_index(bam_path)

    short_regions = [region for region in regions if region.length < 30]
    scannable_regions = [region for region in regions if region.length >= 30]
    present_regions = [region for region in scannable_regions if region.chrom in bam_contigs]
    missing_regions = [region for region in scannable_regions if region.chrom not in bam_contigs]

    blocks = build_scan_blocks(present_regions, merge_gap=merge_gap)
    workers = max(1, min(sample_processes, len(blocks))) if blocks else 1

    LOGGER.info(
        "Sample=%s regions=%d short_regions=%d present_regions=%d missing_contig_regions=%d blocks=%d workers=%d merge_gap=%d",
        sample,
        len(regions),
        len(short_regions),
        len(present_regions),
        len(missing_regions),
        len(blocks),
        workers,
        merge_gap,
    )

    aggregate = {
        "region_reads": defaultdict(int),
        "region_md_reads": defaultdict(int),
        "region_md_aligned_bases": defaultdict(int),
        "region_mismatch_observations": defaultdict(int),
        "region_mismatch_positions": defaultdict(set),
        "blocks": 0,
        "failed_fetch_blocks": 0,
        "fetched_reads": 0,
        "read_region_overlaps": 0,
        "worker_seconds": 0.0,
    }

    if blocks:
        block_chunks = chunk_blocks(blocks, DEFAULT_BLOCKS_PER_TASK)
        tasks = [(bam_path, chunk, min_mapq, ignore_strand) for chunk in block_chunks]

        LOGGER.info(
            "Sample=%s progress enabled blocks_per_task=%d progress_interval=%d sec tasks=%d",
            sample,
            DEFAULT_BLOCKS_PER_TASK,
            DEFAULT_PROGRESS_INTERVAL,
            len(tasks),
        )

        completed_tasks = 0
        completed_blocks = 0
        last_progress_time = time.time()

        if workers == 1:
            for task in tasks:
                result = worker_scan_blocks(task)
                merge_worker_result(aggregate, result)

                completed_tasks += 1
                completed_blocks += int(result["blocks"])

                now = time.time()
                if now - last_progress_time >= DEFAULT_PROGRESS_INTERVAL or completed_tasks == len(tasks):
                    log_sample_progress(
                        sample,
                        completed_blocks,
                        len(blocks),
                        completed_tasks,
                        len(tasks),
                        sample_start_time,
                    )
                    last_progress_time = now
        else:
            context = mp.get_context("spawn")
            with ProcessPoolExecutor(max_workers=workers, mp_context=context) as executor:
                futures = [executor.submit(worker_scan_blocks, task) for task in tasks]

                for future in as_completed(futures):
                    result = future.result()
                    merge_worker_result(aggregate, result)

                    completed_tasks += 1
                    completed_blocks += int(result["blocks"])

                    now = time.time()
                    if now - last_progress_time >= DEFAULT_PROGRESS_INTERVAL or completed_tasks == len(tasks):
                        log_sample_progress(
                            sample,
                            completed_blocks,
                            len(blocks),
                            completed_tasks,
                            len(tasks),
                            sample_start_time,
                        )
                        last_progress_time = now

    sample_seconds = time.time() - sample_start_time
    missing_region_ids = {region.region_id for region in missing_regions}
    short_region_ids = {region.region_id for region in short_regions}

    rows = []
    status_counter = Counter()

    for region in regions:
        reads = int(aggregate["region_reads"].get(region.region_id, 0))
        md_reads = int(aggregate["region_md_reads"].get(region.region_id, 0))
        md_aligned_bases = int(aggregate["region_md_aligned_bases"].get(region.region_id, 0))
        mismatch_observations = int(aggregate["region_mismatch_observations"].get(region.region_id, 0))
        mismatch_positions = sorted(aggregate["region_mismatch_positions"].get(region.region_id, set()))
        unique_mismatch_sites = len(mismatch_positions)

        if region.region_id in short_region_ids:
            status = "SHORT_EXON"
        elif region.region_id in missing_region_ids:
            status = "MISSING_CONTIG"
        elif reads == 0:
            status = "NO_COVERAGE"
        elif md_aligned_bases == 0:
            status = "NO_MD"
        else:
            status = "OK"

        if md_aligned_bases > 0:
            mismatch_rate = mismatch_observations / md_aligned_bases
            mismatch_per_30nt = mismatch_rate * 30.0
            mismatch_rate_text = f"{mismatch_rate:.8f}"
            mismatch_per_30nt_text = f"{mismatch_per_30nt:.6f}"
        else:
            mismatch_rate_text = "NA"
            mismatch_per_30nt_text = "NA"

        status_counter[status] += 1

        rows.append(
            {
                "sample": sample,
                "bam": bam_path,
                "region_id": region.region_id,
                "chrom": region.chrom,
                "strand": region.strand,
                "start": region.start_1based,
                "end": region.end_1based,
                "length": region.length,
                "status": status,
                "reads": reads,
                "md_reads": md_reads,
                "md_aligned_bases": md_aligned_bases,
                "unique_mismatch_sites": unique_mismatch_sites,
                "mismatch_observations": mismatch_observations,
                "mismatch_rate": mismatch_rate_text,
                "mismatch_per_30nt": mismatch_per_30nt_text,
                "mismatch_positions_1based": ",".join(map(str, mismatch_positions)),
                "sample_seconds": f"{sample_seconds:.3f}",
            }
        )

    LOGGER.info(
        "Sample=%s done seconds=%.2f statuses=%s blocks=%d failed_fetch_blocks=%d fetched_reads=%d read_region_overlaps=%d md_aligned_bases=%d unique_mismatch_sites=%d mismatch_observations=%d",
        sample,
        sample_seconds,
        dict(status_counter),
        aggregate["blocks"],
        aggregate["failed_fetch_blocks"],
        aggregate["fetched_reads"],
        aggregate["read_region_overlaps"],
        sum(aggregate["region_md_aligned_bases"].values()),
        sum(len(value) for value in aggregate["region_mismatch_positions"].values()),
        sum(aggregate["region_mismatch_observations"].values()),
    )

    stats = {
        "sample": sample,
        "seconds": sample_seconds,
        "status_counter": dict(status_counter),
        "blocks": aggregate["blocks"],
        "failed_fetch_blocks": aggregate["failed_fetch_blocks"],
        "fetched_reads": aggregate["fetched_reads"],
        "read_region_overlaps": aggregate["read_region_overlaps"],
    }

    return rows, stats


def scan_bams(args) -> None:
    start_time = time.time()

    regions = load_regions(args.regions)
    bam_paths = read_bam_list(args)

    if not bam_paths:
        raise ValueError("No BAM files provided. Use --bam-list or --bams.")

    if args.sample_processes < 1:
        raise ValueError("--sample-processes must be >= 1")

    if args.merge_gap < 0:
        raise ValueError("--merge-gap must be >= 0")

    LOGGER.info("BAM files=%d", len(bam_paths))
    LOGGER.info("Regions to scan=%d", len(regions))
    LOGGER.info(
        "Scan parameters sample_processes=%d merge_gap=%d min_mapq=%d ignore_strand=%s",
        args.sample_processes,
        args.merge_gap,
        args.min_mapq,
        args.ignore_strand,
    )
    LOGGER.info(
        "Internal progress parameters blocks_per_task=%d progress_interval=%d sec",
        DEFAULT_BLOCKS_PER_TASK,
        DEFAULT_PROGRESS_INTERVAL,
    )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    output_path = f"{args.out_prefix}.sample_region_mismatch.tsv"
    fieldnames = [
        "sample",
        "bam",
        "region_id",
        "chrom",
        "strand",
        "start",
        "end",
        "length",
        "status",
        "reads",
        "md_reads",
        "md_aligned_bases",
        "unique_mismatch_sites",
        "mismatch_observations",
        "mismatch_rate",
        "mismatch_per_30nt",
        "mismatch_positions_1based",
        "sample_seconds",
    ]

    LOGGER.info("Writing scan result incrementally: %s", output_path)

    total_status_counter = Counter()

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for index, bam_path in enumerate(bam_paths, start=1):
            LOGGER.info("Starting BAM %d/%d: %s", index, len(bam_paths), bam_path)

            rows, stats = scan_one_bam(
                bam_path=bam_path,
                regions=regions,
                sample_processes=args.sample_processes,
                merge_gap=args.merge_gap,
                min_mapq=args.min_mapq,
                ignore_strand=args.ignore_strand,
            )

            writer.writerows(rows)
            handle.flush()

            total_status_counter.update(stats["status_counter"])

            elapsed = time.time() - start_time
            LOGGER.info(
                "Finished BAM %d/%d sample=%s elapsed_total=%.2f sec",
                index,
                len(bam_paths),
                stats["sample"],
                elapsed,
            )

    LOGGER.info("All scan statuses=%s", dict(total_status_counter))
    LOGGER.info("Scan finished in %.2f sec", time.time() - start_time)


def parse_positions(position_string: str) -> set[int]:
    position_string = (position_string or "").strip()
    if not position_string or position_string == "NA":
        return set()
    return {int(item) for item in position_string.split(",") if item}


def new_merge_accumulator() -> dict:
    return {
        "region_reads": defaultdict(int),
        "region_md_reads": defaultdict(int),
        "region_md_aligned_bases": defaultdict(int),
        "region_observations": defaultdict(int),
        "region_positions": defaultdict(set),
        "region_unique_mismatch_sites": {},
        "region_status_counter": defaultdict(Counter),
        "input_rows": 0,
    }


def merge_scan_results_python(scan_results: list[str]) -> dict:
    aggregate = new_merge_accumulator()

    region_reads = aggregate["region_reads"]
    region_md_reads = aggregate["region_md_reads"]
    region_md_aligned_bases = aggregate["region_md_aligned_bases"]
    region_observations = aggregate["region_observations"]
    region_positions = aggregate["region_positions"]
    region_status_counter = aggregate["region_status_counter"]

    input_rows = 0
    for path in scan_results:
        LOGGER.info("Reading scan result: %s", path)
        if not os.path.exists(path) or os.path.getsize(path) < 1024:
            continue
        for row in iter_tsv(path):
            input_rows += 1
            region_id = row["region_id"]

            region_reads[region_id] += int(row.get("reads", 0) or 0)
            region_md_reads[region_id] += int(row.get("md_reads", 0) or 0)
            region_md_aligned_bases[region_id] += int(row.get("md_aligned_bases", 0) or 0)
            region_observations[region_id] += int(row.get("mismatch_observations", 0) or 0)
            region_positions[region_id].update(parse_positions(row.get("mismatch_positions_1based", "")))
            region_status_counter[region_id][row.get("status", "UNKNOWN")] += 1

    aggregate["input_rows"] = input_rows
    return aggregate


def get_scannable_scan_results(scan_results: list[str]) -> list[str]:
    paths = []
    for path in scan_results:
        LOGGER.info("Reading scan result: %s", path)
        if not os.path.exists(path) or os.path.getsize(path) < 1024:
            continue
        paths.append(path)
    return paths


def merge_scan_results_polars(scan_results: list[str]) -> dict:
    import polars as pl

    aggregate = new_merge_accumulator()

    paths = get_scannable_scan_results(scan_results)
    if not paths:
        return aggregate

    numeric_targets = {
        "reads": aggregate["region_reads"],
        "md_reads": aggregate["region_md_reads"],
        "md_aligned_bases": aggregate["region_md_aligned_bases"],
        "mismatch_observations": aggregate["region_observations"],
    }
    numeric_columns = list(numeric_targets)
    required_columns = [
        "region_id",
        "status",
        "mismatch_positions_1based",
        *numeric_columns,
    ]

    scan = pl.scan_csv(
        paths,
        separator="\t",
        infer_schema=False,
        missing_utf8_is_empty_string=True,
    )
    df = scan.select(required_columns).collect()

    aggregate["input_rows"] = df.height
    if df.is_empty():
        return aggregate

    numeric_df = (
        df.select(["region_id", *numeric_columns])
        .with_columns(
            [
                pl.when(pl.col(column).is_null() | (pl.col(column) == ""))
                .then(pl.lit("0"))
                .otherwise(pl.col(column))
                .cast(pl.Int64)
                .alias(column)
                for column in numeric_columns
            ]
        )
        .group_by("region_id")
        .agg([pl.col(column).sum().alias(column) for column in numeric_columns])
    )

    for row in numeric_df.iter_rows(named=True):
        region_id = row["region_id"]
        for column, target in numeric_targets.items():
            target[region_id] += int(row[column] or 0)

    status_df = (
        df.select(["region_id", "status"])
        .group_by(["region_id", "status"])
        .agg(pl.len().alias("count"))
    )

    for row in status_df.iter_rows(named=True):
        aggregate["region_status_counter"][row["region_id"]][row["status"]] += int(row["count"])

    position_df = (
        df.select(["region_id", "mismatch_positions_1based"])
        .with_columns(
            pl.col("mismatch_positions_1based")
            .fill_null("")
            .str.strip_chars()
            .alias("mismatch_positions_1based")
        )
        .filter(
            (pl.col("mismatch_positions_1based") != "")
            & (pl.col("mismatch_positions_1based") != "NA")
        )
        .with_columns(pl.col("mismatch_positions_1based").str.split(",").alias("position"))
        .explode("position")
        .filter(pl.col("position") != "")
        .with_columns(pl.col("position").cast(pl.Int64).alias("position"))
        .group_by("region_id")
        .agg(pl.col("position").n_unique().alias("unique_mismatch_sites"))
    )

    for row in position_df.iter_rows(named=True):
        aggregate["region_unique_mismatch_sites"][row["region_id"]] = int(row["unique_mismatch_sites"] or 0)

    return aggregate


def merge_scan_results(scan_results: list[str]) -> dict:
    try:
        aggregate = merge_scan_results_polars(scan_results)
        LOGGER.info("Merged scan results with polars")
        return aggregate
    except ImportError:
        LOGGER.info("polars is not available; merging scan results with the streaming Python path")
        return merge_scan_results_python(scan_results)
    except Exception as error:
        LOGGER.warning(
            "polars merge failed; falling back to the streaming Python path: %s",
            error,
        )
        return merge_scan_results_python(scan_results)


def merge_results(args) -> None:
    start_time = time.time()

    if args.qc_count_mode != "read-bases":
        LOGGER.warning(
            "--qc-count-mode=%s is accepted for backward compatibility, but mismatch_per_30nt is now always calculated from read bases: mismatch_observations * 30 / md_aligned_bases",
            args.qc_count_mode,
        )

    transcript_rows = read_tsv(args.transcripts)
    region_rows = read_tsv(args.regions)

    aggregate = merge_scan_results(args.scan_results)

    region_reads = aggregate["region_reads"]
    region_md_reads = aggregate["region_md_reads"]
    region_md_aligned_bases = aggregate["region_md_aligned_bases"]
    region_observations = aggregate["region_observations"]
    region_positions = aggregate["region_positions"]
    region_unique_mismatch_sites = aggregate["region_unique_mismatch_sites"]
    region_status_counter = aggregate["region_status_counter"]

    LOGGER.info("Merged scan rows=%d", aggregate["input_rows"])

    def summarize_region(region_id: str, length: int) -> tuple[str, str, int, int, int, str, str]:
        reads = region_reads[region_id]
        md_aligned_bases = region_md_aligned_bases[region_id]
        unique_sites = region_unique_mismatch_sites.get(region_id)
        if unique_sites is None:
            unique_sites = len(region_positions[region_id])
        observations = region_observations[region_id]

        if length < 30:
            return "SHORT_EXON", "FAIL", unique_sites, observations, md_aligned_bases, "NA", "NA"

        if md_aligned_bases > 0:
            status = "OK"
        elif region_status_counter[region_id]["MISSING_CONTIG"] > 0 and reads == 0:
            status = "MISSING_CONTIG"
        elif reads == 0:
            status = "NO_COVERAGE"
        else:
            status = "NO_MD"

        if status == "OK":
            mismatch_rate_value = observations / md_aligned_bases
            mismatch_per_30nt_value = mismatch_rate_value * 30.0
            mismatch_rate = f"{mismatch_rate_value:.8f}"
            mismatch_per_30nt = f"{mismatch_per_30nt_value:.6f}"
            qc_status = "PASS" if length >= 30 and mismatch_per_30nt_value <= 1.0 else "FAIL"
        else:
            mismatch_rate = "NA"
            mismatch_per_30nt = "NA"
            qc_status = "NA"

        return status, qc_status, unique_sites, observations, md_aligned_bases, mismatch_rate, mismatch_per_30nt

    region_qc_rows = []

    for row in region_rows:
        region_id = row["region_id"]
        length = int(row["length"])

        (
            status,
            qc_status,
            unique_sites,
            observations,
            md_aligned_bases,
            mismatch_rate,
            mismatch_per_30nt,
        ) = summarize_region(region_id, length)

        statuses = region_status_counter[region_id]

        region_qc_rows.append(
            {
                "region_id": region_id,
                "chrom": row["chrom"],
                "strand": row["strand"],
                "start": row["start"],
                "end": row["end"],
                "length": length,
                "transcript_count": row["transcript_count"],
                "total_reads": region_reads[region_id],
                "total_md_reads": region_md_reads[region_id],
                "total_md_aligned_bases": md_aligned_bases,
                "unique_mismatch_sites": unique_sites,
                "mismatch_observations": observations,
                "mismatch_rate": mismatch_rate,
                "mismatch_per_30nt": mismatch_per_30nt,
                "ok_samples": statuses["OK"],
                "no_coverage_samples": statuses["NO_COVERAGE"],
                "no_md_samples": statuses["NO_MD"],
                "missing_contig_samples": statuses["MISSING_CONTIG"],
                "status": status,
                "qc_status": qc_status,
            }
        )

    transcript_qc_rows = []

    for row in transcript_rows:
        region_id = row["region_id"]
        length = int(row["first_exon_length"])

        (
            status,
            qc_status,
            unique_sites,
            observations,
            md_aligned_bases,
            mismatch_rate,
            mismatch_per_30nt,
        ) = summarize_region(region_id, length)

        statuses = region_status_counter[region_id]

        transcript_qc_rows.append(
            {
                "transcript_id": row["transcript_id"],
                "gene_id": row.get("gene_id", "NA"),
                "region_id": region_id,
                "chrom": row["chrom"],
                "strand": row["strand"],
                "first_exon_start": row["first_exon_start"],
                "first_exon_end": row["first_exon_end"],
                "first_exon_length": length,
                "total_reads": region_reads[region_id],
                "total_md_reads": region_md_reads[region_id],
                "total_md_aligned_bases": md_aligned_bases,
                "unique_mismatch_sites": unique_sites,
                "mismatch_observations": observations,
                "mismatch_rate": mismatch_rate,
                "mismatch_per_30nt": mismatch_per_30nt,
                "ok_samples": statuses["OK"],
                "no_coverage_samples": statuses["NO_COVERAGE"],
                "no_md_samples": statuses["NO_MD"],
                "missing_contig_samples": statuses["MISSING_CONTIG"],
                "status": status,
                "qc_status": qc_status,
            }
        )

    write_tsv(
        f"{args.out_prefix}.region_qc.tsv",
        region_qc_rows,
        [
            "region_id",
            "chrom",
            "strand",
            "start",
            "end",
            "length",
            "transcript_count",
            "total_reads",
            "total_md_reads",
            "total_md_aligned_bases",
            "unique_mismatch_sites",
            "mismatch_observations",
            "mismatch_rate",
            "mismatch_per_30nt",
            "ok_samples",
            "no_coverage_samples",
            "no_md_samples",
            "missing_contig_samples",
            "status",
            "qc_status",
        ],
    )

    write_tsv(
        f"{args.out_prefix}.transcript_qc.tsv",
        transcript_qc_rows,
        [
            "transcript_id",
            "gene_id",
            "region_id",
            "chrom",
            "strand",
            "first_exon_start",
            "first_exon_end",
            "first_exon_length",
            "total_reads",
            "total_md_reads",
            "total_md_aligned_bases",
            "unique_mismatch_sites",
            "mismatch_observations",
            "mismatch_rate",
            "mismatch_per_30nt",
            "ok_samples",
            "no_coverage_samples",
            "no_md_samples",
            "missing_contig_samples",
            "status",
            "qc_status",
        ],
    )

    LOGGER.info("Merge finished in %.2f sec", time.time() - start_time)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fast novel first-exon mismatch QC pipeline for LRS transcripts."
    )
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--verbose", action="store_true")

    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--ref-gtf", required=True)
    prepare.add_argument("--lrs-gtf", required=True)
    prepare.add_argument("--sqanti", required=True)
    prepare.add_argument("--out-prefix", required=True)
    prepare.add_argument("--strip-version", action="store_true")
    prepare.add_argument(
        "--compare-to",
        choices=["ref-first-exon", "ref-any-exon-exact"],
        default="ref-first-exon",
        help="For AF events, ref-first-exon is usually the recommended comparison mode.",
    )
    prepare.set_defaults(func=prepare_regions)

    scan = subparsers.add_parser("scan")
    scan.add_argument("--regions", required=True)
    scan.add_argument("--out-prefix", required=True)
    scan.add_argument("--bam-list", default=None)
    scan.add_argument("--bams", nargs="*", default=[])
    scan.add_argument(
        "--sample-processes",
        type=int,
        default=max(1, min(8, os.cpu_count() or 1)),
        help="Worker processes per BAM. Too many workers can slow indexed BAM reads on shared storage.",
    )
    scan.add_argument(
        "--merge-gap",
        type=int,
        default=0,
        help="Merge nearby target regions before BAM fetch. Use 0 for focused first-exon QC.",
    )
    scan.add_argument("--min-mapq", type=int, default=0)
    scan.add_argument("--ignore-strand", action="store_true")
    scan.set_defaults(func=scan_bams)

    merge = subparsers.add_parser("merge")
    merge.add_argument("--transcripts", required=True)
    merge.add_argument("--regions", required=True)
    merge.add_argument("--scan-results", nargs="+", required=True)
    merge.add_argument("--out-prefix", required=True)
    merge.add_argument(
        "--qc-count-mode",
        choices=["read-bases", "unique-sites", "observations"],
        default="read-bases",
        help="Kept for compatibility. Final mismatch_per_30nt is calculated from read bases.",
    )
    merge.set_defaults(func=merge_results)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    setup_logging(args.log_file, args.verbose)

    LOGGER.info("Command: %s", args.command)
    LOGGER.info("Arguments: %s", vars(args))

    args.func(args)


if __name__ == "__main__":
    main()
