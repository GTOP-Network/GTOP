# -*- coding: utf-8 -*-


from __future__ import annotations


import argparse
import csv
import gzip
import logging
import os
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, Iterable, List, Tuple


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

GTF_SOURCE = "GTOP"
DEFAULT_SUPPORT_ATTR = "sources"
DEFAULT_MAX_SPAN = 1_000_000


def format_seconds(seconds: float) -> str:
    return f"{seconds:.2f}s"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge many sample GTF/GFF files by intron chain with multiprocessing."
    )
    parser.add_argument("input_gtfs", help="CSV/TSV config with source name and GTF path.")
    parser.add_argument("output_gtf", help="Merged output GTF path.")
    parser.add_argument(
        "--source-column",
        default=None,
        help="Column name for sample/source label. Auto-detected if omitted.",
    )
    parser.add_argument(
        "--path-column",
        default=None,
        help="Column name for GTF path. Auto-detected if omitted.",
    )
    parser.add_argument(
        "--support-attr",
        default=DEFAULT_SUPPORT_ATTR,
        help="Output GTF attribute name storing supporting source labels.",
    )
    parser.add_argument(
        "--gtf-source",
        default=GTF_SOURCE,
        help="Value written to column 2 (source) in output GTF.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 1) - 1),
        help="Number of worker processes.",
    )
    parser.add_argument(
        "--max-span",
        type=int,
        default=DEFAULT_MAX_SPAN,
        help="Skip merged transcripts longer than this span.",
    )
    return parser.parse_args()


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_attributes(attr_str: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in attr_str.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            value = value.strip().strip('"')
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attrs[key.strip()] = value.strip()
    return attrs


def choose_column(
    fieldnames: Iterable[str],
    explicit_name: str | None,
    candidates: List[str],
    role: str,
) -> str:
    if explicit_name:
        if explicit_name not in fieldnames:
            raise ValueError(f"{role} column '{explicit_name}' not found in config.")
        return explicit_name

    for name in candidates:
        if name in fieldnames:
            return name
    raise ValueError(
        f"Unable to detect {role} column from config. Available columns: {list(fieldnames)}"
    )


def load_input_config(
    path: str,
    source_column: str | None,
    path_column: str | None,
) -> List[Tuple[str, str]]:
    with open(path, "r", newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        except csv.Error:
            dialect = csv.get_dialect("excel")

        reader = csv.DictReader(fh, dialect=dialect)
        if reader.fieldnames is None:
            raise ValueError(f"Invalid config file: {path}")

        src_col = choose_column(
            reader.fieldnames,
            source_column,
            ["source_name", "sample_name", "tool_name", "name", "source"],
            "source",
        )
        gtf_col = choose_column(
            reader.fieldnames,
            path_column,
            ["gtf_path", "path", "gtf", "file"],
            "path",
        )

        items = []
        for row in reader:
            source = (row.get(src_col) or "").strip()
            gtf_path = (row.get(gtf_col) or "").strip()
            if not source or not gtf_path:
                continue
            items.append((source, gtf_path))

    if not items:
        raise ValueError(f"No valid rows found in config: {path}")
    return items


def unique_preserve_order(values: Iterable[str]) -> List[str]:
    seen = set()
    result = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        result.append(value)
    return result


def infer_transcript_ids(feature: str, attrs: Dict[str, str]) -> List[str]:
    if feature == "exon":
        if "transcript_id" in attrs:
            return [attrs["transcript_id"]]
        if "Parent" in attrs:
            return [item.strip() for item in attrs["Parent"].split(",") if item.strip()]
        if "parent" in attrs:
            return [item.strip() for item in attrs["parent"].split(",") if item.strip()]
    return []


def exon_chain_to_ic(exons: List[Tuple[int, int]]) -> str:
    intron_coords: List[str] = []
    for i in range(len(exons) - 1):
        intron_coords.append(str(exons[i][1]))
        intron_coords.append(str(exons[i + 1][0]))
    return "-".join(intron_coords)


def is_valid_exon_chain(exons: List[Tuple[int, int]]) -> bool:
    for i in range(len(exons) - 1):
        if exons[i][1] >= exons[i + 1][0]:
            return False
    return True


def process_gtf_file(task: Tuple[str, str]) -> Tuple[str, Dict[Tuple[str, str, str], Tuple[int, int]]]:
    source_name, gtf_path = task
    start_time = time.perf_counter()

    transcript_exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    transcript_meta: Dict[str, Tuple[str, str]] = {}
    transcript_meta_conflict: set[str] = set()

    exon_count = 0
    duplicate_exon_transcripts = 0
    invalid_chain_transcripts = 0
    singleton_transcripts = 0
    accepted_transcripts = 0

    with open_text(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue

            chrom, _, feature, start, end, _, strand, _, attr_str = cols
            if feature != "exon" or strand not in {"+", "-"}:
                continue

            exon_count += 1

            attrs = parse_attributes(attr_str)
            transcript_ids = infer_transcript_ids(feature, attrs)
            if not transcript_ids:
                continue

            start_i = int(start)
            end_i = int(end)
            if start_i > end_i:
                start_i, end_i = end_i, start_i

            for transcript_id in transcript_ids:
                transcript_exons[transcript_id].append((start_i, end_i))

                curr_meta = (chrom, strand)
                prev_meta = transcript_meta.get(transcript_id)
                if prev_meta is None:
                    transcript_meta[transcript_id] = curr_meta
                elif prev_meta != curr_meta:
                    transcript_meta_conflict.add(transcript_id)

    merged: Dict[Tuple[str, str, str], Tuple[int, int]] = {}

    for transcript_id, exons in transcript_exons.items():
        if transcript_id in transcript_meta_conflict:
            invalid_chain_transcripts += 1
            continue

        if len(exons) < 2:
            singleton_transcripts += 1
            continue

        chrom, strand = transcript_meta[transcript_id]
        sorted_exons = sorted(exons)

        if len(set(sorted_exons)) != len(sorted_exons):
            duplicate_exon_transcripts += 1
            continue

        if not is_valid_exon_chain(sorted_exons):
            invalid_chain_transcripts += 1
            continue

        ic = exon_chain_to_ic(sorted_exons)
        key = (chrom, strand, ic)

        tx_start = sorted_exons[0][0]
        tx_end = sorted_exons[-1][1]

        accepted_transcripts += 1

        prev = merged.get(key)
        if prev is None:
            merged[key] = (tx_start, tx_end)
        else:
            merged[key] = (min(prev[0], tx_start), max(prev[1], tx_end))

    elapsed = time.perf_counter() - start_time

    logger.info(
        "Processed %s in %s: %d exon records, %d transcripts, %d merged intron chains",
        source_name,
        format_seconds(elapsed),
        exon_count,
        len(transcript_exons),
        len(merged),
    )
    logger.info(
        "Filtered %s: accepted=%d, singleton=%d, duplicate_exon=%d, invalid_chain=%d, meta_conflict=%d",
        source_name,
        accepted_transcripts,
        singleton_transcripts,
        duplicate_exon_transcripts,
        invalid_chain_transcripts,
        len(transcript_meta_conflict),
    )

    return source_name, merged


def merge_results(
    tasks: List[Tuple[str, str]],
    workers: int,
) -> Dict[Tuple[str, str, str], Dict[str, object]]:
    start_time = time.perf_counter()
    combined: Dict[Tuple[str, str, str], Dict[str, object]] = {}

    if workers <= 1:
        iterator = map(process_gtf_file, tasks)

        for source_name, merged in iterator:
            for key, (tx_start, tx_end) in merged.items():
                item = combined.get(key)
                if item is None:
                    combined[key] = {
                        "start": tx_start,
                        "end": tx_end,
                        "sources": {source_name},
                    }
                else:
                    item["start"] = min(item["start"], tx_start)
                    item["end"] = max(item["end"], tx_end)
                    item["sources"].add(source_name)

        logger.info(
            "Completed serial merge stage in %s",
            format_seconds(time.perf_counter() - start_time),
        )

    else:
        chunksize = max(1, len(tasks) // max(workers * 4, 1))

        logger.info(
            "Starting parallel merge stage: %d tasks, %d workers, chunksize=%d",
            len(tasks),
            workers,
            chunksize,
        )

        with ProcessPoolExecutor(max_workers=workers) as executor:
            for source_name, merged in executor.map(
                process_gtf_file,
                tasks,
                chunksize=chunksize,
            ):
                for key, (tx_start, tx_end) in merged.items():
                    item = combined.get(key)
                    if item is None:
                        combined[key] = {
                            "start": tx_start,
                            "end": tx_end,
                            "sources": {source_name},
                        }
                    else:
                        item["start"] = min(item["start"], tx_start)
                        item["end"] = max(item["end"], tx_end)
                        item["sources"].add(source_name)

        logger.info(
            "Completed parallel merge stage in %s",
            format_seconds(time.perf_counter() - start_time),
        )

    return combined


def intron_chain_to_exons(tx_start: int, tx_end: int, ic: str) -> List[Tuple[int, int]]:
    coords = [tx_start] + [int(x) for x in ic.split("-")] + [tx_end]
    return [(coords[i], coords[i + 1]) for i in range(0, len(coords) - 1, 2)]


def format_attributes(attrs: Dict[str, str]) -> str:
    return "; ".join(f'{key} "{value}"' for key, value in attrs.items()) + ";"


def get_meta_tsv_path(output_gtf: str) -> str:
    if output_gtf.endswith(".gtf.gz"):
        return output_gtf[:-7] + ".meta.tsv"
    if output_gtf.endswith(".gff.gz"):
        return output_gtf[:-7] + ".meta.tsv"
    if output_gtf.endswith(".gtf"):
        return output_gtf[:-4] + ".meta.tsv"
    if output_gtf.endswith(".gff"):
        return output_gtf[:-4] + ".meta.tsv"
    return output_gtf + ".meta.tsv"


def write_gtf(
    merged: Dict[Tuple[str, str, str], Dict[str, object]],
    output_gtf: str,
    gtf_source: str,
    support_attr: str,
    max_span: int,
    source_names: List[str],
) -> None:
    output_dir = os.path.dirname(output_gtf)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    meta_tsv = get_meta_tsv_path(output_gtf)
    meta_dir = os.path.dirname(meta_tsv)
    if meta_dir:
        os.makedirs(meta_dir, exist_ok=True)

    start_time = time.perf_counter()

    sorted_items = sorted(
        merged.items(),
        key=lambda x: (x[0][0], x[0][1], x[1]["start"], x[1]["end"], x[0][2]),
    )

    written_transcripts = 0
    written_exons = 0

    with open(output_gtf, "w") as out_gtf, open(meta_tsv, "w", newline="") as out_meta:
        meta_writer = csv.writer(out_meta, delimiter="\t", lineterminator="\n")
        meta_writer.writerow(["transcript_id"] + source_names)

        for _, ((chrom, strand, ic), item) in enumerate(sorted_items, start=1):
            tx_start = int(item["start"])
            tx_end = int(item["end"])

            if tx_end - tx_start > max_span:
                continue

            written_transcripts += 1
            transcript_id = f"transcript_{written_transcripts}"

            support_set = set(item["sources"])
            support_values = ",".join(sorted(support_set))

            transcript_attrs = {
                "gene_id": transcript_id,
                "transcript_id": transcript_id,
                support_attr: support_values,
            }

            out_gtf.write(
                "\t".join(
                    [
                        chrom,
                        gtf_source,
                        "transcript",
                        str(tx_start),
                        str(tx_end),
                        ".",
                        strand,
                        ".",
                        format_attributes(transcript_attrs),
                    ]
                )
                + "\n"
            )

            meta_writer.writerow(
                [transcript_id] + ["1" if source_name in support_set else "0" for source_name in source_names]
            )

            exons = intron_chain_to_exons(tx_start, tx_end, ic)

            for exon_number, (start, end) in enumerate(exons, start=1):
                exon_attrs = {
                    "gene_id": transcript_id,
                    "transcript_id": transcript_id,
                    "exon_number": str(exon_number),
                    support_attr: support_values,
                }

                out_gtf.write(
                    "\t".join(
                        [
                            chrom,
                            gtf_source,
                            "exon",
                            str(start),
                            str(end),
                            ".",
                            strand,
                            ".",
                            format_attributes(exon_attrs),
                        ]
                    )
                    + "\n"
                )

                written_exons += 1

    logger.info(
        "Finished writing GTF in %s: %d transcripts, %d exons",
        format_seconds(time.perf_counter() - start_time),
        written_transcripts,
        written_exons,
    )
    logger.info("Wrote transcript support meta table: %s", meta_tsv)


def main() -> None:
    total_start = time.perf_counter()

    args = parse_args()

    config_start = time.perf_counter()
    tasks = load_input_config(args.input_gtfs, args.source_column, args.path_column)
    source_names = unique_preserve_order(source_name for source_name, _ in tasks)

    logger.info(
        "Loaded %d input GTFs from %s in %s",
        len(tasks),
        args.input_gtfs,
        format_seconds(time.perf_counter() - config_start),
    )
    logger.info("Using %d worker processes", args.workers)
    logger.info("Detected %d unique source/tool names", len(source_names))

    merged = merge_results(tasks, args.workers)

    logger.info("Merged into %d unique intron chains", len(merged))

    write_gtf(
        merged=merged,
        output_gtf=args.output_gtf,
        gtf_source=args.gtf_source,
        support_attr=args.support_attr,
        max_span=args.max_span,
        source_names=source_names,
    )

    logger.info("Wrote merged GTF: %s", args.output_gtf)
    logger.info("Total runtime: %s", format_seconds(time.perf_counter() - total_start))


if __name__ == "__main__":
    main()
