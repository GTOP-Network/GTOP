#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@Desc    : Split a GTF/GFF3 annotation file by chromosome and strand.
           When input is GFF3, output chromosome/strand-split files in GTF
           format with gene_id and transcript_id attributes. Optionally filter
           records by transcript ID.

           Output naming:
             + strand -> chr1-1.gtf
             - strand -> chr1-2.gtf
           Records with unknown strand are skipped.
"""

import argparse
import gzip
import os
from typing import Dict, List, Optional, Set, TextIO
from urllib.parse import unquote


PRIMARY_CHROMOSOME_MAP = {
    **{f"chr{i}": f"chr{i}" for i in range(1, 23)},
    **{str(i): f"chr{i}" for i in range(1, 23)},
    "chrx": "chrX",
    "x": "chrX",
    "chry": "chrY",
    "y": "chrY",
    "chrm": "chrM",
    "chrmt": "chrM",
    "m": "chrM",
    "mt": "chrM",
    "nc_000001.11": "chr1",
    "nc_000002.12": "chr2",
    "nc_000003.12": "chr3",
    "nc_000004.12": "chr4",
    "nc_000005.10": "chr5",
    "nc_000006.12": "chr6",
    "nc_000007.14": "chr7",
    "nc_000008.11": "chr8",
    "nc_000009.12": "chr9",
    "nc_000010.11": "chr10",
    "nc_000011.10": "chr11",
    "nc_000012.12": "chr12",
    "nc_000013.11": "chr13",
    "nc_000014.9": "chr14",
    "nc_000015.10": "chr15",
    "nc_000016.10": "chr16",
    "nc_000017.11": "chr17",
    "nc_000018.10": "chr18",
    "nc_000019.10": "chr19",
    "nc_000020.11": "chr20",
    "nc_000021.9": "chr21",
    "nc_000022.11": "chr22",
    "nc_000023.11": "chrX",
    "nc_000024.10": "chrY",
    "nc_012920.1": "chrM",
}


TRANSCRIPT_FEATURE_TYPES = {
    "mrna",
    "transcript",
    "ncrna",
    "lncrna",
    "rrna",
    "trna",
    "snrna",
    "snorna",
    "mirna",
    "primary_transcript",
    "pseudogenic_transcript",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Split a GTF/GFF3 file by chromosome and strand. "
            "When input is GFF3, output split files in GTF format. "
            "Non-primary chromosomes are grouped into chr_other. "
            "Unknown-strand records are skipped."
        )
    )
    parser.add_argument("input_path", help="Input GTF/GFF3 path. Supports .gz files.")
    parser.add_argument("output_dir", help="Output directory.")
    parser.add_argument(
        "-t",
        "--transcript-id-list",
        default=None,
        help=(
            "Optional text file containing transcript IDs, one ID per line. "
            "Only records whose transcript_id is in this file will be written. "
            "Supports .gz files."
        ),
    )
    return parser.parse_args()


def open_maybe_gzip(path: str, mode: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode, encoding="utf-8")


def detect_output_suffix(input_path: str) -> str:
    return ".gtf"


def normalize_chromosome(seqname: str) -> str:
    seqname_key = seqname.strip().lower()
    return PRIMARY_CHROMOSOME_MAP.get(seqname_key, "chr_other")


def normalize_strand(strand: str) -> Optional[str]:
    if strand == "+":
        return "1"
    if strand == "-":
        return "2"
    return None


def is_gff3_input(input_path: str) -> bool:
    with open_maybe_gzip(input_path, "r") as handle:
        for raw_line in handle:
            if raw_line.startswith("##FASTA"):
                break
            if not raw_line.strip():
                continue
            if raw_line.startswith("##gff-version"):
                return True
            if raw_line.startswith("#"):
                continue

            columns = raw_line.rstrip("\n").split("\t")
            if len(columns) < 9:
                continue

            first_attr = columns[8].split(";", 1)[0].strip()
            if not first_attr:
                continue
            return "=" in first_attr

    return False


def format_gtf_value(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')


def load_transcript_id_filter(path: Optional[str]) -> Optional[Set[str]]:
    if not path:
        return None

    transcript_ids: Set[str] = set()
    with open_maybe_gzip(path, "r") as handle:
        for raw_line in handle:
            transcript_id = raw_line.strip()
            if transcript_id:
                transcript_ids.add(transcript_id)

    return transcript_ids


def parse_gff3_attributes(attributes: str) -> Dict[str, str]:
    attributes = attributes.strip()
    if not attributes or attributes == ".":
        return {}

    parsed = {}
    for item in attributes.split(";"):
        item = item.strip()
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            key = key.strip()
            value_parts = [unquote(part.strip()) for part in value.strip().split(",")]
            parsed[key] = ",".join(value_parts)
        else:
            parsed[item.strip()] = ""

    return parsed


def parse_gtf_attributes(attributes: str) -> Dict[str, str]:
    attributes = attributes.strip()
    if not attributes or attributes == ".":
        return {}

    parsed = {}
    for item in attributes.split(";"):
        item = item.strip()
        if not item:
            continue

        if " " not in item:
            parsed[item] = ""
            continue

        key, value = item.split(None, 1)
        value = value.strip()

        if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
            value = value[1:-1]

        parsed[key] = value.replace('\\"', '"').replace("\\\\", "\\")

    return parsed


def get_gtf_transcript_id(raw_line: str) -> str:
    columns = raw_line.rstrip("\n").split("\t")
    if len(columns) < 9:
        raise ValueError(f"Invalid annotation line with fewer than 9 columns: {raw_line.rstrip()}")

    return parse_gtf_attributes(columns[8]).get("transcript_id", "")


def should_keep_transcript_line(raw_line: str, transcript_id_filter: Optional[Set[str]]) -> bool:
    if transcript_id_filter is None:
        return True

    transcript_id = get_gtf_transcript_id(raw_line)
    return transcript_id in transcript_id_filter


def get_first_attr(attributes: Dict[str, str], keys: List[str]) -> str:
    for key in keys:
        value = attributes.get(key, "")
        if value:
            return value
    return ""


def split_parent_ids(parent_value: str) -> List[str]:
    return [parent_id.strip() for parent_id in parent_value.split(",") if parent_id.strip()]


def add_if_present(items: List[tuple], key: str, value: str):
    if value:
        items.append((key, value))


def build_gtf_attributes(
    feature_type: str,
    gff3_attrs: Dict[str, str],
    gene_id_by_feature_id: Dict[str, str],
    transcript_id_by_feature_id: Dict[str, str],
) -> str:
    feature_type_lower = feature_type.lower()
    feature_id = gff3_attrs.get("ID", "")
    parent_ids = split_parent_ids(gff3_attrs.get("Parent", ""))

    gene_id = get_first_attr(gff3_attrs, ["gene_id", "gene", "gene_name"])
    transcript_id = get_first_attr(gff3_attrs, ["transcript_id", "transcript"])

    if feature_type_lower == "gene":
        gene_id = gene_id or feature_id or gff3_attrs.get("Name", "")
        if feature_id and gene_id:
            gene_id_by_feature_id[feature_id] = gene_id

    elif feature_type_lower in TRANSCRIPT_FEATURE_TYPES:
        transcript_id = transcript_id or feature_id or gff3_attrs.get("Name", "")
        if not gene_id and parent_ids:
            gene_id = gene_id_by_feature_id.get(parent_ids[0], parent_ids[0])

        if feature_id:
            if gene_id:
                gene_id_by_feature_id[feature_id] = gene_id
            if transcript_id:
                transcript_id_by_feature_id[feature_id] = transcript_id

    else:
        if parent_ids:
            parent_id = parent_ids[0]
            transcript_id = transcript_id or transcript_id_by_feature_id.get(parent_id, parent_id)
            gene_id = gene_id or gene_id_by_feature_id.get(parent_id, "")

        if not gene_id and transcript_id:
            gene_id = gene_id_by_feature_id.get(transcript_id, "")

        if not gene_id and parent_ids:
            gene_id = gene_id_by_feature_id.get(parent_ids[0], "")

        if feature_id:
            if gene_id:
                gene_id_by_feature_id[feature_id] = gene_id
            if transcript_id:
                transcript_id_by_feature_id[feature_id] = transcript_id

    ordered_attrs = []
    add_if_present(ordered_attrs, "gene_id", gene_id)
    add_if_present(ordered_attrs, "transcript_id", transcript_id)

    for key, value in gff3_attrs.items():
        if key in {"gene_id", "transcript_id"}:
            continue
        ordered_attrs.append((key, value))

    if not ordered_attrs:
        return "."

    return " ".join(f'{key} "{format_gtf_value(value)}";' for key, value in ordered_attrs)


def gff3_line_to_gtf(
    raw_line: str,
    gene_id_by_feature_id: Dict[str, str],
    transcript_id_by_feature_id: Dict[str, str],
) -> str:
    columns = raw_line.rstrip("\n").split("\t")
    if len(columns) < 9:
        raise ValueError(f"Invalid annotation line with fewer than 9 columns: {raw_line.rstrip()}")

    gff3_attrs = parse_gff3_attributes(columns[8])
    columns[8] = build_gtf_attributes(
        columns[2],
        gff3_attrs,
        gene_id_by_feature_id,
        transcript_id_by_feature_id,
    )
    return "\t".join(columns) + "\n"


def get_writer(
    split_key: str,
    output_dir: str,
    output_suffix: str,
    header_lines: List[str],
    writers: Dict[str, TextIO],
) -> TextIO:
    if split_key not in writers:
        output_path = os.path.join(output_dir, f"{split_key}{output_suffix}")
        handle = open(output_path, "w", encoding="utf-8")
        if header_lines:
            handle.writelines(header_lines)
        writers[split_key] = handle
    return writers[split_key]


def manifest_sort_key(split_key: str):
    chrom, _, strand = split_key.partition("-")

    if chrom == "chr_other":
        chrom_rank = (2, 0)
    elif chrom == "chrX":
        chrom_rank = (1, 23)
    elif chrom == "chrY":
        chrom_rank = (1, 24)
    elif chrom == "chrM":
        chrom_rank = (1, 25)
    else:
        try:
            chrom_rank = (1, int(chrom.replace("chr", "")))
        except ValueError:
            chrom_rank = (1, 99)

    try:
        strand_rank = int(strand)
    except ValueError:
        strand_rank = 99

    return (*chrom_rank, strand_rank)


def split_annotation(
    input_path: str,
    output_dir: str,
    transcript_id_list: Optional[str] = None,
) -> Dict[str, str]:
    os.makedirs(output_dir, exist_ok=True)

    output_suffix = detect_output_suffix(input_path)
    convert_gff3_to_gtf = is_gff3_input(input_path)
    transcript_id_filter = load_transcript_id_filter(transcript_id_list)

    writers: Dict[str, TextIO] = {}
    output_paths: Dict[str, str] = {}
    header_lines: List[str] = []
    gene_id_by_feature_id: Dict[str, str] = {}
    transcript_id_by_feature_id: Dict[str, str] = {}

    first_feature_seen = False
    current_split_key: Optional[str] = None

    with open_maybe_gzip(input_path, "r") as reader:
        for raw_line in reader:
            if raw_line.startswith("##FASTA"):
                break

            if not raw_line.strip():
                continue

            if raw_line.startswith("#"):
                if not first_feature_seen:
                    header_lines.append(raw_line)
                elif current_split_key is not None:
                    writer = get_writer(
                        current_split_key,
                        output_dir,
                        output_suffix,
                        header_lines,
                        writers,
                    )
                    writer.write(raw_line)
                continue

            first_feature_seen = True

            columns = raw_line.rstrip("\n").split("\t")
            if len(columns) < 9:
                raise ValueError(f"Invalid annotation line with fewer than 9 columns: {raw_line.rstrip()}")

            seqname = columns[0]
            strand = columns[6]

            strand_suffix = normalize_strand(strand)
            if strand_suffix is None:
                continue

            chrom = normalize_chromosome(seqname)
            split_key = f"{chrom}-{strand_suffix}"

            if convert_gff3_to_gtf:
                output_line = gff3_line_to_gtf(
                    raw_line,
                    gene_id_by_feature_id,
                    transcript_id_by_feature_id,
                )
            else:
                output_line = raw_line

            if not should_keep_transcript_line(output_line, transcript_id_filter):
                continue

            current_split_key = split_key

            writer = get_writer(
                split_key,
                output_dir,
                output_suffix,
                header_lines,
                writers,
            )
            writer.write(output_line)

            output_paths[split_key] = os.path.join(output_dir, f"{split_key}{output_suffix}")

    for writer in writers.values():
        writer.close()

    manifest_path = os.path.join(output_dir, "split_manifest.tsv")
    with open(manifest_path, "w", encoding="utf-8") as handle:
        handle.write("chromosome_strand\tpath\n")
        for split_key in sorted(output_paths.keys(), key=manifest_sort_key):
            handle.write(f"{split_key}\t{output_paths[split_key]}\n")

    return output_paths


def main():
    args = parse_args()
    output_paths = split_annotation(
        args.input_path,
        args.output_dir,
        args.transcript_id_list,
    )

    for split_key in sorted(output_paths.keys(), key=manifest_sort_key):
        print(f"{split_key}\t{output_paths[split_key]}")

    print(f"manifest\t{os.path.join(args.output_dir, 'split_manifest.tsv')}")


if __name__ == "__main__":
    main()
