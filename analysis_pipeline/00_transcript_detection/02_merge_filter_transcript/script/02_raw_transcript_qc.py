#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, Optional, Set, Tuple, Union

import pandas as pd

from split_gtf_by_chr_helper import (
    TRANSCRIPT_FEATURE_TYPES,
    get_first_attr,
    is_gff3_input,
    load_transcript_id_filter,
    open_maybe_gzip,
    parse_gff3_attributes,
    parse_gtf_attributes,
)
from pathlib import Path

EXCLUDED_SIDS = ['GTOP-CA241-5032-LN-YGE8', 'GTOP-BI281-0087-LN-L7VG', 'GTOP-CB271-4155-LN-A5YE']

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
SAMPLE_BASED_DIR = str(Path(PROJECT_DIR) / 'output/assembly/LRS/sample_based')
MAIN_DIR = str(Path(PROJECT_DIR) / 'output/assembly/LRS/tx_merge')

DEFAULT_TOOLS = ("isoquant", "flames", "talon", "isotools", "isoseq", "flair", "bambu")
TOOLS = DEFAULT_TOOLS
N_TASK = int(os.environ.get("RAW_QC_WORKERS", "180"))
OUTPUT_DIR = f"{MAIN_DIR}/raw_transcript_qc"
FAIL_SAMPLE_TABLE = f"{OUTPUT_DIR}/fail_sample_tool.tsv"
RAW_TRANSCRIPT_COUNT_TABLE = f"{OUTPUT_DIR}/raw_transcript_count.tsv"
FILTERED_TRANSCRIPT_COUNT_TABLE = f"{OUTPUT_DIR}/raw_transcript_count.filtered.tsv"
CHECK_MIN_GTF_SIZE_KB = 10

# Set this to a transcript ID list path to apply the same filter to every GTF/GFF.
# Keep it as None to skip this global filter.
TRANSCRIPT_ID_LIST = None

# Set this to True to use per-sample tool filter_tx lists from
# get_tools_transcript_list_path() for the filtered table.
USE_TOOL_FILTER_LIST = True
MP_START_METHOD = "fork"

TranscriptListSpec = Union[
    None,
    str,
    Dict[str, Optional[str]],
    Dict[str, Dict[str, Optional[str]]],
]


def _strip_known_suffix(path: str) -> str:
    base = os.path.basename(path)
    for suffix in (".gtf.gz", ".gff3.gz", ".gff.gz", ".gtf", ".gff3", ".gff"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


def _is_gff3_annotation(annotation_path: str) -> bool:
    if not is_gff3_input(annotation_path):
        return False

    base = annotation_path[:-3] if annotation_path.endswith(".gz") else annotation_path
    lower_base = base.lower()
    if lower_base.endswith(".gff3"):
        return True

    with open_maybe_gzip(annotation_path, "r") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            columns = raw_line.rstrip("\n").split("\t")
            if len(columns) < 9:
                continue

            first_attr = columns[8].split(";", 1)[0].strip()
            return "=" in first_attr

    return False


def get_tools_gtf_path(tools=DEFAULT_TOOLS) -> Dict[str, Dict[str, str]]:
    patterns = {
        'bambu': 'discoveryOnly.ndr.default.gtf',
        'flair': '{sid}.isoforms.gtf',
        'flames': 'isoform_annotated.gff3',
        'isoquant': 'OUT/OUT.transcript_models.gtf',
        'isoseq': '{sid}.collapsed.gff',
        'isotools': '{sid}.isotools.gtf',
        'talon': 'talon_observed*.gtf',
    }
    paths = {tool: {} for tool in tools}
    for sample_dir in sorted(Path(SAMPLE_BASED_DIR).iterdir()):
        sid = sample_dir.name
        if not sample_dir.is_dir() or not sid.startswith('GTOP') or sid in EXCLUDED_SIDS:
            continue
        for tool in tools:
            tool_dir = sample_dir / 'hg38' / tool
            pattern = patterns[tool].format(sid=sid)
            matches = sorted(tool_dir.glob(pattern))
            if len(matches) > 1:
                raise ValueError(f'Multiple {tool} annotations for {sid}: {matches}')
            paths[tool][sid] = str(matches[0] if matches else tool_dir / pattern)
    if not any(paths.values()):
        raise ValueError(f'No samples found in {SAMPLE_BASED_DIR}')
    return paths



def get_tools_transcript_list_path(tools=DEFAULT_TOOLS) -> Dict[str, Dict[str, Optional[str]]]:
    paths = {}
    for tool, samples in get_tools_gtf_path(tools).items():
        paths[tool] = {}
        for sid in samples:
            paths[tool][sid] = (
                str(Path(SAMPLE_BASED_DIR) / sid / 'hg38' / tool / 'filter_tx' / f'{sid}.txt')
                if tool in {'isoseq', 'isotools', 'talon'} else None
            )
    return paths


def _resolve_transcript_list_path(
    transcript_id_list: TranscriptListSpec,
    tool_name: str,
    sample_id: str,
) -> Optional[str]:
    if transcript_id_list is None or isinstance(transcript_id_list, str):
        return transcript_id_list

    tool_value = transcript_id_list.get(tool_name)
    if isinstance(tool_value, dict):
        return tool_value.get(sample_id)
    if tool_value is not None:
        return tool_value

    sample_value = transcript_id_list.get(sample_id)
    if isinstance(sample_value, str) or sample_value is None:
        return sample_value

    return None


def _count_gtf_transcripts(gtf_path: str, transcript_id_filter: Optional[Set[str]]) -> int:
    transcript_ids = set()
    with open_maybe_gzip(gtf_path, "r") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            columns = raw_line.rstrip("\n").split("\t")
            if len(columns) < 9:
                continue
            transcript_id = parse_gtf_attributes(columns[8]).get("transcript_id", "")
            if transcript_id and (transcript_id_filter is None or transcript_id in transcript_id_filter):
                transcript_ids.add(transcript_id)
    return len(transcript_ids)


def _count_gff3_transcripts(gff_path: str, transcript_id_filter: Optional[Set[str]]) -> int:
    transcript_ids = set()
    transcript_id_by_feature_id = {}

    with open_maybe_gzip(gff_path, "r") as handle:
        for raw_line in handle:
            if raw_line.startswith("##FASTA"):
                break
            if not raw_line.strip() or raw_line.startswith("#"):
                continue

            columns = raw_line.rstrip("\n").split("\t")
            if len(columns) < 9:
                continue

            feature_type = columns[2].lower()
            attrs = parse_gff3_attributes(columns[8])
            feature_id = attrs.get("ID", "")
            transcript_id = get_first_attr(attrs, ["transcript_id", "transcript"])

            if feature_type in TRANSCRIPT_FEATURE_TYPES:
                transcript_id = transcript_id or feature_id or attrs.get("Name", "")
            elif attrs.get("Parent"):
                parent_id = attrs["Parent"].split(",")[0].strip()
                transcript_id = transcript_id or transcript_id_by_feature_id.get(parent_id, parent_id)

            if feature_id and transcript_id:
                transcript_id_by_feature_id[feature_id] = transcript_id

            if transcript_id and (transcript_id_filter is None or transcript_id in transcript_id_filter):
                transcript_ids.add(transcript_id)

    return len(transcript_ids)


def count_transcripts_in_annotation(gtf_path: str, transcript_id_list: Optional[str] = None) -> int:
    if not os.path.isfile(gtf_path):
        raise FileNotFoundError(f"Annotation file not found: {gtf_path}")

    transcript_id_filter = load_transcript_id_filter(transcript_id_list)
    if _is_gff3_annotation(gtf_path):
        return _count_gff3_transcripts(gtf_path, transcript_id_filter)
    return _count_gtf_transcripts(gtf_path, transcript_id_filter)


def _count_one_task(task: Tuple[str, str, str, Optional[str]]) -> Tuple[str, str, int]:
    tool_name, sample_id, gtf_path, transcript_list_path = task
    count = count_transcripts_in_annotation(gtf_path, transcript_list_path)
    return tool_name, sample_id, count


def build_transcript_count_table(
    tools_gtf_path: Optional[Dict[str, Dict[str, str]]] = None,
    transcript_id_list: TranscriptListSpec = None,
    max_workers: int = 8,
) -> pd.DataFrame:
    if tools_gtf_path is None:
        tools_gtf_path = get_tools_gtf_path()

    tasks = []
    for tool_name, sample_paths in tools_gtf_path.items():
        for sample_id, gtf_path in sample_paths.items():
            transcript_list_path = _resolve_transcript_list_path(transcript_id_list, tool_name, sample_id)
            tasks.append((tool_name, sample_id, gtf_path, transcript_list_path))

    rows = []
    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=mp.get_context(MP_START_METHOD),
    ) as executor:
        futures = {executor.submit(_count_one_task, task): task for task in tasks}
        for future in as_completed(futures):
            tool_name, sample_id, _, transcript_list_path = futures[future]
            count_tool, count_sample, count = future.result()
            logger.info(
                "Counted %s %s: %s transcripts%s",
                count_tool,
                count_sample,
                count,
                f" after filtering with {transcript_list_path}" if transcript_list_path else "",
            )
            rows.append({"sample_id": sample_id, "tool": tool_name, "transcript_count": count})

    if not rows:
        return pd.DataFrame()

    long_df = pd.DataFrame(rows)
    table = long_df.pivot(index="sample_id", columns="tool", values="transcript_count")
    return table.reindex(sorted(table.index)).reindex(sorted(table.columns), axis=1).astype("Int64")


def write_transcript_count_table(
    output_path: str,
    tools_gtf_path: Optional[Dict[str, Dict[str, str]]] = None,
    transcript_id_list: TranscriptListSpec = None,
    max_workers: int = 8,
) -> pd.DataFrame:
    table = build_transcript_count_table(
        tools_gtf_path=tools_gtf_path,
        transcript_id_list=transcript_id_list,
        max_workers=max_workers,
    )
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    table.to_csv(output_path, sep="\t")
    logger.info("Saved transcript count table: %s", output_path)
    return table


def check(
    tools=TOOLS,
    min_size_kb: int = CHECK_MIN_GTF_SIZE_KB,
    output_path: str = FAIL_SAMPLE_TABLE,
    tools_gtf_path: Optional[Dict[str, Dict[str, str]]] = None,
) -> pd.DataFrame:
    if tools_gtf_path is None:
        tools_gtf_path = get_tools_gtf_path(tools)

    min_size_bytes = min_size_kb * 1024
    rows = []
    for tool_name, sample_paths in tools_gtf_path.items():
        for sample_id, gtf_path in sample_paths.items():
            exists = os.path.isfile(gtf_path)
            size_bytes = os.path.getsize(gtf_path) if exists else 0
            if not exists or size_bytes < min_size_bytes:
                rows.append(
                    {
                        "sample_id": sample_id,
                        "tool": tool_name,
                        "status": "missing" if not exists else f"<{min_size_kb}KB",
                        "size_kb": round(size_bytes / 1024, 2),
                        "gtf_path": gtf_path,
                    }
                )

    failed_table = pd.DataFrame(
        rows,
        columns=["sample_id", "tool", "status", "size_kb", "gtf_path"],
    )
    if output_path:
        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        failed_table.to_csv(output_path, sep="\t", index=False)
        logger.info("Saved failed sample/tool table: %s", output_path)

    if failed_table.empty:
        logger.info("No failed GTF/GFF found. Criterion: missing or smaller than %sKB.", min_size_kb)
    else:
        failed_table = failed_table.sort_values(["sample_id", "tool"]).reset_index(drop=True)
        print(failed_table.to_csv(sep="\t", index=False), end="")
        logger.info(
            "Found %s failed GTF/GFF records. Criterion: missing or smaller than %sKB.",
            failed_table.shape[0],
            min_size_kb,
        )

    return failed_table


def run():
    tools_gtf_path = get_tools_gtf_path(TOOLS)
    failed = check(tools_gtf_path=tools_gtf_path)
    if (failed['status'] == 'missing').any():
        raise FileNotFoundError(f'Missing annotations; see {FAIL_SAMPLE_TABLE}')

    raw_table = write_transcript_count_table(
        output_path=RAW_TRANSCRIPT_COUNT_TABLE,
        tools_gtf_path=tools_gtf_path,
        transcript_id_list=None,
        max_workers=N_TASK,
    )
    logger.info(
        "Finished raw table. Matrix shape: %s samples x %s tools",
        raw_table.shape[0],
        raw_table.shape[1],
    )

    transcript_id_list: TranscriptListSpec = TRANSCRIPT_ID_LIST
    if USE_TOOL_FILTER_LIST:
        transcript_id_list = get_tools_transcript_list_path(TOOLS)

    if transcript_id_list is None:
        logger.info("Skip filtered table because TRANSCRIPT_ID_LIST is None and USE_TOOL_FILTER_LIST is False.")
        return

    filtered_table = write_transcript_count_table(
        output_path=FILTERED_TRANSCRIPT_COUNT_TABLE,
        tools_gtf_path=tools_gtf_path,
        transcript_id_list=transcript_id_list,
        max_workers=N_TASK,
    )
    logger.info(
        "Finished filtered table. Matrix shape: %s samples x %s tools",
        filtered_table.shape[0],
        filtered_table.shape[1],
    )


if __name__ == "__main__":
    run()
