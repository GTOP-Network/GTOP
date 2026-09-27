#!/usr/bin/env python3

import os
from pathlib import Path

CONFIG_EXCLUDED_SIDS = ['GTOP-CA241-5032-LN-YGE8', 'GTOP-BI281-0087-LN-L7VG', 'GTOP-CB271-4155-LN-A5YE']

import argparse
import re
import shlex
import sys
from collections import defaultdict


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

OUTPUT_DIR = os.path.join(os.getcwd(), "junction_srs_hpc")
PARTITION = os.environ.get("SLURM_PARTITION", os.environ.get("PARTITION", ""))
CPU_PER_NODE = int(os.environ.get("CPU_PER_NODE", "4"))
NT_PER_TASK = int(os.environ.get("NT_PER_TASK", str(CPU_PER_NODE)))
SAMPLES_PER_JOB = int(os.environ.get("SAMPLES_PER_JOB", "20"))

excluded_sids = CONFIG_EXCLUDED_SIDS


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


def load_excluded_sids(path):
    sids = set(excluded_sids)

    if not path:
        return sids

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                sids.add(line.split()[0])

    return sids


def discover_sj_files(sj_dir, suffix, recursive):
    sj_files = []

    if recursive:
        walker = os.walk(sj_dir)
    else:
        walker = [(sj_dir, [], os.listdir(sj_dir))]

    for root, _, files in walker:
        for name in files:
            if name.endswith(suffix) or name.endswith(suffix + ".gz"):
                path = os.path.join(root, name)
                sid = infer_sample_id(path)
                sj_files.append((sid, path))

    return sorted(sj_files, key=lambda x: x[0])


def chunked(items, size):
    for i in range(0, len(items), size):
        yield items[i : i + size]


def summarize_transcript_qc(args):
    merged_path = args.merged_file
    if not os.path.exists(merged_path):
        raise FileNotFoundError(merged_path)

    transcript_stats = defaultdict(lambda: {"junctions": 0, "pass_junctions": 0, "min_support": None})

    with open(merged_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: idx for idx, name in enumerate(header)}

        required = ["transcript_id", "transcript_junction_id", "samples_ge_5"]
        missing = [x for x in required if x not in col]
        if missing:
            raise ValueError(f"missing columns in {merged_path}: {missing}")

        for line in fh:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            tid = fields[col["transcript_id"]]
            support_samples = int(fields[col["samples_ge_5"]])

            item = transcript_stats[tid]
            item["junctions"] += 1
            if support_samples >= args.qc_min_samples:
                item["pass_junctions"] += 1

            if item["min_support"] is None or support_samples < item["min_support"]:
                item["min_support"] = support_samples

    if not transcript_stats:
        raise RuntimeError(f"no transcript rows found in {merged_path}")

    out_path = args.qc_out or f"{merged_path}.transcript_qc.tsv"
    passed_count = 0
    total_count = 0

    with open(out_path, "w") as out:
        out.write(
            "transcript_id\tjunction_count\tpass_junction_count\tmin_samples_ge_5\tqc_pass\n"
        )
        for tid in sorted(transcript_stats):
            item = transcript_stats[tid]
            total_count += 1
            qc_pass = item["junctions"] > 0 and item["junctions"] == item["pass_junctions"]
            if qc_pass:
                passed_count += 1
            out.write(
                f"{tid}\t{item['junctions']}\t{item['pass_junctions']}\t{item['min_support']}\t"
                f"{'PASS' if qc_pass else 'FAIL'}\n"
            )

    ratio = passed_count / total_count if total_count else 0.0

    print(f"[DONE] write: {out_path}", file=sys.stderr)
    print(
        f"[QC] transcripts_pass={passed_count} transcripts_total={total_count} pass_ratio={ratio:.4f}",
        file=sys.stderr,
    )
    print(
        f"QC summary: {passed_count}/{total_count} transcripts passed "
        f"({ratio:.2%}); rule=all junctions require at least {args.qc_min_samples} samples with reads>=5"
    )


def submit_junction_srs_tasks(args):
    step = args.job_name_prefix
    root_dir = args.root_dir
    index_dir = args.index_dir
    unique_junctions = os.path.join(index_dir, "unique_junctions.tsv")
    transcript_junctions = os.path.join(index_dir, "transcript_junctions.tsv")
    junction_py = args.junction_py or os.path.join(CURRENT_DIR, "junction_srs.py")
    work_dir = args.work_dir or CURRENT_DIR

    if not os.path.exists(unique_junctions):
        raise FileNotFoundError(unique_junctions)
    if not os.path.exists(transcript_junctions):
        raise FileNotFoundError(transcript_junctions)
    if not os.path.exists(junction_py):
        raise FileNotFoundError(junction_py)

    os.makedirs(root_dir, exist_ok=True)

    excluded = load_excluded_sids(args.exclude_sids)
    sj_files = discover_sj_files(args.sj_dir, args.suffix, args.recursive)
    sj_files = [(sid, path) for sid, path in sj_files if sid not in excluded]

    if args.max_samples > 0:
        sj_files = sj_files[: args.max_samples]

    if not sj_files:
        raise RuntimeError("no SJ.out.tab files found")

    all_count_list = os.path.join(root_dir, "all_count_files.list")

    with open(all_count_list, "w") as out:
        for sid, _ in sj_files:
            batch_placeholder = "BATCH_PLACEHOLDER"
            # rewritten below with real batch dirs
            _ = batch_placeholder

    submitted_cmds = []
    expected_count_paths = []

    for batch_no, batch in enumerate(chunked(sj_files, args.samples_per_job), 1):
        batch_id = f"batch_{batch_no:04d}"
        node_dir = os.path.join(root_dir, batch_id)
        counts_dir = os.path.join(node_dir, "counts")
        os.makedirs(counts_dir, exist_ok=True)

        sj_list_path = os.path.join(node_dir, "sj_list.tsv")
        with open(sj_list_path, "w") as out:
            for sid, path in batch:
                out.write(f"{sid}\t{path}\n")
                expected_count_paths.append((sid, os.path.join(counts_dir, f"{safe_name(sid)}.junction_counts.tsv")))

        hpc_log_path = os.path.join(node_dir, "job_log")
        hpc_job_shell = os.path.join(node_dir, "job_shell.sh")

        sbatch_resource = f"#SBATCH -N 1 -n {args.cpu_per_node}"
        if args.partition:
            sbatch_resource = f"#SBATCH -p {args.partition} -N 1 -n {args.cpu_per_node}"

        count_cmd = (
            f"python {shlex.quote(junction_py)} count "
            f"--junctions {shlex.quote(unique_junctions)} "
            f"--sj-list {shlex.quote(sj_list_path)} "
            f"--out-dir {shlex.quote(counts_dir)} "
            f"--sample-processes {args.sample_processes} "
            f"--count-mode unique"
        )

        shell_lines = [
            "#!/bin/bash",
            f"#SBATCH -J {step}.{batch_id}",
            f"#SBATCH -o {hpc_log_path}.out",
            f"#SBATCH -e {hpc_log_path}.err",
            sbatch_resource,
            "set -euo pipefail",
            f"mkdir -p {shlex.quote(counts_dir)}",
            f"cd {shlex.quote(work_dir)}",
            count_cmd,
        ]

        with open(hpc_job_shell, "w") as out:
            for line in shell_lines:
                out.write(line + "\n")

        cmd = f"sbatch {shlex.quote(hpc_job_shell)}"
        print(cmd)

        if not args.dry_run:
            __import__('subprocess').run(cmd, shell=True, executable='/bin/bash', check=True)

        submitted_cmds.append(cmd)

    with open(all_count_list, "w") as out:
        for sid, count_path in expected_count_paths:
            out.write(f"{sid}\t{count_path}\n")

    merge_dir = os.path.join(root_dir, "merged")
    os.makedirs(merge_dir, exist_ok=True)
    merge_out_prefix = os.path.join(merge_dir, "srs")

    merge_cmd = (
        f"python {shlex.quote(junction_py)} merge "
        f"--transcript-junctions {shlex.quote(transcript_junctions)} "
        f"--count-list {shlex.quote(all_count_list)} "
        f"--out-prefix {shlex.quote(merge_out_prefix)} "
        f"--skip-missing"
    )

    print(f"submit {len(submitted_cmds)} tasks")
    print(f"sample count: {len(sj_files)}")
    print(f"count list: {all_count_list}")
    print("run merge after all jobs finished:")
    print(merge_cmd)


def main():
    parser = argparse.ArgumentParser(description="Submit STAR SJ.out.tab junction SRS counting jobs to SLURM.")
    parser.add_argument(
        "--merged-file",
        help="merged transcript junction summary file, e.g. srs.transcript_junction_summary.tsv; when set, run transcript QC summary only",
    )
    parser.add_argument("--qc-min-samples", type=int, default=3, help="minimum number of samples with reads>=5 per junction")
    parser.add_argument("--qc-out", help="output transcript QC summary TSV")
    parser.add_argument("--sj-dir", help="directory containing STAR *.SJ.out.tab files")
    parser.add_argument("--index-dir", help="directory from junction_srs.py build-index")
    parser.add_argument("--root-dir", default=OUTPUT_DIR, help="output root for HPC jobs")
    parser.add_argument("--junction-py", default=os.path.join(CURRENT_DIR, "junction_srs.py"))
    parser.add_argument("--work-dir", default=CURRENT_DIR)
    parser.add_argument("--suffix", default=".SJ.out.tab")
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--exclude-sids", help="one sample id per line")
    parser.add_argument("--partition", default=PARTITION)
    parser.add_argument("--cpu-per-node", type=int, default=CPU_PER_NODE)
    parser.add_argument("--sample-processes", type=int, default=NT_PER_TASK)
    parser.add_argument("--samples-per-job", type=int, default=SAMPLES_PER_JOB)
    parser.add_argument("--job-name-prefix", default="srsj")
    parser.add_argument("--max-samples", type=int, default=0, help="debug only; 0 means all samples")
    parser.add_argument("--dry-run", action="store_true", help="write job scripts but do not sbatch")
    args = parser.parse_args()

    if args.merged_file:
        summarize_transcript_qc(args)
        return

    if not args.sj_dir:
        parser.error("--sj-dir is required unless --merged-file is set")
    if not args.index_dir:
        parser.error("--index-dir is required unless --merged-file is set")

    submit_junction_srs_tasks(args)


if __name__ == "__main__":
    main()
