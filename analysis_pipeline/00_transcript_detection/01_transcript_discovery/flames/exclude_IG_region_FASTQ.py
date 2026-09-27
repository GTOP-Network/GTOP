# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/8/18 17:57
@Email   : xuechao@szbl.ac.cn
@Desc    : Export non-IG long reads from one or more pbmm2 BAM files as FASTQ.GZ.
"""


import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_FLAMES_FASTQ_DIR = os.environ.get('FLAMES_FASTQ_DIR', str(Path(CONFIG_PROJECT_DIR) / 'input/flames_non_ig'))

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import argparse
import csv
import shutil
import subprocess
import sys
import traceback


IG_REGIONS = {
    "chr14": [(103000000, 107500000)],  # IGH
    "chr2": [(88000000, 90000000)],  # IGK
    "chr22": [(22000000, 23500000)],  # IGL
}

DEFAULT_OUTPUT_ROOT = Path(CONFIG_FLAMES_FASTQ_DIR)
DEFAULT_PIGZ = 'pigz'
DEFAULT_SAMTOOLS = 'samtools'
DEFAULT_SORT = "sort"
DEFAULT_AWK = "awk"
GZIP_COMPRESSLEVEL = 1


class SampleJob:
    def __init__(
        self,
        sample_id,
        input_bam,
        output_fastq,
        bam_threads,
        gzip_threads,
        samtools,
        pigz,
        sort=DEFAULT_SORT,
        awk=DEFAULT_AWK,
        force=False,
    ):
        self.sample_id = sample_id
        self.input_bam = input_bam
        self.output_fastq = output_fastq
        self.bam_threads = bam_threads
        self.gzip_threads = gzip_threads
        self.samtools = samtools
        self.pigz = pigz
        self.sort = sort
        self.awk = awk
        self.force = force

    @property
    def output_dir(self):
        return self.output_fastq.parent


def log(message, log_file=None, stream=sys.stdout):
    line = "[{}] {}".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), message)
    print(line, file=stream)
    stream.flush()
    if log_file:
        log_dir = os.path.dirname(log_file)
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)
        with open(log_file, "a") as handle:
            handle.write(line + "\n")


def has_bam_index(input_bam):
    input_bam = Path(input_bam)
    return (
        Path(str(input_bam) + ".bai").is_file()
        or input_bam.with_suffix(".bai").is_file()
        or Path(str(input_bam) + ".csi").is_file()
    )


def resolve_executable(executable):
    if Path(executable).is_file():
        return executable
    resolved = shutil.which(executable)
    if resolved:
        return resolved
    raise RuntimeError("executable not found: {}".format(executable))


def samtools_reference_names(input_bam, samtools):
    completed = subprocess.run(
        [samtools, "view", "-H", str(input_bam)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "samtools view -H failed with exit_code={}:\n{}".format(completed.returncode, completed.stderr)
        )

    reference_names = set()
    for line in completed.stdout.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        for field in line.split("\t"):
            if field.startswith("SN:"):
                reference_names.add(field[3:])
                break
    return reference_names


def samtools_ig_region_args(reference_names):
    region_args = []
    for chrom, regions in sorted(IG_REGIONS.items()):
        if chrom in reference_names:
            reference_name = chrom
        elif chrom.startswith("chr") and chrom[3:] in reference_names:
            reference_name = chrom[3:]
        elif "chr" + chrom in reference_names:
            reference_name = "chr" + chrom
        else:
            continue

        for region_start, region_end in regions:
            region_args.append("{}:{}-{}".format(reference_name, region_start + 1, region_end))
    return region_args


def create_excluded_read_name_file(job, excluded_names):
    if not has_bam_index(job.input_bam):
        raise RuntimeError("BAM index is required for samtools region fetch: {}".format(job.input_bam))

    samtools = resolve_executable(job.samtools)
    sort = resolve_executable(job.sort)
    region_args = samtools_ig_region_args(samtools_reference_names(job.input_bam, samtools))
    if not region_args:
        log("[WARN] {}: no IG reference names found in BAM header".format(job.sample_id), stream=sys.stderr)
        excluded_names.write_text("")
        return 0

    view_command = [samtools, "view", "-@", str(job.bam_threads), str(job.input_bam)] + region_args
    sort_command = [sort, "-u"]
    log("[CMD] {} | cut -f1 | {}".format(" ".join(view_command), " ".join(sort_command)))

    with open(str(excluded_names), "wb") as output_handle:
        view_process = subprocess.Popen(view_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cut_process = subprocess.Popen(["cut", "-f", "1"], stdin=view_process.stdout, stdout=subprocess.PIPE)
        sort_process = subprocess.Popen(
            sort_command, stdin=cut_process.stdout, stdout=output_handle, stderr=subprocess.PIPE
        )
        view_process.stdout.close()
        cut_process.stdout.close()

        view_stderr = view_process.stderr.read().decode("utf-8", errors="replace")
        sort_stderr = sort_process.stderr.read().decode("utf-8", errors="replace")
        view_exit = view_process.wait()
        cut_exit = cut_process.wait()
        sort_exit = sort_process.wait()

    if view_exit != 0 or cut_exit != 0 or sort_exit != 0:
        raise RuntimeError(
            "failed to create excluded read names: samtools={}, cut={}, sort={}\n{}\n{}".format(
                view_exit, cut_exit, sort_exit, view_stderr, sort_stderr
            )
        )

    with open(str(excluded_names), "rb") as input_handle:
        return sum(1 for _ in input_handle)


def run_samtools_fastq_pipeline(job, excluded_names, tmp_fastq):
    awk_program = r'''
NR == FNR {
    drop[$1] = 1
    next
}
{
    read_name = $1
    if ((read_name in drop) || (read_name in seen)) {
        next
    }
    seen[read_name] = 1
    seq = $10
    if (seq == "*" || seq == "") {
        next
    }
    qual = $11
    if (qual == "*" || length(qual) != length(seq)) {
        qual = ""
        for (i = 1; i <= length(seq); i++) {
            qual = qual "I"
        }
    }
    print "@" read_name "\n" seq "\n+\n" qual
}
'''
    samtools = resolve_executable(job.samtools)
    awk = resolve_executable(job.awk)
    pigz = resolve_executable(job.pigz)
    view_command = [samtools, "view", "-@", str(job.bam_threads), str(job.input_bam)]
    awk_command = [awk, awk_program, str(excluded_names), "-"]
    pigz_command = [pigz, "-p", str(job.gzip_threads), "-{}".format(GZIP_COMPRESSLEVEL), "-c"]
    log("[CMD] {} | {} <excluded_names> - | {}".format(" ".join(view_command), awk, " ".join(pigz_command)))

    with open(str(tmp_fastq), "wb") as output_handle:
        view_process = subprocess.Popen(view_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        awk_process = subprocess.Popen(
            awk_command, stdin=view_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        pigz_process = subprocess.Popen(
            pigz_command, stdin=awk_process.stdout, stdout=output_handle, stderr=subprocess.PIPE
        )
        view_process.stdout.close()
        awk_process.stdout.close()

        view_stderr = view_process.stderr.read().decode("utf-8", errors="replace")
        awk_stderr = awk_process.stderr.read().decode("utf-8", errors="replace")
        pigz_stderr = pigz_process.stderr.read().decode("utf-8", errors="replace")
        view_exit = view_process.wait()
        awk_exit = awk_process.wait()
        pigz_exit = pigz_process.wait()

    if view_exit != 0 or awk_exit != 0 or pigz_exit != 0:
        raise RuntimeError(
            "FASTQ pipeline failed: samtools={}, awk={}, pigz={}\n{}\n{}\n{}".format(
                view_exit, awk_exit, pigz_exit, view_stderr, awk_stderr, pigz_stderr
            )
        )


def run_one_sample(job):
    """Run one sample: BAM -> non-IG FASTQ.GZ."""
    tmp_fastq = job.output_fastq.with_suffix(job.output_fastq.suffix + ".tmp")
    excluded_names = job.output_dir / "{}.IG_region.read_names.tmp".format(job.sample_id)

    try:
        job.output_dir.mkdir(parents=True, exist_ok=True)

        if job.output_fastq.exists() and not job.force:
            return job.sample_id, True, "{}: output exists, skipped: {}".format(job.sample_id, job.output_fastq)

        if tmp_fastq.exists():
            tmp_fastq.unlink()
        if excluded_names.exists():
            excluded_names.unlink()

        log("[START] {}: input={}, output={}".format(job.sample_id, job.input_bam, job.output_fastq))
        excluded_read_count = create_excluded_read_name_file(job, excluded_names)
        run_samtools_fastq_pipeline(job, excluded_names, tmp_fastq)
        os.replace(str(tmp_fastq), str(job.output_fastq))

        if excluded_names.exists():
            excluded_names.unlink()

        message = "{}: excluded_reads={}, output={}, bytes={}".format(
            job.sample_id,
            excluded_read_count,
            job.output_fastq,
            job.output_fastq.stat().st_size,
        )
        return job.sample_id, True, message
    except Exception:
        if tmp_fastq.exists():
            tmp_fastq.unlink()
        if excluded_names.exists():
            excluded_names.unlink()
        return job.sample_id, False, "{} failed:\n{}".format(job.sample_id, traceback.format_exc())


def load_jobs_from_conf(conf_csv, output_root, nt_per_task, samtools, pigz, force=False, sort=DEFAULT_SORT, awk=DEFAULT_AWK):
    jobs = []
    with open(conf_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sample_id = row["sample_id"]
            input_bam = Path(row["input_bam"])
            output_fastq = row.get("output_fastq", "").strip()
            if output_fastq:
                output_fastq = Path(output_fastq)
            else:
                output_fastq = Path(output_root) / sample_id / "{}.exclude_IG_region.fastq.gz".format(sample_id)

            jobs.append(
                SampleJob(
                    sample_id=sample_id,
                    input_bam=input_bam,
                    output_fastq=output_fastq,
                    bam_threads=max(1, int(nt_per_task)),
                    gzip_threads=max(1, int(nt_per_task)),
                    samtools=row.get("samtools", "").strip() or samtools,
                    pigz=row.get("pigz", "").strip() or pigz,
                    sort=row.get("sort", "").strip() or sort,
                    awk=row.get("awk", "").strip() or awk,
                    force=force,
                )
            )
    return jobs


def run_samples_from_conf(conf_csv, output_root, log_name, n_task, nt_per_task, samtools, pigz, force=False):
    """Read sample config and run multiple samples in parallel."""
    jobs = load_jobs_from_conf(conf_csv, output_root, nt_per_task, samtools, pigz, force=force)
    if not jobs:
        log("[ERROR] no sample found in config: {}".format(conf_csv), log_file=log_name, stream=sys.stderr)
        return 1

    worker_count = max(1, min(int(n_task), len(jobs)))
    log("[INFO] samples={}, workers={}, nt_per_task={}".format(len(jobs), worker_count, nt_per_task), log_file=log_name)
    failed = []

    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = [executor.submit(run_one_sample, job) for job in jobs]
        for future in as_completed(futures):
            sample_id, ok, message = future.result()
            if ok:
                log("[OK] {}".format(message), log_file=log_name)
            else:
                failed.append(sample_id)
                log("[ERROR] {}".format(message), log_file=log_name, stream=sys.stderr)

    if failed:
        log("[ERROR] failed samples: {}".format(",".join(sorted(failed))), log_file=log_name, stream=sys.stderr)
        return 1
    log("[INFO] all samples finished", log_file=log_name)
    return 0


def check_outputs_from_conf(conf_csv, output_root, log_name, nt_per_task, samtools, pigz):
    jobs = load_jobs_from_conf(conf_csv, output_root, nt_per_task, samtools, pigz)
    missing = []
    empty = []
    for job in jobs:
        if not job.output_fastq.exists():
            missing.append(job.sample_id)
            log("[MISSING] {}: {}".format(job.sample_id, job.output_fastq), log_file=log_name, stream=sys.stderr)
        elif job.output_fastq.stat().st_size <= 0:
            empty.append(job.sample_id)
            log("[EMPTY] {}: {}".format(job.sample_id, job.output_fastq), log_file=log_name, stream=sys.stderr)
        else:
            log("[EXISTS] {}: {}".format(job.sample_id, job.output_fastq), log_file=log_name)

    log(
        "[INFO] check summary: expected={}, missing={}, empty={}".format(len(jobs), len(missing), len(empty)),
        log_file=log_name,
    )
    return 0 if not missing and not empty else 1


def parse_args():
    parser = argparse.ArgumentParser(description="Export non-IG long reads from BAM to FASTQ.GZ.")
    parser.add_argument("step", choices=("run", "check"))
    parser.add_argument("conf_csv")
    parser.add_argument("output_root")
    parser.add_argument("log_name")
    parser.add_argument("n_task", type=int)
    parser.add_argument("nt_per_task", type=int)
    parser.add_argument("--samtools", default=DEFAULT_SAMTOOLS)
    parser.add_argument("--pigz", default=DEFAULT_PIGZ)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.step == "run":
        return run_samples_from_conf(
            conf_csv=args.conf_csv,
            output_root=args.output_root,
            log_name=args.log_name,
            n_task=args.n_task,
            nt_per_task=args.nt_per_task,
            samtools=args.samtools,
            pigz=args.pigz,
            force=args.force,
        )
    if args.step == "check":
        return check_outputs_from_conf(
            conf_csv=args.conf_csv,
            output_root=args.output_root,
            log_name=args.log_name,
            nt_per_task=args.nt_per_task,
            samtools=args.samtools,
            pigz=args.pigz,
        )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
