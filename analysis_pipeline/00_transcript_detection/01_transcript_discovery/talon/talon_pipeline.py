# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2026/09/10
@Email   : xuechao@szbl.ac.cn
@Desc    : Execute TALON transcript discovery tasks.
"""

import glob
import os
import shlex
import shutil
import sys
from pathlib import Path

import pandas as pd

from datetime import datetime

def log(msg, log_file=None):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line)
    if log_file:
        with open(log_file, 'a') as f:
            f.write(line + '\n')

import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_commands_threadpool(cmd_list, main_log, max_workers=4, main_log_prefix='batch'):
    if main_log:
        Path(main_log).parent.mkdir(parents=True, exist_ok=True)
    def worker(commands):
        for item in commands:
            log_path = Path(item['log_path']).resolve()
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with log_path.open('w') as handle:
                subprocess.run(item['cmd'], shell=True, executable='/bin/bash',
                               check=True, stdout=handle, stderr=subprocess.STDOUT,
                               cwd=log_path.parent)
    failures = []
    with ThreadPoolExecutor(max_workers=max(1, int(max_workers))) as pool:
        futures = [pool.submit(worker, commands) for commands in cmd_list]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as error:
                failures.append(str(error))
    if main_log:
        with open(main_log, 'a') as handle:
            handle.write(f'{main_log_prefix}: tasks={len(cmd_list)}, failures={len(failures)}\n')
            handle.writelines(error + '\n' for error in failures)
    if failures:
        raise RuntimeError(f'{len(failures)} task(s) failed; see {main_log}')

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

ARGS = sys.argv[1:]
STEP = ARGS[0] if ARGS else None
if STEP == 'tx_abundance_one':
    CONF_CSV = None
    LOG_NAME = None
    N_TASK = None
    NT_PER_TASK = None
else:
    CONF_CSV = ARGS[1]
    LOG_NAME = ARGS[2]
    N_TASK = int(ARGS[3])
    NT_PER_TASK = int(ARGS[4])

REQUIRED_OUTPUT_NAMES = [
    'talon.db',
    'TALON_AND_GTF_SUCCESS.done',
    'tx_abundance/read_cov.tsv.gz',
]
REQUIRED_OUTPUT_GLOBS = [
    'talon_observed*.gtf',
]
TRANSFER_OUTPUT_NAMES = [
    'tx_abundance/read_cov.tsv.gz',
    'run_talon.log',
    'TALON_AND_GTF_SUCCESS.done',
]
TRANSFER_OUTPUT_GLOBS = [
    'talon_observed*.gtf',
]

IG_REGIONS = {
    'chr14': [(103000000, 107500000)],
    'chr2': [(88000000, 90000000)],
    'chr22': [(22000000, 23500000)],
}


def _sample_label(row):
    sample_id = row.get('sample_id', row.get('mapped_bam_path', 'NA'))
    ref_genome_tag = row.get('ref_genome_tag', '')
    return f'{sample_id}:{ref_genome_tag}' if ref_genome_tag else sample_id


def _is_nonempty_file(file_path):
    return os.path.exists(file_path) and os.path.getsize(file_path) > 0


def _quote_cmd(parts):
    return ' '.join(shlex.quote(str(x)) for x in parts)


def _shell_quote(value):
    return shlex.quote(str(value))


def _sql_quote(value):
    return "'" + str(value).replace("'", "''") + "'"


def _write_shell_script(script_path, script_text):
    os.makedirs(os.path.dirname(script_path), exist_ok=True)
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(script_text.strip() + '\n')


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


def tx_abundance_one(tx_gtf, input_path, output_path):
    out_dir = os.path.dirname(output_path)
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    df = pd.read_csv(input_path, sep='\t')
    if 'annot_transcript_id' not in df.columns:
        raise Exception(f'annot_transcript_id column is required: {input_path}')
    count_series = df['annot_transcript_id'].value_counts()
    cov_df = pd.DataFrame({
        'transcript_id': count_series.index,
        'read_count': count_series.values,
    })
    ig_transcript_ids = get_ig_transcript_ids(tx_gtf)
    cov_df['in_IG_region'] = cov_df['transcript_id'].isin(ig_transcript_ids).astype('int8')
    cov_df.to_csv(output_path, sep='\t', index=False, encoding='utf-8')


def _deduplicate_prepare_jobs(df):
    required_cols = ['ref_gtf_path', 'clean_db_path', 'genome_build', 'annot_name']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise Exception(f'{required_cols} columns are required for prepare step; missing={missing_cols}')
    return df.drop_duplicates(subset=['clean_db_path']).copy()


def _check_prepare_sample(clean_db_path):
    missing_files = []
    if not _is_nonempty_file(clean_db_path):
        missing_files.append(clean_db_path)
    success_file = f'{clean_db_path}.prepare_success.done'
    if not _is_nonempty_file(success_file):
        missing_files.append(success_file)
    return len(missing_files) == 0, missing_files


def _build_integrity_check_cmd(db_path, required=True):
    db_q = _shell_quote(db_path)
    missing_cmd = (
        f'echo "ERROR: DB missing or empty: {db_q}" && exit 1'
        if required else
        f'echo "DB missing or empty: {db_q}" && exit 1'
    )
    return (
        f'if [ ! -s {db_q} ]; then {missing_cmd}; fi && '
        f'if command -v sqlite3 >/dev/null 2>&1; then '
        f'result=$(sqlite3 {db_q} "PRAGMA integrity_check;" || true); '
        f'if [ "$result" != "ok" ]; then '
        f'echo "ERROR: SQLite integrity_check failed for {db_q}: $result"; exit 1; '
        f'fi; '
        f'else echo "WARNING: sqlite3 not found; only non-empty DB check was done."; fi'
    )


def _build_annotation_check_cmd(db_path, annot_name):
    db_q = _shell_quote(db_path)
    annot_sql = _sql_quote(annot_name)
    return (
        f'if command -v sqlite3 >/dev/null 2>&1; then '
        f'table_n=$(sqlite3 {db_q} "SELECT COUNT(*) FROM sqlite_master WHERE type = \'table\' AND name = \'gene_annotations\';"); '
        f'if [ "$table_n" != "1" ]; then '
        f'echo "ERROR: TALON DB missing gene_annotations table: {db_q}"; exit 1; '
        f'fi; '
        f'annot_n=$(sqlite3 {db_q} "SELECT COUNT(*) FROM gene_annotations WHERE annot_name = {annot_sql};"); '
        f'if [ "$annot_n" = "0" ]; then '
        f'echo "ERROR: TALON annotation name not found after initialization: {annot_name}"; '
        f'echo "Available annotation names:"; '
        f'sqlite3 {db_q} "SELECT DISTINCT annot_name FROM gene_annotations ORDER BY annot_name;"; '
        f'exit 1; '
        f'fi; '
        f'fi'
    )


def prepare(conf_path, n_task, log_name):
    parallel_cmds = []
    df = _deduplicate_prepare_jobs(pd.read_csv(conf_path))
    for _, row in df.iterrows():
        ref_gtf_path = row['ref_gtf_path']
        clean_db_path = row['clean_db_path']
        genome_build = row['genome_build']
        annot_name = row['annot_name']

        out_dir = os.path.dirname(clean_db_path)
        os.makedirs(out_dir, exist_ok=True)
        db_prefix = clean_db_path[:-3] if clean_db_path.endswith('.db') else clean_db_path
        prepare_prefix = os.path.basename(db_prefix)
        log_path = f'{out_dir}/{prepare_prefix}.prepare_talon_database.log'
        script_path = f'{out_dir}/{prepare_prefix}.prepare_talon_database.sh'
        talon_gtf_path = f'{out_dir}/{prepare_prefix}.talon_init.gtf'
        success_file = f'{clean_db_path}.prepare_success.done'
        script = f'''
            #!/bin/bash
            set -euo pipefail
            test -s {_shell_quote(ref_gtf_path)}
            mkdir -p {_shell_quote(out_dir)}

            if [ -s {_shell_quote(clean_db_path)} ] && [ -s {_shell_quote(success_file)} ]; then
                if command -v sqlite3 >/dev/null 2>&1; then
                    result=$(sqlite3 {_shell_quote(clean_db_path)} "PRAGMA integrity_check;" || true)
                    if [ "$result" = "ok" ]; then
                        annot_n=$(sqlite3 {_shell_quote(clean_db_path)} "SELECT COUNT(*) FROM gene_annotations WHERE annot_name = {_sql_quote(annot_name)};" || echo 0)
                        if [ "$annot_n" != "0" ]; then
                            echo "Existing TALON DB is valid. Skip initialization: {_shell_quote(clean_db_path)}"
                            exit 0
                        fi
                        echo "Existing TALON DB has no annotation {_shell_quote(annot_name)}; rebuilding."
                    fi
                    echo "Existing TALON DB failed integrity_check: $result"
                else
                    echo "sqlite3 not found; reuse non-empty TALON DB: {_shell_quote(clean_db_path)}"
                    exit 0
                fi
            fi

            if [ -e {_shell_quote(clean_db_path)} ]; then
                bad_db={_shell_quote(clean_db_path)}.bad_$(date +%Y%m%d%H%M%S)
                mv {_shell_quote(clean_db_path)} "$bad_db"
            fi
            rm -f {_shell_quote(success_file)} {_shell_quote(clean_db_path)} {_shell_quote(clean_db_path)}-wal {_shell_quote(clean_db_path)}-shm {_shell_quote(clean_db_path)}-journal

            python {_shell_quote(f'{CURRENT_DIR}/prepare_talon_gtf.py')} \\
                {_shell_quote(ref_gtf_path)} \\
                {_shell_quote(talon_gtf_path)}
            test -s {_shell_quote(talon_gtf_path)}

            talon_initialize_database \\
                --f {_shell_quote(talon_gtf_path)} \\
                --g {_shell_quote(genome_build)} \\
                --a {_shell_quote(annot_name)} \\
                --o {_shell_quote(db_prefix)}
            {_build_integrity_check_cmd(clean_db_path)}
            {_build_annotation_check_cmd(clean_db_path, annot_name)}
            {{
                echo "TALON clean database initialized successfully."
                echo "Date: $(date)"
                echo "GTF: {ref_gtf_path}"
                echo "TALON init GTF: {talon_gtf_path}"
                echo "Genome build: {genome_build}"
                echo "Annotation: {annot_name}"
                echo "DB: {clean_db_path}"
            }} > {_shell_quote(success_file)}
        '''
        _write_shell_script(script_path, script)
        cmd = _quote_cmd(['bash', script_path])
        parallel_cmds.append([{'cmd': cmd, 'log_path': log_path}])

    print(f'prepare_n={len(parallel_cmds)}')
    if not parallel_cmds:
        return
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task)


def prepare_check(conf_path):
    df = _deduplicate_prepare_jobs(pd.read_csv(conf_path))
    fail_ids = []
    succ_n = 0
    for _, row in df.iterrows():
        sample_label = _sample_label(row)
        clean_db_path = row['clean_db_path']
        is_success, missing_files = _check_prepare_sample(clean_db_path)
        if is_success:
            succ_n += 1
        else:
            fail_ids.append(sample_label)
            print(f'[WARN] prepare sample={sample_label} missing files: {missing_files}')
    print(f'prepare_succ_n={succ_n}; prepare_fail_n={len(fail_ids)}')
    print(f'prepare_failed sample_hap: {fail_ids}')


def talon(conf_path, n_task, log_name, nt_per_task):
    parallel_cmds = []
    df = pd.read_csv(conf_path)
    for _, row in df.iterrows():
        sample_id = row['sample_id']
        mapped_bam_path = row['mapped_bam_path']
        out_dir = row['out_dir']
        ref_genome_path = row['ref_genome_path']
        clean_db_path = row['clean_db_path']
        genome_build = row['genome_build']
        annot_name = row['annot_name']
        sample_desc = row.get('ref_genome_tag', sample_id)
        platform = row.get('platform', 'PacBio')

        os.makedirs(out_dir, exist_ok=True)
        log_path = f'{out_dir}/run_talon.log'
        raw_sam = f'{out_dir}/{sample_id}.flnc_minimap2.sam'
        labeled_prefix = f'{out_dir}/talon'
        labeled_sam = f'{labeled_prefix}_labeled.sam'
        sample_config = f'{out_dir}/talon_config.csv'
        sample_db = f'{out_dir}/talon.db'
        talon_out_prefix = f'{out_dir}/talon'
        talon_read_annot = f'{talon_out_prefix}_talon_read_annot.tsv'
        gtf_out_prefix = f'{out_dir}/talon_observed'
        tx_abundance = f'{out_dir}/tx_abundance/read_cov.tsv.gz'
        success_file = f'{out_dir}/TALON_AND_GTF_SUCCESS.done'
        script_path = f'{out_dir}/run_talon.sh'

        script = f'''
            #!/bin/bash
            set -euo pipefail
            mkdir -p {_shell_quote(out_dir)}
            short_tmp_root="${{TALON_TMP_ROOT:-${{TMPDIR:-/tmp}}}}"
            tmp_dir="${{short_tmp_root}}/talon_${{SLURM_JOB_ID:-$$}}_${{SLURM_PROCID:-0}}"
            mkdir -p "$tmp_dir"
            export TMPDIR="$tmp_dir"
            export SQLITE_TMPDIR="$tmp_dir"

            cleanup_on_exit() {{
                status=$?
                if [ "$status" -eq 0 ]; then
                    rm -rf "$tmp_dir"
                else
                    echo "ERROR: TALON run failed with status $status"
                    echo "Temporary directory kept for debugging: $tmp_dir"
                fi
                exit "$status"
            }}
            trap cleanup_on_exit EXIT

            if [ -s {_shell_quote(success_file)} ] &&
               [ -s {_shell_quote(sample_db)} ] &&
               [ -s {_shell_quote(tx_abundance)} ] &&
               ls {_shell_quote(gtf_out_prefix)}*.gtf >/dev/null 2>&1; then
                echo "Skip completed TALON sample: {_shell_quote(sample_id)}"
                {_build_integrity_check_cmd(sample_db)}
                exit 0
            fi

            test -s {_shell_quote(mapped_bam_path)}
            test -s {_shell_quote(ref_genome_path)}
            test -s {_shell_quote(clean_db_path)}
            {_build_integrity_check_cmd(clean_db_path)}

            if ! samtools quickcheck -v {_shell_quote(mapped_bam_path)}; then
                echo "ERROR: input BAM failed samtools quickcheck: {_shell_quote(mapped_bam_path)}"
                exit 1
            fi

            rm -f {_shell_quote(sample_db)} {_shell_quote(sample_db)}-wal {_shell_quote(sample_db)}-shm {_shell_quote(sample_db)}-journal
            cp -av {_shell_quote(clean_db_path)} {_shell_quote(sample_db)}
            {_build_integrity_check_cmd(sample_db)}

            ref_fai={_shell_quote(ref_genome_path)}.fai
            if [ ! -s "$ref_fai" ]; then
                samtools faidx {_shell_quote(ref_genome_path)}
            fi
            echo "Using TALON SAM cleanup version: lossless_sam_filter_v6"

            header_refs="$tmp_dir/header_refs.txt"
            original_header="$tmp_dir/original_header.sam"
            filtered_sam="$tmp_dir/filtered.sam"
            filter_stats="$tmp_dir/filter_stats.txt"
            filtered_bam="$tmp_dir/filtered.bam"
            md_bam="$tmp_dir/md.bam"
            validate_bam="$tmp_dir/validate_talon_input.bam"

            cut -f1 "$ref_fai" > "$header_refs"
            test -s "$header_refs"
            samtools view -H {_shell_quote(mapped_bam_path)} > "$original_header"
            test -s "$original_header"

            # Build a canonical SAM header from the target FASTA. Some BAMs have
            # reference names in their binary dictionary that are absent from
            # the textual @SQ header; retaining that header makes the SAM invalid.
            # Rebuilding @SQ from the FAI guarantees that every retained RNAME is
            # declared, while preserving @HD, @RG, @PG and @CO metadata.
            {{
                awk '/^@HD/{{print}}' "$original_header"
                awk 'BEGIN{{OFS="\\t"}} {{print "@SQ", "SN:"$1, "LN:"$2}}' "$ref_fai"
                awk '!/^@HD/ && !/^@SQ/{{print}}' "$original_header"
                samtools view -@ {int(nt_per_task)} {_shell_quote(mapped_bam_path)} | \\
                    awk -F '\\t' -v refs="$header_refs" -v stats="$filter_stats" 'BEGIN{{while((getline<refs)>0) valid[$1]=1}} {{total++}} !($3 in valid){{dropped_rname++; invalid_rname[$3]++; next}} ($7!="=" && $7!="*" && !($7 in valid)){{dropped_rnext++; invalid_rnext[$7]++; next}} {{print; kept++}} END{{printf "input_alignment_n=%d; kept_alignment_n=%d; dropped_invalid_rname_n=%d; dropped_invalid_rnext_n=%d\\n", total, kept, dropped_rname, dropped_rnext > stats; for(ref in invalid_rname) printf "RNAME\\t%s\\t%d\\n", ref, invalid_rname[ref] >> stats; for(ref in invalid_rnext) printf "RNEXT\\t%s\\t%d\\n", ref, invalid_rnext[ref] >> stats}}'
            }} > "$filtered_sam"
            test -s "$filtered_sam"
            test -s "$filter_stats"
            head -n 11 "$filter_stats"

            kept_alignment_n=$(awk -F'[=;]' 'NR==1{{print $4}}' "$filter_stats")
            if [ "${{kept_alignment_n:-0}}" -eq 0 ]; then
                echo "ERROR: no alignments remain after filtering against FASTA references"
                exit 1
            fi

            invalid_ref_n=$(awk -v refs="$header_refs" 'BEGIN{{while((getline<refs)>0) valid[$1]=1}} !/^@/ && !($3 in valid){{if(n<10) print $3 > "/dev/stderr"; n++}} END{{print n+0}}' "$filtered_sam")
            if [ "$invalid_ref_n" != "0" ]; then
                echo "ERROR: filtered TALON SAM still contains invalid RNAME records: $invalid_ref_n"
                exit 1
            fi
            if ! samtools view -@ {int(nt_per_task)} -bS "$filtered_sam" > "$validate_bam"; then
                echo "ERROR: filtered SAM validation failed; rerunning single-threaded parser for an exact line number"
                samtools view -S "$filtered_sam" >/dev/null || true
                exit 1
            fi

            set +o pipefail
            has_md=$(awk '!/^@/{{for(i=12;i<=NF;i++) if($i~/^MD:Z:/){{print "yes"; exit}} seen++; if(seen>1000) exit}}' "$filtered_sam")
            set -o pipefail

            rm -f {_shell_quote(raw_sam)}
            if [ "$has_md" = "yes" ]; then
                mv "$filtered_sam" {_shell_quote(raw_sam)}
            else
                mv "$validate_bam" "$filtered_bam"
                samtools calmd -@ {int(nt_per_task)} -b "$filtered_bam" {_shell_quote(ref_genome_path)} > "$md_bam"
                samtools view -@ {int(nt_per_task)} -h "$md_bam" > {_shell_quote(raw_sam)}
            fi
            test -s {_shell_quote(raw_sam)}
            samtools view -@ {int(nt_per_task)} -bS {_shell_quote(raw_sam)} > "$validate_bam"
            rm -f "$validate_bam"

            cd {_shell_quote(out_dir)}
            rm -f {_shell_quote(labeled_sam)} {_shell_quote(labeled_prefix + '_read_labels.tsv')}
            talon_label_reads \\
                --f {_shell_quote(raw_sam)} \\
                --g {_shell_quote(ref_genome_path)} \\
                --t {int(nt_per_task)} \\
                --o {_shell_quote(labeled_prefix)}
            rm -rf tmp_label_reads
            rm -f {_shell_quote(raw_sam)}

            printf '%s,%s,%s,%s\\n' \\
                {_shell_quote(sample_id)} \\
                {_shell_quote(sample_desc)} \\
                {_shell_quote(platform)} \\
                {_shell_quote(labeled_sam)} > {_shell_quote(sample_config)}

            talon \\
                --f {_shell_quote(sample_config)} \\
                --db {_shell_quote(sample_db)} \\
                --build {_shell_quote(genome_build)} \\
                --threads {int(nt_per_task)} \\
                --tmpDir "$tmp_dir/" \\
                --o {_shell_quote(talon_out_prefix)}
            rm -f {_shell_quote(labeled_sam)}

            {_build_integrity_check_cmd(sample_db)}
            rm -f {_shell_quote(gtf_out_prefix)}*.gtf
            talon_create_GTF \\
                --db {_shell_quote(sample_db)} \\
                --build {_shell_quote(genome_build)} \\
                --annot {_shell_quote(annot_name)} \\
                --observed \\
                --o {_shell_quote(gtf_out_prefix)}

            ls {_shell_quote(gtf_out_prefix)}*.gtf >/dev/null 2>&1
            test -s {_shell_quote(talon_read_annot)}
            tx_gtf=$(ls {_shell_quote(gtf_out_prefix)}*.gtf | head -n 1)
            python {_shell_quote(os.path.abspath(__file__))} tx_abundance_one \\
                "$tx_gtf" \\
                {_shell_quote(talon_read_annot)} \\
                {_shell_quote(tx_abundance)}
            test -s {_shell_quote(tx_abundance)}

            {{
                echo "TALON and GTF generation completed successfully."
                echo "Date: $(date)"
                echo "Sample: {sample_id}"
                echo "BAM: {mapped_bam_path}"
                echo "Sample DB: {sample_db}"
                echo "Transcript abundance: {tx_abundance}"
            }} > {_shell_quote(success_file)}
        '''
        _write_shell_script(script_path, script)
        cmd = _quote_cmd(['bash', script_path])
        parallel_cmds.append([{'cmd': cmd, 'log_path': log_path}])
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task)


def _has_required_glob(out_dir, pattern):
    return any(_is_nonempty_file(path) for path in glob.glob(f'{out_dir}/{pattern}'))


def _build_required_files(out_dir):
    return [f'{out_dir}/{name}' for name in REQUIRED_OUTPUT_NAMES]


def _build_transfer_files(out_dir):
    files = []
    for name in TRANSFER_OUTPUT_NAMES:
        file_path = f'{out_dir}/{name}'
        if os.path.exists(file_path):
            files.append(file_path)
    for pattern in TRANSFER_OUTPUT_GLOBS:
        files.extend(glob.glob(f'{out_dir}/{pattern}'))
    return sorted(set(files))


def _build_final_required_files(final_out_dir):
    return [f'{final_out_dir}/{name}' for name in TRANSFER_OUTPUT_NAMES]


def _check_sample(out_dir):
    missing_files = []
    for file_path in _build_required_files(out_dir):
        if not _is_nonempty_file(file_path):
            missing_files.append(file_path)
    for pattern in REQUIRED_OUTPUT_GLOBS:
        if not _has_required_glob(out_dir, pattern):
            missing_files.append(f'{out_dir}/{pattern}')
    return len(missing_files) == 0, missing_files


def _check_final_transfer(final_out_dir):
    missing_files = []
    for file_path in _build_final_required_files(final_out_dir):
        if not _is_nonempty_file(file_path):
            missing_files.append(file_path)
    for pattern in TRANSFER_OUTPUT_GLOBS:
        if not _has_required_glob(final_out_dir, pattern):
            missing_files.append(f'{final_out_dir}/{pattern}')
    return len(missing_files) == 0, missing_files


def _check_transfer_source(out_dir):
    return _check_final_transfer(out_dir)


def _is_relative_to(child_path, parent_path):
    try:
        child_path.relative_to(parent_path)
        return True
    except ValueError:
        return False


def _validate_transfer_paths(out_dir, final_out_dir, run_out_main_dir=None):
    out_path = Path(out_dir).resolve()
    final_path = Path(final_out_dir).resolve()
    if out_path == final_path:
        raise Exception(f'out_dir and final_out_dir are the same: {out_path}')
    if _is_relative_to(final_path, out_path):
        raise Exception(f'final_out_dir is inside out_dir; refuse to delete: {final_path}')
    if run_out_main_dir is not None:
        run_out_path = Path(run_out_main_dir).resolve()
        if not _is_relative_to(out_path, run_out_path):
            raise Exception(f'out_dir is not inside run_out_main_dir; refuse to delete: {out_path}')


def _transfer_one_row(row, main_log=None):
    sample_label = _sample_label(row)
    out_dir = row['out_dir']
    final_out_dir = row['final_out_dir']
    if os.path.realpath(out_dir) == os.path.realpath(final_out_dir):
        return sample_label, True
    run_out_main_dir = row.get('run_out_main_dir')
    is_success, missing_files = _check_sample(out_dir)
    if not is_success:
        log(f'[WARN] Skip transfer sample={sample_label}; missing files: {missing_files}', main_log)
        return sample_label, False

    is_transfer_ready, transfer_missing_files = _check_transfer_source(out_dir)
    if not is_transfer_ready:
        log(
            f'[WARN] Skip transfer sample={sample_label}; '
            f'missing transfer files: {transfer_missing_files}',
            main_log,
        )
        return sample_label, False

    _validate_transfer_paths(out_dir, final_out_dir, run_out_main_dir)
    os.makedirs(final_out_dir, exist_ok=True)
    copied_files = []
    for source_file in _build_transfer_files(out_dir):
        rel_path = os.path.relpath(source_file, out_dir)
        target_file = f'{final_out_dir}/{rel_path}'
        os.makedirs(os.path.dirname(target_file), exist_ok=True)
        shutil.copy2(source_file, target_file)
        copied_files.append(target_file)

    removed_files = []
    is_transferred, final_missing_files = _check_final_transfer(final_out_dir)
    if not is_transferred:
        log(
            f'[WARN] Skip remove sample={sample_label}; '
            f'final missing files after copy: {final_missing_files}',
            main_log,
        )
        return sample_label, False

    # Keep the working output and logs for reproducibility.
    log(
        f'[INFO] Transferred sample={sample_label}; copied_files={len(copied_files)}; '
        f'pruned_final_files={len(removed_files)}; retained_work_dir={out_dir}',
        main_log,
    )
    return sample_label, True


def transfer_one(conf_path, row_idx, log_name):
    df = pd.read_csv(conf_path)
    row = df.iloc[int(row_idx)]
    _, is_transferred = _transfer_one_row(row, log_name)
    if not is_transferred:
        sys.exit(1)


def transfer(conf_path, n_task=1, log_name=None):
    df = pd.read_csv(conf_path)
    if 'final_out_dir' not in df.columns:
        raise Exception('final_out_dir column is required for transfer step')

    n_task = max(1, int(n_task))
    parallel_cmds = []
    for i, row in df.iterrows():
        out_dir = row['out_dir']
        sample_id = row.get('sample_id', os.path.basename(out_dir))
        cmd = _quote_cmd([
            sys.executable,
            os.path.abspath(__file__),
            'transfer_one',
            conf_path,
            log_name,
            1,
            1,
            i,
        ])
        log_path = f'{out_dir}/{sample_id}.run_transfer.log'
        parallel_cmds.append([{'cmd': cmd, 'log_path': log_path}])
    run_commands_threadpool(parallel_cmds, main_log=log_name, max_workers=n_task, main_log_prefix='transfer')


if __name__ == '__main__':
    if STEP == 'tx_abundance_one':
        tx_abundance_one(tx_gtf=ARGS[1], input_path=ARGS[2], output_path=ARGS[3])
    if STEP == 'prepare':
        prepare(conf_path=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME)
    if STEP == 'prepare_check':
        prepare_check(CONF_CSV)
    if STEP == 'run':
        talon(conf_path=CONF_CSV, n_task=N_TASK, log_name=LOG_NAME, nt_per_task=NT_PER_TASK)
    if STEP == 'transfer':
        transfer(CONF_CSV, n_task=N_TASK, log_name=LOG_NAME)
    if STEP == 'transfer_one':
        transfer_one(CONF_CSV, row_idx=ARGS[5], log_name=LOG_NAME)
