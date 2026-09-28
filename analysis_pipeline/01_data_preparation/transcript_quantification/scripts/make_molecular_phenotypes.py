import logging
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

_WORKER_FEATURE_POS = None
_WORKER_GENE_MAP = None
_WORKER_TU_DIR = None
_WORKER_MIN_MEDIAN_TPM = None
_WORKER_MIN_SAMPLE_SIZE = None
_WORKER_CHR_KEEP = None


def _clean_expression_matrix(df, matrix_name, tissue=None):
    label = f"{tissue} {matrix_name}" if tissue is not None else matrix_name
    df = df.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
    missing_n = int(df.isna().sum().sum())
    if missing_n:
        logging.warning(f"{label} has {missing_n} missing/non-finite values; fill with 0")
        df = df.fillna(0)
    negative_n = int((df < 0).sum().sum())
    if negative_n:
        logging.warning(f"{label} has {negative_n} negative values; set to 0")
        df = df.clip(lower=0)
    return df


def _collapse_duplicate_samples(df, matrix_name, tissue):
    duplicated = df.columns[df.columns.duplicated()].unique()
    if len(duplicated) == 0:
        return df
    logging.warning(
        f"{tissue} {matrix_name} has {len(duplicated)} duplicated sample IDs after renaming; "
        "collapse with mean"
    )
    return df.T.groupby(level=0, sort=False).mean().T


def _init_worker(feature_pos, gene_map, tu_dir, min_median_tpm, min_sample_size, chr_keep):
    global _WORKER_FEATURE_POS, _WORKER_GENE_MAP, _WORKER_TU_DIR
    global _WORKER_MIN_MEDIAN_TPM, _WORKER_MIN_SAMPLE_SIZE, _WORKER_CHR_KEEP
    _WORKER_FEATURE_POS = feature_pos
    _WORKER_GENE_MAP = gene_map
    _WORKER_TU_DIR = tu_dir
    _WORKER_MIN_MEDIAN_TPM = min_median_tpm
    _WORKER_MIN_SAMPLE_SIZE = min_sample_size
    _WORKER_CHR_KEEP = chr_keep


def _normalize_tissue_name(tissue):
    tissue = str(tissue).strip().replace(" ", "_")
    return "_".join(word.capitalize() for word in tissue.split("_") if word)


def _parse_tissue_from_expr_path(expr_path):
    return _normalize_tissue_name(os.path.basename(expr_path).split(".")[0])


def _read_tissue_tpm_file(tpm_path):
    df = pd.read_csv(tpm_path, sep="\t", low_memory=False)
    if "ID" in df.columns:
        id_col_i = df.columns.get_loc("ID")
        expr_df = df.iloc[:, id_col_i + 1 :].copy()
        expr_df.index = df["ID"].astype(str)
        return expr_df
    return pd.read_csv(tpm_path, index_col=0, sep="\t")


def _read_gtf(gtf_file):
    return pd.read_csv(
        gtf_file,
        sep="\t",
        comment="#",
        header=None,
        names=["#chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"],
        dtype={"#chr": str, "source": str, "feature": str, "score": str,
               "strand": str, "frame": str, "attribute": str},
        low_memory=False,
    )


def _extract_attr(series, attr_name):
    return series.astype(str).str.extract(r'%s "([^"]+)"' % attr_name, expand=False)


def _infer_gene_tss_from_transcripts(transcripts):
    gene_rows = []
    for gene_id, gdf in transcripts.dropna(subset=["gene_id"]).groupby("gene_id"):
        strand = gdf["strand"].dropna().iloc[0] if not gdf["strand"].dropna().empty else "+"
        chrom = gdf["#chr"].dropna().iloc[0] if not gdf["#chr"].dropna().empty else None
        tss = gdf["end"].max() if strand == "-" else gdf["start"].min()
        gene_rows.append({"gene_id": gene_id, "#chr": chrom, "start": tss - 1, "end": tss})
    return pd.DataFrame(gene_rows)


def extract_transcripts_from_gtf(gtf_file):
    gtf = _read_gtf(gtf_file)
    transcripts = gtf[gtf["feature"] == "transcript"].copy()
    transcripts["transcript_id"] = _extract_attr(transcripts["attribute"], "transcript_id")
    transcripts["gene_id"] = _extract_attr(transcripts["attribute"], "gene_id")
    df = transcripts.drop_duplicates(subset="transcript_id")

    gene_df = gtf[gtf["feature"] == "gene"].copy()
    if not gene_df.empty:
        gene_df["gene_id"] = _extract_attr(gene_df["attribute"], "gene_id")
        gene_df["TSS"] = gene_df.apply(
            lambda row: row["start"] if row["strand"] == "+" else row["end"], axis=1
        )
        gene_tss = dict(zip(gene_df["gene_id"], gene_df["TSS"]))
    else:
        gene_pos = _infer_gene_tss_from_transcripts(transcripts)
        gene_tss = dict(zip(gene_pos["gene_id"], gene_pos["end"]))

    df = df.loc[df["gene_id"].isin(gene_tss)].copy()
    df["end"] = df["gene_id"].map(gene_tss)
    df["start"] = df["end"] - 1
    return df.set_index("transcript_id")[["#chr", "start", "end"]]


def preprocess_expression_to_bed(df, feature_pos, missing_threshold=0.25):
    def qqnorm(values):
        n = len(values)
        a = 3.0 / 8.0 if n <= 10 else 0.5
        return norm.ppf((rankdata(values) - a) / (n + 1.0 - 2.0 * a))

    df = df.loc[df.isna().mean(axis=1) <= missing_threshold]
    df = df.apply(lambda row: row.fillna(row.mean()), axis=1)
    df = df.loc[df.std(axis=1) > 0]
    df_qn = pd.DataFrame(
        data=[qqnorm(row.values) for _, row in df.iterrows()],
        index=df.index,
        columns=df.columns,
    ).round(6)
    keep_features = df_qn.index.intersection(feature_pos.index)
    df_qn = df_qn.loc[keep_features]
    bed_df = feature_pos.loc[keep_features].copy().assign(ID=keep_features)
    return pd.concat([bed_df, df_qn], axis=1)


def _process_tissue_file_task(task):
    tissue, tpm_path = task
    out_bed = f"{_WORKER_TU_DIR}/{tissue}.tu.bed.gz"
    if os.path.isfile(out_bed):
        return tissue, None

    sdf = _read_tissue_tpm_file(tpm_path)
    if sdf.shape[1] < _WORKER_MIN_SAMPLE_SIZE:
        return tissue, None
    sdf.columns = sdf.columns.map(lambda sample: "-".join(str(sample).split("-")[:2]))
    sdf = _clean_expression_matrix(sdf, "TPM", tissue=tissue)
    sdf = _collapse_duplicate_samples(sdf, "TPM", tissue)
    if sdf.shape[1] < _WORKER_MIN_SAMPLE_SIZE:
        return tissue, None

    sdf = sdf.loc[sdf.median(axis=1) > _WORKER_MIN_MEDIAN_TPM]
    expr_cols = sdf.columns
    sdf["gene_id"] = sdf.index.map(_WORKER_GENE_MAP)
    df_gene = sdf.dropna(subset=["gene_id"])
    gene_sum = df_gene.groupby("gene_id")[expr_cols].transform("sum")
    usage = df_gene[expr_cols] / gene_sum[expr_cols]
    usage.index = df_gene.index
    s_bed = preprocess_expression_to_bed(usage, _WORKER_FEATURE_POS)
    s_bed["ID"] = s_bed["ID"].apply(lambda tx: f"{tx}_{_WORKER_GENE_MAP[tx]}")
    before_chr_filter = s_bed.shape[0]
    s_bed = s_bed.loc[s_bed["#chr"].isin(_WORKER_CHR_KEEP), :]
    s_bed.to_csv(out_bed, index=False, sep="\t")
    return tissue, (s_bed.shape[0], sdf.shape[0], before_chr_filter)


def generate_tissue_tu_bed_by_tissue_based_expr_file(n_thread=33):
    min_median_tpm = 0.1
    min_sample_size = 10
    project_dir = os.environ.get("GTOP_PROJECT_DIR")
    if not project_dir:
        raise EnvironmentError("Set GTOP_PROJECT_DIR before running this script")
    release_dir = f"{project_dir}/release"
    tissue_tpm_dir = f"{project_dir}/output/Salmon_RSEM_Tier3_average_TPM"
    ref_gtf_prefix = f"{release_dir}/gtf/GTOP_novel-GENCODE_v47"
    tu_dir = f"{release_dir}/molec_pheno/tu/averaged"

    os.makedirs(tu_dir, exist_ok=True)
    gene_map_df = pd.read_csv(f"{ref_gtf_prefix}.gene_transcript_map.txt", sep="\t")
    gene_map = dict(zip(gene_map_df["transcript_id"], gene_map_df["gene_id"]))
    feature_pos = extract_transcripts_from_gtf(f"{ref_gtf_prefix}.gtf")
    chr_keep = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

    task_items = []
    for filename in sorted(os.listdir(tissue_tpm_dir)):
        tpm_file = os.path.join(tissue_tpm_dir, filename)
        if filename.startswith(".") or not os.path.isfile(tpm_file):
            continue
        tissue = _parse_tissue_from_expr_path(tpm_file)
        if not os.path.isfile(f"{tu_dir}/{tissue}.tu.bed.gz"):
            task_items.append((tissue, tpm_file))

    if not task_items:
        logging.info("no tissue files require processing")
        return

    max_workers = min(len(task_items), n_thread, os.cpu_count() or n_thread)
    worker_args = (feature_pos, gene_map, tu_dir, min_median_tpm, min_sample_size, chr_keep)
    if "fork" in mp.get_all_start_methods():
        context = mp.get_context("fork")
        _init_worker(*worker_args)
        executor_kwargs = {"max_workers": max_workers, "mp_context": context}
    else:
        context = mp.get_context()
        executor_kwargs = {
            "max_workers": max_workers,
            "mp_context": context,
            "initializer": _init_worker,
            "initargs": worker_args,
        }

    with ProcessPoolExecutor(**executor_kwargs) as executor:
        for tissue, result in executor.map(_process_tissue_file_task, task_items):
            if result is None:
                continue
            remain_n, raw_n, before_chr_filter = result
            logging.info(
                f"{tissue}: retained {remain_n} features "
                f"(expression-filtered={raw_n}; before chromosome filter={before_chr_filter})"
            )


if __name__ == "__main__":
    generate_tissue_tu_bed_by_tissue_based_expr_file()
