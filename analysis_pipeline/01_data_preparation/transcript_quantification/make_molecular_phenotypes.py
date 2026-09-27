# -*- coding: utf-8 -*-
"""
@Author  : Chao Xue
@Time    : 2025/10/6 17:36
@Desc    : Description
"""
import logging
import os

from scipy.stats import norm, rankdata


import pandas as pd
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

PROJ_DIR = os.environ.get("PROJECT_ROOT")

def extract_transcripts_from_gtf(gtf_file):
    gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None,
                      names=['#chr','source','feature','start','end','score','strand','frame','attribute'])
    transcripts = gtf[gtf['feature'] == 'transcript'].copy()
    transcripts['transcript_id'] = transcripts['attribute'].str.extract('transcript_id "([^"]+)"')
    transcripts['gene_id'] = transcripts['attribute'].str.extract('gene_id "([^"]+)"')
    df = transcripts.drop_duplicates(subset='transcript_id')
    # gene-level chr, TSS-1, TSS
    gene_df = gtf[gtf['feature'] == 'gene'].copy()
    gene_df['gene_id'] = gene_df['attribute'].str.extract('gene_id "([^"]+)"')
    gene_df['TSS']=gene_df.apply(lambda row: row['start'] if row['strand'] == '+' else row['end'],axis=1)
    gene_tss=dict(zip(gene_df['gene_id'], gene_df['TSS']))
    # assign gene tss to isoform
    df['end']=df['gene_id'].apply(lambda x: gene_tss[x])
    df['start']=df['end']-1
    df = df.set_index('transcript_id')[['#chr','start','end']]
    return df

def preprocess_expression_to_bed(df, feature_pos, missing_threshold=0.25):
    def qqnorm(x):
        n = len(x)
        a = 3.0 / 8.0 if n <= 10 else 0.5
        return norm.ppf((rankdata(x) - a) / (n + 1.0 - 2.0 * a))
    df = df.loc[df.isna().mean(axis=1) <= missing_threshold]
    df = df.apply(lambda x: x.fillna(x.mean()), axis=1)
    df = df.loc[df.std(axis=1) > 0]
    df_qn = pd.DataFrame(
        data=[qqnorm(row.values) for _, row in df.iterrows()],
        index=df.index,
        columns=df.columns
    )
    df_qn=df_qn.round(6)
    bed_df = feature_pos.loc[df_qn.index].copy()
    bed_df = bed_df.assign(ID=df_qn.index)
    bed_df = pd.concat([bed_df, df_qn], axis=1)
    return bed_df


def generate_tissue_tu_bed():
    min_median_tpm=0.1
    min_sample_size=10
    rna_meta_path=f'{PROJ_DIR}/raw_data/GMTiP/meta/RNA/tissue_code.csv'
    main_dir=f'{PROJ_DIR}/project/GMTiP-RNA/20260131'
    tpm_path=f'{main_dir}/output/SRS/quantification/HPC/output/combine/salmon/quant/enhanced.transcript.tpm.salmon.tsv'
    ref_gtf_prefix=f'{main_dir}/output/LRS/isoform_discovery/merged/enhanced_gtf/GTOP_novel-GENCODE_v47'
    tu_dir=f'{main_dir}/release/molecular_phenotype/tu_bed'
    os.makedirs(tu_dir,exist_ok=True)
    logging.info('start')
    gene_map_df=pd.read_csv(f'{ref_gtf_prefix}.gene_transcript_map.txt', sep='\t')
    gene_map=dict(zip(gene_map_df['transcript_id'], gene_map_df['gene_id']))
    logging.info(f'load {len(gene_map)} transcripts from map file.')
    df=pd.read_csv(tpm_path,index_col=0,sep='\t')
    rna_meta_df=pd.read_csv(rna_meta_path,dtype=str)
    code_tissue=dict(zip(rna_meta_df['Tissue_Code'], rna_meta_df['Tissue']))
    ref_gtf=f'{ref_gtf_prefix}.gtf'
    feature_pos=extract_transcripts_from_gtf(ref_gtf)
    logging.info(f'load {len(feature_pos)} transcripts location annotation from gtf')
    tissue_columns={}
    for t in df.columns:
        tc=str(t).split('-')[2]
        tissue=code_tissue[tc]
        if not tissue in tissue_columns:
            tissue_columns[tissue]=[]
        tissue_columns[tissue].append(t)
    for tissue,cols in tissue_columns.items():
        logging.info(f'start {tissue} with {len(cols)} samples')
        out_bed=f'{tu_dir}/{tissue}.tu.bed.gz'
        if os.path.isfile(out_bed):
            continue
        if len(cols)<min_sample_size:
            continue
        sdf = df.loc[:, cols].copy()
        sdf.columns = sdf.columns.map(lambda x: "-".join(x.split("-")[:2]))
        # filter out isoforms with median TPM < 0.1.
        median_expr = sdf.median(axis=1)
        sdf = sdf.loc[median_expr > min_median_tpm]
        # calculate transcript usage.
        expr_cols = sdf.columns
        sdf['gene_id'] = sdf.index.map(gene_map)
        df_gene = sdf.dropna(subset=['gene_id'])
        gene_sum = df_gene.groupby('gene_id')[expr_cols].transform('sum')
        usage = df_gene[expr_cols] / gene_sum[expr_cols]
        usage.index = df_gene.index
        s_bed = preprocess_expression_to_bed(usage, feature_pos, missing_threshold=0.25)
        s_bed['ID']=s_bed['ID'].apply(lambda x:'_'.join([x,gene_map[x]]))
        raw_gene_n=s_bed.shape[0]
        # retain only chr1-22,X
        s_bed=s_bed.loc[s_bed['#chr'].isin([f'chr{i}' for i in range(1,23)]+['chrX']),:]
        s_bed.to_csv(out_bed, index=False, sep='\t')
        logging.info(f'remain {s_bed.shape[0]} features (raw: {sdf.shape[0]}; all chr gene: {raw_gene_n}) in {tissue}')
    logging.info(f'done')

if __name__ == '__main__':
    generate_tissue_tu_bed()