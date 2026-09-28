# -*- coding: utf-8 -*-

import os
from pathlib import Path

CONFIG_PROJECT_DIR = os.environ.get('GTOP_PROJECT_DIR', str(Path.cwd() / 'gtop_run'))
CONFIG_REF_DIR = os.environ.get('GTOP_REF_DIR', str(Path(CONFIG_PROJECT_DIR) / 'reference'))
CONFIG_TISSUE_META = os.environ.get('TISSUE_META', str(Path(CONFIG_PROJECT_DIR) / 'input/tissue_code.csv'))
CONFIG_GENCODE_PROTEINS = os.environ.get('GENCODE_PROTEINS', str(Path(CONFIG_REF_DIR) / 'gencode.v47.pc_translations.fa.gz'))
CONFIG_MS_RAW_DIR = os.environ.get('MS_RAW_DIR', str(Path(CONFIG_PROJECT_DIR) / 'input/ms_raw'))

import gzip
import logging
import re
import shutil

import numpy as np
import pandas as pd


PROJ_DIR = os.environ.get("PROJECT_ROOT")

Tissue_meta_path = CONFIG_TISSUE_META

Low_than_1k_sample_IDs = ['GTOP-BF221-0378-LN-2ZKY', 'GTOP-BG171-0378-LN-TDV5', 'GTOP-CF242-0378-LN-F6UE', 'GTOP-BL021-0378-LN-RZ42', 'GTOP-AJ221-1053-LN-DP1V', 'GTOP-AJ221-4019-LN-D5AL', 'GTOP-CF241-1236-LN-NJ1W', 'GTOP-BA131-2099-LN-H19X', 'GTOP-CF242-2099-LN-6Y37', 'GTOP-CA081-1392-LN-20Q9', 'GTOP-CB161-1392-LN-I067', 'GTOP-CB271-1584-LN-N1C6', 'GTOP-CG041-5156-LN-75CF']

class GeneAnnotation:
    def __init__(self):
        self.df=None
        pass

    def load_df(self, gene_annot_prefix, remove_ver=True):
        path=f'{gene_annot_prefix}.gene_transcript_map.txt'
        df=pd.read_csv(path,sep='\t',dtype=str)
        if remove_ver:
            for gid in ['gene_id','transcript_id']:
                df[gid]=df[gid].apply(lambda x:x.split('.')[0] if x.startswith('ENS') else x)
        self.df=df
        return self
        pass
    def get_gene_id_name_map(self):
        return dict(zip(self.df['gene_id'],self.df['gene_name']))

    def get_gene_isoforms_map(self,gene_key='gene_id'):
        gene_isoform_map = self.df.groupby(gene_key)['transcript_id'].apply(list).to_dict()
        return gene_isoform_map

    def get_isoform_gene_map(self,gene_key='gene_id',no_version=False):
        self.df['transcript_id_no_version']=self.df['transcript_id'].map(lambda x:x.split('.')[0] if x.startswith('ENS') else x)
        if no_version:
            gene_isoform_map = dict(zip(self.df['transcript_id_no_version'], self.df[gene_key]))
        else:
            gene_isoform_map = dict(zip(self.df['transcript_id'], self.df[gene_key]))
        return gene_isoform_map

    def get_isoform_pos(self):
        self.df['pos_str']=self.df['chr']+':'+self.df['start']+'-'+self.df['end']
        return dict(zip(self.df['transcript_id'],self.df['pos_str']))

    def get_isoform_type_map(self):
        return dict(zip(self.df['transcript_id'],self.df['transcript_type']))

    def get_gene_type_map(self):
        return dict(zip(self.df['gene_id'],self.df['gene_type']))

class GTOP_Tissue_Meta:
    def __init__(self):
        self.df = None
        self.tissue_name_meta_map = None
        self.tissue_code_meta_map = None
        self.tissue_meta_path=Tissue_meta_path
        self.load_meta()

    def load_meta(self):
        meta_df = pd.read_csv(self.tissue_meta_path, dtype=str)
        self.df=meta_df
        tissue_code_meta_map={}
        for k in meta_df.columns:
            tissue_code_meta_map[k]=dict(zip(meta_df['Tissue_Code'],meta_df[k]))
        tissue_name_meta_map={}
        for k in meta_df.columns:
            tissue_name_meta_map[k]=dict(zip(meta_df['Tissue'],meta_df[k]))
        self.tissue_code_meta_map=tissue_code_meta_map
        self.tissue_name_meta_map=tissue_name_meta_map

    def get_meta_by_tissue_code(self,tissue_code,key='Tissue'):
        return self.tissue_code_meta_map[key][tissue_code]

    def get_meta_by_tissue_name(self,tissue_name,key='Tissue_Color_Code'):
        return self.tissue_name_meta_map[key][tissue_name]

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

result_version='20260815'

def __write_faa(faa_path,id_set,out_f,tid_gene=None,transcript_index=0,gene_index=1):
    if faa_path.endswith('.gz'):
        open_func = gzip.open
        open_kwargs = {'mode': 'rt', 'encoding': 'utf-8'}
    else:
        open_func = open
        open_kwargs = {'mode': 'r', 'encoding': 'utf-8'}
    k=0
    current_id=None
    current_seq = []
    with open_func(faa_path, **open_kwargs) as in_f:
        for line in in_f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id in id_set:
                    gene_name=current_gene_id
                    if tid_gene:
                        gene_name = tid_gene[current_id]
                    out_f.write(f'>{current_id} {gene_name}\n')
                    out_f.write(''.join(current_seq) + '\n')
                    k+=1
                header = line.lstrip('>')
                parse_arr=[part.strip() for part in re.split(r'\s+|\|', header) if part]
                current_id = parse_arr[transcript_index]
                current_gene_id = parse_arr[gene_index] if len(parse_arr) > gene_index else ''
                current_seq = []
            else:
                current_seq.append(line.strip('*').strip())
        if current_id in id_set:
            gene_name = current_gene_id
            if tid_gene:
                gene_name = tid_gene[current_id]
            out_f.write(f'>{current_id} {gene_name}\n')
            out_f.write(''.join(current_seq) + '\n')
            k+=1
    return k


def make_enhanced_cds_fa():
    gencode_v47_cds_protein_fa_path=CONFIG_GENCODE_PROTEINS
    main_dir=f'{CONFIG_PROJECT_DIR}/release'
    gtop_cds_protein_fa_path=f'{main_dir}/LRS_assembly/sqanti3/filter.clean.faa'
    enhanced_gtf_prefix=f'{main_dir}/gtf/GTOP_novel-GENCODE_v47'
    output_faa=f'{main_dir}/cds/global/GTOP_novel-GENCODE_v47.faa'

    os.makedirs(os.path.dirname(output_faa),exist_ok=True)
    ga=GeneAnnotation().load_df(enhanced_gtf_prefix,remove_ver=False)
    tid_gene=ga.get_isoform_gene_map(gene_key='gene_id')
    tid_type=ga.get_isoform_type_map()
    print(f'load {len(tid_gene)} isoforms; novel {len([x for x in tid_gene.keys() if not x.startswith("ENS")])}.')
    id_set=set([t for t in tid_gene.keys() if tid_type[t] == 'protein_coding'])
    all_k=0
    with open(output_faa, 'w', encoding='utf-8') as out_f:
        # write gencode isoform
        logger.info(f'start merge gencode')
        k=__write_faa(gencode_v47_cds_protein_fa_path, id_set, out_f, tid_gene, transcript_index=1)
        logger.info(f'merge gencode {k}')
        all_k+=k
        # write novel isoform
        k=__write_faa(gtop_cds_protein_fa_path,id_set,out_f,tid_gene,transcript_index=0)
        logger.info(f'merge novel {k}')
        all_k+=k
        logger.info(f'done {all_k} proteins; save to {output_faa}')


def tissue_based_protein_cds_fa(min_tpm=1,min_sample=None):
    # if min_sample is None ,use median to determine.
    main_dir=f'{CONFIG_PROJECT_DIR}/release'
    isoform_expr_path = f'{main_dir}/LRS_quant/GTOP/GTOP.transcript.tpm.flair.tsv'
    enhanced_faa_path=f'{main_dir}/cds/global/GTOP_novel-GENCODE_v47.faa'
    faa_dir=f'{main_dir}/cds/tissue_based'

    os.makedirs(f'{faa_dir}',exist_ok=True)
    df=pd.read_csv(f'{isoform_expr_path}',index_col=0,sep='\t')
    # filter out 13 samples with averaged read length <1000bp.
    df=df.loc[:,[c for c in df.columns if c not in Low_than_1k_sample_IDs]]
    logger.info(f'load {df.shape[1]} samples')
    gtm=GTOP_Tissue_Meta()
    tissue_id={}
    for x in df.columns:
        tis_code=x.split('-')[2]
        tis=gtm.get_meta_by_tissue_code(tis_code,'Tissue')
        if tis not in tissue_id:
            tissue_id[tis]=[]
        tissue_id[tis].append(x)
    for tis,sid in tissue_id.items():
        sdf=df.loc[:,sid]
        if min_sample is None:
            median_expr=sdf.median(axis=1,skipna=True)
            target_isoforms=set(median_expr.loc[median_expr > min_tpm].index.tolist())
        else:
            target_isoforms = set(sdf.index[(sdf > min_tpm).sum(axis=1) >= min_sample])

        logger.info(f'{tis} expressed {len(target_isoforms)} protein-coding genes')
        correct_tis=re.sub(r'\s+',"_",tis)
        faa_output=f'{faa_dir}/{correct_tis}.faa'
        with open(faa_output, 'w', encoding='utf-8') as out_f:
            # write gencode isoform
            logger.info(f'start write')
            k = __write_faa(enhanced_faa_path, target_isoforms, out_f, tid_gene=None, transcript_index=0, gene_index=1)
            logger.info(f'done {len(target_isoforms)} isoforms, {k} proteins; save to {faa_output}')

def prepare_diann_conf():
    faa_dir=f'{CONFIG_PROJECT_DIR}/release/cds/tissue_based'
    workdir = f'{CONFIG_PROJECT_DIR}/output/MS/diann/input'
    os.makedirs(workdir,exist_ok=True)
    with open(f"{workdir}/tissue_fa_list.txt",'w') as bw:
        for f in sorted(os.listdir(faa_dir)):
            if not f.endswith(".faa"):
                continue
            bw.write(f+'\n')
        pass

    # copy input/tissue_faa
    des_faa_dir=f'{workdir}/tissue_faa'
    os.makedirs(des_faa_dir,exist_ok=True)
    for f in sorted(os.listdir(faa_dir)):
        if not f.endswith(".faa"):
            continue
        shutil.copy(f'{faa_dir}/{f}',des_faa_dir)

    # copy MS raw file.
    # read tissue code
    gtm=GTOP_Tissue_Meta()
    raw_ms_dir=CONFIG_MS_RAW_DIR
    n = 0
    for sid in os.listdir(raw_ms_dir):
        tissue=gtm.get_meta_by_tissue_code(sid.split('-')[1],'Tissue')
        if sid in ['CA241-5032']:
            continue
        n+=1
        des_ms_dir=f'{workdir}/rawfiles/tissues/{tissue}'
        os.makedirs(des_ms_dir,exist_ok=True)
        if os.path.exists(f'{des_ms_dir}/{sid}.mzML'):
            print(f'skip {tissue}: {sid}')
            continue
        shutil.copy(f'{raw_ms_dir}/{sid}/{sid}.mzML',des_ms_dir)
        print(f'copy {tissue}: {sid}')
    print(f'copied {n} files.')


if __name__ == '__main__':
    ## make enhanced CDS by merging predicted CDS of novel tx and GENCODE v47 CDS.
    make_enhanced_cds_fa()
    ## tpm>5 in >=1 samples.
    tissue_based_protein_cds_fa(min_tpm=5, min_sample=1)
    ## prepare diann config
    prepare_diann_conf()
