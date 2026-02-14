# inspire by https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/ase_aggregate_by_individual.py
# 
import numpy as np
import scipy.stats
import pandas as pd
import argparse
import pyBigWig
import os
import subprocess
import io
import gzip
import pickle
import tempfile
import pybedtools
import qtl.stats


parser = argparse.ArgumentParser(description='ASE')
parser.add_argument('read_count_file_list', help='Read count file list (one per sample); [sample_id, tissue_site_detail, file_path]')
parser.add_argument('het_vcf')
parser.add_argument('gene_bed', help='Gene annotation BED file (replaces vep_dict)') #edit by wyn
parser.add_argument('--simulation_bias_file', help='?')
parser.add_argument('mappability_bigwig', help='Mappability track in bigWig format')
parser.add_argument('--tissue_abbreviations', help='File mapping tissue_site_detail to abbreviation')
parser.add_argument('lamp_values', help='Table with foreign allele frequency per individual')
parser.add_argument('individual_id', help='individual_id')
parser.add_argument('--coverage_cutoff', default=8, type=int, help='')
parser.add_argument('--other_ratio_cutoff', default=0.05, type=float, help='')
parser.add_argument('--mono_cutoff', default=0.01, type=float, help='')
parser.add_argument('-o', '--output_dir', default='.')
args = parser.parse_args()


print('Parsing inputs')
readcount_file_df = pd.read_csv(args.read_count_file_list, sep='\t', index_col=0)



print('Loading gene annotations')
gene_bed_df = pd.read_csv(args.gene_bed, sep='\t', header=None, 
                         usecols=[0, 1, 2, 3],  # chr, start, end, gene_id
                         names=['chr', 'start', 'end', 'gene_id'],
                         dtype={'chr': str, 'start': int, 'end': int, 'gene_id': str})


print('Parsing read count files')
readcount_df_list = []
for i,rfile in enumerate(readcount_file_df['ase_readcount_file']):
    readcount_df = pd.read_csv(rfile, sep='\t', index_col=2)
    readcount_df = readcount_df[['contig', 'position', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'otherBases']]
    readcount_df = readcount_df.rename(columns={'contig':'chr', 'position':'coord', 'refAllele':'ref', 'altAllele':'alt',
         'refCount':'refcount', 'altCount':'altcount', 'totalCount':'totalcount', 'otherBases':'othercount'})
    readcount_df = readcount_df[readcount_df['totalcount']>=args.coverage_cutoff]
    
    readcount_df['refratio'] = readcount_df['refcount']/readcount_df['totalcount']
    readcount_df['otherratio'] = readcount_df['othercount'] / (readcount_df['othercount'] + readcount_df['totalcount'])
    readcount_df['otherflag'] = (readcount_df['otherratio']>=args.other_ratio_cutoff)*1
    readcount_df['allcount'] = readcount_df['totalcount'] + readcount_df['othercount']
    sample_id = readcount_file_df.index[i]
    readcount_df['sampid'] = sample_id
    readcount_df['subjid'] = '-'.join(sample_id.split('-')[:2])
    readcount_df['tissue'] = readcount_file_df.loc[sample_id, 'tissue_site_detail']
    readcount_df['tissueabrv'] = readcount_file_df.loc[sample_id, 'tissue_site_detail']
    readcount_df['covflag'] = 0  # covflag is never 1, since filtered above (coverage_cutoff)
    
    readcount_df_list.append(readcount_df)


print('Loading VCF')
vcf_df = pd.read_csv(args.het_vcf, sep='\t', comment='#', header=None,
    names=['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'genotype'], dtype=str,
    usecols=['chr', 'pos', 'id', 'info','format', 'genotype'], index_col=2)

vcf_snp_id_df = pd.DataFrame(index=vcf_df.index, columns=['chr', 'coord', 'genotype', 'ensg', 'vtype', 'mapbias', 'mapflag', 'monoflag', 'mono_refcount', 'mono_altcount', 'mono_totalcount', 'mono_othercount'])
vcf_snp_id_df[['chr', 'coord']] = vcf_df[['chr', 'pos']]
vcf_snp_id_df['genotype'] = vcf_df['format']+';'+vcf_df['genotype']

print('filtering variant of ase res in vcf')
filtered_list = []
for df in readcount_df_list:
    filtered_df = df[df.index.isin(vcf_snp_id_df.index)]
    if not filtered_df.empty:
        filtered_list.append(filtered_df)

readcount_df_list = filtered_list


print('Annotating genes with BED file (replacing VEP)')

# vcf2bed
with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_bed:
    for var_id, row in vcf_df.iterrows():
        bed_line = f"{row['chr']}\t{int(row['pos'])-1}\t{row['pos']}\t{var_id}\n"
        tmp_bed.write(bed_line)
    vcf_bed_file = tmp_bed.name

try:
    vcf_bed = pybedtools.BedTool(vcf_bed_file)
    gene_bed = pybedtools.BedTool(args.gene_bed)
    

    intersected = vcf_bed.intersect(gene_bed, loj=True, wa=True, wb=True)
    

    var_to_genes = {}
    for feature in intersected:
        var_id = feature[3]  
        gene_id = feature[7] if feature[7] != '.' else 'NA' 
        
        
        if var_id not in var_to_genes:
            var_to_genes[var_id] = set()
    
        var_to_genes[var_id].add(gene_id)
    

    ensg = []
    vtype = []
    
    for i in vcf_df.index:
        if i in var_to_genes and var_to_genes[i]:
            gene_list = var_to_genes[i]
            ensg.append(','.join(gene_list))
            vtype.append('genic_region')
        else:
            ensg.append('NA')
            vtype.append('NA')
    
    vcf_snp_id_df['ensg'] = ensg
    vcf_snp_id_df['vtype'] = vtype
    
finally:

    os.unlink(vcf_bed_file)


print('Adding mappability')
mp = []
bw = pyBigWig.open(args.mappability_bigwig)
for c,p in zip(vcf_df['chr'], vcf_df['pos']):
    mp.append((bw.stats(c, int(p)-1, int(p), exact=True)[0]!=1) * 1)  # BED coordinates, 0-indexed; input must be int (not numpy)
bw.close()


vcf_snp_id_df['mapbias'] = 0
vcf_snp_id_df['mapflag'] = mp
vcf_snp_id_df['monoflag'] = 0
vcf_snp_id_df['mono_refcount'] = 0
vcf_snp_id_df['mono_altcount'] = 0
vcf_snp_id_df['mono_totalcount'] = 0
vcf_snp_id_df['mono_othercount'] = 0
for readcount_df in readcount_df_list:
     # combine read counts for each variant
     vcf_snp_id_df.loc[readcount_df.index, 'mono_refcount'] += readcount_df['refcount']
     vcf_snp_id_df.loc[readcount_df.index, 'mono_altcount'] += readcount_df['altcount']
     vcf_snp_id_df.loc[readcount_df.index, 'mono_totalcount'] += readcount_df['totalcount']
     vcf_snp_id_df.loc[readcount_df.index, 'mono_othercount'] += readcount_df['othercount']


print('Calculating statistics')
lamp = pd.read_csv(args.lamp_values, sep='\t', index_col=0).squeeze('columns').median()
ref = vcf_snp_id_df['mono_refcount']
tot = vcf_snp_id_df['mono_totalcount']
monop_list = scipy.stats.binom.cdf(tot-ref, tot, 1-lamp) + scipy.stats.binom.cdf(ref, tot, 1-lamp)  # monoallelic_p
monop_adj_list = qtl.stats.padjust_bh(monop_list)
vcf_snp_id_df['monoflag'] = (monop_adj_list > args.mono_cutoff) * 1

indiv_cov75_counts = []
for readcount_df in readcount_df_list:
    readcount_df['GENOTYPE_WARNING'] = vcf_snp_id_df.loc[readcount_df.index, 'monoflag']
    idx = (vcf_snp_id_df.loc[readcount_df.index, ['monoflag', 'mapbias', 'mapflag']].sum(axis=1)==0) & (readcount_df['otherflag']==0)
    indiv_cov75_counts.extend(list(readcount_df.loc[idx, 'totalcount']))
cov75 = np.percentile(indiv_cov75_counts, 75)


print('Calculating bias')
genomewide_bias = [0.0, 0.0, 0]
for readcount_df in readcount_df_list:
    idx = (readcount_df[['covflag', 'otherflag']].sum(axis=1) + vcf_snp_id_df.loc[readcount_df.index, ['mapbias', 'mapflag', 'monoflag']].sum(axis=1)) == 0
    refcountcov = readcount_df.loc[idx, 'refcount']
    altcountcov = readcount_df.loc[idx, 'altcount']
    totcountcov = refcountcov + altcountcov
    
    bias_keys = readcount_df.loc[idx, 'ref']+'/'+readcount_df.loc[idx, 'alt']
    
    idx2 = (refcountcov+altcountcov) > cov75
    refcountcov[idx2] = cov75*(refcountcov[idx2]/totcountcov[idx2])
    altcountcov[idx2] = cov75 - refcountcov[idx2]
    totcountcov[idx2] = cov75
    
    genomewide_bias[0] += refcountcov.sum()
    genomewide_bias[1] += totcountcov.sum()
    genomewide_bias[2] += refcountcov.shape[0]

genomewide_bias_value = float(genomewide_bias[0]) / genomewide_bias[1]


print('Calculating binomial tests, adjusted p-values')
for readcount_df in readcount_df_list:
    readcount_df['binom_p'] = [
         scipy.stats.binomtest(i, j, genomewide_bias_value, alternative='two-sided').pvalue 
         for i, j in zip(readcount_df['refcount'], readcount_df['totalcount'])
    ]
    readcount_df['nullratio'] = genomewide_bias_value
    idx = (readcount_df[['covflag', 'otherflag']].sum(axis=1) + vcf_snp_id_df.loc[readcount_df.index, ['mapbias', 'mapflag', 'monoflag']].sum(axis=1))==0
    readcount_df.loc[idx, 'binom_p_adj'] = qtl.stats.padjust_bh(readcount_df.loc[idx, 'binom_p'])
    readcount_df.loc[~idx, 'binom_p_adj'] = 'NA'


print('Writing output')
with gzip.open(os.path.join(args.output_dir, args.individual_id+'.ase_table.tsv.gz'), 'wt') as f:
    f.write('\t'.join([
        'CHR',
        'POS',
        'VARIANT_ID',
        'REF_ALLELE',
        'ALT_ALLELE',
        'SAMPLE_ID',
        'SUBJECT_ID',
        'TISSUE_ID',
        'REF_COUNT',
        'ALT_COUNT',
        'TOTAL_COUNT',
        'REF_RATIO',
        'OTHER_ALLELE_COUNT',
        'NULL_RATIO',
        'BINOM_P',
        'BINOM_P_ADJUSTED',
        'MAMBA_POST_SINGLETIS',
        'MAMBA_POST_MULTITIS',
        'GENOTYPE',
        'VARIANT_ANNOTATION',
        'GENE_ID',
        'LOW_MAPABILITY',
        'MAPPING_BIAS_SIM',
        'GENOTYPE_WARNING'])+'\n')
    
    merged_df = []
    for readcount_df in readcount_df_list:
        readcount_df['id'] = readcount_df.index
        readcount_df['blank'] = 'NA'
        out_df = readcount_df[['chr', 'coord', 'id', 'ref', 'alt', 'sampid', 'subjid', 'tissueabrv', 'refcount', 'altcount', 'totalcount', 'refratio', 'othercount', 'nullratio',
            'binom_p', 'binom_p_adj', 'blank', 'blank']]
        merged_df.append(pd.concat([out_df, vcf_snp_id_df.loc[readcount_df.index, ['genotype', 'vtype', 'ensg', 'mapflag', 'mapbias']], readcount_df['GENOTYPE_WARNING']], axis=1))
    merged_df = pd.concat(merged_df, axis=0)
    merged_df = merged_df.sort_values(['chr', 'coord', 'tissueabrv'])
    merged_df.to_csv(f, sep='\t', index=False, header=False, float_format='%.6g')

print('Done')
