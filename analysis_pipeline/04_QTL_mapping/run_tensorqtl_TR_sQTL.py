import pandas as pd
import numpy as np
import tensorqtl
from tensorqtl import genotypeio, cis, trans
import argparse

# Create ArgumentParser object
parser = argparse.ArgumentParser(description='Function: Perform tensorQTL permutation, conditional, and trans analysis for specified tissue')

# Add arguments
parser.add_argument('-t', '--tissue', dest="tissue", type=str, required=True, help='Name of the tissue to analyze')
parser.add_argument('-g', '--geno', dest="genotype", type=str, required=True, help='Path to genotype file')
parser.add_argument('-p', '--pheno', dest="phenotype", type=str, required=True, help='Path to phenotype file in BED format')
parser.add_argument('-c', '--cov', dest="covariate", type=str, required=True, help='Path to covariates file')
parser.add_argument('-P', '--permutation', dest="permutation_output", type=str, required=True, help='Output path for permutation analysis results')
parser.add_argument('-N', '--nominal', dest="nominal_output_dir", type=str, required=True, help='Output directory for nominal analysis results')
parser.add_argument('-I', '--independent', dest="independent_output", type=str, required=True, help='Output path for independent analysis results')
parser.add_argument('-T', '--trans', dest="trans_output", type=str, required=True, help='Output path for trans analysis results')
parser.add_argument('-G', '--group', dest="cluster_group", type=str, required=True, help='Path to splicing group file')

# Parse arguments
args = parser.parse_args()

# Print received arguments for verification
print(f"Tissue: {args.tissue}")
print(f"Genotype file: {args.genotype}")
print(f"Phenotype file: {args.phenotype}")
print(f"Covariate file: {args.covariate}")
print(f"group file: {args.cluster_group}")

# Set random seed for reproducibility
seed = 12345

# Load genotype data
# Note: na_values='.' ensures proper handling of missing values
genotype_df = pd.read_table(args.genotype, sep='\t', index_col=0, na_values='.')
#group_df = pd.read_table(args.cluster_group, sep='\t', index_col=0, na_values='.')
#group_s = pd.Series(group_df.iloc[:,1].values, index=group_df.iloc[:,0])
group_df = pd.read_table(args.cluster_group, sep='\t', header=None, na_values='.', names=['phenotype_id', 'group'])
group_s = pd.Series(group_df['group'].values, index=group_df['phenotype_id'])

# Load phenotype data (BED format)
pheno_df, pheno_pos_df = tensorqtl.read_phenotype_bed(args.phenotype)

# Get sample intersection between genotype and phenotype data
geno_sample = genotype_df.columns
pheno_sample = pheno_df.columns
common_samples = list(set(geno_sample).intersection(set(pheno_sample)))
print(f"Common samples for {args.tissue}: {len(common_samples)}")

# Filter data to only include common samples
pheno_df = pheno_df[common_samples]
genotype_df = genotype_df[common_samples]
#group_df = group_df[common_samples]

# Load and prepare covariates data
covariates_df = pd.read_csv(args.covariate, sep='\t', index_col=0)[common_samples].T

print(f"Genotype matrix shape: {genotype_df.shape}")
print(f"Phenotype matrix shape: {pheno_df.shape}")
print(f"Covariates matrix shape: {covariates_df.shape}")

# Create variant information dataframe
# Split variant IDs into components (chromosome, position, etc.)
variant_df_empty = pd.DataFrame({"ID": genotype_df.index}, index=genotype_df.index)
variant_df = variant_df_empty["ID"].str.split('_', expand=True)
variant_df.columns = ["chrom", "pos", "end", "RU"]

# Keep only chromosome and position columns
variant_df.drop(columns=["RU", "end"], inplace=True)

# Ensure proper data types
variant_df = variant_df.astype({"chrom": str, "pos": int})

# Perform nominal QTL analysis
# This identifies all variant-phenotype pairs within cis-window
cis.map_nominal(genotype_df, variant_df, pheno_df, pheno_pos_df,covariates_df=covariates_df,output_dir=args.nominal_output_dir,
               write_top=True,  # Write top associations per gene
               prefix=args.tissue,
               window=1000000)  # 1Mb cis-window

# Perform permutation analysis
# This calculates empirical p-values by permuting phenotypes
cis_df = cis.map_cis(genotype_df, variant_df, pheno_df, pheno_pos_df,
                    covariates_df=covariates_df,group_s=group_s,
                    window=1000000,  # 1Mb cis-window
                    seed=seed)  # For reproducibility

# Clean data by removing NA p-values
print(f"Before removing NA values: {cis_df.shape}")
cis_df = cis_df.dropna(subset=['pval_beta'])
print(f"After removing NA values: {cis_df.shape}")

# Calculate q-values (FDR-adjusted p-values)
tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85)

# Sort by q-value and save results
cis_df = cis_df.sort_values("qval")
cis_df.to_csv(args.permutation_output, index=True, sep='\t')

# Perform conditional analysis
# Identifies independent signals after conditioning on top variant
try:
    indep_df = cis.map_independent(genotype_df, variant_df, cis_df, pheno_df,
                                 pheno_pos_df, covariates_df=covariates_df,
                                 seed=seed)
    indep_df.to_csv(args.independent_output, index=True, sep='\t')
except Exception as e:
    print(f"Conditional analysis failed: {str(e)}")
    print("Possible reason: No significant phenotypes at FDR ≤ 0.05")

# Perform trans-QTL analysis
# Identifies associations between variants and distant phenotypes
trans_df = trans.map_trans(genotype_df, pheno_df, covariates_df,
                          return_sparse=True,  # Efficient memory usage
                          pval_threshold=1e-5,  # Significance threshold
                          maf_threshold=0,  # No MAF filtering
                          batch_size=20000)  # Processing batch size

# Filter out cis-associations from trans results
print(f"Trans-QTL associations before cis-filtering: {trans_df.shape}")
trans_df = trans.filter_cis(trans_df, pheno_pos_df, variant_df, window=5000000)
print(f"Trans-QTL associations after cis-filtering: {trans_df.shape}")

# Save trans-QTL results
trans_df.to_csv(args.trans_output, index=False, sep='\t')
