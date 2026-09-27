#!/usr/bin/env bash
# MAJIQTL analysis for GTOP. Edit the settings below before running:
# bash run_MAJIQTL.sh > majiqtl.log 2>&1
# Requirements: the study's MAJIQTL Python environment, bgzip, tabix, awk,
# sort, cut and gzip. Activate the environment before running this script.
# No scheduler submission, chromosome splitting, or input-check functions.
set -euo pipefail

# Software and input paths. Use the MAJIQTL version used in the study.
majiqtl_dir="/path/to/majiqtl"
weight_model="$majiqtl_dir/data/example/input/gp_weight_model.pkl"
phenotype_dir="/path/to/original_phenotypes"
covariate_dir="/path/to/covariates_for_sQTL/juQTL"
snp_vcf="/path/to/SNP.genotypes.vcf.gz"
sv_vcf="/path/to/SV.genotypes.vcf.gz"
tr_vcf="/path/to/TR.genotypes.vcf.gz"

# Input VCFs must already be filtered, bgzip-compressed and indexed.
# Supply all desired chromosomes in each VCF; do not split by chromosome.
# Prepare TR VCF upstream: the original convert_tr_to_vcf.py is not included.
# Phenotype, genotype and covariate sample IDs must correspond.
tissues="Adipose Adrenal_Gland Gallbladder Liver Muscle Pancreas_Body Pancreas_Head Pancreas_Tail Skin Spleen Whole_Blood"
gtypes="SNP TR SV"
covsets="PMI_RIN base"
threads=32
sv_threads=5
window=1000000
outdir="results_majiqtl"

# Optional extraction of nominal associations for TensorQTL significant genes.
# Disabled by default because it requires independently computed TensorQTL results.
extract_tensorqtl=false
tensor_covsets="PMI_RIN"
tensor_base="/path/to/tensorqtl/output"
tr_tensor_base="/path/to/TR_tensorqtl/output"
# Original workflow used the penultimate column as the significance measure.
# Set this to the actual q-value column name (e.g. qval) to select by header.
tensor_qvalue_column="penultimate"
fdr=0.05

# Restrict BLAS threads to avoid oversubscription with MAJIQTL workers.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

pre_phenotype() {
    # Input: <phenotype_dir>/<tissue>/splicing.phenotype.sorted.bed
    # Columns: chromosome, start, end, phenotype ID, then sample values.
    # Study ID convention: the last colon-separated field is the gene ID;
    # the last underscore-separated token of field 4 specifies the strand.
    # Preserve input coordinate order for tabix indexing.
    for tis in $tissues; do
        awk 'BEGIN {FS=OFS="\t"}
            NR==1 {
                printf "#Chr\t%s\t%s\tpid\tgid\tstrand", $2, $3
                for (i=5; i<=NF; i++) printf "\t%s", $i
                printf "\n"
                next
            }
            {
                n=split($4, a, ":")
                m=split(a[4], b, "_")
                printf "%s\t%s\t%s\t%s\t%s\t%s", $1,$2,$3,$4,a[n],b[m]
                for (i=5; i<=NF; i++) printf "\t%s", $i
                printf "\n"
            }' "$phenotype_dir/$tis/splicing.phenotype.sorted.bed" |
            bgzip -c > "$outdir/phenotype/$tis.phenotype.bed.gz"
        tabix -p bed "$outdir/phenotype/$tis.phenotype.bed.gz"
    done
}

run_nominal() {
    # Run once for the entire VCF, using PSI and a 1 Mb gene-based cis window.
    # -CHR is intentionally omitted. The installed MAJIQTL build must support
    # genome-wide execution without this option; verify with its --help.
    python "$majiqtl_dir/src/sqtl/majiqtl_sqtl.py" \
        --vcf "$gvcf" \
        --bed "$bed" \
        --cov "$cov" \
        --out "$prefix" \
        --mode psi \
        --window_pos gene \
        --window "$window" \
        -J "$jobs"
}

sort_nominal() {

    mkdir "$run_dir/unsorted_nominal"
    for nominal_file in "$prefix".[0-9]*.tmp; do
        case "$nominal_file" in
            *.sgene.tmp|*.ori.tmp) continue ;;
        esac
        cp "$nominal_file" "$run_dir/unsorted_nominal/"
        awk '{
            n=split($1,a,":")
            gene=a[n]
            sub(/\..*/,"",gene)
            print gene "\t" NR "\t" $0
        }' "$nominal_file" |
            LC_ALL=C sort -k1,1 -k2,2n |
            cut -f3- > "$nominal_file.sorted"
        mv "$nominal_file.sorted" "$nominal_file"
    done
}

run_sgene() {
    # Use the same VCF, phenotype, prefix, cis window and worker count as
    # the nominal run. All nominal files are ready before this step starts.
    python "$majiqtl_dir/src/sqtl/majiqtl_sgene.py" \
        "$gvcf" "$bed" "$prefix" "$window" "$jobs" "$weight_model"
}

merge_nominal() {
    # Combine worker-level nominal results; there are no chromosome files.
    # Each association remains one row; concatenation does not change tests.
    for nominal_file in "$prefix".[0-9]*.tmp; do
        case "$nominal_file" in
            *.sgene.tmp|*.ori.tmp) continue ;;
        esac
        cat "$nominal_file"
    done | gzip -c > "$run_dir/$tis.all.tmp.gz"
}



# Use a new output directory for each run. mkdir stops on an existing directory.
# This script does not modify input files or the original cluster scripts.
mkdir "$outdir"
mkdir "$outdir/phenotype"
pre_phenotype

# Process covariate sets, variant types and tissues sequentially.
for covset in $covsets; do
    for gtype in $gtypes; do
        jobs="$threads"
        if [ "$gtype" = "SNP" ]; then
            gvcf="$snp_vcf"
        elif [ "$gtype" = "SV" ]; then
            gvcf="$sv_vcf"
            jobs="$sv_threads"
        else
            gvcf="$tr_vcf"
        fi
        for tis in $tissues; do
            echo "Running MAJIQTL: $covset / $gtype / $tis"
            bed="$outdir/phenotype/$tis.phenotype.bed.gz"
            cov="$covariate_dir/covariates_${covset}/$tis.GPC0.covariates.txt"
            run_dir="$outdir/$covset/$gtype/$tis"
            mkdir -p "$run_dir"
            prefix="$run_dir/$tis"
            run_nominal
            sort_nominal
            run_sgene
            merge_nominal

        done
    done
done

# Keep MAJIQTL's native sGene outputs in each run directory.
# The original cross-chromosome R summary depended on bin/merged.cov.r,
# which was not supplied; no replacement statistical aggregation is assumed.
echo "Completed: $outdir"
