
main(){
    qtl_mapping
    mash
}

qtl_mapping(){
    cd /path/to/dir/

    tail -n+2 /path/to/input/gtop_gtex_tissues.txt | cut -f 2 | uniq | parallel -j 2 /path/to/anaconda3/envs/work/bin/python bin/scripts/cis_map.py -i input/gtexv8_eqtl_cellproportion -o input/gtexv8_eqtl_cellproportion -g /path/to/src/data/GTEx/v8/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased -w 1000000 -t {}

    tail -n+2 /path/to/input/gtop_gtex_tissues.txt | cut -f 2 | uniq | parallel -j 2 /path/to/anaconda3/envs/work/bin/python /path/to/bin/qtl_process/calculate_qvalues.py -i input/gtexv8_eqtl_cellproportion/03_permutations/tmp_permutation -o input/gtexv8_eqtl_cellproportion/03_permutations -t {}

    tail -n+2 path/to/input/gtop_gtex_tissues.txt | cut -f 2 | uniq | parallel -j 2 path/to/anaconda3/envs/work/bin/python /path/to/bin/qtl_process/cis_nominal.py -i input/gtexv8_eqtl_cellproportion/05_nominal/tmp_nominal_parquet -o input/gtexv8_eqtl_cellproportion/05_nominal -t {}
}

mash(){
    # GTEx data


    QTL_type=ct_eqtl
    QTL_type=noct_eqtl

    dir_permutation=output/data/mash/${QTL_type}/permutation
    dir_nominal=output/data/mash/${QTL_type}/nominal
    dir_output_Strong=output/data/mash/${QTL_type}/Strong
    dir_output_Random=output/data/mash/${QTL_type}/Random
    dir_slurm=output/data/mash/${QTL_type}/slurm

    /path/to/bin/MashR/1-1Prepare-MashR-StrongPairs.sh ${dir_permutation} ${dir_output_Strong}
    /path/to/bin/MashR/1-2Prepare-MashR-StrongPairs.sh ${dir_permutation} ${dir_nominal} ${dir_output_Strong}
    /path/to/bin/MashR/1-3Prepare-MashR-StrongPairs.sh ${dir_permutation} ${dir_output_Strong}
    /path/to/bin/MashR/2-1Prepare-MashR-RandomPairs.sh ${dir_nominal} ${dir_output_Random}

    mkdir -p ${dir_slurm}/log
    rm ${dir_slurm}/${QTL_type}_*
    perm_files=(`find ${dir_permutation} -name "*.txt"`)
    echo ${#perm_files[@]}
    nFile=${#perm_files[@]}
    NAMEs=(${perm_files[@]/*\//})
    NAMEs=(${NAMEs[@]/.txt/})

    for i in `seq 0 $[nFile-1]`
    do
        name=${NAMEs[i]}
        nominal_files=(`find ${dir_nominal} -name "${name}.txt.gz"`)
        echo "/bin/sh

#SBATCH --job-name=${QTL_type}_${name}
#SBATCH --partition=Compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --error=${dir_slurm}/log/${QTL_type}_${name}.err
#SBATCH --output=${dir_slurm}/log/${QTL_type}_${name}.out

# if [[ ! -f ${dir_output_Random}/${name}.nominal_pairs.extracted_pairs.txt.gz ]];then
echo ${dir_nominal}/${name}.txt.gz > ${dir_output_Random}/${name}.nominal_files.txt
python /path/to/bin/MashR/extract_pairs_tjy.py ${dir_output_Random}/${name}.nominal_files.txt ${dir_output_Random}/nominal_pairs.combined_signifpairs.txt.gz ${name}.nominal_pairs -o ${dir_output_Random}
# rm -f ${dir_output_Random}/${name}.nominal_files.txt

"   > ${dir_slurm}/${QTL_type}_${name}.slurm
        sed -i "s/\/bin\/sh/\#\!\/bin\/sh/" ${dir_slurm}/${QTL_type}_${name}.slurm
    done

    ls ${dir_slurm}/${QTL_type}_*.slurm | parallel -j 3 bash {}

    ls output/data/mash/${QTL_type}/Random/*.nominal_pairs.extracted_pairs.txt.gz | wc -l

    /path/to/bin/MashR/2-3Prepare-MashR-RandomPairs.sh ${dir_output_Random}

    subset_size=1000000
    Rscript /path/to/bin/MashR/MashR-random_subset.R ${dir_output_Random} ${subset_size}

    ll -h output/data/mash/${QTL_type}/Random/*.MashR_input.txt.gz
    ll -h output/data/mash/${QTL_type}/Random/MashR.random_subset_1000000.RDS

    Rscript /path/to/bin/MashR/run_MashR.R ${dir_output_Strong}/strong_pairs.MashR_input.txt.gz ${dir_output_Random}/MashR.random_subset_${subset_size}.RDS 0 ./output/data/mash/${QTL_type}/top_pairs

    ## finemapping result
    QTL_type=noct_eqtl

    dir_finemapping=output/data/mash/${QTL_type}/finemapping
    dir_nominal=output/data/mash/${QTL_type}/nominal
    dir_output_finemapping=output/data/mash/${QTL_type}/Strong2

    /path/to/bin/MashR/3-1Prepare-MashR-finemappingPairs.sh ${dir_finemapping} ${dir_output_finemapping}
    /path/to/bin/MashR/3-2Prepare-MashR-finemappingPairs.sh ${dir_finemapping} ${dir_nominal} ${dir_output_finemapping}
    /path/to/bin/MashR/3-3Prepare-MashR-finemappingPairs.sh ${dir_finemapping} ${dir_output_finemapping}
}

main
