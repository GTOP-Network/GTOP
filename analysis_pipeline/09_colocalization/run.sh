
main(){
    finemapping_for_GWAS_1KGP_EAS

    GWAS_QTL_pair

    finemapping_for_QTL_GTOP

    coloc_for_GTOP_prior_combination
    
    coloc_for_GTOP_LD_combination

    merge_coloc_res
}


finemapping_for_GWAS_1KGP_EAS(){

    mkdir -p input/slurm/fmGWAS_1KGPEAS/log
   
    tail -n+2 input/EAS_GWAS_noMHC_all_sentinal_loci.txt | while read line
    do
        GWAS_name=`echo ${line} | cut -f 1 -d " "`
        CHR_name=`echo ${line} | cut -f 2 -d " "`
        sentinel_rsid=`echo ${line} | cut -f 5 -d " "`
        loci_start=`echo ${line} | cut -f 11 -d " "`
        loci_end=`echo ${line} | cut -f 12 -d " "`
        task=${GWAS_name}_${sentinel_rsid}

        files=(input/finemapping_GWAS_1KGP_EAS/${GWAS_name}_${sentinel_rsid}*.RData)

        if [[ ! -f "${files[0]}" ]]; then
        
        echo $task
        echo "#!/bin/bash

#SBATCH --job-name=$task
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --error=input/slurm/fmGWAS_1KGPEAS/log/${task}.err
#SBATCH --output=input/slurm/fmGWAS_1KGPEAS/log/${task}.out

Rscript a1.finemapping_for_GWAS.R ${GWAS_name} ${CHR_name} ${sentinel_rsid} ${loci_start} ${loci_end}

" > input/slurm/fmGWAS_1KGPEAS/${task}.slurm
        sbatch input/slurm/fmGWAS_1KGPEAS/${task}.slurm
    
        process=$(squeue|wc -l)
        while((process >= 200))
        do
            echo "Current jobs' count is larger than 200"
            echo "Wait for another 60 s"
            sleep 60
            process=$(squeue|wc -l)
        done
        fi
    done
}


GWAS_QTL_pair(){
    mkdir -p input/slurm/GWAS_QTL_overlap/log

    declare -A TISSUES

    TISSUES[snv_eqtl]="Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood"
    TISSUES[snv_juqtl]="Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood"
    TISSUES[snv_tuqtl]="Adipose,Adrenal_Gland,Gallbladder,Liver,Muscle,Pancreas_Body,Pancreas_Head,Pancreas_Tail,Skin,Spleen,Whole_Blood"
    TISSUES[gtexv8_eqtl]="Adipose_Visceral_Omentum,Adrenal_Gland,Liver,Muscle_Skeletal,Pancreas,Skin_Not_Sun_Exposed_Suprapubic,Spleen,Whole_Blood"
    TISSUES[jctf_eqtl]="Whole_Blood"
    TISSUES[mage_eqtl]="LCL"


    for xQTL_type in snv_eqtl snv_juqtl snv_tuqtl gtexv8_eqtl jctf_eqtl mage_eqtl
    do
        for tissue in ${TISSUES[$xQTL_type]//,/ }
        do
            task="${xQTL_type}_${tissue}"

            echo "${task}"

            cat > input/slurm/GWAS_QTL_overlap/${task}.slurm <<EOF
#!/bin/bash

#SBATCH --job-name=${task}
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --error=input/slurm/GWAS_QTL_overlap/log/${task}.err
#SBATCH --output=input/slurm/GWAS_QTL_overlap/log/${task}.out

Rscript a2.prepare_genes.R ${xQTL_type} ${tissue}

EOF

            process=$(squeue -u "$USER" -h | wc -l)
            while (( process >= 200 ))
            do
                echo "Current job count: ${process} >= 200"
                echo "Wait for another 20 s"
                sleep 20

                process=$(squeue -u "$USER" -h | wc -l)
            done

            sbatch input/slurm/GWAS_QTL_overlap/${task}.slurm

        done
    done
}


## QTL: finemapping 
finemapping_for_QTL_GTOP(){

    mkdir input/slurm/run_data
    awk 'FNR == 1 && NR != 1 {next} {print}' input/GWAS_QTL_overlap/{snv_eqtl,snv_juqtl,snv_tuqtl}/*.txt > input/slurm/gtop_gwas_qtl_overlap_all.txt
    tail -n+2 input/slurm/gtop_gwas_qtl_overlap_all.txt | cut -f 1-3 | sort -k2,2 -k3,4 | uniq > input/slurm/run_data/gtop_prepare_genes.txt
    rm -r input/slurm/run_data/gtop_prepare_genes
    mkdir -p input/slurm/run_data/gtop_prepare_genes
    split -l 300 input/slurm/run_data/gtop_prepare_genes.txt input/slurm/run_data/gtop_prepare_genes/task_

    tasknum=task1
    mkdir -p input/slurm/fmQTL_GTOP/${tasknum}/log 

    for f in `ls input/slurm/run_data/gtop_prepare_genes/task_* `
    do
        task=`basename $f`
        echo $task
        echo "#!/bin/bash

#SBATCH --job-name=${tasknum}_$task
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --error=input/slurm/fmQTL_GTOP/${tasknum}/log/${task}.err
#SBATCH --output=input/slurm/fmQTL_GTOP/${tasknum}/log/${task}.out

cat $f | while read line
do
    QTLTYPE=\`echo \$line | awk '{print \$1}'\`
    TISSUE=\`echo \$line | awk '{print \$2}'\`
    EGENE=\`echo \$line | awk '{print \$3}'\`
    
    files=(input/finemapping_QTL_GTOP/\${QTLTYPE}/\${TISSUE}/\${EGENE}.*.RData)

    if [[ ! -f \"\${files[0]}\" ]]; then
        Rscript a3.finemapping_for_QTL_GTOP_LD.R \${QTLTYPE} \${TISSUE} \${EGENE}
    fi
done

"   >   input/slurm/fmQTL_GTOP/${tasknum}/${task}.slurm

        sbatch input/slurm/fmQTL_GTOP/${tasknum}/${task}.slurm

        process=$(squeue|wc -l)
        while((process >= 200))
        do
            echo "Current jobs' count is larger than 200"
            echo "Wait for another 60 s"
            sleep 20
            process=$(squeue|wc -l)
        done

    done
}



coloc_for_GTOP_prior_combination(){

    tail -n+2 input/slurm/gtop_gwas_qtl_overlap_all.txt | grep -v "xqtl_type" > input/slurm/run_data/coloc_for_gtop.txt
    wc -l input/slurm/run_data/coloc_for_gtop.txt
    rm -r input/slurm/run_data/coloc_for_gtop
    mkdir -p input/slurm/run_data/coloc_for_gtop
    split -l 800 input/slurm/run_data/coloc_for_gtop.txt input/slurm/run_data/coloc_for_gtop/task_

    rm -r input/slurm/coloc

    tail -n+2 input/prior_setting.txt | awk '{OFS="_";print $2,$3,$4}' | while read line
    do
        mkdir -p input/slurm/coloc/prior_setting/${line}/{log,output}

        for f in `ls input/slurm/run_data/coloc_for_gtop/task_* `
        do
            task=`basename $f`
            echo $task
            echo "#!/bin/bash

#SBATCH --job-name=prior_${task}_${line}
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --error=input/slurm/coloc/prior_setting/${line}/log/${task}.err
#SBATCH --output=input/slurm/coloc/prior_setting/${line}/log/${task}.out

Rscript a4.coloc.R $f finemapping_GWAS_1KGP_EAS finemapping_QTL_GTOP ${line//_/ } 3 input/slurm/coloc/prior_setting/${line}/output/${task}.out
"   >       input/slurm/coloc/prior_setting/${line}/${task}.slurm

            sbatch input/slurm/coloc/prior_setting/${line}/${task}.slurm

            process=$(squeue|wc -l)
            while((process >= 400))
            do
                echo "Current jobs' count is larger than 400"
                echo "Wait for another 60 s"
                sleep 20
                process=$(squeue|wc -l)
            done

        done
    done
}


coloc_for_GTOP_LD_combination(){
    for gwas_ld in GTOP
    do
        for qtl_ld in 1KGP_EAS GTOP
        do
            mkdir -p input/slurm/coloc/ld_combination/GTOP_xqtl_${gwas_ld}_${qtl_ld}/{log,output}

            for f in `ls input/slurm/run_data/coloc_for_gtop/task_* `
            do
                task=`basename $f`
                echo $task
                echo "#!/bin/bash

#SBATCH --job-name=ld_${gwas_ld}_${qtl_ld}_${task}
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --error=input/slurm/coloc/ld_combination/GTOP_xqtl_${gwas_ld}_${qtl_ld}/log/${task}.err
#SBATCH --output=input/slurm/coloc/ld_combination/GTOP_xqtl_${gwas_ld}_${qtl_ld}/log/${task}.out

Rscript a4.coloc.R $f finemapping_GWAS_${gwas_ld} finemapping_QTL_${qtl_ld} 1e-04 1e-04 5e-06 3 input/slurm/coloc/ld_combination/GTOP_xqtl_${gwas_ld}_${qtl_ld}/output/${task}.out

"   >           input/slurm/coloc/ld_combination/GTOP_xqtl_${gwas_ld}_${qtl_ld}/${task}.slurm

                sbatch input/slurm/coloc/ld_combination/GTOP_xqtl_${gwas_ld}_${qtl_ld}/${task}.slurm

                process=$(squeue|wc -l)
                while((process >= 400))
                do
                    echo "Current jobs' count is larger than 400"
                    echo "Wait for another 60 s"
                    sleep 20
                    process=$(squeue|wc -l)
                done

            done
        done
    done
}


merge_coloc_res(){
    # for i in ld_combination prior_setting; do ls input/slurm/coloc/${i} | while read line;do  ls input/slurm/coloc/${i}/${line}/output | wc -l;done;done
    # for i in naive_coloc; do ls input/slurm/coloc/${i} | while read line;do  ls input/slurm/coloc/${i}/${line}/output | wc -l;done;done

    for source_type in gtop
    do
        mkdir -p input/slurm/merge_res/${source_type}/log

        echo "#!/bin/bash

#SBATCH --job-name=${source_type}
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --error=input/slurm/merge_res/${source_type}/log/info.err
#SBATCH --output=input/slurm/merge_res/${source_type}/log/info.out

find input/slurm/coloc/ld_combination -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | \
xargs -I {} Rscript a5.merge_result.R input/slurm/coloc/ld_combination/{}/output output/${source_type}/ld_combination/{}

find input/slurm/coloc/prior_setting -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | \
xargs -I {} Rscript a5.merge_result.R input/slurm/coloc/prior_setting/{}/output output/${source_type}/prior_setting/{}

" > input/slurm/merge_res/${source_type}/info.slurm

    sbatch input/slurm/merge_res/${source_type}/info.slurm

    done
}

main
