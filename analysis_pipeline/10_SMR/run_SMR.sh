main(){
	#make_esd
	#make_flist
	#make_besd_file
	run_SMR

}




make_esd(){
for xQTL_type in "tuQTL" "eQTL" "sQTL"
do
echo "#!/bin/bash
#SBATCH --job-name='make_esd_${xQTL_type}'
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=25
#SBATCH --error=make_esd_${xQTL_type}.err
#SBATCH --output=make_esd_${xQTL_type}.out

##################################

echo \"process will start at : \"

date

module load R/4.1.2-anaconda3

Rscript /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/bin/00_sort_esd.R ${xQTL_type}

echo \"process will end at : \"

date
" > slurm/make_esd_${xQTL_type}.slurm

done



}

make_flist(){
for xQTL_type in "eQTL" "tuQTL" "sQTL"
do

module load R/4.1.2-anaconda3
realpath /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/input/esd_${xQTL_type}/*/*.esd > /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/input/${xQTL_type}.flist 
Rscript /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/bin/01_make_flist.R ${xQTL_type}

done
}


make_besd_file(){
for tiss in `ls /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/input/esd_eQTL`
do

	for xQTL_type in "tuQTL" "sQTL" "eQTL"
do
smr="/flashfs1/scratch.global/hchen/biosoft/SMR_HEIDI_analysis/smr_Linux"
my_flist=/flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/input/flist_${xQTL_type}/${xQTL_type}_${tiss}_sorted.flist

echo "#!/bin/bash
#SBATCH --job-name='make_esd_${xQTL_type}_${tiss}'
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --error=make_besd_${xQTL_type}_${tiss}.err
#SBATCH --output=make_besd_${xQTL_type}_${tiss}.out 
    
##################################
    
echo \"process will start at : \"

date

mkdir -p /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/input/mybesd_${xQTL_type}/
module load R/4.1.2-anaconda3
${smr} --eqtl-flist ${my_flist} --make-besd --out /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/input/mybesd_${xQTL_type}/${tiss}
echo \"process will end at : \"

date
" > slurm/make_besd/make_besd_${xQTL_type}_${tiss}.slurm
	done
done

}

run_SMR(){

for xQTL_type in "eQTL" "tuQTL" "sQTL"
do
for GWAS in `ls /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/input/GWASs/*/*ma`
do
GWAS_name=`basename ${GWAS%.ma}`
smr="/flashfs1/scratch.global/hchen/biosoft/SMR_HEIDI_analysis/smr_Linux"
mkdir -p /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/slurm/run_smr

	for tiss in `less /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/input/tissue.list`
	do
	my_besd=/flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/input/mybesd_${xQTL_type}/${tiss}
	smr_out=/flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/output/${xQTL_type}/${GWAS_name}/${GWAS_name}_${tiss}.smr

if [ ! -f "${smr_out}" ]; then

echo "#!/bin/bash
#SBATCH --job-name='run_smr_${xQTL_type}_${tiss}_${GWAS_name}'
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --error=run_smr_${xQTL_type}_${tiss}_${GWAS_name}.err
#SBATCH --output=run_smr_${xQTL_type}_${tiss}_${GWAS_name}.out 
    
##################################
    
echo \"process will start at : \"

date

mkdir -p /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/output/${xQTL_type}/${GWAS_name}
module load R/4.1.2-anaconda3
EAS_plink=/flashfs1/scratch.global/hchen/00data/g1000_eas
my_smr=/flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/output/${xQTL_type}/${GWAS_name}/${GWAS_name}_${tiss}
${smr} --bfile \${EAS_plink} \
	--gwas-summary ${GWAS} \
	--beqtl-summary ${my_besd} \
	--out \${my_smr} \
	--thread-num 10 \
	--diff-freq 0.2 \
	--diff-freq-prop 0.05 \
	--heidi-mtd 1 \
	--peqtl-smr 5e-8 \
	--ld-upper-limit 0.9 \
	--ld-lower-limit 0.05 \
	--peqtl-heidi 1.57e-3 \
	--heidi-min-m 3 --heidi-max-m 20 --cis-wind 2000

${smr} --bfile \${EAS_plink} \
        --gwas-summary ${GWAS} \
        --beqtl-summary ${my_besd} \
        --out \${my_smr}.multi \
        --thread-num 10 \
        --diff-freq 0.2 \
        --diff-freq-prop 0.05 \
        --heidi-mtd 1 \
        --peqtl-smr 5e-8 \
        --ld-upper-limit 0.9 \
        --ld-lower-limit 0.05 \
        --peqtl-heidi 1.57e-3 \
        --heidi-min-m 3 \
	--heidi-max-m 20 \
	--cis-wind 2000 \
	--smr-multi --set-wind 500 --ld-multi-snp 0.1
 
echo \"process will end at : \"
date
" > /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/slurm/run_smr/run_smr_${xQTL_type}_${GWAS_name}_${tiss}.slurm
cd /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/slurm/run_smr/log
sbatch /flashfs1/scratch.global/hchen/2025-12-31-GTOP/2025-12-31-SMR/slurm/run_smr/run_smr_${xQTL_type}_${GWAS_name}_${tiss}.slurm

MAX_JOBS=150
USER=hchen
while true; do
    RUNNING_JOBS=$(squeue -u $USER -t R,PD -h |grep cu-1 | wc -l)
    
    if [ "$RUNNING_JOBS" -lt "$MAX_JOBS" ]; then
        echo "Current jobs: $RUNNING_JOBS，no more than $MAX_JOBS，submission ongoing..."
        break
    else
        echo "Current jobs: $RUNNING_JOBS, reaching limits  $MAX_JOBS waiting 2min!"
        sleep 120
    fi
done

else
	echo "${smr_out} file exists!"

fi

		done

	done

done

}

main
