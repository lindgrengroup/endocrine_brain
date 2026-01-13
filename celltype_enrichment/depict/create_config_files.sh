#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J create-batch-infert

#SBATCH --output /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/logs/create-batch-infert-%j.out 

#SBATCH --error /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/logs/create-batch-infert-%j.err 
#SBATCH -p short 

#SBATCH --cpus-per-task=1

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/config_files

while IFS=$'\t' read -r TRAIT MEASURE ANALYSIS	
	do
	echo $TRAIT $MEASURE $ANALYSIS
        scp HT_volume_comb.cfg ${TRAIT}_${MEASURE}_${ANALYSIS}.cfg
        sed -i "s/HT/${TRAIT}/g" ${TRAIT}_${MEASURE}_${ANALYSIS}.cfg
        sed -i "s/volume/${MEASURE}/g" ${TRAIT}_${MEASURE}_${ANALYSIS}.cfg
        sed -i "s/comb/${ANALYSIS}/g" ${TRAIT}_${MEASURE}_${ANALYSIS}.cfg
        sed -i "s/gwas_updated/gwas_update_oct2023/g" ${TRAIT}_${MEASURE}_${ANALYSIS}.cfg
    done < /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/inputfile.txt



echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0