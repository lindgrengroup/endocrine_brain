#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J create-batch

#SBATCH --output /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/logs/create_batch-%j.out 

#SBATCH --error /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/logs/create_batch-%j.err
#SBATCH -p short 

#SBATCH --cpus-per-task=1

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren/ferreira/PROJECTS/MRI/scripts/depict_analysis

while IFS=$'\t' read -r TRAIT MEASURE ANALYSIS	
	do
	echo $TRAIT $MEASURE $ANALYSIS
        scp run_depict_HT_volume_comb.sh run_depict_${TRAIT}_${MEASURE}_${ANALYSIS}.sh
        sed -i "s/HT_volume_comb/${TRAIT}_${MEASURE}_${ANALYSIS}/g" run_depict_${TRAIT}_${MEASURE}_${ANALYSIS}.sh
    done < /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/inputfile.txt



echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0