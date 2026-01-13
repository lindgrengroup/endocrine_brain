#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J run-plots

#SBATCH --output /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/logs/run_plots-%j.out 

#SBATCH --error /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/logs/run_plots-%j.err 
#SBATCH -p short

#SBATCH --cpus-per-task=3

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################
 
module purge 

module load R/3.6.2-foss-2019b

cd /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/tissue_plots

while IFS=$'\t' read -r TRAIT MEASURE ANALYSIS	
	do
	echo $TRAIT $MEASURE $ANALYSIS
    Rscript /well/lindgren/ferreira/PROJECTS/MRI/scripts/depict_analysis/tissue_enrich_plot.R /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/results/${TRAIT}_${MEASURE}_${ANALYSIS}_tissueenrichment.txt ${TRAIT}_${MEASURE}_${ANALYSIS}
done < /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/inputfile.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
