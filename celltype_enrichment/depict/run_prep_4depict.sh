#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J prep4depict

#SBATCH --output /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/logs/prep4depict-%j.out 

#SBATCH --error /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/logs/prep4depict-%j.err 
#SBATCH -p short 

#SBATCH --cpus-per-task=1



echo `date`: Executing job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################
module purge 

module load R/3.6.2-foss-2019b


while IFS=$'\t' read -r TRAIT MEASURE ANALYSIS	
	do
	echo $TRAIT $MEASURE $ANALYSIS
	Rscript /well/lindgren/ferreira/PROJECTS/MRI/scripts/depict_analysis/prep_4depict.R $TRAIT $MEASURE $ANALYSIS  > /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/logs/${TRAIT}_${MEASURE}_${ANALYSIS}_prep4depict.log
    done < /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/inputfile.txt


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

