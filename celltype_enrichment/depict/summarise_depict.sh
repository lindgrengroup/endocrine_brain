#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J depict_summ

#SBATCH --output /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/logs/depict_summ-%j.out 

#SBATCH --error /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/logs/depict_summ-%j.err 
#SBATCH -p short

#SBATCH --cpus-per-task=1

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################
 
module purge 

module load R/3.6.2-foss-2019b

cd /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/results

while IFS=$'\t' read -r TRAIT MEASURE ANALYSIS	
	do
	echo $TRAIT $MEASURE $ANALYSIS
    awk  -F'\t' '{if($1 == "MeSH term" || ($5 < 0.05 && ($6 == "<0.01" || $6 == "<0.05"))){print $0}}' /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/results/${TRAIT}_${MEASURE}_${ANALYSIS}_tissueenrichment.txt > /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/results/${TRAIT}_${MEASURE}_${ANALYSIS}_tissueenrichment_signif.txt
    awk  -F'\t' '{if($1 == "Locus" || ($7 < 0.05 && ($10 == "<0.01" || $10 == "<0.05"))){print $0}}' /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/results/${TRAIT}_${MEASURE}_${ANALYSIS}_geneprioritization.txt > /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/results/${TRAIT}_${MEASURE}_${ANALYSIS}_geneprioritization_signif.txt
    awk  -F'\t' '{if($1 == "Original gene set ID" || ($3 < 0.05 && ($4 == "<0.01" || $4 == "<0.05"))){print $0}}' /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/results/${TRAIT}_${MEASURE}_${ANALYSIS}_genesetenrichment.txt > /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/results/${TRAIT}_${MEASURE}_${ANALYSIS}_genesetenrichment_signif.txt
    done < /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/inputfile.txt


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0