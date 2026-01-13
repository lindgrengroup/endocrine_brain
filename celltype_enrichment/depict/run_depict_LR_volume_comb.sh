#!/bin/bash 

#SBATCH -A lindgren.prj 
#SBATCH -J LR_volume_comb_run-depict

#SBATCH --output /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/logs/LR_volume_comb_run_depict-%j.out 

#SBATCH --error /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/logs/LR_volume_comb_run_depict-%j.err 
#SBATCH -p short 

#SBATCH --cpus-per-task=1

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module purge 

module load Python/2.7.18-GCCcore-10.2.0

source /well/lindgren/ferreira/RESOURCES/venvs/DEPICT/bin/activate

cd /well/lindgren/ferreira/PROJECTS/MRI


        /well/lindgren/ferreira/RESOURCES/DEPICT/src/python/depict.py depict_analysis/config_files/LR_volume_comb.cfg



echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
