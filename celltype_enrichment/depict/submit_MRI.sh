#!/bin/bash

while IFS=$'\t' read -r TRAIT MEASURE ANALYSIS	
	do
	echo $TRAIT $MEASURE $ANALYSIS
	sbatch -p short run_depict_${TRAIT}_${MEASURE}_${ANALYSIS}.sh
    done < /well/lindgren/ferreira/PROJECTS/MRI/depict_analysis/inputfile.txt
