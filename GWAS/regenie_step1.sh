#######################
###RUN REGENIE STEP1###
#######################

run_regenie_step_one="regenie --step 1\
 --lowmem --out ${pheno_name} --bed ukb22418_c1_22X_v2_merged\
 --phenoFile ${pheno_file} \
 --covarFile ${covar_file}\
 --extract genotype_array_snps_qc_pass.snplist\
 --bsize 1000 --qt --loocv --gz --threads 32"

dx run swiss-army-knife -iin="${project_out}/genotype_process/ukb22418_c1_22X_v2_merged.bed" \
   -iin="${project_out}/genotype_process/ukb22418_c1_22X_v2_merged.bim" \
   -iin="${project_out}/genotype_process/ukb22418_c1_22X_v2_merged.fam"\
   -iin="${pheno_input}" \
   -iin="${covar_input}" \
   -iin="${snplist_qc}" \	
   -icmd="${run_regenie_step_one}" --tag="Step1" --instance-type "mem3_ssd1_v2_x32"\
   --destination="${project_out}/genotype_process/${pheno_group}/regenie_step1/" --brief
