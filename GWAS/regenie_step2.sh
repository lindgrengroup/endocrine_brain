genotype_file_dir="Bulk/Imputation/UKB imputation from genotype/"

#loco_input=$(dx ls "${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_*.loco.gz" | awk '{ sub("^", "-iin=\"${project_out}/genotype_process/${pheno_group}/regenie_step1/", $0); printf("%s\" ", $0) }' | tr -d '\n')

for chr in {1..22}; do
  run_regenie_cmd="regenie --step 2 --bgen ukb22828_c${chr}_b0_v3.bgen --out assoc.c${chr}\
    --sample ukb22828_c${chr}_b0_v3.sample \
    --phenoFile ${pheno_file} --covarFile ${covar_file}\
    --qt \
    --pred ${pheno_name}_pred.list --bsize 200\
    --threads 32 --gz"

###############################
###ORIGINAL BRAIN PHENOTYPES###
###############################

 dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.bgen" \
   -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.sample" -iin="${pheno_input}" -iin="${covar_input}" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_pred.list" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_1.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_2.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_3.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_4.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_5.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_6.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_7.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_8.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_9.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_10.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_11.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_12.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_13.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_14.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_15.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_16.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_17.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_18.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_19.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_20.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_21.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_22.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_23.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_24.loco.gz" \
   -icmd=${run_regenie_cmd} --tag="Step2" --instance-type "mem3_ssd1_v2_x32"\
   --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/" --brief -y
done

dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.bgen" \
   -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.sample" -iin="${pheno_input}" -iin="${covar_input}" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_pred.list" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_1.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_2.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_3.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_4.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_5.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_6.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_7.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_8.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_9.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_10.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_11.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_12.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_13.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_14.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_15.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_16.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_17.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_18.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_19.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_20.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_21.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_22.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_23.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_24.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_25.loco.gz" \
   -icmd=${run_regenie_cmd} --tag="Step2" --instance-type "mem3_ssd1_v2_x32"\
   --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/" --brief -y
done



#################
###HEIGHT GWAS###
#################

dx run  swiss-army-knife -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.bgen" \
   -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.sample" -iin="${pheno_input}" -iin="${covar_input}" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_pred.list" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_1.loco.gz" \
   -icmd=${run_regenie_cmd} --tag="Step2" --instance-type "mem3_ssd1_v2_x32"\
   --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/" --brief -y
done

###########################
###OB LR MEAN PHENOTYPES###
###########################

 dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.bgen" \
   -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.sample" -iin="${pheno_input}" -iin="${covar_input}" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_pred.list" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_1.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_2.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_3.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_4.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_5.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_6.loco.gz" \
   -icmd=${run_regenie_cmd} --tag="Step2" --instance-type "mem3_ssd1_v2_x32"\
   --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/" --brief -y
done

#################################
###GM 510 INTENSITY PHENOTYPES###
#################################


 dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.bgen" \
   -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.sample" -iin="${pheno_input}" -iin="${covar_input}" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_pred.list" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_1.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_2.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_3.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_4.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_5.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_6.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_7.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_8.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_9.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_10.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_11.loco.gz" \
   -icmd=${run_regenie_cmd} --tag="Step2" --instance-type "mem3_ssd1_v2_x32"\
   --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/" --brief -y
done



############################
############################
###CONVERT CHR X TO PLINK###
############################
############################

run_convert="plink2 --bgen ukb22828_cX_b0_v3.bgen ref-first --sample ukb22828_cX_b0_v3.sample\
	--make-bed --out ukb22828_cX_b0_v3_imputed_plink"

dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_cX_b0_v3.bgen" \
	-iin="${genotype_file_dir}/ukb22828_cX_b0_v3.sample" \
	-icmd=${run_convert} --tag="convert" --instance-type "mem3_ssd1_v2_x32"\
	 --destination="${project_out}/genotype_process/" --brief -y

#################################
#################################
###RUN REGENIE STEP 2 ON CHR X###
#################################
#################################

  run_regenie_X="regenie --step 2 --bed ukb22828_cX_b0_v3_imputed_plink --out assoc.cX\
    --phenoFile ${pheno_file} --covarFile ${covar_file}\
    --qt \
    --pred ${pheno_name}_pred.list --bsize 200\
    --threads 32 --gz"

###############################
###ORIGINAL BRAIN PHENOTYPES###
###############################

dx run swiss-army-knife -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.bed" \
	-iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.bim" \
	-iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.fam" \
	 -iin="${pheno_input}" -iin="${covar_input}" \
	   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_pred.list" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_1.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_2.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_3.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_4.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_5.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_6.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_7.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_8.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_9.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_10.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_11.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_12.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_13.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_14.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_15.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_16.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_17.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_18.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_19.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_20.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_21.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_22.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_23.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_24.loco.gz" \
   -icmd=${run_regenie_X} --tag="Step2X" --instance-type "mem3_ssd1_v2_x32"\
   --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/" --brief -y


#################################
###NORMALISED BRAIN PHENOTYPES###
#################################

dx run swiss-army-knife -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.bed" \
        -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.bim" \
        -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.fam" \
         -iin="${pheno_input}" -iin="${covar_input}" \
           -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_pred.list" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_1.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_2.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_3.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_4.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_5.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_6.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_7.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_8.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_9.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_10.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_11.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_12.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_13.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_14.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_15.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_16.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_17.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_18.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_19.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_20.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_21.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_22.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_23.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_24.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_25.loco.gz" \
   -icmd=${run_regenie_X} --tag="Step2X" --instance-type "mem3_ssd1_v2_x32"\
   --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/" --brief -y


##################
###HEIGHT GWAS ###
##################

dx run swiss-army-knife -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.bed" \
        -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.bim" \
        -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.fam" \
         -iin="${pheno_input}" -iin="${covar_input}" \
           -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_pred.list" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_1.loco.gz" \
   -icmd=${run_regenie_X} --tag="Step2X" --instance-type "mem3_ssd1_v2_x32"\
   --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/" --brief -y

###########################
###OB LR MEAN PHENOTYPES###
###########################

dx run swiss-army-knife -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.bed" \
        -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.bim" \
        -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.fam" \
         -iin="${pheno_input}" -iin="${covar_input}" \
           -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_pred.list" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_1.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_2.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_3.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_4.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_5.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_6.loco.gz" \
   -icmd=${run_regenie_X} --tag="Step2X" --instance-type "mem3_ssd1_v2_x32"\
   --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/" --brief -y

#################################
###GM 510 INTENSITY PHENOTYPES###
#################################

dx run swiss-army-knife -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.bed" \
        -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.bim" \
        -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.fam" \
         -iin="${pheno_input}" -iin="${covar_input}" \
           -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_pred.list" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_1.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_2.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_3.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_4.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_5.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_6.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_7.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_8.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_9.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_10.loco.gz" \
   -iin="${project_out}/genotype_process/${pheno_group}/regenie_step1/${pheno_name}_11.loco.gz" \
   -icmd=${run_regenie_X} --tag="Step2X" --instance-type "mem3_ssd1_v2_x32"\
   --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/" --brief -y

