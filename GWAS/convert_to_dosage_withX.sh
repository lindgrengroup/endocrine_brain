###############################
###PREP FILES FOR CONVERSION###
###############################

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
cojo_result=assoc.$item.regenie.merged.withX.filtered.merged.txt
pheno_name=assoc.$item.regenie.merged.withX.filtered.merged
run_dosage_prep="Rscript dosage_prep_tmp.R --g ukb22828_c1_b0_v3.sample --p ${pheno_file} --c ${cojo_result} --s ${pheno_name}_incl_samples --r ${pheno_name}_incl_rsids"

dx run swiss-army-knife -iin="${script_loc}/dosage_prep_tmp.R" \
		-iin="${genotype_file_dir}/ukb22828_c1_b0_v3.sample" \
		-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/merged/${cojo_result}" \
		-iin="${pheno_input}" \
                -icmd="${run_dosage_prep}" --tag="dosage_prep" --instance-type "mem3_ssd1_v2_x32" \
                --destination="${project_out}/genotype_process/${pheno_group}/dosage/"  --brief -y
done

############################
### RUN CONVERT TO BIMBAM###
############################
pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
pheno_name=assoc.$item.regenie.merged.withX.filtered.merged

for chr in {1..22}; do
	run_convert="qctool -g ukb22828_c${chr}_b0_v3.bgen \
	-og ${pheno_name}_chr${chr}.dosage \
	-s ukb22828_c${chr}_b0_v3.sample \
	-incl-samples ${pheno_name}_incl_samples \
	-incl-rsids ${pheno_name}_incl_rsids"

dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.bgen" \
	-iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.sample" \
	-iin="data/brain_all/genotype_process/${pheno_group}/dosage/${pheno_name}_incl_samples" \
	-iin="data/brain_all/genotype_process/${pheno_group}/dosage/${pheno_name}_incl_rsids" \ 
	-icmd="${run_convert}" --tag="convert" --instance-type "mem3_ssd1_v2_x32"\
	--destination="data/brain_all/genotype_process/${pheno_group}/dosage/soft_calls/" --brief -y
done
done
done

        run_convert="qctool -g ukb22828_cX_b0_v3.bgen \
        -og ${pheno_name}_chr23.dosage \
        -s ukb22828_cX_b0_v3.sample \
        -incl-samples ${pheno_name}_incl_samples \
        -incl-rsids ${pheno_name}_incl_rsids"

dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_cX_b0_v3.bgen" \
        -iin="${genotype_file_dir}/ukb22828_cX_b0_v3.sample" \
        -iin="data/brain_all/genotype_process/${pheno_group}/dosage/${pheno_name}_incl_samples" \
        -iin="data/brain_all/genotype_process/${pheno_group}/dosage/${pheno_name}_incl_rsids" \
        -icmd="${run_convert}" --tag="convert" --instance-type "mem3_ssd1_v2_x32"\
        --destination="data/brain_all/genotype_process/${pheno_group}/dosage/soft_calls/" --brief -y


################################
###RUN DOSAGE POST FORMATTING###
################################
pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
pheno_name=assoc.$item.regenie.merged.withX.filtered.merged

run_dosage_format="Rscript dosage_format.R --i /mnt/project/data/brain_all/genotype_process/${pheno_group}/dosage/soft_calls/ --n ${pheno_name}_ --d  ${pheno_name}_dosage --m ${pheno_name}_metadata"

dx run swiss-army-knife -iin="${script_loc}/dosage_format.R" \
                -icmd="${run_dosage_format}" --tag="dosage_format" --instance-type "mem3_ssd1_v2_x32" \
                --destination="${project_out}/genotype_process/${pheno_group}/dosage/soft_calls/merged/"  --brief -y
done

run_dosage_format="Rscript dosage_format_withX.R --i /mnt/project/data/brain_all/genotype_process/${pheno_group}/dosage/soft_calls/ --n ${pheno_name}_ --d  ${pheno_name}_dosage --m ${pheno_name}_metadata"

dx run swiss-army-knife -iin="${script_loc}/dosage_format_withX.R" \
                -icmd="${run_dosage_format}" --tag="dosage_format" --instance-type "mem3_ssd1_v2_x32" \
                --destination="${project_out}/genotype_process/${pheno_group}/dosage/soft_calls/merged/"  --brief -y


###########################
###CONVERT TO HARD CALL ###
###########################

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
cojo_result=assoc.$item.regenie.merged.withX.filtered.merged.txt
pheno_name=assoc.$item.regenie.merged.withX.filtered.merged

for chr in {1..22}; do
                run_hard_call_p1="plink2 --bgen ukb22828_c${chr}_b0_v3.bgen ref-first\
                        --sample ukb22828_c${chr}_b0_v3.sample --keep ${pheno_file} --autosome\
                        --extract ${pheno_name}_incl_rsids \
                        --make-pgen erase-dosage \
                        --out ukb22828_c${chr}_b0_v3_hard_call_pgen

		plink2 --pfile ukb22828_c${chr}_b0_v3_hard_call_pgen\
                        --keep ${pheno_file} --autosome\
                        --extract ${pheno_name}_incl_rsids \
                        --hard-call-threshold 0.1 \
                        --export A \
                        --out ${pheno_name}_ukb22828_c${chr}_b0_v3_hard_call

		rm ukb22828_c${chr}_b0_v3_hard_call_pgen.*"

                dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.bgen" \
                        -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.sample"\
                        -iin="${pheno_input}"\
                        -iin="data/brain_all/genotype_process/${pheno_group}/dosage/${pheno_name}_incl_rsids" \
                        -icmd="${run_hard_call_p1}" --tag="hard_call" --instance-type "mem3_ssd1_v2_x32" \
                --destination="${project_out}/genotype_process/${pheno_group}/dosage/hard_calls/"  --brief -y
done
done

                run_hard_call_p1="plink2 --bgen ukb22828_cX_b0_v3.bgen ref-first\
                        --sample ukb22828_cX_b0_v3.sample --keep ${pheno_file}\
                        --extract ${pheno_name}_incl_rsids \
                        --make-pgen erase-dosage \
                        --out ukb22828_cX_b0_v3_hard_call_pgen

                plink2 --pfile ukb22828_cX_b0_v3_hard_call_pgen\
                        --keep ${pheno_file} \
                        --extract ${pheno_name}_incl_rsids \
                        --hard-call-threshold 0.1 \
                        --export A \
                        --out ${pheno_name}_ukb22828_c23_b0_v3_hard_call

                rm ukb22828_cX_b0_v3_hard_call_pgen.*"

 dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_cX_b0_v3.bgen" \
                        -iin="${genotype_file_dir}/ukb22828_cX_b0_v3.sample"\
                        -iin="${pheno_input}"\
                        -iin="data/brain_all/genotype_process/${pheno_group}/dosage/${pheno_name}_incl_rsids" \
                        -icmd="${run_hard_call_p1}" --tag="hard_call" --instance-type "mem3_ssd1_v2_x32" \
                --destination="${project_out}/genotype_process/${pheno_group}/dosage/hard_calls/"  --brief -y

#####################################
###HARD RUN DOSAGE POST FORMATTING###
#####################################

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
pheno_name=assoc.$item.regenie.merged.withX.filtered.merged

run_dosage_format="Rscript hard_dosage_format.R --i /mnt/project/data/brain_all/genotype_process/${pheno_group}/dosage/hard_calls/ --n ${pheno_name} --d  ${pheno_name}_dosage"

dx run swiss-army-knife -iin="${script_loc}/hard_dosage_format.R" \
                -icmd="${run_dosage_format}" --tag="dosage_format" --instance-type "mem3_ssd1_v2_x32" \
                --destination="${project_out}/genotype_process/${pheno_group}/dosage/hard_calls/merged/"  --brief -y

done

		
####################
###ENCRYPT DOSAGE###
####################

dx run app-cloud_workstation --ssh
unset container-GQ7b248JBJ6G4xxqq83XBZF9
dx cd project-GGqJK68J6fpyp6qQGJxq4853:
dx download  project-GGqJK68J6fpyp6qQGJxq4853:/data/brain_all/genotype_process/all/dosage/hard_calls/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged_dosage
gpg -c assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged_dosage
dx upload assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged_dosage.gpg --path "project-GGqJK68J6fpyp6qQGJxq4853:/data/brain_all/genotype_process/all/dosage/hard_calls/merged/encrypted/"

