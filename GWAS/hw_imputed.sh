#################################################
###RUN HW VIA PLINK ON IMPUTED DATA FOR FILTER###
#################################################

genotype_file_dir="Bulk/Imputation/UKB imputation from genotype/"

for chr in {1..22} 
do

	run_hw="plink2 --bgen ukb22828_c${chr}_b0_v3.bgen ref-first\
		--sample ukb22828_c${chr}_b0_v3.sample --keep ${pheno_file} --autosome\
		--hwe 1e-15 --write-snplist allow-dups --write-samples \
		--no-id-header --out ukb22828_c${chr}_b0_v3_hw_snp_pass"

	dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.bgen" \
	-iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.sample" -iin="${pheno_input}" \
  	-icmd=${run_hw} --tag="hw" --instance-type "mem3_ssd1_v2_x32"\
  	--destination="${project_out}/genotype_process/${pheno_group}" --brief -y
done

merge_cmd='out_file="'"ukb22828_c1_22_b0_v3_hw_snp_pass.snplist"'"
files='"/mnt/project/${project_out}/genotype_process/${pheno_group}/ukb22828_c*_b0_v3_hw_snp_pass.snplist"'
for file in $files
do
	cat $file >> "${out_file}"
done'

 dx run swiss-army-knife -iin="${project_out}/genotype_process/${pheno_group}/ukb22828_c1_b0_v3_hw_snp_pass.snplist" \
                -icmd="${merge_cmd}" --tag="merge" --instance-type "mem1_ssd1_v2_x16"\
                --destination="${project_out}/genotype_process/${pheno_group}/" --brief --yes

###############
###RUN FOR X###
###############

run_hw="plink2 --bfile ukb22828_cX_b0_v3_imputed_plink\
       --keep ${pheno_file} \
       --hwe 1e-15 --write-snplist allow-dups --write-samples\
	--no-id-header --out ukb22828_cX_b0_v3_hw_snp_pass"

dx run swiss-army-knife -iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.bed" \
	-iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.bim" \
	-iin="${project_out}/genotype_process/ukb22828_cX_b0_v3_imputed_plink.fam" \
        -iin="${pheno_input}" \
        -icmd=${run_hw} --tag="hw" --instance-type "mem3_ssd1_v2_x32"\
        --destination="${project_out}/genotype_process/${pheno_group}" --brief -y







