##################################
###CONVERT BGEN TO PLINK FORMAT###
##################################

for chr in {{1..22},X}; do
	run_convert="plink2 --bgen ukb22828_c${chr}_b0_v3.bgen ref-first\
		--sample ukb22828_c${chr}_b0_v3.sample --keep brain_morph_threshold_resid_norm_euro_all.phe --autosome\
		--make-bed --out ukb22828_c${chr}_b0_v3_cojo_format"
	
	dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.bgen"\
	-iin="${genotype_file_dir}/ukb22828_c${chr}_b0_v3.sample" -iin="${pheno_input}"\
	-icmd=${run_convert} --tag="convert"  --instance-type "mem3_ssd1_v2_x32" \ 
	--destination="${project_out}/genotype_process/plink_format/" --brief -y
done

run_convert_x="plink2 --bgen ukb22828_cX_b0_v3.bgen ref-first\
                --sample ukb22828_cX_b0_v3.sample --keep brain_morph_threshold_resid_norm_euro_all.phe \
                --make-bed --out ukb22828_cX_b0_v3_cojo_format"

        dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_cX_b0_v3.bgen"\
        -iin="${genotype_file_dir}/ukb22828_cX_b0_v3.sample" -iin="${pheno_input}"\
        -icmd=${run_convert_x} --tag="convert"  --instance-type "mem3_ssd1_v2_x32" \
        --destination="${project_out}/genotype_process/plink_format/" --brief -y


#############################################
###RUN FORMATTING OF GWAS RESULTS FOR COJO###
#############################################

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=assoc.$item.regenie.merged.withX.filtered.txt
out_name=$(basename "$file" .txt).cojo_format.ma
run_format="Rscript gwas_cojo_format.R --i $file --o ${out_name}"

dx run swiss-army-knife -iin="${script_loc}/gwas_cojo_format.R" \
	-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/${file}"\
	-icmd="${run_format}" --tag="format" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/" --brief -y
done

#####################################
###RUN COJO ACROSS EACH CHROMOSOME###
#####################################

#Install software
run_get_cojo="wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip ;\
        unzip gcta-1.94.1-linux-kernel-3-x86_64.zip"
dx run swiss-army-knife  -iin="${genotype_file_dir}/ukb22828_c22_b0_v3.bgen" \
                -icmd="${run_get_cojo}" --tag="get_cojo" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${software_loc}/" --brief --yes

#run through all the phenotype files and per chromosome
pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=assoc.$item.regenie.merged.withX.filtered.cojo_format.ma
	out_file=$(basename "$file" .cojo_format.ma)
		run_cojo="chmod +x gcta-1.94.1; \
			  ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $file \
			  --cojo-slct --thread-num 32 --out ${out_file}_cX_cojo_result"
		
		dx run swiss-army-knife \
				-iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \ 
		               	-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
				-iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
				-iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
				-iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
				-icmd="${run_cojo}" --tag="cojo" --instance-type "mem3_ssd1_v2_x32"\
		                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/" --brief --yes
done

######################################
###MERGE X WITH AUTOSOME COJO FILES###
######################################

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
	merge_cmd='out_name="'"assoc.$item.regenie.merged.filtered.merged.withX.txt"'"
	cp '"/mnt/project/${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/assoc.$item.regenie.merged.withX.filtered_cX_cojo_result.jma.cojo"' .
	cp '"/mnt/project/${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/merged/assoc.$item.regenie.merged.filtered.merged.txt"' .
	echo -e "Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp\tn\tfreq_geno\tbJ\tbJ_se\tpJ\tLD_r" > $out_name
	files=('"assoc.$item.regenie.merged.withX.filtered_cX_cojo_result.jma.cojo"' '"assoc.$item.regenie.merged.filtered.merged.txt"')
	for f in ${files[@]}
	do
		tail -n+2 $f | tr " " "\t" >> "${out_name}"
	done
	rm '"assoc.$item.regenie.merged.withX.filtered_cX_cojo_result.jma.cojo"'
	rm '"assoc.$item.regenie.merged.filtered.merged.txt"''
	dx run swiss-army-knife -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/assoc.$item.regenie.merged.withX.filtered_cX_cojo_result.log" \
		-icmd="${merge_cmd}" --tag="merge" --instance-type "mem1_ssd1_v2_x16" \
		--destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/merged/" --brief --yes
done

#######################
###DOWNLOAD ALL COJO###
#######################

mkdir 
dx download "${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/merged/*"  
