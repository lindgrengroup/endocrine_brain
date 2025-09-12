#############################
###RUN GWAS FILTER RSCRIPT###
#############################

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=assoc.cX_$item.txt
out_name=$(basename "$file" .txt).filtered.txt
run_filter="Rscript gwas_filter.R --i $file --o ${out_name} --hw ukb22828_cX_b0_v3_hw_snp_pass.snplist --maf 0.01 --info 0.4"

dx run swiss-army-knife -iin="${script_loc}/gwas_filter.R" \
	-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/${file}" \
	-iin="${project_out}/genotype_process/${pheno_group}/ukb22828_cX_b0_v3_hw_snp_pass.snplist" \
	-icmd="${run_filter}" --tag="Filter" --instance-type "mem3_ssd1_v2_x32" \
	--destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/" --brief -y
done

##########################
###RUN GWAS PLOT SCRIPT###
##########################

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=assoc.cX_$item.filtered.txt
out_man=$(basename "$file" .txt).manhattan.png
out_qq=$(basename "$file" .txt).qqplot.png
run_plot="Rscript gwas_plot.R --i $file --man ${out_man} --qq ${out_qq}"

dx run swiss-army-knife -iin="${script_loc}/gwas_plot.R" \
	-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/${file}" \
	-icmd="${run_plot}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
	--destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/plots/" --brief -y
done

###########################################
###RUN PLOT FOR MERGED AUTOMOSOMES AND X###
###########################################

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=assoc.$item.regenie.merged.withX.filtered.txt
out_man=$(basename "$file" .txt).manhattan.png
out_qq=$(basename "$file" .txt).qqplot.png
run_plot="Rscript gwas_plot.R --i $file --man ${out_man} --qq ${out_qq}"

dx run swiss-army-knife -iin="${script_loc}/gwas_plot.R" \
        -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/${file}" \
        -icmd="${run_plot}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/plots/" --brief -y
done



########################
###RUN GWAS SEX PLOTS###
########################

# ONLY NEEDS TO BE RUN ONCE, NOT PER SEX GROUPING

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
name=assoc.$item.regenie.merged.filtered.txt
	out_miami=$(basename "$name" .txt).miami.png

	run_miami="Rscript gwas_miami.R --d /mnt/project/${project_out}/genotype_process/female_only/regenie_step2/filtered/$name \
	--h /mnt/project/${project_out}/genotype_process/male_only/regenie_step2/filtered/$name --miami ${out_miami}"

	dx run swiss-army-knife -iin="${script_loc}/gwas_miami.R" \
		-icmd="${run_miami}" --tag="miami" --instance-type "mem3_ssd1_v2_x32" \
 		--destination="${project_out}/genotype_process/plots/"  --brief -y
done

#RUN LEFT RIGHT MIAMI PLOTS
run_lr_volume="Rscript gwas_miami.R --d /mnt/project/${project_out}/genotype_process/all/regenie_step2/filtered/assoc.resid_OB_left.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.txt \
        --h /mnt/project/${project_out}/genotype_process/all/regenie_step2/filtered/assoc.resid_OB_right.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.txt --miami assoc.resid_OB_LR_comp.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.miami.png"
dx run swiss-army-knife -iin="${script_loc}/gwas_miami.R" \
                -icmd="${run_lr_volume}" --tag="miami" --instance-type "mem3_ssd1_v2_x32" \
                --destination="${project_out}/genotype_process/plots/"  --brief -y

