#################################
###RUN GWAS PLOT SCRIPT - GENE###
#################################

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list=("${pheno_list[@]//[-=]/_}")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
pheno_list_b=("${pheno_list_a[@]//./_}")
for item in ${pheno_list_b[@]}; do
allfile=$item.gene.EUR.merged.all.txt
allout_man=$(basename "$allfile" .txt).manhattan.png
allout_qq=$(basename "$allfile" .txt).qqplot.png
ffile=$item.gene.EUR.merged.female.txt
fout_man=$(basename "$ffile" .txt).manhattan.png
fout_qq=$(basename "$ffile" .txt).qqplot.png
mfile=$item.gene.EUR.merged.male.txt
mout_man=$(basename "$mfile" .txt).manhattan.png
mout_qq=$(basename "$mfile" .txt).qqplot.png
all_run_plot="Rscript exome_plot_v7_SKAT.R --i $allfile --man ${allout_man} --qq ${allout_qq} --gwsl 6.7E-7"
female_run_plot="Rscript exome_plot_v7_SKAT.R --i $ffile --man ${fout_man} --qq ${fout_qq} --gwsl 6.7E-7"
male_run_plot="Rscript exome_plot_v7_SKAT.R --i $mfile --man ${mout_man} --qq ${mout_qq} --gwsl 6.7E-7"

dx run swiss-army-knife -iin="${script_loc}/exome_plot_v7_SKAT.R" \
	-iin="${output_location}/step2/merged/gene/${allfile}" \
	-icmd="${all_run_plot}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
	--destination="${output_location}/step2/merged/gene/skat/plots/" --brief -y

dx run swiss-army-knife -iin="${script_loc}/exome_plot_v7_SKAT.R" \
        -iin="${output_location}/step2/merged/gene/${ffile}" \
        -icmd="${female_run_plot}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${output_location}/step2/merged/gene/skat/plots/" --brief -y

dx run swiss-army-knife -iin="${script_loc}/exome_plot_v7_SKAT.R" \
        -iin="${output_location}/step2/merged/gene/${mfile}" \
        -icmd="${male_run_plot}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${output_location}/step2/merged/gene/skat/plots/" --brief -y

done

####################################
###RUN GWAS PLOT SCRIPT - SINGVAR###
####################################

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list=("${pheno_list[@]//[-=]/_}")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
pheno_list_b=("${pheno_list_a[@]//./_}")
for item in ${pheno_list_b[@]}; do
allfile=$item.EUR.merged.all.txt.singleAssoc.txt
allout_man=$(basename "$allfile" .txt).manhattan.png
allout_qq=$(basename "$allfile" .txt).qqplot.png
ffile=$item.EUR.merged.female.txt.singleAssoc.txt
fout_man=$(basename "$ffile" .txt).manhattan.png
fout_qq=$(basename "$ffile" .txt).qqplot.png
mfile=$item.EUR.merged.male.txt.singleAssoc.txt
mout_man=$(basename "$mfile" .txt).manhattan.png
mout_qq=$(basename "$mfile" .txt).qqplot.png
all_run_plot="Rscript exome_plot_sing_var_v7.R --i $allfile --man ${allout_man} --qq ${allout_qq} --gwsl 8E-9"
female_run_plot="Rscript exome_plot_sing_var_v7.R --i $ffile --man ${fout_man} --qq ${fout_qq} --gwsl 8E-9"
male_run_plot="Rscript exome_plot_sing_var_v7.R --i $mfile --man ${mout_man} --qq ${mout_qq} --gwsl 8E-9"

dx run swiss-army-knife -iin="${script_loc}/exome_plot_sing_var_v7.R" \
        -iin="${output_location}/step2/merged/singleAssoc/${allfile}" \
        -icmd="${all_run_plot}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${output_location}/step2/merged/singleAssoc/plots/" --brief -y

dx run swiss-army-knife -iin="${script_loc}/exome_plot_sing_var_v7.R" \
        -iin="${output_location}/step2/merged/singleAssoc/${ffile}" \
        -icmd="${female_run_plot}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${output_location}/step2/merged/singleAssoc/plots/" --brief -y

dx run swiss-army-knife -iin="${script_loc}/exome_plot_sing_var_v7.R" \
        -iin="${output_location}/step2/merged/singleAssoc/${mfile}" \
        -icmd="${male_run_plot}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${output_location}/step2/merged/singleAssoc/plots/" --brief -y

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

