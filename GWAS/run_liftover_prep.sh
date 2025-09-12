pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")

pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.txt
out_name=$item.liftover_format_hg19.bed

liftover_format="Rscript liftover_prep_parse.R --i $file --o ${out_name}"
dx run swiss-army-knife -iin="${script_loc}/liftover_prep_parse.R" \
        -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/${file}" \
        -icmd="${liftover_format}" --tag="liftover_format" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/liftover/input/" --brief -y
done


######################
###INSTALL LIFTOVER###
######################

run_get_liftover="wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ "
dx run swiss-army-knife  -iin="${genotype_file_dir}/ukb22828_c22_b0_v3.bgen" \
                -icmd="${run_get_liftover}" --tag="get_liftover" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${software_loc}/" --brief --yes

run_get_chain="wget https://hgdownload2.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"

dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_c22_b0_v3.bgen" \
		-icmd="${run_get_chain}" --tag="get_chain" --instance-type "mem3_ssd1_v2_x32"\
		--destination="data/resources/" --brief --yes


###################
###TEST LIFTOVER###
###################

test_liftover="liftOver -h"

dx run swiss-army-knife  -iin="${genotype_file_dir}/ukb22828_c22_b0_v3.bgen" -icmd=${test_liftover} --tag="test_liftover" --instance-type "mem3_ssd1_v2_x96"\
                --destination="/mnt/project/" -iimage="quay.io/biocontainers/ucsc-liftover:357--h446ed27_4" --brief --yes


##################
###RUN LIFTOVER###
##################

pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")

pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.liftover_format_hg19.bed
output=$item.liftover_markers_hg38
liftover_command="liftOver ${file} hg19ToHg38.over.chain.gz ${output}.bed ${output}_unlifted.bed"

dx run swiss-army-knife -iin="/data/resources/hg19ToHg38.over.chain.gz" \
        -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/liftover/input/${file}" \
        -icmd="${liftover_command}" --tag="liftover_command" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/liftover/output/markers/" \
	-iimage="quay.io/biocontainers/ucsc-liftover:357--h446ed27_4" --brief -y
done

##################
###FORMAT FILES###
##################

pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")

pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.txt
liftover_markers=$item.liftover_markers_hg38.bed
output=$item.hg38.txt

liftover_post="Rscript liftover_post_format.R --g $file --l ${liftover_markers} --o $output"
dx run swiss-army-knife -iin="${script_loc}/liftover_post_format.R" \
        -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/${file}" \
	-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/liftover/output/markers/${liftover_markers}" \
        -icmd="${liftover_post}" --tag="liftover_post" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/liftover/output/formatted/" --brief -y
done
