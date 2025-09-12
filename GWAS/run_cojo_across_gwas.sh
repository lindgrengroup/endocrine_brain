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
file=assoc.$item.regenie.merged.filtered.txt
out_name=$(basename "$file" .txt).cojo_format.ma
run_format="Rscript gwas_cojo_format.R --i $file --o ${out_name}"

dx run swiss-army-knife -iin="${script_loc}/gwas_cojo_format.R" \
	-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/${file}"\
	-icmd="${run_format}" --tag="format" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/" --brief -y
done

###############################
###ESTABLISH COJO COND LISTS###
###############################

pheno_list=$(dx ls ${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/merged/)
delete=(mega_volume_results_table.csv)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done

pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
out_name=$(basename $item .txt).merged.cond.snplist

cond_prep="Rscript cojo_cond_pre_format.R --i $item --o ${out_name}"
dx run swiss-army-knife -iin="${script_loc}/cojo_cond_pre_format.R" \
        -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/merged/${item}" \
        -icmd="${cond_prep}" --tag="cond_prep" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/" --brief -y
done

#####################################
###RUN COJO ACROSS EACH CHROMOSOME###
#####################################

#Condition on PG vol sex-combined
pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
	out_file=${item}_PG_vol
	for chr in {1..22}; do
		run_cojo="chmod +x gcta-1.94.1; \
			  ./gcta-1.94.1 --bfile ukb22828_c${chr}_b0_v3_cojo_format --chr $chr --cojo-file $file \
			  --cojo-cond assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_c${chr}_cojo_result"
		
		dx run swiss-army-knife \
				-iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \ 
		               	-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
				-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist" \
				-iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bim"\
				-iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bed"\
				-iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.fam"\
				-icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
		                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
	done
done


pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_PG_vol
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $file \
                          --cojo-cond assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_cX_cojo_result"
		dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \ 
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done

#Condition on HT vol sex-combined
pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_HT_vol
        for chr in {1..22}; do
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_c${chr}_b0_v3_cojo_format --chr $chr --cojo-file $file \
                          --cojo-cond assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX.merged.cond.snplist --thread-num 32 --out ${out_file}_c${chr}_cojo_result"

                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done
done

pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_HT_vol
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $file \
                          --cojo-cond assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX.merged.cond.snplist --thread-num 32 --out ${out_file}_cX_cojo_result"
                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done

#Condition on LR vol sex-combined
pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_LR_vol
        for chr in {1..22}; do
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_c${chr}_b0_v3_cojo_format --chr $chr --cojo-file $file \
                          --cojo-cond assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_c${chr}_cojo_result"

                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done
done

pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_LR_vol
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $file \
                          --cojo-cond assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_cX_cojo_result"
                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done

#Condition on HTGM vol sex-combined
pheno_list=("assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_HTGM_vol
        for chr in {1..22}; do
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_c${chr}_b0_v3_cojo_format --chr $chr --cojo-file $file \
                          --cojo-cond assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_c${chr}_cojo_result"

                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done
done

pheno_list=("assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_HTGM_vol
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $file \
                          --cojo-cond assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_cX_cojo_result"
                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done


#Condition on PG vol female only
pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_PG_vol
        for chr in {1..22}; do
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_c${chr}_b0_v3_cojo_format --chr $chr --cojo-file $file \
                          --cojo-cond assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_c${chr}_cojo_result"

                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done
done


pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_PG_vol
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $file \
                          --cojo-cond assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_cX_cojo_result"
                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done

#Condition on HT vol female-only
pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_HT_vol
        for chr in {1..22}; do
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_c${chr}_b0_v3_cojo_format --chr $chr --cojo-file $file \
                          --cojo-cond assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_c${chr}_cojo_result"

                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done
done

pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_HT_vol
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $file \
                          --cojo-cond assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_cX_cojo_result"
                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done

#Condition on LR vol female-only
pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_LR_vol
        for chr in {1..22}; do
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_c${chr}_b0_v3_cojo_format --chr $chr --cojo-file $file \
                          --cojo-cond assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_c${chr}_cojo_result"

                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done
done

pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_LR_vol
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $file \
                          --cojo-cond assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_cX_cojo_result"
                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done

#Condition on HTGM vol female-only
pheno_list=("assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_HTGM_vol
        for chr in {1..22}; do
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_c${chr}_b0_v3_cojo_format --chr $chr --cojo-file $file \
                          --cojo-cond assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_c${chr}_cojo_result"

                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done
done

pheno_list=("assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_HTGM_vol
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $file \
                          --cojo-cond assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_cX_cojo_result"
                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done

#Condition on PG vol male-only
pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_PG_vol
        for chr in {1..22}; do
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_c${chr}_b0_v3_cojo_format --chr $chr --cojo-file $file \
                          --cojo-cond assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_c${chr}_cojo_result"

                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done
done

pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_PG_vol
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $file \
                          --cojo-cond assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_cX_cojo_result"
                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done


#Condition on HT vol male-only
pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_HT_vol
        for chr in {1..22}; do
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_c${chr}_b0_v3_cojo_format --chr $chr --cojo-file $file \
                          --cojo-cond assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_c${chr}_cojo_result"

                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done
done

pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_HT_vol
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $file \
                          --cojo-cond assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_cX_cojo_result"
                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done


#Condition on LR vol male-only
pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_LR_vol
        for chr in {1..22}; do
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_c${chr}_b0_v3_cojo_format --chr $chr --cojo-file $file \
                          --cojo-cond assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_c${chr}_cojo_result"

                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done
done

pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_LR_vol
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $file \
                          --cojo-cond assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_cX_cojo_result"
                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done


#Condition on HTGM vol male-only
pheno_list=("assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_HTGM_vol
        for chr in {1..22}; do
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_c${chr}_b0_v3_cojo_format --chr $chr --cojo-file $file \
                          --cojo-cond assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_c${chr}_cojo_result"

                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done
done

pheno_list=("assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.cojo_format.ma
        out_file=${item}_HTGM_vol
                run_cojo="chmod +x gcta-1.94.1; \
                          ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $file \
                          --cojo-cond assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist --thread-num 32 --out ${out_file}_cX_cojo_result"
                dx run swiss-army-knife \
                                -iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/${file}" \
                                -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_format/cojo_cond/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.merged.cond.snplist" \
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
                                -iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
                                -icmd="${run_cojo}" --tag="cojo_cond" --instance-type "mem3_ssd1_v2_x32"\
                                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/" --brief --yes
        done




####################################
###MERGE BY CHROMOSOME COJO FILES###
####################################
#Some will fail as the loop includes self correlated files which don't happen

pheno_list=("assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))

cond_list=("HT" "PG" "LR" "HTGM")
cond_list_a=($(echo ${cond_list} | tr " " "\n"))

for item in ${pheno_list_a[@]}; do
	for cond in ${cond_list_a[@]}; do
	merge_cmd='out_name="'"${item}.${cond}_vol.cond.merged.txt"'"
	cp '"/mnt/project/${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/${item}_${cond}_vol_c*_cojo_result.cma.cojo"' .
	#files=$(dx ls '${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/${item}_${cond}_vol_c*_cojo_result.cma.cojo')
	echo -e "Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp\tn\tfreq_geno\tbC\tbC_se\tpC" > $out_name
	files="./*cma.cojo"
	for f in $files
	do
		tail -n+2 $f | tr " " "\t" >> "${out_name}"
	done
	rm *cma.cojo
	'
	dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_c22_b0_v3.bgen" \
		-icmd="${merge_cmd}" --tag="merge" --instance-type "mem1_ssd1_v2_x16" \
		--destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/merged/" --brief --yes
done
done


##############################
###FILTER MERGED COJO FILES###
##############################

files=$(dx ls ${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/merged/)
files_a=($(echo ${files} | tr " " "\n"))
for f in ${files_a[@]};
do
	out_name=$(basename $f .txt)_filtered.txt
	filter_cmd="Rscript cojo_cond_filter.R --i $f --o ${out_name}"
	
	dx run swiss-army-knife -iin="${script_loc}/cojo_cond_filter.R" \
		-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/merged/$f" \
		-icmd="${filter_cmd}" --tag="filted_cojo" --instance-type "mem1_ssd1_v2_x16" \
		--destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/cojo_cond/merged/filtered/" --brief --yes
done
