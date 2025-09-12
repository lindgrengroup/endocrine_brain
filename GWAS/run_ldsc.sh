####################################
###RUN LDSC FORMAT ALL PHENOTYPES###
####################################

pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
	"assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")

pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.txt
out_name=$item.ldsc_format.txt
ldsc_format="Rscript gwas_ldsc_format.R --i $file --o ${out_name}"
dx run swiss-army-knife -iin="${script_loc}/gwas_ldsc_format.R" \
        -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/${file}" \
        -icmd="${ldsc_format}" --tag="ldsc_format" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_format/update_231125/" --brief -y
done


################
###TEST LDSC ###
################

test_ldsc="ls ;\
	/ldsc/ldsc.py -h ;\
	/ldsc/munge_sumstats.py -h"

dx run swiss-army-knife  -iin="${genotype_file_dir}/ukb22828_c22_b0_v3.bgen" -icmd=${test_ldsc} --tag="test_ldsc" --instance-type "mem3_ssd1_v2_x96"\
                --destination="/mnt/project/" -iimage="zijingliu/ldsc:latest" --brief --yes


###############################
###DOWNLOAD REFERENCE PANELS###
###############################

get_refs="wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2 ;\
                tar -jxvf eur_w_ld_chr.tar.bz2"

dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_c22_b0_v3.bgen" -icmd=${get_refs} --tag="refs"  --instance-type "mem3_ssd1_v2_x32"\
        --destination="software/" --brief --yes

get_refs2="wget https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/w_hm3.snplist"

dx run swiss-army-knife -iin="${genotype_file_dir}/ukb22828_c22_b0_v3.bgen" -icmd=${get_refs2} --tag="refs"  --instance-type "mem3_ssd1_v2_x32"\
        --destination="software/" --brief --yes


################################
###RUN SUMSTATS SUMSTATS STEP###
################################

pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")

pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=$item.ldsc_format.txt

out_name=${item}_sumstats
run_sumstats="/ldsc/munge_sumstats.py \
	--sumstats $file \
	--out ${out_name}\
	--merge-alleles /mnt/project/software/w_hm3.snplist "

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_format/update_231125/$file" \
                -icmd="${run_sumstats}" --tag="sumstats" --instance-type "mem3_ssd1_v2_x32"\
		--destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/" -iimage="zijingliu/ldsc:latest" --brief --yes
done

########################
###RUN MAIN LDSC STEP###
########################

pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")

pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
file=${item}_sumstats.sumstats.gz

out_name=assoc.${item}_ldsc
run_ldsc="/ldsc/ldsc.py \
	--h2 $file \
	--ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
	--w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
	--out ${out_name}"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/$file" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/" -iimage="zijingliu/ldsc:latest" --brief --yes
done

 
#########################
###GENETIC CORRELATION###
#########################

###RUN GENETIC CORRELATION###
###ALL HORMONES###

###TESTOSTERONE
run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Testosterone_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_test_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Testosterone_sex_comb_EUR_for_ldsc.sumstats.gz" \
		-icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_megre_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Testosterone_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_test_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Testosterone_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Testosterone_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_test_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Testosterone_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Testosterone_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_test_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Testosterone_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes


###PROGESTERONE

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Progesterone_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_prog_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Progesterone_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Progesterone_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_prog_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Progesterone_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Progesterone_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_prog_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Progesterone_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Progesterone_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_prog_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Progesterone_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes


###LEUTENISING HORMONE

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,LH_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_LH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/LH_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,LH_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_LH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/LH_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,LH_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_LH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/LH_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,LH_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_LH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/LH_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes


###FOLICLE STIMULATING HORMONE

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,FSH_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_fsh_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/FSH_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,FSH_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_fsh_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/FSH_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,FSH_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_fsh_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/FSH_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,FSH_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_FSH_ldsc"
        
dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/FSH_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

###OESTRODIOL

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Oestradiol_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_oest_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Oestradiol_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Oestradiol_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_oest_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Oestradiol_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Oestradiol_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_oest_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Oestradiol_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Oestradiol_sex_comb_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_oest_ldsc"
        
dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Oestradiol_sex_comb_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

###female_infertility1

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis1_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_fem_infert1_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
		-iin="data/hormone_fertility_update_231125/female_infertility_analysis1_EUR_for_ldsc.sumstats.gz" \
		-icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
		--destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis1_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_fem_infert1_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis1_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis1_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_fem_infert1_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis1_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis1_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_fem_infert1_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis1_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

###female_infertility2

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis2_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_fem_infert2_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis2_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis2_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_fem_infert2_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis2_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis2_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_fem_infert2_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis2_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis2_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_fem_infert2_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis2_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes


###female_infertility3

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis3_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_fem_infert3_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis3_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis3_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_fem_infert3_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis3_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis3_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_fem_infert3_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis3_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis3_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_fem_infert3_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis3_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

###female_infertility4

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis4_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_fem_infert4_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis4_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis4_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_fem_infert4_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis4_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis4_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_fem_infert4_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis4_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis4_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_fem_infert4_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis4_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

###female_infertility5

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis5_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_fem_infert5_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis5_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis5_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_fem_infert5_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis5_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis5_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_fem_infert5_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis5_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,female_infertility_analysis5_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_fem_infert5_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/female_infertility_analysis5_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/female_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

###male_infertility

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,male_infertility_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_male_infert_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/male_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/male_infertility_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/male_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,male_infertility_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_male_infert_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/male_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/male_infertility_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/male_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,male_infertility_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_male_infert_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/male_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/male_infertility_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/male_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,male_infertility_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_male_infert_ldsc"

dx run swiss-army-knife -iin="${project_out}/genotype_process/male_only/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/male_infertility_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/male_only/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes


###FEMALE HORMONES###

###TESTOSTERONE
run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Testosterone_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_test_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Testosterone_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Testosterone_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_test_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Testosterone_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Testosterone_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_test_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Testosterone_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Testosterone_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_test_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Testosterone_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes


###PROGESTERONE

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Progesterone_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_prog_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Progesterone_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Progesterone_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_prog_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Progesterone_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Progesterone_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_prog_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Progesterone_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Progesterone_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_prog_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Progesterone_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes


###LEUTENISING HORMONE

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,LH_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_LH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/LH_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,LH_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_LH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/LH_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,LH_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_LH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/LH_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,LH_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_LH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/LH_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

###FOLICLE STIMULATING HORMONE

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,FSH_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_fsh_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/FSH_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,FSH_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_fsh_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/FSH_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,FSH_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_fsh_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/FSH_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,FSH_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_FSH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/FSH_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

###OESTRODIOL

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Oestradiol_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_oest_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Oestradiol_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Oestradiol_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_oest_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Oestradiol_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Oestradiol_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_oest_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Oestradiol_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Oestradiol_F_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_oest_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Oestradiol_F_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

###MALE HORMONES

###TESTOSTERONE
run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Testosterone_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_test_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Testosterone_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Testosterone_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_test_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Testosterone_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Testosterone_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_test_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Testosterone_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Testosterone_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_test_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Testosterone_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

###LEUTENISING HORMONE

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,LH_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_LH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/LH_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,LH_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_LH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/LH_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,LH_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_LH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/LH_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,LH_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_LH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/LH_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

###FOLICLE STIMULATING HORMONE

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,FSH_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_fsh_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/FSH_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,FSH_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_fsh_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/FSH_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,FSH_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_fsh_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/FSH_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,FSH_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_FSH_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/FSH_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

###OESTRODIOL

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Oestradiol_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out PGvol_oest_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Oestradiol_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Oestradiol_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTvol_oest_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Oestradiol_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Oestradiol_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out OBLRvol_oest_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Oestradiol_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

run_ldsc="/ldsc/ldsc.py \
        --rg assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz,Oestradiol_M_EUR_for_ldsc.sumstats.gz \
        --ref-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --w-ld-chr /mnt/project/software/eur_w_ld_chr/ \
        --out HTGMvol_oest_ldsc"

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/sumstats_merge_alleles/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sumstats.sumstats.gz" \
                -iin="data/hormone_fertility_update_231125/Oestradiol_M_EUR_for_ldsc.sumstats.gz" \
                -icmd="${run_ldsc}" --tag="ldsc" --instance-type "mem3_ssd1_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/gc/" -iimage="zijingliu/ldsc:latest" --brief --yes

