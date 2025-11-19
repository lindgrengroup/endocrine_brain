##############
###GET MTAG###
##############

get_mtag="git clone https://github.com/omeed-maghzian/mtag.git"

dx run swiss-army-knife  -iin="${genotype_file_dir}/ukb22828_c22_b0_v3.bgen" -icmd=${get_mtag} --tag="get_mtag" --instance-type "mem3_ssd1_v2_x96"\
                --destination="/mnt/project/software/mtag/" --brief --yes

################
###TEST LDSC ###
################

test_mtag="python /mnt/project/mnt/project/software/mtag/mtag.py -h"

mnt/project/software/mtag/mtag.py -h 

dx run swiss-army-knife  -iin="${genotype_file_dir}/ukb22828_c22_b0_v3.bgen" -icmd=${test_mtag} --tag="test_mtag" --instance-type "mem3_ssd1_v2_x96"\
                --destination="/mnt/project/"  -iimage="lifebitai/mtag:latest" --brief --yes



##############
###RUN MTAG###
##############

run_mtag="python /mnt/project/mnt/project/software/mtag/mtag.py \
	--sumstats "assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.ldsc_format.txt","assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.ldsc_format.txt","assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.ldsc_format.txt","assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.ldsc_format.txt" \
	--snp_name SNP \
	--a1_name A1 --a2_name A2 \
	--eaf_name freq \
	--z_name Z \
	--beta_name b --se_name se \
	--n_name N \
	--p_name p \
	--out "./HT_HTGM_PG_OB_mtag_output" \
	--n_min 0.0 "

dx run swiss-army-knife  -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/mtag_format/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.mtag_format.txt" \
	-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/mtag_format/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.mtag_format.txt" \
	-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/mtag_format/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.mtag_format.txt" \
	-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/mtag_format/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.mtag_format.txt" \
                -icmd="${run_mtag}" --tag="mtag" --instance-type "mem3_ssd1_v2_x32"\
		--destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/mtag_results/" -iimage="lifebitai/mtag:latest" --brief --yes



