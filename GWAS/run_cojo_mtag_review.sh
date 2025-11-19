
genotype_file_dir="Bulk/Imputation/UKB imputation from genotype"
script_loc="scripts"
software_loc="software"

############################
###ALL BRAIN PHENOS - M&F###
############################

pheno_input="Phenotypes_MRI/regenie_input/phenotype_files/brain_morph_threshold_resid_norm_euro_all.phe"
project_out="data/brain_all"
pheno_file="brain_morph_threshold_resid_norm_euro_all.phe"
pheno_group="all"
pheno_name="brain_morph_threshold_resid_norm_euro_all"
covar_file="brain_morph_threshold_resid_norm_euro_all.cov"
covar_input="Phenotypes_MRI/regenie_input/covar_files/brain_morph_threshold_resid_norm_euro_all.cov"
snplist_qc="${project_out}/genotype_process/${pheno_group}/genotype_array_snps_qc_pass.snplist"

pheno_list=("HT_mtag_cojo_input.ma" "PG_mtag_cojo_input.ma" \
	"OB_mtag_cojo_input.ma" "HT_GM_mtag_cojo_input.ma")

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
pheno_list=("HT_mtag_cojo_input.ma" "PG_mtag_cojo_input.ma" \
	"OB_mtag_cojo_input.ma" "HT_GM_mtag_cojo_input.ma")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))

for item in ${pheno_list_a[@]}; do
	out_file=$(basename "$item" _cojo_input.ma)
	for chr in {1..22}; do
		run_cojo="chmod +x gcta-1.94.1; \
			  ./gcta-1.94.1 --bfile ukb22828_c${chr}_b0_v3_cojo_format --chr $chr --cojo-file $item \
			  --cojo-slct --thread-num 32 --out ${out_file}_c${chr}_cojo_result"
		
		dx run swiss-army-knife \
				-iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
		        -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/mtag_results/cojo/input/${item}" \
				-iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bim"\
				-iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.bed"\
				-iin="${project_out}/genotype_process/plink_format/ukb22828_c${chr}_b0_v3_cojo_format.fam"\
				-icmd="${run_cojo}" --tag="cojo" --instance-type "mem3_ssd1_v2_x32"\
		                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/mtag_results/cojo/output/" --brief --yes
	done
done

#run through all the phenotype files for chromosome X
pheno_list=("HT_mtag_cojo_input.ma" "PG_mtag_cojo_input.ma" \
	"OB_mtag_cojo_input.ma" "HT_GM_mtag_cojo_input.ma")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
	out_file=$(basename "$item" _cojo_input.ma)
		run_cojo="chmod +x gcta-1.94.1; \
			  ./gcta-1.94.1 --bfile ukb22828_cX_b0_v3_cojo_format --chr 23 --cojo-file $item \
			  --cojo-slct --thread-num 32 --out ${out_file}_cX_cojo_result"
		
		dx run swiss-army-knife \
				-iin="software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1" \
		               	-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/mtag_results/cojo/input/${item}" \
				-iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim"\
				-iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bed"\
				-iin="${project_out}/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.fam"\
				-icmd="${run_cojo}" --tag="cojo" --instance-type "mem3_ssd1_v2_x32"\
		                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/mtag_results/cojo/output/" --brief --yes
done


##########################################
###MERGE BY CHROMOSOME COJO FILES INC X###
##########################################

pheno_list=("HT_mtag" "PG_mtag" \
	"OB_mtag" "HT_GM_mtag")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}; do
        merge_cmd='out_name="'"${item}.merged.txt"'"
        cp '"/mnt/project/${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/mtag_results/cojo/output/${item}_c*_cojo_result.jma.cojo"' .
        echo -e "Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp\tn\tfreq_geno\tbJ\tbJ_se\tpJ\tLD_r" > $out_name
        files="./*jma.cojo"
        for f in $files
        do
                tail -n+2 $f | tr " " "\t" >> "${out_name}"
        done
        rm *jma.cojo'
        dx run swiss-army-knife -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/mtag_results/cojo/output/${item}_c1_cojo_result.log" \
                -icmd="${merge_cmd}" --tag="merge" --instance-type "mem1_ssd1_v2_x16" \
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/mtag_results/cojo/output/merged/" --brief --yes
done

#######################
###DOWNLOAD ALL COJO###
#######################

mkdir 
dx download "${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/cojo_results/merged/*"  
