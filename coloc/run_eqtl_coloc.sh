###############
###LOCATIONS###
###############

tissue_loc="data/resources/gtex"

#####################################
###RUN COLOC ACROSS ALL CELL TYPES###
#####################################

tissue_files=$(dx ls ${tissue_loc})
tissue_files_a=($(echo ${tissue_files} | tr " " "\n"))

pheno_list=("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered" "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" \
        "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered" "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))

for item in ${pheno_list_a[@]}; do

	for tissue in ${tissue_files_a[@]}; do
	tissue_name=$(basename "$tissue" .allpairs.txt.gz)
	gwas_file=${item}.hg38.txt
	gene_file=$(basename "$item" .regenie.merged.withX.filtered)_TSS_1Mb_genes.csv
	output_file=$(basename "$item" .regenie.merged.withX.filtered).${tissue_name}.coloc.txt
	output_sig=$(basename "$item" .regenie.merged.withX.filtered).${tissue_name}.coloc.sig.txt 	

	run_coloc="Rscript coloc_gtex_parse.R --i ${gwas_file} --g ${gene_file} --t ${tissue} --o ${output_file} --s ${output_sig}"
	dx run swiss-army-knife -iin="${script_loc}/coloc_gtex_parse.R" \
		-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/liftover/output/formatted/${gwas_file}" \
		-iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/coloc/genes/${gene_file}" \
		-iin="${tissue_loc}/${tissue}" \
		-icmd="${run_coloc}" --tag="eqtl" --instance-type "mem3_ssd1_v2_x32" \
		--destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/coloc/gtex_results/" \
		 -iimage="naotokubota/coloc-locuscomparer:1.0" --brief -y
	done
done



