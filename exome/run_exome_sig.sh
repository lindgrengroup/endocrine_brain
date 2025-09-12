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
allout=$(basename "$allfile" .txt).sig.csv
ffile=$item.gene.EUR.merged.female.txt
fout=$(basename "$ffile" .txt).sig.csv
mfile=$item.gene.EUR.merged.male.txt
mout=$(basename "$mfile" .txt).sig.csv
all_run_sig="Rscript exome_sig_gene_v7_skat.R --i $allfile --out ${allout} --gwsl 6.7E-7"
female_run_sig="Rscript exome_sig_gene_v7_skat.R --i $ffile --out ${fout} --gwsl 6.7E-7"
male_run_sig="Rscript exome_sig_gene_v7_skat.R --i $mfile --out ${mout} --gwsl 6.7E-7"

dx run swiss-army-knife -iin="${script_loc}/exome_sig_gene_v7_skat.R" \
	-iin="${output_location}/step2/merged/gene/${allfile}" \
	-icmd="${all_run_sig}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
	--destination="${output_location}/step2/merged/gene/skat/sig/" --brief -y

dx run swiss-army-knife -iin="${script_loc}/exome_sig_gene_v7_skat.R" \
        -iin="${output_location}/step2/merged/gene/${ffile}" \
        -icmd="${female_run_sig}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${output_location}/step2/merged/gene/skat/sig/" --brief -y

dx run swiss-army-knife -iin="${script_loc}/exome_sig_gene_v7_skat.R" \
        -iin="${output_location}/step2/merged/gene/${mfile}" \
        -icmd="${male_run_sig}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${output_location}/step2/merged/gene/skat/sig/" --brief -y

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
allout=$(basename "$allfile" .txt).sig.csv
ffile=$item.EUR.merged.female.txt.singleAssoc.txt
fout=$(basename "$ffile" .txt).sig.csv
mfile=$item.EUR.merged.male.txt.singleAssoc.txt
mout=$(basename "$mfile" .txt).sig.csv
all_run_sig="Rscript exome_sig_sing_var_v7.R --i $allfile --out ${allout} --gwsl 8E-9"
female_run_sig="Rscript exome_sig_sing_var_v7.R --i $ffile --out ${fout} --gwsl 8E-9"
male_run_sig="Rscript exome_sig_sing_var_v7.R --i $mfile --out ${mout} --gwsl 8E-9"

dx run swiss-army-knife -iin="${script_loc}/exome_sig_sing_var_v7.R" \
        -iin="${output_location}/step2/merged/singleAssoc/${allfile}" \
        -icmd="${all_run_sig}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${output_location}/step2/merged/singleAssoc/sig/" --brief -y

dx run swiss-army-knife -iin="${script_loc}/exome_sig_sing_var_v7.R" \
        -iin="${output_location}/step2/merged/singleAssoc/${ffile}" \
        -icmd="${female_run_sig}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${output_location}/step2/merged/singleAssoc/sig/" --brief -y

dx run swiss-army-knife -iin="${script_loc}/exome_sig_sing_var_v7.R" \
        -iin="${output_location}/step2/merged/singleAssoc/${mfile}" \
        -icmd="${male_run_sig}" --tag="Plot" --instance-type "mem3_ssd1_v2_x32" \
        --destination="${output_location}/step2/merged/singleAssoc/sig/" --brief -y

done

