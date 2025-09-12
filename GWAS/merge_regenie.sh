######################
###MEGRE AUTOSOMSES###
######################

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
	pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}
do
	merge_cmd='out_file="'"assoc.$item.regenie.merged.txt"'"
	cp '"/mnt/project/${project_out}/genotype_process/${pheno_group}/regenie_step2/*$item.regenie.gz"' . 
	gunzip *.regenie.gz 
	echo -e "CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tINFO\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" > $out_file
	files="./*.regenie"
	for f in $files
	do
		tail -n+2 $f | tr " " "\t" >> "${out_file}"
	done
	rm *.regenie'
	dx run swiss-army-knife -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/assoc.c1_${item}.regenie.gz" \
   		-icmd="${merge_cmd}" --tag="merge" --instance-type "mem1_ssd1_v2_x16"\
   		--destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/" --brief --yes

done

#############################
###MERGE WITH X CHROMOSOME###
#############################

###THIS NEEDS TO BE DONE AFTER THE UNZIP AND THE ADDING OF THE INFO AND FILTERING###

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}
do
        merge_cmd='out_file="'"assoc.$item.regenie.merged.withX.filtered.txt"'"
        cp '"/mnt/project/${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/assoc.$item.regenie.merged.filtered.txt"' .
	cp '"/mnt/project/${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/assoc.cX_$item.filtered.txt"' .
        echo -e "CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tINFO\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" > $out_file
	#EDIT FROM HERE#
        files=('"assoc.$item.regenie.merged.filtered.txt"' '"assoc.cX_$item.filtered.txt"')
        for f in ${files[@]}
	do
                tail -n+2 $f | tr " " "\t" >> "${out_file}"
        done
        rm '"assoc.cX_$item.filtered.txt"'
	rm '"assoc.$item.regenie.merged.filtered.txt"''
        dx run swiss-army-knife -iin="${project_out}/genotype_process/ukb22418_c1_22X_v2_merged.bed" \
                -icmd="${merge_cmd}" --tag="merge" --instance-type "mem1_ssd1_v2_x16"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/filtered/" --brief --yes

done

#############################
###UNZIP CHR X FOR TESTING###
#############################

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}
do
     unzipx_cmd='cp '"/mnt/project/${project_out}/genotype_process/${pheno_group}/regenie_step2/assoc.cX_$item.regenie.gz"' .
	gunzip *.regenie.gz'
        dx run swiss-army-knife -iin="${project_out}/genotype_process/ukb22418_c1_22X_v2_merged.bed" \
                -icmd="${unzipx_cmd}" --tag="merge" --instance-type "mem1_ssd1_v2_x16"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/" --brief --yes

done

###################################
###MERGE X WITH INFO INFORMATION###
###################################

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
for item in ${pheno_list_a[@]}
do
    input=assoc.cX_${item}.regenie
    output=assoc.cX_${item}.txt
    x_info="Rscript xchrom_info_merge.R --i $input --o $output"
        dx run swiss-army-knife -iin="${project_out}/genotype_process/${pheno_group}/regenie_step2/assoc.cX_${item}.regenie" \
		-iin="${script_loc}/xchrom_info_merge.R" \
                -icmd="${x_info}" --tag="x_info" --instance-type "mem2_ssd2_v2_x32"\
                --destination="${project_out}/genotype_process/${pheno_group}/regenie_step2/" --brief --yes

done

