#########################
###MERGE BY CHROMOSOME###
#########################

###GENE WISE###

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list=("${pheno_list[@]//[-=]/_}")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
pheno_list_b=("${pheno_list_a[@]//./_}")
for item in ${pheno_list_b[@]}
do
	merge_cmd='out_all="'"$item.gene.EUR.merged.all.txt"'"
	out_female="'"$item.gene.EUR.merged.female.txt"'"
	out_male="'"$item.gene.EUR.merged.male.txt"'"
	cp '"/mnt/project/${output_location}/step2/*${item}_EUR_all.txt.gz"' .
	gunzip *.txt.gz
	echo -e "Region\tGroup\tmax_MAF\tPvalue\tPvalue_Burden\tPvalue_SKAT\tBETA_Burden\tSE_Burden\tMAC\tNumber_rare\tNumber_ultra_rare" > $out_all
	files="./*_EUR_all.txt"
	for f in $files
	do
	tail -n+2 $f | tr " " "\t" >> "${out_all}"
	done
	rm *EUR_all.txt
	cp '"/mnt/project/${output_location}/step2/*${item}_EUR_female.txt.gz"' .
	gunzip *.txt.gz
        echo -e "Region\tGroup\tmax_MAF\tPvalue\tPvalue_Burden\tPvalue_SKAT\tBETA_Burden\tSE_Burden\tMAC\tNumber_rare\tNumber_ultra_rare" > $out_female
        files="./*EUR_female.txt"
        for f in $files
        do
        tail -n+2 $f | tr " " "\t" >> "${out_female}"
        done
        rm *EUR_female.txt
        cp '"/mnt/project/${output_location}/step2/*${item}_EUR_male.txt.gz"' .
	gunzip *.txt.gz
        echo -e "Region\tGroup\tmax_MAF\tPvalue\tPvalue_Burden\tPvalue_SKAT\tBETA_Burden\tSE_Burden\tMAC\tNumber_rare\tNumber_ultra_rare" > $out_male
        files="./*EUR_male.txt"
        for f in $files
        do
        tail -n+2 $f | tr " " "\t" >> "${out_male}"
        done 
        rm *EUR_male.txt'
	dx run swiss-army-knife -iin="${pheno_input}" \
                -icmd="${merge_cmd}" --tag="merge" --instance-type "mem3_ssd1_v2_x16"\
                --destination="${output_location}/step2/merged/gene" --brief --yes
done


###SINGLE VARAIANT###

pheno_list=$(dx head -n 1 ${pheno_input})
delete=(FID IID)
for del in ${delete[@]}
do
        pheno_list=("${pheno_list[@]/$del}")
done
pheno_list=("${pheno_list[@]//[-=]/_}")
pheno_list_a=($(echo ${pheno_list} | tr " " "\n"))
pheno_list_b=("${pheno_list_a[@]//./_}")
for item in ${pheno_list_b[@]}
do
        merge_cmd='out_all="'"$item.EUR.merged.all.txt.singleAssoc.txt"'"
        out_female="'"$item.EUR.merged.female.txt.singleAssoc.txt"'"
        out_male="'"$item.EUR.merged.male.txt.singleAssoc.txt"'"
        cp '"/mnt/project/${output_location}/step2/*${item}_EUR_all.txt.singleAssoc.txt.gz"' .
        gunzip *.txt.gz
	echo -e "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\tMissingRate\tBETA\tSE\tTstat\tvar\tp.value\tN" > $out_all
        files="./*EUR_all.txt.singleAssoc.txt"
        for f in $files
        do
        tail -n+2 $f | tr " " "\t" >> "${out_all}"
        done
        rm *EUR_all.txt.singleAssoc.txt
        cp '"/mnt/project/${output_location}/step2/*${item}_EUR_female.txt.singleAssoc.txt.gz"' .
        gunzip *.txt.gz
	echo -e "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\tMissingRate\tBETA\tSE\tTstat\tvar\tp.value\tN" > $out_female
        files="./*EUR_female.txt.singleAssoc.txt"
        for f in $files
        do
        tail -n+2 $f | tr " " "\t" >> "${out_female}"
        done
        rm *EUR_female.txt.singleAssoc.txt
        cp '"/mnt/project/${output_location}/step2/*${item}_EUR_male.txt.singleAssoc.txt.gz"' .
        gunzip *.txt.gz
	echo -e "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\tMissingRate\tBETA\tSE\tTstat\tvar\tp.value\tN" > $out_male
        files="./*EUR_male.txt.singleAssoc.txt"
        for f in $files
        do
        tail -n+2 $f | tr " " "\t" >> "${out_male}"
        done
        rm *EUR_male.txt.singleAssoc.txt'
        dx run swiss-army-knife -iin="${pheno_input}" \
                -icmd="${merge_cmd}" --tag="merge" --instance-type "mem3_ssd1_v2_x16"\
                --destination="${output_location}/step2/merged/singleAssoc" --brief --yes
done

