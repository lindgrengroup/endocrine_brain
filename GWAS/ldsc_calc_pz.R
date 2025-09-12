###############
###LIBRARIES###
###############

library(data.table)

##########
###DATA###
##########

ldsc_all <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/ldsc_all.csv", data.table=FALSE)

##############
###ANALYSIS###
##############

ldsc_all$V7 <- NULL

ldsc_all$intZ <- (ldsc_all$ldscIntercept-1)/ldsc_all$ldscinterceptSE
ldsc_all$intP <- pnorm(ldsc_all$intZ, lower.tail = FALSE)

write.csv(ldsc_all, "ldsc_all_zp.csv", row.names=FALSE, quote=FALSE)

system("dx upload ldsc_all_zp.csv --path data/brain_all/genotype_process/all/regenie_step2/filtered/ldsc_results/update_231125/ldsc_merge_alleles/")


