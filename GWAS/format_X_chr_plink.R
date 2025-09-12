###############
###LIBRARIES###
###############

library(data.table)

##########
###DATA###
##########

X_chr <- fread("/mnt/project/data/brain_all/genotype_process/plink_format/ukb22828_cX_b0_v3_cojo_format.bim", data.table=FALSE)

##############
###ANALYSIS###
##############

X_chr$V1 <- 23

write.table(tmp, "ukb22828_cX_b0_v3_cojo_format.bim", col.names=FALSE, sep="\t", row.names=FALSE, quote=FALSE)
system("dx upload ukb22* --path data/brain_all/genotype_process/plink_format/")