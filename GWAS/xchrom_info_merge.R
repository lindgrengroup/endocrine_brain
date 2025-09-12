###############
###LIBRARIES###
###############

install.packages("data.table")
library(data.table)
install.packages("docopt")
library(docopt)
install.packages("dplyr")
library("dplyr")

########################
###PARSING PARAMETERS###
########################

doc <- "Usage:
    xchrom_info_merge.R --i <input> --o <output>
Arguments:
    --i Path to x chromosome regenie output [default: NULL]
    --o File to be written with INFO added and aligned with other regenie results [default: NULL]
"

opts <- docopt(doc)

##########
###DATA###
##########

x_mfi <- fread("/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_cX_b0_v3.mfi.txt", data.table=FALSE)
gwas <- fread(opts$i, data.table=FALSE)

##############
###ANALYSIS###
##############

x_mfi <- rename(x_mfi, alt_id=V1, ID=V2, POS=V3, ALLELE0=V4, ALLELE1=V5, MAF=V6, ALTAL=V7, INFO=V8)
x_mfi_select <- x_mfi[,c("ID", "ALLELE0", "ALLELE1", "INFO")]

merged_data <- left_join(gwas, x_mfi_select, by=c("ID", "ALLELE0", "ALLELE1"))
merged_data_order <- relocate(merged_data, INFO, .after=A1FREQ)

############
###OUTPUT###
############
write.table(merged_data_order, opts$o, sep=" ", row.names=FALSE)

