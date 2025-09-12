###############
###LIBRARIES###
###############

install.packages("data.table")
library(data.table)
install.packages("docopt")
library(docopt)

########################
###PARSING PARAMETERS###
########################

doc <- "Usage:
    dosage_prep.R --g <genotype_sample> --p <phenotype_file> --c <cojo_result> --s <sample_out> --r <rsids_out>
Arguments:
    --g Path to genotype sample file [default: NULL]
    --p Path to phenotype file from GWAS [default: NULL]
    --c Path to cojo filtered file of rsids [default: NULL]
    --s Path to write the sample file [default: NULL]
    --r Path to write rsid file [default: NULL]
"

opts <- docopt(doc)


##########
###DATA###
##########

genotype_sample <- fread(opts$g, data.table=FALSE)
phenotype_file <- fread(opts$p, data.table=FALSE)
gwas_cojo <- fread(opts$c, data.table=FALSE)

##############
###ANALYSIS###
##############

genotype_sample_subset <- subset(genotype_sample, ID_1 %in% c(phenotype_file$FID, 0))

rsids <- gwas_cojo$SNP

write.table(genotype_sample_subset, opts$s, sep=" ", row.names=FALSE, quote=FALSE)
write.table(rsids, opts$r, sep="\n", row.names=FALSE, col.names=FALSE, quote=FALSE)


