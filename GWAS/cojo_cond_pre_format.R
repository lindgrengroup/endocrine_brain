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
    cojo_cond_pre_forma.R --i <input> --o <output> 
Arguments:
    --i Path to input cojo file from gwas sumstats [default: NULL]
    --o Path to output list of snps for use in cojo_cond analysis [default: NULL]
"

opts <- docopt(doc)

##########
###DATA###
##########

cojo <- fread(opts$i, data.table=FALSE)
snps <- cojo[,"SNP", drop=FALSE]

write.table(snps, opts$o, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
