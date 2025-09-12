#!/usr/bin/env Rscript
#############
###LIBRARY###
#############

install.packages("docopt", repos = "http://cran.us.r-project.org")
install.packages("data.table", repos = "http://cran.us.r-project.org")
library(docopt)
library(data.table)

###############
###FUNCTIONS###
###############

biallelic <- function(o_bgenie){
  multi_rsID <- o_bgenie$ID[duplicated(o_bgenie$ID)]
  tmp <- o_bgenie[!(o_bgenie$ID %in% multi_rsID),]
  multi_pos <- tmp$GENPOS[duplicated(tmp$GENPOS)]
  tmp <- tmp[!(tmp$ID %in% multi_pos),]
  return(tmp)
}

########################
###PARSING PARAMETERS###
########################

doc <- "Usage:
    gwas_filter.R --i <input> --o <output> --hw <hwe_list> [--maf <maf> --info <info>]

Arguments:
    --i Path to GWAS results to be filtered [default: NULL]
    --o File (with full path) in which to write the filtered file [default: NULL]
    --hw Path to HWE filtered SNPs for imputed data [default: NULL]
Options:
    --maf Minimum maf value on which to filter [default: 0.01]
    --info Minimum info score on which to filter [default: 0.4]
"

opts <- docopt(doc)

#print(doc)
##########
###DATA###
##########

gwas <- fread(opts$i, data.table=FALSE)
hw_list <- fread(opts$hw, data.table=FALSE, header=FALSE)

##############
###ANALYSIS###
##############

opts$maf <- as.numeric(opts$maf)
opts$info <- as.numeric(opts$info)

# Remove variants with AF of 1/0
gwas_f <- gwas[!gwas$A1FREQ %in% c(0,1),]

# Filter based on MAF
gwas_f  <- gwas_f[gwas_f$A1FREQ > opts$maf & (1 - gwas_f$A1FREQ) > opts$maf,]

# Filter on info criterion
gwas_f <- gwas_f[gwas_f$INFO > opts$info,]

# Filter based on HWE
gwas_f <- subset(gwas_f, ID %in% hw_list$V1)

# Filter for biallelic
gwas_f <- biallelic(gwas_f)

write.table(gwas_f, opts$o, sep=" ", row.names=FALSE)

