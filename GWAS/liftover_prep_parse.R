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
    liftover_prep.R --i <input> --o <output> 
Arguments:
    --i Path to input GWAS sumstats [default: NULL]
    --o Path to output input bed file for liftover [default: NULL]
"

opts <- docopt(doc)

##############
###ANALYSIS###
##############

data_in <- fread(opts$i, data.table=FALSE)
data_in$CHROM <- ifelse(data_in$CHROM == 23, "X", data_in$CHROM)
data_in$chr_format <- paste("chr", data_in$CHROM, sep="")
data_in$start <- as.integer(data_in$GENPOS -1)
data_in$end <- as.integer(data_in$GENPOS)
data_in_bed <- data_in[,c("chr_format", "start", "end", "ID")]
write.table(data_in_bed,opts$o, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
