install.packages("data.table")
library(data.table)
install.packages("docopt")
library(docopt)

########################
###PARSING PARAMETERS###
########################

doc <- "Usage:
    liftover_post_format.R --g <gwas> --l <liftover> --o <output> 
Arguments:
    --g Path to input GWAS sumstats [default: NULL]
    --l Path to liftover mapping file [default: NULL]
    --o Path to output merged sumstats liftover file [default: NULL]
"

opts <- docopt(doc)

##############
###ANALYSIS###
##############

gwas <- fread(opts$g, data.table=FALSE)
gwas$hg19pos <- paste("chr", gwas$CHROM, ":", gwas$GENPOS, sep="")

liftover_map <- fread(opts$l, data.table=FALSE)
liftover_map$hg38pos <- paste(liftover_map$V1, ":", liftover_map$V3, sep="")

gwas_liftover <- merge(gwas, liftover_map[,c("V4", "hg38pos")], by.x="ID", by.y="V4")

write.table(gwas_liftover, opts$o, row.names=FALSE, quote=FALSE, sep="\t")
