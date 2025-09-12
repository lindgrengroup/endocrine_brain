###############
###LIBRARIES###
###############

install.packages("data.table")
install.packages("docopt")
library(data.table)
library(docopt)

########################
###PARSING PARAMETERS###
########################

doc <- "Usage:
    gwas_cojo_format.R --i <input> --o <output>
Arguments:
    --i Path to filtered GWAS results to be formatted for cojo analysis [default: NULL]
    --o Path to output file [default: NULL]
"

opts <- docopt(doc)

##########
###DATA###
##########

d <- fread(opts$i, data.table=FALSE)

##############
###ANALYSIS###
##############

d$p <- 10^(-d$LOG10P)
d_select <- d[,c("ID", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "p", "N")]

names(d_select)[names(d_select)=="ID"] <- "SNP"
names(d_select)[names(d_select)=="ALLELE1"] <- "A1"
names(d_select)[names(d_select)=="ALLELE0"] <- "A2"
names(d_select)[names(d_select)=="A1FREQ"] <- "freq"
names(d_select)[names(d_select)=="BETA"] <- "b"
names(d_select)[names(d_select)=="SE"] <- "se"

write.table(d_select, opts$o, quote=FALSE, row.names=FALSE)


