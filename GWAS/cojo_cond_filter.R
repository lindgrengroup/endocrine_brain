###############
###LIBRARIES###
###############

install.packages("docopt", repos = "http://cran.us.r-project.org")
install.packages("data.table", repos = "http://cran.us.r-project.org")
library(docopt)
library(data.table)

##########
###DATA###
##########

doc <- "Usage:
    cojo_cond_filter.R --i <input> --o <output>

Arguments:
    --i Path to cojo cond results to be filtered [default: NULL]
    --o File (with full path) in which to write the filtered file [default: NULL]
"

opts <- docopt(doc)

##############
###ANALYSIS###
##############

cojo_cond <- fread(opts$i, data.table=FALSE)
cojo_cond_sig <- subset(cojo_cond, pC < 0.05 & p < 5e-8)

write.table(cojo_cond_sig, opts$o, sep="\t", quote=FALSE, row.names=FALSE)