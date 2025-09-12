###############
###LIBRARIES###
###############

install.packages("purrr")
library(purrr)
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
    hard_dosage_format.R --i <input> --n <name> --d <dosage>
Arguments:
    --i Path to output of qctool dosage conversion [default: NULL]
    --n Name of dosage set [default: NULL]
    --d Path to output dosage file [default: NULL]
"

opts <- docopt(doc)

##########
###DATA###
##########

print(opts$i)
print(opts$n)
print((paste(opts$i, opts$n, "*.raw", sep="")))
dataFiles <- lapply(Sys.glob(paste(opts$i, opts$n, "*.raw", sep="")), function(x) fread(x, header=TRUE, sep="\t", stringsAsFactors=FALSE, data.table=FALSE))


##############
###ANALYSIS###
##############


dataFiles_clean <- dataFiles[sapply(dataFiles, function(x) dim(x)[1]) > 0]

dataFilesTable <- dataFiles_clean %>% reduce(full_join, by=c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))

write.table(dataFilesTable, opts$d, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
