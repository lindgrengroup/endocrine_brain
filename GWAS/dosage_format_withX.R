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
    dosage_format.R --i <input> --n <name> --d <dosage> --m <metadata>
Arguments:
    --i Path to output of qctool dosage conversion [default: NULL]
    --n Name of dosage set [default: NULL]
    --d Path to output dosage file [default: NULL]
    --m Path to output meta data file [default: NULL] 
"

opts <- docopt(doc)

##########
###DATA###
##########

print(opts$i)
print(opts$n)
print((paste(opts$i, opts$n, "*.dosage", sep="")))
dataFiles <- lapply(Sys.glob(paste(opts$i, opts$n, "*.dosage", sep="")), function(x) fread(x, header=TRUE, sep=" ", stringsAsFactors=FALSE, data.table=FALSE))


##############
###ANALYSIS###
##############


dataFiles_clean <- dataFiles[sapply(dataFiles, function(x) dim(x)[1]) > 0]
dataFiles_clean <- lapply(dataFiles_clean, function(x){
    x$chromosome <- as.character(x$chromosome)
    return(x)
})
                           


dataFilesTable <- do.call(bind_rows, dataFiles_clean)
print(dataFilesTable[1:3,1:3])
print(dim(dataFilesTable))
dataFilesTable <- dataFilesTable %>% distinct(rsid, .keep_all=TRUE)
print(dim(dataFilesTable))
rownames(dataFilesTable) <- dataFilesTable$rsid
dataFilesTable_t <- t(dataFilesTable)
dataFilesTable_t <- as.data.frame(dataFilesTable_t)
meta_data <- dataFilesTable_t[which(row.names(dataFilesTable_t) %in% c("chromosome", "SNPID", "rsid", "position", "alleleA", "alleleB")),]
dosage_data <- dataFilesTable_t[which(!row.names(dataFilesTable_t) %in% c("chromosome", "SNPID", "rsid", "position", "alleleA", "alleleB")),]

write.table(dosage_data, opts$d, col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(meta_data, opts$m, col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")