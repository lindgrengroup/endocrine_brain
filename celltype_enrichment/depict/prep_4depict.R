
## add chr:pos IDs




args = commandArgs(trailingOnly=TRUE)

TRAIT <- args[1]
MEASURE <- args[2]
ANALYSIS <- args[3]
setwd(paste("/well/lindgren/ferreira/PROJECTS/MRI/gwas_updated/", ANALYSIS, sep=""))

INPUT_FILE <- paste("assoc.resid_", TRAIT,".", MEASURE, ".threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt",sep="")
data <- read.table(INPUT_FILE, header=T, stringsAsFactors=F, sep="\t")

data$P = 10^(-data$LOG10P)
data = data[, which(colnames(data) %in% c("ID", "CHROM", "GENPOS", "P"))]


OUTPUT_FILENAME <- paste("assoc.resid_", TRAIT,".", MEASURE, ".", ANALYSIS, ".threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_4DEPICT.txt",sep="")

write.table(data, OUTPUT_FILENAME, quote=F, row.names=F, col.names=T, sep="\t")
