#############
###LIBRARY###
#############

library(data.table)

##########
###DATA###
##########

HT <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
PG <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
OB <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)


##############
###ANALYSIS###
##############

HT$p <- 10^(-HT$LOG10P)
HT_select <- HT[,c("ID", "CHROM", "GENPOS", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "p", "N")]

PG$p <- 10^(-PG$LOG10P)
PG_select <- PG[,c("ID", "CHROM", "GENPOS", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "p", "N")]

OB$p <- 10^(-OB$LOG10P)
OB_select <- OB[,c("ID", "CHROM", "GENPOS", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "p", "N")]

HT_GM$p <- 10^(-HT_GM$LOG10P)
HT_GM_select <- HT_GM[,c("ID", "CHROM", "GENPOS", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "p", "N")]

names(HT_select)[names(HT_select)=="ID"] <- "SNP"
names(HT_select)[names(HT_select)=="ALLELE1"] <- "A1"
names(HT_select)[names(HT_select)=="ALLELE0"] <- "A2"
names(HT_select)[names(HT_select)=="A1FREQ"] <- "freq"
names(HT_select)[names(HT_select)=="BETA"] <- "b"
names(HT_select)[names(HT_select)=="SE"] <- "se"

names(PG_select)[names(PG_select)=="ID"] <- "SNP"
names(PG_select)[names(PG_select)=="ALLELE1"] <- "A1"
names(PG_select)[names(PG_select)=="ALLELE0"] <- "A2"
names(PG_select)[names(PG_select)=="A1FREQ"] <- "freq"
names(PG_select)[names(PG_select)=="BETA"] <- "b"
names(PG_select)[names(PG_select)=="SE"] <- "se"

names(OB_select)[names(OB_select)=="ID"] <- "SNP"
names(OB_select)[names(OB_select)=="ALLELE1"] <- "A1"
names(OB_select)[names(OB_select)=="ALLELE0"] <- "A2"
names(OB_select)[names(OB_select)=="A1FREQ"] <- "freq"
names(OB_select)[names(OB_select)=="BETA"] <- "b"
names(OB_select)[names(OB_select)=="SE"] <- "se"

names(HT_GM_select)[names(HT_GM_select)=="ID"] <- "SNP"
names(HT_GM_select)[names(HT_GM_select)=="ALLELE1"] <- "A1"
names(HT_GM_select)[names(HT_GM_select)=="ALLELE0"] <- "A2"
names(HT_GM_select)[names(HT_GM_select)=="A1FREQ"] <- "freq"
names(HT_GM_select)[names(HT_GM_select)=="BETA"] <- "b"
names(HT_GM_select)[names(HT_GM_select)=="SE"] <- "se"

HT_select$Z <- HT_select$b/HT_select$se
PG_select$Z <- PG_select$b/PG_select$se
OB_select$Z <- OB_select$b/OB_select$se
HT_GM_select$Z <- HT_GM_select$b/HT_GM_select$se

write.table(HT_select,"assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.mtag_format.txt",quote=FALSE, row.names=FALSE)
write.table(PG_select, "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.mtag_format.txt",quote=FALSE, row.names=FALSE)
write.table(OB_select, "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.mtag_format.txt",quote=FALSE, row.names=FALSE)
write.table(HT_GM_select, "assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.mtag_format.txt",quote=FALSE, row.names=FALSE)

system("dx upload *mtag_format.txt --path data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_format/")
