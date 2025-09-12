###############
###LIBRARIES###
###############

library(data.table)

##########
###DATA###
##########

HT_all <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_female <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_male <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

HT_GM_all <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM_female <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM_male <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

PG_all <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
PG_female <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
PG_male <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

OB_all <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
OB_female <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
OB_male <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

##############
###ANALYSIS###
##############

HT_all$z_score <- HT_all$BETA/HT_all$SE
HT_female$z_score <- HT_female$BETA/HT_female$SE
HT_male$z_score <- HT_male$BETA/HT_male$SE

HT_all$pval <- 10^(-HT_all$LOG10P)
HT_female$pval <- 10^(-HT_female$LOG10P)
HT_male$pval <- 10^(-HT_male$LOG10P)

HT_GM_all$z_score <- HT_GM_all$BETA/HT_GM_all$SE
HT_GM_female$z_score <- HT_GM_female$BETA/HT_GM_female$SE
HT_GM_male$z_score <- HT_GM_male$BETA/HT_GM_male$SE

HT_GM_all$pval <- 10^(-HT_GM_all$LOG10P)
HT_GM_female$pval <- 10^(-HT_GM_female$LOG10P)
HT_GM_male$pval <- 10^(-HT_GM_male$LOG10P)

PG_all$z_score <- PG_all$BETA/PG_all$SE
PG_female$z_score <- PG_female$BETA/PG_female$SE
PG_male$z_score <- PG_male$BETA/PG_male$SE

PG_all$pval <- 10^(-PG_all$LOG10P)
PG_female$pval <- 10^(-PG_female$LOG10P)
PG_male$pval <- 10^(-PG_male$LOG10P)

OB_all$z_score <- OB_all$BETA/OB_all$SE
OB_female$z_score <- OB_female$BETA/OB_female$SE
OB_male$z_score <- OB_male$BETA/OB_male$SE

OB_all$pval <- 10^(-OB_all$LOG10P)
OB_female$pval <- 10^(-OB_female$LOG10P)
OB_male$pval <- 10^(-OB_male$LOG10P)

colnames(HT_all)[which(names(HT_all) == "ALLELE0")] <- "A2"
colnames(HT_all)[which(names(HT_all) == "ALLELE1")] <- "A1"
colnames(HT_female)[which(names(HT_female) == "ALLELE0")] <- "A2"
colnames(HT_female)[which(names(HT_female) == "ALLELE1")] <- "A1"
colnames(HT_male)[which(names(HT_male) == "ALLELE0")] <- "A2"
colnames(HT_male)[which(names(HT_male) == "ALLELE1")] <- "A1"
 
colnames(HT_GM_all)[which(names(HT_GM_all) == "ALLELE0")] <- "A2"
colnames(HT_GM_all)[which(names(HT_GM_all) == "ALLELE1")] <- "A1"
colnames(HT_GM_female)[which(names(HT_GM_female) == "ALLELE0")] <- "A2"
colnames(HT_GM_female)[which(names(HT_GM_female) == "ALLELE1")] <- "A1"
colnames(HT_GM_male)[which(names(HT_GM_male) == "ALLELE0")] <- "A2"
colnames(HT_GM_male)[which(names(HT_GM_male) == "ALLELE1")] <- "A1"

colnames(PG_all)[which(names(PG_all) == "ALLELE0")] <- "A2"
colnames(PG_all)[which(names(PG_all) == "ALLELE1")] <- "A1"
colnames(PG_female)[which(names(PG_female) == "ALLELE0")] <- "A2"
colnames(PG_female)[which(names(PG_female) == "ALLELE1")] <- "A1"
colnames(PG_male)[which(names(PG_male) == "ALLELE0")] <- "A2"
colnames(PG_male)[which(names(PG_male) == "ALLELE1")] <- "A1"

colnames(OB_all)[which(names(OB_all) == "ALLELE0")] <- "A2"
colnames(OB_all)[which(names(OB_all) == "ALLELE1")] <- "A1"
colnames(OB_female)[which(names(OB_female) == "ALLELE0")] <- "A2"
colnames(OB_female)[which(names(OB_female) == "ALLELE1")] <- "A1"
colnames(OB_male)[which(names(OB_male) == "ALLELE0")] <- "A2"
colnames(OB_male)[which(names(OB_male) == "ALLELE1")] <- "A1"


############
###OUTPUT###
############

write.table(HT_all, "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sex_combined.txt", quote=FALSE, row.names=FALSE)
write.table(HT_female, "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_female.txt", quote=FALSE, row.names=FALSE)
write.table(HT_male, "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_male.txt", quote=FALSE, row.names=FALSE) 

write.table(HT_GM_all, "assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_sex_combined.txt", quote=FALSE, row.names=FALSE)
write.table(HT_GM_female, "assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_female.txt", quote=FALSE, row.names=FALSE)
write.table(HT_GM_male, "assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered_male.txt", quote=FALSE, row.names=FALSE)

write.table(PG_all, "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sex_combined.txt", quote=FALSE, row.names=FALSE)
write.table(PG_female, "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_female.txt", quote=FALSE, row.names=FALSE)
write.table(PG_male, "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_male.txt", quote=FALSE, row.names=FALSE)

write.table(OB_all, "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_sex_combined.txt", quote=FALSE, row.names=FALSE)
write.table(OB_female, "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_female.txt", quote=FALSE, row.names=FALSE)
write.table(OB_male, "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered_male.txt", quote=FALSE, row.names=FALSE)

system("dx upload *.txt --path data/brain_all/genotype_process/for_saskia/withpval/")
system("dx upload *.txt --path data/brain_all/genotype_process/for_saskia/withpval/allele_label/")
system("dx upload *.R --path scripts/")







