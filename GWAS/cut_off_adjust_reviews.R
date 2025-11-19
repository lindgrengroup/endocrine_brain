###############
###LIBRARIES###
###############

library(data.table)

##########
###DATA###
##########

system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/meta_table_cojo_cond_edit.csv")

HT <- fread("assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
PG <- fread("assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
OB <- fread("assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM <- fread("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

meta_cojo_cond <- fread("meta_table_cojo_cond_edit.csv", data.table=FALSE)

##############

cut_off <- -log10(5e-8)
cut_off_multi <- -log10((5e-8)/4)

meta_cojo_cond$log10pall <- apply(meta_cojo_cond[,c("LOG10P_HT", "LOG10P_HT_GM", "LOG10P_PG", "LOG10P_OB")], 1, max)

meta_cojo_cond_lost <- subset(meta_cojo_cond, log10pall < cut_off_multi)
