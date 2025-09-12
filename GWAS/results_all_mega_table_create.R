###############
###LIBRARIES###
###############

library(data.table)
library(tidyverse)
library(dplyr)

##########
###DATA###
##########

HT_vol_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX.txt", data.table=FALSE)
HT_GM_vol_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)
PG_vol_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table = FALSE)
OB_vol_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table = FALSE)

HT_vol_all <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table = FALSE)
HT_GM_vol_all <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
PG_vol_all <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
OB_vol_all <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table = FALSE)


##############
###ANALYSIS###
##############

all_snps <- c(HT_vol_cojo$SNP, HT_GM_vol_cojo$SNP, PG_vol_cojo$SNP, OB_vol_cojo$SNP)
all_snps_unique <- unique(all_snps)

HT_vol_all_snps <- subset(HT_vol_all, ID %in% all_snps_unique)
HT_GM_vol_all_snps <- subset(HT_GM_vol_all, ID %in% all_snps_unique)
PG_vol_all_snps <- subset(PG_vol_all, ID %in% all_snps_unique)
OB_vol_all_snps <- subset(OB_vol_all, ID %in% all_snps_unique)

mega_table <- full_join(HT_vol_all_snps, HT_GM_vol_all_snps, by=c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "EXTRA"), suffix=c("", "_HT_GM")) %>%
  full_join(PG_vol_all_snps, by=c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "EXTRA"), suffix=c("", "_PG")) %>%
  full_join(OB_vol_all_snps, by=c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "EXTRA"), suffix=c("_HT", "_OB"))

mega_table$P_HT <- 10^(-mega_table$LOG10P_HT)
mega_table$P_PG <- 10^(-mega_table$LOG10P_PG)
mega_table$P_OB <- 10^(-mega_table$LOG10P_OB)
mega_table$P_HT_GM <- 10^(-mega_table$LOG10P_HT_GM)

mega_table_format <- mega_table[,c("ID", "CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "EXTRA", "BETA_HT", 
                                   "BETA_HT_GM", "BETA_PG", "BETA_OB", "SE_HT", "SE_HT_GM", "SE_PG", "SE_OB", "P_HT", "P_HT_GM", "P_PG", "P_OB")]

mega_table_format_order <- mega_table_format[order(mega_table_format$CHROM),]

write.csv(mega_table_format_order, "mega_volume_results_table_paper.csv", quote=FALSE, row.names=FALSE)

HT_vol_all_set <- subset(HT_vol_all, ID %in% HT_vol_cojo$SNP)
HT_GM_vol_all_set <- subset(HT_GM_vol_all, ID %in% HT_GM_vol_cojo$SNP)
PG_vol_all_set <- subset(PG_vol_all, ID %in% PG_vol_cojo$SNP)
OB_vol_all_set <- subset(OB_vol_all, ID %in% OB_vol_cojo$SNP)

HT_vol_all_set$Phenotype <- "HT"
HT_GM_vol_all_set$Phenotype <- "HT-GM"
PG_vol_all_set$Phenotype <- "PG"
OB_vol_all_set$Phenotype <- "OB"

HT_vol_all_set$P <- 10^(-HT_vol_all_set$LOG10P)
HT_GM_vol_all_set$P <- 10^(-HT_GM_vol_all_set$LOG10P)
PG_vol_all_set$P <- 10^(-PG_vol_all_set$LOG10P)
OB_vol_all_set$P <- 10^(-OB_vol_all_set$LOG10P)

HT_vol_all_set_order <- HT_vol_all_set[order(HT_vol_all_set$CHROM, HT_vol_all_set$GENPOS),]
HT_GM_vol_all_set_order <- HT_GM_vol_all_set[order(HT_GM_vol_all_set$CHROM, HT_GM_vol_all_set$GENPOS),]
PG_vol_all_set_order <- PG_vol_all_set[order(PG_vol_all_set$CHROM, PG_vol_all_set$GENPOS),]
OB_vol_all_set_order <- OB_vol_all_set[order(OB_vol_all_set$CHROM, OB_vol_all_set$GENPOS),]

mega_table_assoc_prep <- do.call(rbind, list(HT_vol_all_set_order, PG_vol_all_set_order, OB_vol_all_set_order, HT_GM_vol_all_set_order))
mega_table_assoc_prep_select <- mega_table_assoc_prep[,c("ID", "CHROM", "Phenotype", "P")]
write.csv(mega_table_assoc_prep_select, "mega_table_assoc_prep_select_paper.csv", quote=FALSE, row.names=FALSE)
