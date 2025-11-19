##############
###LIBRARRY###
##############

library(data.table)
library(dplyr)

##########
###DATA###
##########

HT <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/HT_HTGM_PG_OB_mtag_output_trait_2.txt", data.table = FALSE)
PG <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/HT_HTGM_PG_OB_mtag_output_trait_4.txt", data.table = FALSE)
OB <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/HT_HTGM_PG_OB_mtag_output_trait_3.txt", data.table = FALSE)
HT_GM <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/HT_HTGM_PG_OB_mtag_output_trait_1.txt", data.table = FALSE)

HT_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/cojo/output/merged/HT_mtag.merged.txt", data.table=FALSE)
PG_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/cojo/output/merged/PG_mtag.merged.txt", data.table=FALSE)
OB_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/cojo/output/merged/OB_mtag.merged.txt", data.table=FALSE)
HT_GM_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/cojo/output/merged/HT_GM_mtag.merged.txt", data.table=FALSE)

##############
###ANALYSIS###
##############

HT_select <- HT[,c("SNP", "mtag_beta", "mtag_se", "mtag_pval")]
PG_select <- PG[,c("SNP", "mtag_beta", "mtag_se", "mtag_pval")]
OB_select <- OB[,c("SNP", "mtag_beta", "mtag_se", "mtag_pval")]
HT_GM_select <- HT_GM[,c("SNP", "mtag_beta", "mtag_se", "mtag_pval")]
central_info <- HT[,c("SNP", "CHR", "BP", "A1", "A2", "FRQ")]

colnames(HT_select) <- c("SNP", "HT_beta", "HT_se", "HT_pval")
colnames(PG_select) <- c("SNP", "PG_beta", "PG_se", "PG_pval")
colnames(OB_select) <- c("SNP", "OB_beta", "OB_se", "OB_pval")
colnames(HT_GM_select) <- c("SNP", "HT_GM_beta", "HT_GM_se", "HT_GM_pval")

all_mtag <- left_join(central_info, HT_select, by='SNP') %>%
  left_join(., PG_select, by='SNP') %>% 
  left_join(., OB_select, by='SNP') %>%
  left_join(., HT_GM_select, by='SNP')

write.csv(all_mtag, "all_mtag_merged.csv", quote=FALSE, row.names=FALSE)
system("dx upload all_mtag_merged.csv --path data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/")

all_mtag_sig <- subset(all_mtag, HT_pval < 5E-8 | PG_pval < 5E-8 | OB_pval < 5E-8 | HT_GM_pval < 5E-8)
all_mtag_sig$multiple_testing <- 0
all_mtag_sig[which(all_mtag_sig$HT_pval < 1.25E-8 |all_mtag_sig$PG_pval < 1.25E-8 | all_mtag_sig$OB_pval < 1.25E-8 | all_mtag_sig$HT_GM_pval < 1.25E-8),"multiple_testing"] <- 1
write.csv(all_mtag_sig, "all_mtag_merged_sig.csv", quote=FALSE, row.names=FALSE)
system("dx upload all_mtag_merged_sig.csv --path data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/")

all_mtag_cojo <- do.call(rbind, list(HT_cojo, PG_cojo, OB_cojo, HT_GM_cojo))
all_mtag_sig_cojo <- subset(all_mtag, SNP %in% all_mtag_cojo$SNP)
