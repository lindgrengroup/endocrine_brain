###############
###LIBRARIES###
###############

library(data.table)

########################
###PARSING PARAMETERS###
########################

##########
###DATA###
##########

HT <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/HT_HTGM_PG_OB_mtag_output_trait_2.txt", data.table = FALSE)
PG <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/HT_HTGM_PG_OB_mtag_output_trait_4.txt", data.table = FALSE)
OB <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/HT_HTGM_PG_OB_mtag_output_trait_3.txt", data.table = FALSE)
HT_GM <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/HT_HTGM_PG_OB_mtag_output_trait_1.txt", data.table = FALSE)

##############
###ANALYSIS###
##############

HT_select <- HT[,c("SNP", "A1", "A2", "FRQ", "mtag_beta", "mtag_se", "mtag_pval", "N")]
names(HT_select)[names(HT_select)=="FRQ"] <- "freq"
names(HT_select)[names(HT_select)=="mtag_beta"] <- "b"
names(HT_select)[names(HT_select)=="mtag_se"] <- "se"

PG_select <- PG[,c("SNP", "A1", "A2", "FRQ", "mtag_beta", "mtag_se", "mtag_pval", "N")]
names(PG_select)[names(PG_select)=="FRQ"] <- "freq"
names(PG_select)[names(PG_select)=="mtag_beta"] <- "b"
names(PG_select)[names(PG_select)=="mtag_se"] <- "se"

OB_select <- OB[,c("SNP", "A1", "A2", "FRQ", "mtag_beta", "mtag_se", "mtag_pval", "N")]
names(OB_select)[names(OB_select)=="FRQ"] <- "freq"
names(OB_select)[names(OB_select)=="mtag_beta"] <- "b"
names(OB_select)[names(OB_select)=="mtag_se"] <- "se"

HT_GM_select <- HT_GM[,c("SNP", "A1", "A2", "FRQ", "mtag_beta", "mtag_se", "mtag_pval", "N")]
names(HT_GM_select)[names(HT_GM_select)=="FRQ"] <- "freq"
names(HT_GM_select)[names(HT_GM_select)=="mtag_beta"] <- "b"
names(HT_GM_select)[names(HT_GM_select)=="mtag_se"] <- "se"

write.table(HT_select, "HT_mtag_cojo_input.ma", quote=FALSE, row.names=FALSE)
write.table(PG_select, "PG_mtag_cojo_input.ma", quote=FALSE, row.names=FALSE)
write.table(OB_select, "OB_mtag_cojo_input.ma", quote=FALSE, row.names=FALSE)
write.table(HT_GM_select, "HT_GM_mtag_cojo_input.ma", quote=FALSE, row.names=FALSE)

system("dx upload *.ma --path data/brain_all/genotype_process/all/regenie_step2/filtered/mtag_results/cojo/input/")



