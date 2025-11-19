###############
###LIBRARIES###
###############

library(data.table)
library(mgcv)
library(dplyr)
library(tidyr)
library(purrr)

###############
###FUNCTIONS###
###############

invnorm = function(x) {
  qnorm((rank(x, na.last="keep", ties.method="random") - 0.5) / sum(!is.na(x)));
}

##########
###DATA###
##########

confounds <- fread("/mnt/project/data/replication/replication_confounds_13j_2025-05-15.csv", data.table=FALSE)
PCs <- fread("/mnt/project/Phenotypes_MRI/genetic_pcs21_participant.csv", data.table=FALSE)
genetic_sex <- fread("/mnt/project/Phenotypes_MRI/genetic_sex_participant.csv", data.table=FALSE)
euro_ids <- fread("/mnt/project/final_EUR_list.tsv", data.table=FALSE)
fam_file <- fread("/mnt/project/QC_byUKB/ukb11867_cal_chr1_v2_s488363.fam", data.table=FALSE)
pheno_file <- fread("/mnt/project/data/replication/replication_allIDPs_v6_2025-05-15.csv", data.table=FALSE)
bridge <- fread("/mnt/project/bridge_11867_8107.csv", data.table=FALSE)
og_pheno <- fread("/mnt/project/Phenotypes_MRI/allIDPs_v5_2023-07-18_norm.csv", data.table=FALSE)

##############
###ANALYSIS###
##############

#select only the 0.3 threshold values
pheno_select <- select(pheno_file, matches("HT-volume-threshold=0.3-warpResolution=2mm|PG-volume-threshold=0.3-warpResolution=2mm|OB_left-volume-threshold=0.3-warpResolution=2mm|OB_right-volume-threshold=0.3-warpResolution=2mm|HT-SumGM-threshold=0.3-warpResolution=2mm|ID"))


#Make all ID variables characters
pheno_select$subjectID <- as.character(pheno_select$subjectID)
pheno_select$subjectID <- substring(pheno_select$subjectID, 2)
confounds$subject_ID <- as.character(confounds$subject_ID)
confounds$subject_ID <- substring(confounds$subject_ID, 2)
PCs$eid <- as.character(PCs$eid)
euro_ids$V1 <- as.character(euro_ids$V1)
genetic_sex$`Participant ID` <- as.character(genetic_sex$`Participant ID`)
fam_file$V1 <- as.character(fam_file$V1)
bridge$eid_11867 <- as.character(bridge$eid_11867)
bridge$eid_8107 <- as.character(bridge$eid_8107)
og_pheno$subjectID <- as.character(og_pheno$subjectID)
og_pheno$subjectID <- substring(og_pheno$subjectID, 2)

pheno_select <- left_join(pheno_select, bridge, by=c("subjectID"="eid_8107"))

#select European population population
pheno_select <- subset(pheno_select, eid_11867 %in% euro_ids$V1)

#select non-replicates
pheno_select <- subset(pheno_select, ! subjectID %in% og_pheno$subjectID)

#remove all duplicates
pheno_select <- pheno_select %>% group_by(subjectID) %>% 
  filter(n()==1) %>% 
  ungroup() %>% 
  as.data.frame()


pheno_select$`LR-volume-threshold=0.3-warpResolution=2mm` <- rowMeans(pheno_select[,c("OB_left-volume-threshold=0.3-warpResolution=2mm", "OB_right-volume-threshold=0.3-warpResolution=2mm")], na.rm=TRUE)

pheno_select_sex <- left_join(pheno_select, genetic_sex, by=c("eid_11867"="Participant ID"))
pheno_select_f <- subset(pheno_select_sex, `Genetic sex`=="Female")
pheno_select_m <- subset(pheno_select_sex, `Genetic sex`=="Male")
pheno_select_f$`Genetic sex` <- NULL
pheno_select_m$`Genetic sex` <- NULL

#stretch the data into long format
pheno_select_f_long <- pheno_select_f %>% gather("measure", "value", -c(subjectID,eid_11867))
pheno_select_m_long <- pheno_select_m %>% gather("measure", "value", -c(subjectID,eid_11867))

#Join the phenotypic and smith variables
pheno_select_f_long <- pheno_select_f_long %>% left_join(confounds, by=c("subjectID" = "subject_ID"))
pheno_select_m_long <- pheno_select_m_long %>% left_join(confounds, by=c("subjectID" = "subject_ID"))

#Calculate residuals from models, regressing out Smith covariates
pheno_select_f_long_resid <- pheno_select_f_long %>% group_by(measure) %>% nest() %>% 
        mutate(.,
                data_clean=map(data, ~ .x %>% filter(if_all(everything(),~!is.na(.)))),
                resid=map(data_clean, ~ resid(gam(value ~ s(scan_date1) + s(scan_time1) +
                                                                 Site_1_vs_2 + Site_1_vs_3 + Site_1_vs_4 + Age_Site_1 + Age_Site_2 + Age_Site_3 + Age_Site_4 +
                                                                 Sex_1_Site_1 + Sex_1_Site_2 + Sex_1_Site_3 + Sex_1_Site_4 + HeadSize_Site_1 +
                                                                 HeadSize_Site_2 + HeadSize_Site_3 + HeadSize_Site_4 + TablePos_Table_Site_1 + TablePos_Table_Site_2 +
                                                                 TablePos_Table_Site_3 + TablePos_Table_Site_4 + TablePos_COG_X_Site_1 + TablePos_COG_Y_Site_1 +
                                                                 TablePos_COG_Z_Site_1 + TablePos_COG_X_Site_2 + TablePos_COG_Y_Site_2 +
                                                                 TablePos_COG_Z_Site_2 + TablePos_COG_X_Site_3 + TablePos_COG_Y_Site_3 +
                                                                 TablePos_COG_Z_Site_3 + TablePos_COG_X_Site_4 + TablePos_COG_Y_Site_4 + TablePos_COG_Z_Site_4 +
                                                                 HeadMotion_mean_rfMRI_rel_Site_1 +
                                                                 HeadMotion_mean_rfMRI_rel_Site_2 + HeadMotion_mean_rfMRI_rel_Site_3 + HeadMotion_mean_rfMRI_rel_Site_4 +
                                                                 HeadMotion_mean_tfMRI_rel_Site_1 + HeadMotion_mean_tfMRI_rel_Site_2 +
                                                                 HeadMotion_mean_tfMRI_rel_Site_3 + HeadMotion_mean_tfMRI_rel_Site_4, data=., na.action=na.omit))))

pheno_select_m_long_resid <- pheno_select_m_long %>% group_by(measure) %>% nest() %>% 
  mutate(.,
         data_clean=map(data, ~ .x %>% filter(if_all(everything(),~!is.na(.)))),
         resid=map(data_clean, ~ resid(gam(value ~ s(scan_date1) + s(scan_time1) +
                                             Site_1_vs_2 + Site_1_vs_3 + Site_1_vs_4 + Age_Site_1 + Age_Site_2 + Age_Site_3 + Age_Site_4 +
                                             Sex_1_Site_1 + Sex_1_Site_2 + Sex_1_Site_3 + Sex_1_Site_4 + HeadSize_Site_1 +
                                             HeadSize_Site_2 + HeadSize_Site_3 + HeadSize_Site_4 + TablePos_Table_Site_1 + TablePos_Table_Site_2 +
                                             TablePos_Table_Site_3 + TablePos_Table_Site_4 + TablePos_COG_X_Site_1 + TablePos_COG_Y_Site_1 +
                                             TablePos_COG_Z_Site_1 + TablePos_COG_X_Site_2 + TablePos_COG_Y_Site_2 +
                                             TablePos_COG_Z_Site_2 + TablePos_COG_X_Site_3 + TablePos_COG_Y_Site_3 +
                                             TablePos_COG_Z_Site_3 + TablePos_COG_X_Site_4 + TablePos_COG_Y_Site_4 + TablePos_COG_Z_Site_4 +
                                             HeadMotion_mean_rfMRI_rel_Site_1 +
                                             HeadMotion_mean_rfMRI_rel_Site_2 + HeadMotion_mean_rfMRI_rel_Site_3 + HeadMotion_mean_rfMRI_rel_Site_4 +
                                             HeadMotion_mean_tfMRI_rel_Site_1 + HeadMotion_mean_tfMRI_rel_Site_2 +
                                             HeadMotion_mean_tfMRI_rel_Site_3 + HeadMotion_mean_tfMRI_rel_Site_4, data=., na.action=na.omit))))

#Reformat into dataframe by unnesting etc.
pheno_select_f_long_resid <- pheno_select_f_long_resid %>% mutate(data_resid=map2(data_clean,resid,~.x %>% mutate(resid=.y)))
pheno_select_f_wide_resid <- pheno_select_f_long_resid %>% select(-c(data, resid, data_clean)) %>% unnest(data_resid)

pheno_select_f_wide_resid_df <- pheno_select_f_wide_resid %>% select(-c(value)) %>%
  pivot_wider(names_from=measure, values_from=resid, names_prefix="replication_resid_") %>% as.data.frame()

pheno_select_m_long_resid <- pheno_select_m_long_resid %>% mutate(data_resid=map2(data_clean,resid,~.x %>% mutate(resid=.y)))
pheno_select_m_wide_resid <- pheno_select_m_long_resid %>% select(-c(data, resid, data_clean)) %>% unnest(data_resid)

pheno_select_m_wide_resid_df <- pheno_select_m_wide_resid %>% select(-c(value)) %>%
  pivot_wider(names_from=measure, values_from=resid, names_prefix="replication_resid_") %>% as.data.frame()


#Inverse normalise each value
pheno_select_f_wide_resid_norm <- pheno_select_f_wide_resid_df %>% mutate(across(colnames(pheno_select_f_wide_resid_df)[grep("*2mm", colnames(pheno_select_f_wide_resid_df))],invnorm, .names = "{col}_norm"))
pheno_select_m_wide_resid_norm <- pheno_select_m_wide_resid_df %>% mutate(across(colnames(pheno_select_m_wide_resid_df)[grep("*2mm", colnames(pheno_select_m_wide_resid_df))],invnorm, .names = "{col}_norm"))

pheno_select_f_wide_resid_norm_PCs <- pheno_select_f_wide_resid_norm %>% left_join(PCs, by=c("eid_11867" = "eid"))
pheno_select_m_wide_resid_norm_PCs <- pheno_select_m_wide_resid_norm %>% left_join(PCs, by=c("eid_11867" = "eid"))



replication_resid_norm_PCs_euro_f <- pheno_select_f_wide_resid_norm_PCs[,grep("*2mm_norm|eid_11867", colnames(pheno_select_f_wide_resid_norm_PCs))]
replication_resid_norm_PCs_euro_f_format <- right_join(fam_file, replication_resid_norm_PCs_euro_f, by=c("V1" = "eid_11867"))
replication_resid_norm_PCs_euro_f_format <- select(replication_resid_norm_PCs_euro_f_format, !c("V3", "V4", "V5", "V6"))
replication_resid_norm_PCs_euro_f_format <- replication_resid_norm_PCs_euro_f_format %>% rename(FID=V1) %>% rename(IID=V2)
write.table(replication_resid_norm_PCs_euro_f_format, "replication_resid_norm_euro_female.phe", quote = FALSE, row.names=FALSE)
replication_resid_norm_PCs_euro_f_covars <- pheno_select_f_wide_resid_norm_PCs[,grep("p22009|eid_11867", colnames(pheno_select_f_wide_resid_norm_PCs))]
replication_resid_norm_PCs_euro_f_covars <- right_join(fam_file, replication_resid_norm_PCs_euro_f_covars, by=c("V1" = "eid_11867"))
replication_resid_norm_PCs_euro_f_covars <- select(replication_resid_norm_PCs_euro_f_covars, !c("V3", "V4", "V5", "V6"))
replication_resid_norm_PCs_euro_f_covars <- replication_resid_norm_PCs_euro_f_covars %>% rename(FID=V1) %>% rename(IID=V2)
write.table(replication_resid_norm_PCs_euro_f_covars, "replication_resid_norm_euro_female.cov", quote = FALSE, row.names=FALSE)

replication_resid_norm_PCs_euro_m <- pheno_select_m_wide_resid_norm_PCs[,grep("*2mm_norm|eid_11867", colnames(pheno_select_m_wide_resid_norm_PCs))]
replication_resid_norm_PCs_euro_m_format <- right_join(fam_file, replication_resid_norm_PCs_euro_m, by=c("V1" = "eid_11867"))
replication_resid_norm_PCs_euro_m_format <- select(replication_resid_norm_PCs_euro_m_format, !c("V3", "V4", "V5", "V6"))
replication_resid_norm_PCs_euro_m_format <- replication_resid_norm_PCs_euro_m_format %>% rename(FID=V1) %>% rename(IID=V2)
write.table(replication_resid_norm_PCs_euro_m_format, "replication_resid_norm_euro_male.phe", quote = FALSE, row.names = FALSE)
replication_resid_norm_PCs_euro_m_covars <- pheno_select_m_wide_resid_norm_PCs[,grep("p22009|eid_11867", colnames(pheno_select_m_wide_resid_norm_PCs))]
replication_resid_norm_PCs_euro_m_covars <- right_join(fam_file, replication_resid_norm_PCs_euro_m_covars, by=c("V1" = "eid_11867"))
replication_resid_norm_PCs_euro_m_covars <- select(replication_resid_norm_PCs_euro_m_covars, !c("V3", "V4", "V5", "V6"))
replication_resid_norm_PCs_euro_m_covars <- replication_resid_norm_PCs_euro_m_covars %>% rename(FID=V1) %>% rename(IID=V2)
write.table(replication_resid_norm_PCs_euro_m_covars, "replication_resid_norm_euro_male.cov", quote = FALSE, row.names = FALSE)


replication_resid_norm_PCs_euro_all_LR_format <- bind_rows(replication_resid_norm_PCs_euro_f_format, replication_resid_norm_PCs_euro_m_format) 
write.table(replication_resid_norm_PCs_euro_all_LR_format, "replication_resid_norm_euro_all_GM_510_03.phe", quote=FALSE, row.names=FALSE)
replication_resid_norm_PCs_euro_all_LR_covars <- bind_rows(replication_resid_norm_PCs_euro_f_covars, replication_resid_norm_PCs_euro_m_covars) 
write.table(replication_resid_norm_PCs_euro_all_LR_covars, "replication_resid_norm_euro_all_GM_510_03.cov", quote=FALSE, row.names=FALSE)


system("dx upload *.phe --path Phenotypes_MRI/regenie_input/phenotype_files/")
system("dx upload *.cov --path Phenotypes_MRI/regenie_input/covar_files/")










