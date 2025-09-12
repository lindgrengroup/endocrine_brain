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

#brain_morph <- fread("/mnt/project/Phenotypes_MRI/allIDPs_v1_2021-11-18_IDs11867.csv", data.table=FALSE)
smith_confounds <- fread("/mnt/project/Phenotypes_MRI/confounds_12eBETA_2022-02-11.csv", data.table=FALSE)
PCs <- fread("/mnt/project/Phenotypes_MRI/genetic_pcs21_participant.csv", data.table=FALSE)
genetic_sex <- fread("/mnt/project/Phenotypes_MRI/genetic_sex_participant.csv", data.table=FALSE)
euro_ids <- fread("/mnt/project/final_EUR_list.tsv", data.table=FALSE)
fam_file <- fread("/mnt/project/QC_byUKB/ukb11867_cal_chr1_v2_s488363.fam", data.table=FALSE)
#intensity_norm <- fread("/mnt/project/Phenotypes_MRI/allIDPs_v3_2023-03-06.csv", data.table=FALSE)
pheno_file <- fread("/mnt/project/Phenotypes_MRI/allIDPs_v5_2023-07-18_norm.csv", data.table=FALSE)
bridge <- fread("/mnt/project/bridge_11867_8107.csv", data.table=FALSE)

##############
###ANALYSIS###
##############

#select only the 0.3 threshold values
brain_morph_threshold <- select(pheno_file, matches("HT-SumGM-threshold=0.3-warpResolution=2mm|ID"))



#Make all ID variables characters
brain_morph_threshold$subjectID <- as.character(brain_morph_threshold$subjectID)
brain_morph_threshold$subjectID <- substring(brain_morph_threshold$subjectID, 2)
smith_confounds$subject_ID <- as.character(smith_confounds$subject_ID)
smith_confounds$subject_IDs_internal <- as.character(smith_confounds$subject_IDs_internal)
PCs$eid <- as.character(PCs$eid)
euro_ids$V1 <- as.character(euro_ids$V1)
genetic_sex$`Participant ID` <- as.character(genetic_sex$`Participant ID`)
fam_file$V1 <- as.character(fam_file$V1)
bridge$eid_11867 <- as.character(bridge$eid_11867)
bridge$eid_8107 <- as.character(bridge$eid_8107)

brain_morph_threshold <- left_join(brain_morph_threshold, bridge, by=c("subjectID"="eid_8107"))

#select European population population
brain_morph_threshold <- subset(brain_morph_threshold, eid_11867 %in% euro_ids$V1)

brain_morph_threshold_sex <- left_join(brain_morph_threshold, genetic_sex, by=c("eid_11867"="Participant ID"))
brain_morph_threshold_f <- subset(brain_morph_threshold_sex, `Genetic sex`=="Female")
brain_morph_threshold_m <- subset(brain_morph_threshold_sex, `Genetic sex`=="Male")
brain_morph_threshold_f$`Genetic sex` <- NULL
brain_morph_threshold_m$`Genetic sex` <- NULL


#brain_morph_threshold_f_OB_LR <- select(brain_morph_threshold_f_OB, matches("LR|subjectID|ID11867|ID8107"))
#brain_morph_threshold_m_OB_LR <- select(brain_morph_threshold_m_OB, matches("LR|subjectID|ID11867|ID8107"))

#stretch the data into long format
brain_morph_threshold_f_long <- brain_morph_threshold_f %>% gather("measure", "value", -c(subjectID,eid_11867,subjectID_internal))
brain_morph_threshold_m_long <- brain_morph_threshold_m %>% gather("measure", "value", -c(subjectID,eid_11867, subjectID_internal))

#Join the phenotypic and smith variables
brain_morph_threshold_f_long <- brain_morph_threshold_f_long %>% left_join(smith_confounds, by=c("subjectID" = "subject_ID"))
brain_morph_threshold_m_long <- brain_morph_threshold_m_long %>% left_join(smith_confounds, by=c("subjectID" = "subject_ID"))

#Calculate residuals from models, regressing out Smith covariates
brain_morph_threshold_f_long_resid <- brain_morph_threshold_f_long %>% group_by(measure) %>% nest() %>% 
        mutate(.,
                data_clean=map(data, ~ .x %>% filter(if_all(everything(),~!is.na(.)))),
                resid=map(data_clean, ~ resid(gam(value ~ s(scan_date1) + s(scan_time1) +
                                                                 Site_1_vs_2 + Site_1_vs_3 + Age_Site_1 + Age_Site_2 + Age_Site_3 +
                                                                 Sex_1_Site_1 + Sex_1_Site_2 + Sex_1_Site_3 + HeadSize_Site_1 +
                                                                 HeadSize_Site_2 + HeadSize_Site_3 + TablePos_Table_Site_1 + TablePos_Table_Site_2 +
                                                                 TablePos_Table_Site_3 + TablePos_COG_X_Site_1 + TablePos_COG_Y_Site_1 +
                                                                 TablePos_COG_Z_Site_1 + TablePos_COG_X_Site_2 + TablePos_COG_Y_Site_2 +
                                                                 TablePos_COG_Z_Site_2 + TablePos_COG_X_Site_3 + TablePos_COG_Y_Site_3 +
                                                                 TablePos_COG_Z_Site_3 + HeadMotion_mean_rfMRI_rel_Site_1 +
                                                                 HeadMotion_mean_rfMRI_rel_Site_2 + HeadMotion_mean_rfMRI_rel_Site_3 +
                                                                 HeadMotion_mean_tfMRI_rel_Site_1 + HeadMotion_mean_tfMRI_rel_Site_2 +
                                                                 HeadMotion_mean_tfMRI_rel_Site_3, data=., na.action=na.omit))))

brain_morph_threshold_m_long_resid <- brain_morph_threshold_m_long %>% group_by(measure) %>% nest() %>% 
  mutate(.,
         data_clean=map(data, ~ .x %>% filter(if_all(everything(),~!is.na(.)))),
         resid=map(data_clean, ~ resid(gam(value ~ s(scan_date1) + s(scan_time1) +
                                             Site_1_vs_2 + Site_1_vs_3 + Age_Site_1 + Age_Site_2 + Age_Site_3 +
                                             Sex_1_Site_1 + Sex_1_Site_2 + Sex_1_Site_3 + HeadSize_Site_1 +
                                             HeadSize_Site_2 + HeadSize_Site_3 + TablePos_Table_Site_1 + TablePos_Table_Site_2 +
                                             TablePos_Table_Site_3 + TablePos_COG_X_Site_1 + TablePos_COG_Y_Site_1 +
                                             TablePos_COG_Z_Site_1 + TablePos_COG_X_Site_2 + TablePos_COG_Y_Site_2 +
                                             TablePos_COG_Z_Site_2 + TablePos_COG_X_Site_3 + TablePos_COG_Y_Site_3 +
                                             TablePos_COG_Z_Site_3 + HeadMotion_mean_rfMRI_rel_Site_1 +
                                             HeadMotion_mean_rfMRI_rel_Site_2 + HeadMotion_mean_rfMRI_rel_Site_3 +
                                             HeadMotion_mean_tfMRI_rel_Site_1 + HeadMotion_mean_tfMRI_rel_Site_2 +
                                             HeadMotion_mean_tfMRI_rel_Site_3, data=., na.action=na.omit))))

#Reformat into dataframe by unnesting etc.
brain_morph_threshold_f_long_resid <- brain_morph_threshold_f_long_resid %>% mutate(data_resid=map2(data_clean,resid,~.x %>% mutate(resid=.y)))
brain_morph_threshold_f_wide_resid <- brain_morph_threshold_f_long_resid %>% select(-c(data, resid, data_clean)) %>% unnest(data_resid)

brain_morph_threshold_f_wide_resid_df <- brain_morph_threshold_f_wide_resid %>% select(-c(value)) %>%
  pivot_wider(names_from=measure, values_from=resid, names_prefix="resid_") %>% as.data.frame()

brain_morph_threshold_m_long_resid <- brain_morph_threshold_m_long_resid %>% mutate(data_resid=map2(data_clean,resid,~.x %>% mutate(resid=.y)))
brain_morph_threshold_m_wide_resid <- brain_morph_threshold_m_long_resid %>% select(-c(data, resid, data_clean)) %>% unnest(data_resid)

brain_morph_threshold_m_wide_resid_df <- brain_morph_threshold_m_wide_resid %>% select(-c(value)) %>%
  pivot_wider(names_from=measure, values_from=resid, names_prefix="resid_") %>% as.data.frame()


#Inverse normalise each value
brain_morph_threshold_f_wide_resid_norm <- brain_morph_threshold_f_wide_resid_df %>% mutate(across(colnames(brain_morph_threshold_f_wide_resid_df)[grep("*2mm-norm|SumGM", colnames(brain_morph_threshold_f_wide_resid_df))],invnorm, .names = "{col}_norm"))
brain_morph_threshold_m_wide_resid_norm <- brain_morph_threshold_m_wide_resid_df %>% mutate(across(colnames(brain_morph_threshold_m_wide_resid_df)[grep("*2mm-norm|SumGM", colnames(brain_morph_threshold_m_wide_resid_df))],invnorm, .names = "{col}_norm"))

brain_morph_threshold_f_wide_resid_norm_PCs <- brain_morph_threshold_f_wide_resid_norm %>% left_join(PCs, by=c("eid_11867" = "eid"))
brain_morph_threshold_m_wide_resid_norm_PCs <- brain_morph_threshold_m_wide_resid_norm %>% left_join(PCs, by=c("eid_11867" = "eid"))



brain_morph_threshold_resid_norm_PCs_euro_f <- brain_morph_threshold_f_wide_resid_norm_PCs[,grep("resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm|*2mm-norm_norm|eid_11867", colnames(brain_morph_threshold_f_wide_resid_norm_PCs))]
brain_morph_threshold_resid_norm_PCs_euro_f_format <- right_join(fam_file, brain_morph_threshold_resid_norm_PCs_euro_f, by=c("V1" = "eid_11867"))
brain_morph_threshold_resid_norm_PCs_euro_f_format <- select(brain_morph_threshold_resid_norm_PCs_euro_f_format, !c("V3", "V4", "V5", "V6"))
brain_morph_threshold_resid_norm_PCs_euro_f_format <- brain_morph_threshold_resid_norm_PCs_euro_f_format %>% rename(FID=V1) %>% rename(IID=V2)
write.table(brain_morph_threshold_resid_norm_PCs_euro_f_format, "brain_morph_threshold_resid_norm_euro_female_GM_510_03.phe", quote = FALSE, row.names=FALSE)
brain_morph_threshold_resid_norm_PCs_euro_f_covars <- brain_morph_threshold_f_wide_resid_norm_PCs[,grep("p22009|eid_11867", colnames(brain_morph_threshold_f_wide_resid_norm_PCs))]
brain_morph_threshold_resid_norm_PCs_euro_f_covars <- right_join(fam_file, brain_morph_threshold_resid_norm_PCs_euro_f_covars, by=c("V1" = "eid_11867"))
brain_morph_threshold_resid_norm_PCs_euro_f_covars <- select(brain_morph_threshold_resid_norm_PCs_euro_f_covars, !c("V3", "V4", "V5", "V6"))
brain_morph_threshold_resid_norm_PCs_euro_f_covars <- brain_morph_threshold_resid_norm_PCs_euro_f_covars %>% rename(FID=V1) %>% rename(IID=V2)
write.table(brain_morph_threshold_resid_norm_PCs_euro_f_covars, "brain_morph_threshold_resid_norm_euro_female_GM_510_03.cov", quote = FALSE, row.names=FALSE)

brain_morph_threshold_resid_norm_PCs_euro_m <- brain_morph_threshold_m_wide_resid_norm_PCs[,grep("resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm|*2mm-norm_norm|eid_11867", colnames(brain_morph_threshold_m_wide_resid_norm_PCs))]
brain_morph_threshold_resid_norm_PCs_euro_m_format <- right_join(fam_file, brain_morph_threshold_resid_norm_PCs_euro_m, by=c("V1" = "eid_11867"))
brain_morph_threshold_resid_norm_PCs_euro_m_format <- select(brain_morph_threshold_resid_norm_PCs_euro_m_format, !c("V3", "V4", "V5", "V6"))
brain_morph_threshold_resid_norm_PCs_euro_m_format <- brain_morph_threshold_resid_norm_PCs_euro_m_format %>% rename(FID=V1) %>% rename(IID=V2)
write.table(brain_morph_threshold_resid_norm_PCs_euro_m_format, "brain_morph_threshold_resid_norm_euro_male_GM_510_03.phe", quote = FALSE, row.names = FALSE)
brain_morph_threshold_resid_norm_PCs_euro_m_covars <- brain_morph_threshold_m_wide_resid_norm_PCs[,grep("p22009|eid_11867", colnames(brain_morph_threshold_m_wide_resid_norm_PCs))]
brain_morph_threshold_resid_norm_PCs_euro_m_covars <- right_join(fam_file, brain_morph_threshold_resid_norm_PCs_euro_m_covars, by=c("V1" = "eid_11867"))
brain_morph_threshold_resid_norm_PCs_euro_m_covars <- select(brain_morph_threshold_resid_norm_PCs_euro_m_covars, !c("V3", "V4", "V5", "V6"))
brain_morph_threshold_resid_norm_PCs_euro_m_covars <- brain_morph_threshold_resid_norm_PCs_euro_m_covars %>% rename(FID=V1) %>% rename(IID=V2)
write.table(brain_morph_threshold_resid_norm_PCs_euro_m_covars, "brain_morph_threshold_resid_norm_euro_male_GM_510_03.cov", quote = FALSE, row.names = FALSE)


brain_morph_threshold_resid_norm_PCs_euro_all_LR_format <- bind_rows(brain_morph_threshold_resid_norm_PCs_euro_f_format, brain_morph_threshold_resid_norm_PCs_euro_m_format) 
write.table(brain_morph_threshold_resid_norm_PCs_euro_all_LR_format, "brain_morph_threshold_resid_norm_euro_all_GM_510_03.phe", quote=FALSE, row.names=FALSE)
brain_morph_threshold_resid_norm_PCs_euro_all_LR_covars <- bind_rows(brain_morph_threshold_resid_norm_PCs_euro_f_covars, brain_morph_threshold_resid_norm_PCs_euro_m_covars) 
write.table(brain_morph_threshold_resid_norm_PCs_euro_all_LR_covars, "brain_morph_threshold_resid_norm_euro_all_GM_510_03.cov", quote=FALSE, row.names=FALSE)


system("dx upload *.phe --path Phenotypes_MRI/regenie_input/phenotype_files/")
system("dx upload *.cov --path Phenotypes_MRI/regenie_input/covar_files/")










