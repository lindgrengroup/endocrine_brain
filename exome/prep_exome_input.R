###############
###LIBRARIES###
###############

library(data.table)
library(tidyverse)

##########
###DATA###
##########

pheno_vols_all <- fread("/mnt/project/Phenotypes_MRI/regenie_input/phenotype_files/brain_morph_threshold_resid_norm_euro_all.phe", data.table=FALSE)
pheno_OB_LR_vol_all <- fread("/mnt/project/Phenotypes_MRI/regenie_input/phenotype_files/brain_morph_threshold_resid_norm_euro_all_OB_LR.phe", data.table=FALSE)
pheno_HT_WM_vol_all <- fread("/mnt/project/Phenotypes_MRI/regenie_input/phenotype_files/brain_morph_threshold_resid_norm_euro_all_GM_510_03.phe", data.table=FALSE)
covar_vols_all <- fread("/mnt/project/Phenotypes_MRI/regenie_input/covar_files/brain_morph_threshold_resid_norm_euro_all.cov", data.table=FALSE)
covar_OB_LR_vol_all <- fread("/mnt/project/Phenotypes_MRI/regenie_input/covar_files/brain_morph_threshold_resid_norm_euro_all_OB_LR.cov", data.table=FALSE)
covar_HT_WM_vol_all <- fread("/mnt/project/Phenotypes_MRI/regenie_input/covar_files/brain_morph_threshold_resid_norm_euro_all_GM_510_03.cov", data.table=FALSE)

pheno_vols_female <- fread("/mnt/project/Phenotypes_MRI/regenie_input/phenotype_files/brain_morph_threshold_resid_norm_euro_female.phe", data.table=FALSE)
pheno_OB_LR_vol_female <- fread("/mnt/project/Phenotypes_MRI/regenie_input/phenotype_files/brain_morph_threshold_resid_norm_euro_female_OB_LR.phe", data.table=FALSE)
pheno_HT_WM_vol_female <- fread("/mnt/project/Phenotypes_MRI/regenie_input/phenotype_files/brain_morph_threshold_resid_norm_euro_female_GM_510_03.phe", data.table=FALSE)
covar_vols_female <- fread("/mnt/project/Phenotypes_MRI/regenie_input/covar_files/brain_morph_threshold_resid_norm_euro_female.cov", data.table=FALSE)
covar_OB_LR_vol_female <- fread("/mnt/project/Phenotypes_MRI/regenie_input/covar_files/brain_morph_threshold_resid_norm_euro_female_OB_LR.cov", data.table=FALSE)
covar_HT_WM_vol_female <- fread("/mnt/project/Phenotypes_MRI/regenie_input/covar_files/brain_morph_threshold_resid_norm_euro_female_GM_510_03.cov", data.table=FALSE)

pheno_vols_male <- fread("/mnt/project/Phenotypes_MRI/regenie_input/phenotype_files/brain_morph_threshold_resid_norm_euro_male.phe", data.table=FALSE)
pheno_OB_LR_vol_male <- fread("/mnt/project/Phenotypes_MRI/regenie_input/phenotype_files/brain_morph_threshold_resid_norm_euro_male_OB_LR.phe", data.table=FALSE)
pheno_HT_WM_vol_male <- fread("/mnt/project/Phenotypes_MRI/regenie_input/phenotype_files/brain_morph_threshold_resid_norm_euro_male_GM_510_03.phe", data.table=FALSE)
covar_vols_male <- fread("/mnt/project/Phenotypes_MRI/regenie_input/covar_files/brain_morph_threshold_resid_norm_euro_male.cov", data.table=FALSE)
covar_OB_LR_vol_male <- fread("/mnt/project/Phenotypes_MRI/regenie_input/covar_files/brain_morph_threshold_resid_norm_euro_male_OB_LR.cov", data.table=FALSE)
covar_HT_WM_vol_male <- fread("/mnt/project/Phenotypes_MRI/regenie_input/covar_files/brain_morph_threshold_resid_norm_euro_male_GM_510_03.cov", data.table=FALSE)

##############
###ANALYSIS###
##############

pheno_vols_HT_PG_all <- pheno_vols_all[,c("FID", "IID", "resid_HT.volume.threshold.0.3.warpResolution.2mm_norm", "resid_PG.volume.threshold.0.3.warpResolution.2mm_norm")]
pheno_vols_OB_LR_all <- pheno_OB_LR_vol_all[,c("FID", "IID", "resid_LR.volume.threshold.0.3.warpResolution.2mm_norm")]
phenos_list_all <- list(pheno_vols_HT_PG_all, pheno_vols_OB_LR_all, pheno_HT_WM_vol_all)
phenos_merge_all <- phenos_list_all %>% reduce(inner_join, by=c("FID", "IID"))

pheno_vols_HT_PG_female <- pheno_vols_female[,c("FID", "IID", "resid_HT.volume.threshold.0.3.warpResolution.2mm_norm", "resid_PG.volume.threshold.0.3.warpResolution.2mm_norm")]
pheno_vols_OB_LR_female <- pheno_OB_LR_vol_female[,c("FID", "IID", "resid_LR.volume.threshold.0.3.warpResolution.2mm_norm")]
phenos_list_female <- list(pheno_vols_HT_PG_female, pheno_vols_OB_LR_female, pheno_HT_WM_vol_female)
phenos_merge_female <- phenos_list_female %>% reduce(inner_join, by=c("FID", "IID"))

pheno_vols_HT_PG_male <- pheno_vols_male[,c("FID", "IID", "resid_HT.volume.threshold.0.3.warpResolution.2mm_norm", "resid_PG.volume.threshold.0.3.warpResolution.2mm_norm")]
pheno_vols_OB_LR_male <- pheno_OB_LR_vol_male[,c("FID", "IID", "resid_LR.volume.threshold.0.3.warpResolution.2mm_norm")]
phenos_list_male <- list(pheno_vols_HT_PG_male, pheno_vols_OB_LR_male, pheno_HT_WM_vol_male)
phenos_merge_male <- phenos_list_male %>% reduce(inner_join, by=c("FID", "IID"))

identical(covar_vols_all, covar_OB_LR_vol_all) #TRUE
identical(covar_vols_all, covar_HT_WM_vol_all) #TRUE

identical(covar_vols_female, covar_OB_LR_vol_female) #TRUE
identical(covar_vols_female, covar_HT_WM_vol_female) #TRUE

identical(covar_vols_male, covar_OB_LR_vol_male) #TRUE
identical(covar_vols_male, covar_HT_WM_vol_male) #TRUE

write.table(phenos_merge_all, "phenos_complete_vols_exome_0923_all.phe", quote=FALSE, row.names=FALSE)
write.table(phenos_merge_female, "phenos_complete_vols_exome_0923_female.phe", quote=FALSE, row.names=FALSE)
write.table(phenos_merge_male, "phenos_complete_vols_exome_0923_male.phe", quote=FALSE, row.names=FALSE)

write.table(covar_vols_all, "covars_complete_vols_exome_0923_all.cov", quote=FALSE, row.names=FALSE)
write.table(covar_vols_female, "covars_complete_vols_exome_0923_female.cov", quote=FALSE, row.names=FALSE)
write.table(covar_vols_male, "covars_complete_vols_exome_0923_male.cov", quote=FALSE, row.names=FALSE)

system("dx upload *.phe --path Phenotypes_MRI/regenie_input/phenotype_files/")
system("dx upload *.cov --path Phenotypes_MRI/regenie_input/covar_files/")
