###############
###LIBRARIES###
###############

install.packages("devtools")
install.packages(c("dplyr","tidyr","ggplot2","MASS","meta","ggrepel","DT", "rlang"))
devtools::install_github("PheWAS/PheWAS")
install.packages("tidyverse")
devtools::install_github("nstrayer/phewas_helper")
library(PheWAS)
library(data.table)
library(stringi)
library(tidyverse)
library(phewasHelper)

##########
###DATA###
##########

all_icd <- fread("/mnt/project/data/ICD10_complete_code_hesin_diag.csv", data.table = FALSE)

hard_PG_genotype <- fread("/mnt/project/data/brain_all/genotype_process/all/dosage/hard_calls/merged/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged_dosage", data.table=FALSE)
hard_HT_genotype <- fread("/mnt/project/data/brain_all/genotype_process/all/dosage/hard_calls/merged/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX_dosage", data.table=FALSE)
hard_LR_genotype <- fread("/mnt/project/data/brain_all/genotype_process/all/dosage/hard_calls/merged/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged_dosage", data.table=FALSE)
hard_HT_GM_genotype <- fread("/mnt/project/data/brain_all/genotype_process/all/dosage/hard_calls/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged_dosage", data.table = FALSE)

covars <- fread("/mnt/project/Phenotypes_MRI/regenie_input/covar_files/brain_morph_threshold_resid_norm_euro_all.cov", data.table=FALSE)
smith_confounds <- fread("/mnt/project/Phenotypes_MRI/confounds_12eBETA_2022-02-11.csv", data.table=FALSE)
bridge <- fread("/mnt/project/bridge_11867_8107.csv", data.table=FALSE)

##############
###ANALYSIS###
##############

###FORMAT ICD10 FILE
all_icd$vocabulary_id <- "ICD10CM"
all_icd$id <- all_icd$eid
all_icd$code <- all_icd$diag_icd10
stri_sub(all_icd$code, 4, 3) <- "."
all_icd$code <- gsub('^\\.|\\.$', '', all_icd$code)

all_icd <- all_icd %>% group_by(eid) %>% 
  mutate(count=ave(eid, code, FUN = length))

all_icd_distinct <- distinct_at(all_icd, vars(-dnx_hesin_diag_id))

all_icd_final <- all_icd_distinct[,c("id", "vocabulary_id", "code", "count")]

my_phenotypes <- createPhenotypes(all_icd_final, min.code.count=1)


###FORMAT GENOTYPE FILES
hard_PG_genotype$id <- hard_PG_genotype$IID
hard_PG_genotype[,c("FID","IID", "PAT", "MAT", "SEX", "PHENOTYPE")] <- NULL
hard_HT_genotype$id <- hard_HT_genotype$IID
hard_HT_genotype[,c("FID","IID", "PAT", "MAT", "SEX", "PHENOTYPE")] <- NULL
hard_LR_genotype$id <- hard_LR_genotype$IID
hard_LR_genotype[,c("FID","IID", "PAT", "MAT", "SEX", "PHENOTYPE")] <- NULL
hard_HT_GM_genotype$id <- hard_HT_GM_genotype$IID
hard_HT_GM_genotype[,c("FID","IID", "PAT", "MAT", "SEX", "PHENOTYPE")] <- NULL

###FORMAT COVAR FILES
smith_confounds_small <- smith_confounds[,c("subject_ID", "Age_Site_1", "Age_Site_2", "Age_Site_3", "Sex_1_Site_1", "Sex_1_Site_2", "Sex_1_Site_3")]
smith_confounds_small_bridge <- left_join(smith_confounds_small, bridge, by=c("subject_ID"="eid_8107"))

all_covars <- merge(smith_confounds_small_bridge, covars, by.x="eid_11867", by.y="IID")
all_covars$id <- all_covars$eid_11867
all_covars_format <- select(all_covars,-c("eid_11867", "subject_ID", "FID"))

###RUN PHEWAS
PG_results <- phewas(my_phenotypes, hard_PG_genotype, covariates=all_covars_format, cores=16, significance.threshold = c("fdr"))
saveRDS(PG_results, "PG_phewas_results.rds")
HT_results <- phewas(my_phenotypes, hard_HT_genotype, covariates=all_covars_format, cores=16, significance.threshold = c("fdr"))
LR_results <- phewas(my_phenotypes, hard_LR_genotype, covariates=all_covars_format, cores=16, significance.threshold = c("fdr"))
HT_GM_results <- phewas(my_phenotypes, hard_HT_GM_genotype, covariates=all_covars_format, cores=16, significance.threshold=c("fdr"))
saveRDS(HT_GM_results, "HT_GM_phewas_results.rds")

###READ PHEWAS 
PG_results <- readRDS("/mnt/project/data/PG_phewas_results_raw.rds")
HT_results <- readRDS("/mnt/project/data/HT_phewas_results_raw.rds")
LR_results <- readRDS("/mnt/project/data/LR_phewas_results_raw.rds")
HT_GM_results <- readRDS("/mnt/project/data/HT_GM_phewas_results.rds")

###ANNOTATE
PG_results_ann=addPhecodeInfo(PG_results)
HT_results_ann=addPhecodeInfo(HT_results)
LR_results_ann=addPhecodeInfo(LR_results)
HT_GM_results_ann=addPhecodeInfo(HT_GM_results)

###LIST SIGNIFICANT RESULTS
PG_results_ann$fdr_adj_p <- p.adjust(PG_results_ann$p, method="fdr", n=length(PG_results_ann$phenotype))
PG_results_ann$fdr_adj_sig <- PG_results_ann$fdr_adj_p <0.05
PG_sig <- subset(PG_results_ann, fdr_adj_sig==TRUE)
PG_clean <- subset(PG_results_ann, p != "NA")
PG_clean_order <- PG_clean[order(PG_clean$p),]
PG_clean_order <- PG_clean_order[,c("snp", "phenotype", "description", "group", "beta", "SE", "OR", "p", "fdr_adj_p", "n_cases", "n_controls")]
write.table(PG_clean_order, "PG_phewas_clean_order.txt", quote=FALSE, sep="\t", row.names=FALSE)
PG_clean_order_top <- head(PG_clean_order, n=20)
write.csv(PG_clean_order_top, "PG_phewas_top20.csv")

HT_results_ann$fdr_adj_p <- p.adjust(HT_results_ann$p, method="fdr", n=length(HT_results_ann$phenotype))
HT_results_ann$fdr_adj_sig <- HT_results_ann$fdr_adj_p <0.05
HT_sig <- subset(HT_results_ann, fdr_adj_sig==TRUE)
HT_clean <- subset(HT_results_ann, p != "NA")
HT_clean_order <- HT_clean[order(HT_clean$p),]
HT_clean_order <- HT_clean_order[,c("snp", "phenotype", "description", "group", "beta", "SE", "OR", "p", "fdr_adj_p", "n_cases", "n_controls")]
write.table(HT_clean_order, "HT_phewas_clean_order.txt", quote=FALSE, sep="\t", row.names=FALSE)
HT_clean_order_top <- head(HT_clean_order, n=20)
write.csv(HT_clean_order_top, "HT_phewas_top20.csv")

LR_results_ann$fdr_adj_p <- p.adjust(LR_results_ann$p, method="fdr", n=length(LR_results_ann$phenotype))
LR_results_ann$fdr_adj_sig <- LR_results_ann$fdr_adj_p <0.05
LR_sig <- subset(LR_results_ann, fdr_adj_sig==TRUE)
LR_clean <- subset(LR_results_ann, p!="NA")
LR_clean_order <- LR_clean[order(LR_clean$p),]
LR_clean_order <- LR_clean_order[,c("snp", "phenotype", "description", "group", "beta", "SE", "OR", "p", "fdr_adj_p", "n_cases", "n_controls")]
write.table(LR_clean_order, "LR_phewas_clean_order.txt", quote=FALSE, sep="\t", row.names=FALSE)
LR_clean_order_top <- head(LR_clean_order, n=20)
write.csv(LR_clean_order_top, "LR_phewas_top20.csv")


HT_GM_results_ann$fdr_adj_p <- p.adjust(HT_GM_results_ann$p, method="fdr", n=length(HT_GM_results_ann$phenotype))
HT_GM_results_ann$fdr_adj_sig <- HT_GM_results_ann$fdr_adj_p <0.05
HT_GM_sig <- subset(HT_GM_results_ann, fdr_adj_sig==TRUE)
HT_GM_clean <- subset(HT_GM_results_ann, p!="NA")
HT_GM_clean_order <- HT_GM_clean[order(HT_GM_clean$p),]
HT_GM_clean_order <- HT_GM_clean_order[,c("snp", "phenotype", "description", "group", "beta", "SE", "OR", "p", "fdr_adj_p", "n_cases", "n_controls")]
write.table(HT_GM_clean_order, "HT_GM_phewas_clean_order.txt", quote=FALSE, sep="\t", row.names=FALSE)
HT_GM_clean_order_top <- head(HT_GM_clean_order, n=20)
write.csv(HT_GM_clean_order_top, "HT_GM_phewas_tops20.csv")


system("dx upload *phewas_clean_order.txt --path data/")

###PLOT RESULTS
#PG_plot <- phewasManhattan(PG_clean, OR.direction = T)

PG_clean_order_group <- PG_clean[order(PG_clean$group),]
PG_clean_order_group$index <- 1:dim(PG_clean_order_group)[1]
pdf("PG_phewas_manhattan.pdf", width=8, height=4.5)
ggplot(PG_clean_order_group, aes(x=index, y=-log(p))) + geom_point(aes(col=group)) + theme_classic() + 
  theme(axis.text.x = element_blank(), axis.ticks=element_blank()) + 
  labs(color="Category", x="Phenotypes", y="log(p-value)") +
  scale_color_phecode() +
  theme_phewas()
dev.off()

HT_clean_order_group <- HT_clean[order(HT_clean$group),]
HT_clean_order_group$index <- 1:dim(HT_clean_order_group)[1]
pdf("HT_phewas_manhattan.pdf", width=8, height=4.5)
ggplot(HT_clean_order_group, aes(x=index, y=-log(p))) + geom_point(aes(col=group)) + theme_classic() + 
  theme(axis.text.x = element_blank(), axis.ticks=element_blank()) + 
  labs(color="Category", x="Phenotypes", y="log(p-value)") +
  scale_color_phecode() +
  theme_phewas()
dev.off()

LR_clean_order_group <- LR_clean[order(LR_clean$group),]
LR_clean_order_group$index <- 1:dim(LR_clean_order_group)[1]
pdf("LR_phewas_manhattan.pdf", width=8, height=4.5)
ggplot(LR_clean_order_group, aes(x=index, y=-log(p))) + geom_point(aes(col=group)) + theme_classic() + 
  theme(axis.text.x = element_blank(), axis.ticks=element_blank()) + 
  labs(color="Category", x="Phenotypes", y="log(p-value)") +
  scale_color_phecode() +
  theme_phewas()
dev.off()

HT_GM_clean_order_group <- HT_GM_clean[order(HT_GM_clean$group),]
HT_GM_clean_order_group$index <- 1:dim(HT_GM_clean_order_group)[1]
pdf("HT_GM_phewas_manhattan.pdf", width=8, height=4.5)
ggplot(HT_GM_clean_order_group, aes(x=index, y=-log(p))) + geom_point(aes(col=group)) + theme_classic() + 
  theme(axis.text.x = element_blank(), axis.ticks=element_blank()) + 
  labs(color="Category", x="Phenotypes", y="log(p-value)") +
  scale_color_phecode() +
  theme_phewas()
dev.off()






