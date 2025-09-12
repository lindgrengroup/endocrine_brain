###############
###LIBRARIES###
###############

devtools::install_github("boxiangliu/locuscomparer")
library(locuscomparer)
BiocManager::install("snpStats")
if(!require("remotes"))
  install.packages("remotes") # if necessary
library(remotes)
install_github("chr1swallace/coloc@main",build_vignettes=TRUE)
library(coloc)
library(data.table)
install.packages("R.utils")
library(R.utils)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(dplyr)
library(tidyr)

###############
###FUNCTIONS###
###############

coloc_run <- function(brain_snps, brain, brain_hormone){
  coloc_results <- lapply(brain_snps, function(i){
  lead_chr <- brain$hg38chr_brain[which(brain$ID_brain==i)]
  lead_pos <- brain$hg38_genpos_brain[which(brain$ID_brain==i)]
  brain_hormone_subset <- subset(brain_hormone, hg38chr_brain==lead_chr & hg38_genpos_brain>lead_pos-500000 & hg38_genpos_brain<lead_pos+500000)
  brain_hormone_subset <- brain_hormone_subset[!(duplicated(brain_hormone_subset$hg38pos_brain) | duplicated(brain_hormone_subset$hg38pos_brain, fromLast = TRUE)), ]
  coloc_df <- coloc.abf(dataset1=list(snp=brain_hormone_subset$hg38pos_brain, pvalues=brain_hormone_subset$pval_brain, beta=brain_hormone_subset$BETA_brain, varbeta=brain_hormone_subset$varbeta_brain, MAF=brain_hormone_subset$A1FREQ_brain, N=brain_hormone_subset$N_brain, type="quant"),
                        dataset2 = list(snp=brain_hormone_subset$hg38pos_brain, pvalues=brain_hormone_subset$PVALUE_hormone, beta=brain_hormone_subset$BETA_hormone, varbeta=brain_hormone_subset$varbeta_hormone, MAF=brain_hormone_subset$MAF_hormone, N=brain_hormone_subset$N_hormone, type="quant"),
                        p1 = 1e-04, p2 = 1e-04, p12 = 1e-06)
  coloc_df_res <- as.list(coloc_df$summary)
  coloc_df_res$SNP <- i
  return(coloc_df_res)
})
coloc_results_table <- do.call(rbind, coloc_results)
coloc_results_df <- as.data.frame(coloc_results_table)
coloc_results_df <- unnest(coloc_results_df, cols = c(nsnps, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf, SNP))
return(coloc_results_df)
}

coloc_infert_run <- function(brain_snps, brain, brain_infert){
  coloc_results <- lapply(brain_snps, function(i){
  lead_chr <- brain$hg38chr_brain[which(brain$ID_brain==i)]
  lead_pos <- brain$hg38_genpos_brain[which(brain$ID_brain==i)]
  brain_infert_subset <- subset(brain_infert, hg38chr_brain==lead_chr & hg38_genpos_brain>lead_pos-500000 & hg38_genpos_brain<lead_pos+500000)
  brain_infert_subset <- brain_infert_subset[!(duplicated(brain_infert_subset$hg38pos_brain) | duplicated(brain_infert_subset$hg38pos_brain, fromLast = TRUE)), ]
  brain_infert_subset$N_infert <- brain_infert_subset$N_CASES_infert + brain_infert_subset$N_CONTROLS_infert
  brain_infert_subset$case_prop_infert <- brain_infert_subset$N_CASES_infert/brain_infert_subset$N_infert
  brain_infert_subset <- subset(brain_infert_subset, MAF_infert>0 & MAF_infert<1)
  coloc_df <- coloc.abf(dataset1=list(snp=brain_infert_subset$hg38pos_brain, pvalues=brain_infert_subset$pval_brain, beta=brain_infert_subset$BETA_brain, varbeta=brain_infert_subset$varbeta_brain, MAF=brain_infert_subset$MAF_brain, N=brain_infert_subset$N_brain, type="quant"),
                        dataset2 = list(snp=brain_infert_subset$hg38pos_brain, pvalues=brain_infert_subset$PVALUE_infert, beta=brain_infert_subset$BETA_infert, varbeta=brain_infert_subset$varbeta_infert, MAF=brain_infert_subset$MAF_infert, N=brain_infert_subset$N_infert, s=brain_infert_subset$case_prop_infert, type="cc"),
                        p1 = 1e-04, p2 = 1e-04, p12 = 1e-06)
  coloc_df_res <- as.list(coloc_df$summary)
  coloc_df_res$SNP <- i
  return(coloc_df_res)
})
coloc_results_table <- do.call(rbind, coloc_results)
coloc_results_df <- as.data.frame(coloc_results_table)
coloc_results_df <- unnest(coloc_results_df, cols = c(nsnps, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf, SNP))
return(coloc_results_df)
}

##########
###DATA###
##########

PG_vol <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/liftover/output/formatted/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.hg38.txt", data.table=FALSE)
PG_vol_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
HT_vol <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/liftover/output/formatted/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.hg38.txt", data.table=FALSE)
HT_vol_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX.txt", data.table=FALSE)
LR_vol <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/liftover/output/formatted/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.hg38.txt", data.table=FALSE)
LR_vol_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
HT_GM_vol <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/liftover/output/formatted/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.hg38.txt", data.table=FALSE)
HT_GM_vol_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)

PG_vol_F <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/liftover/output/formatted/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.hg38.txt", data.table=FALSE)
PG_vol_cojo_F <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
HT_vol_F <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/liftover/output/formatted/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.hg38.txt", data.table=FALSE)
HT_vol_cojo_F <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
LR_vol_F <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/liftover/output/formatted/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.hg38.txt", data.table=FALSE)
LR_vol_cojo_F <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)
HT_GM_vol_F <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/liftover/output/formatted/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.hg38.txt", data.table=FALSE)
HT_GM_vol_cojo_F <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)

PG_vol_M <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/liftover/output/formatted/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.hg38.txt", data.table=FALSE)
PG_vol_cojo_M <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
HT_vol_M <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/liftover/output/formatted/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.hg38.txt", data.table=FALSE)
HT_vol_cojo_M <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
LR_vol_M <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/liftover/output/formatted/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.hg38.txt", data.table=FALSE)
LR_vol_cojo_M <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)
HT_GM_vol_M <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/liftover/output/formatted/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.hg38.txt", data.table=FALSE)
HT_GM_vol_cojo_M <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)

test_all_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/Testosterone_sex_comb_EUR_filtered.txt.gz", data.table=FALSE)
fsh_all_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/FSH_sex_comb_EUR_filtered.txt.gz", data.table=FALSE)
lh_all_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/LH_sex_comb_EUR_filtered.txt.gz", data.table=FALSE)
prog_all_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/Progesterone_sex_comb_EUR_filtered.txt.gz", data.table=FALSE)
oest_all_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/Oestradiol_sex_comb_EUR_filtered.txt.gz", data.table=FALSE)

test_female_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/Testosterone_F_EUR_filtered.txt.gz", data.table=FALSE)
fsh_female_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/FSH_F_EUR_filtered.txt.gz", data.table=FALSE)
lh_female_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/LH_F_EUR_filtered.txt.gz", data.table=FALSE)
prog_female_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/Progesterone_F_EUR_filtered.txt.gz", data.table=FALSE)
oest_female_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/Oestradiol_F_EUR_filtered.txt.gz", data.table=FALSE)

test_male_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/Testosterone_M_EUR_filtered.txt.gz", data.table=FALSE)
fsh_male_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/FSH_M_EUR_filtered.txt.gz", data.table=FALSE)
lh_male_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/LH_M_EUR_filtered.txt.gz", data.table=FALSE)
oest_male_eur <- fread("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/Oestradiol_M_EUR_filtered.txt.gz", data.table=FALSE)

female_infert1 <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/female_infertility_analysis1_eur_with_rsids.txt.gz"), header=TRUE)
female_infert2 <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/female_infertility_analysis2_eur_with_rsids.txt.gz"), header=TRUE)
female_infert3 <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/female_infertility_analysis3_eur_with_rsids.txt.gz"), header=TRUE)
female_infert4 <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/female_infertility_analysis4_eur_with_rsids.txt.gz"), header=TRUE)
female_infert5 <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/female_infertility_analysis5_eur_with_rsids.txt.gz"), header=TRUE)
male_infert <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/male_infertility_eur_with_rsids.txt.gz"), header=TRUE)


#chain <- import.chain("hg38ToHg19.over.chain")
##############
###ANALYSIS###
##############

PG_vol$hg38chr <- gsub(".*chr(.+):.*", "\\1", PG_vol$hg38pos)
PG_vol$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", PG_vol$hg38pos)))
PG_vol$pval <- 10^(-PG_vol$LOG10P)
PG_vol$varbeta <- (PG_vol$SE)^2
PG_vol$MAF <- NA
PG_vol$MAF[which(PG_vol$A1FREQ < 0.5)] <- PG_vol$A1FREQ[which(PG_vol$A1FREQ < 0.5)]
PG_vol$MAF[which(PG_vol$A1FREQ >= 0.5)] <- 1-PG_vol$A1FREQ[which(PG_vol$A1FREQ >= 0.5)]
colnames(PG_vol) <- paste(colnames(PG_vol), "brain", sep="_")
PG_vol$hg38pos_brain <- gsub("chrX", "chr23", PG_vol$hg38pos_brain)

HT_vol$hg38chr <- gsub(".*chr(.+):.*", "\\1", HT_vol$hg38pos)
HT_vol$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", HT_vol$hg38pos)))
HT_vol$pval <- 10^(-HT_vol$LOG10P)
HT_vol$varbeta <- (HT_vol$SE)^2
HT_vol$MAF <- NA
HT_vol$MAF[which(HT_vol$A1FREQ < 0.5)] <- HT_vol$A1FREQ[which(HT_vol$A1FREQ < 0.5)]
HT_vol$MAF[which(HT_vol$A1FREQ >= 0.5)] <- 1-HT_vol$A1FREQ[which(HT_vol$A1FREQ >= 0.5)]
colnames(HT_vol) <- paste(colnames(HT_vol), "brain", sep="_")
HT_vol$hg38pos_brain <- gsub("chrX", "chr23", HT_vol$hg38pos_brain)

LR_vol$hg38chr <- gsub(".*chr(.+):.*", "\\1", LR_vol$hg38pos)
LR_vol$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", LR_vol$hg38pos)))
LR_vol$pval <- 10^(-LR_vol$LOG10P)
LR_vol$varbeta <- (LR_vol$SE)^2
LR_vol$MAF <- NA
LR_vol$MAF[which(LR_vol$A1FREQ < 0.5)] <- LR_vol$A1FREQ[which(LR_vol$A1FREQ < 0.5)]
LR_vol$MAF[which(LR_vol$A1FREQ >= 0.5)] <- 1-LR_vol$A1FREQ[which(LR_vol$A1FREQ >= 0.5)]
colnames(LR_vol) <- paste(colnames(LR_vol), "brain", sep="_")
LR_vol$hg38pos_brain <- gsub("chrX", "chr23", LR_vol$hg38pos_brain)

HT_GM_vol$hg38chr <- gsub(".*chr(.+):.*", "\\1", HT_GM_vol$hg38pos)
HT_GM_vol$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", HT_GM_vol$hg38pos)))
HT_GM_vol$pval <- 10^(-HT_GM_vol$LOG10P)
HT_GM_vol$varbeta <- (HT_GM_vol$SE)^2
HT_GM_vol$MAF <- NA
HT_GM_vol$MAF[which(HT_GM_vol$A1FREQ < 0.5)] <- HT_GM_vol$A1FREQ[which(HT_GM_vol$A1FREQ < 0.5)]
HT_GM_vol$MAF[which(HT_GM_vol$A1FREQ >= 0.5)] <- 1-HT_GM_vol$A1FREQ[which(HT_GM_vol$A1FREQ >= 0.5)]
colnames(HT_GM_vol) <- paste(colnames(HT_GM_vol), "brain", sep="_")
HT_GM_vol$hg38pos_brain <- gsub("chrX", "chr23", HT_GM_vol$hg38pos_brain)

PG_vol_F$hg38chr <- gsub(".*chr(.+):.*", "\\1", PG_vol_F$hg38pos)
PG_vol_F$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", PG_vol_F$hg38pos)))
PG_vol_F$pval <- 10^(-PG_vol_F$LOG10P)
PG_vol_F$varbeta <- (PG_vol_F$SE)^2
PG_vol_F$MAF <- NA
PG_vol_F$MAF[which(PG_vol_F$A1FREQ < 0.5)] <- PG_vol_F$A1FREQ[which(PG_vol_F$A1FREQ < 0.5)]
PG_vol_F$MAF[which(PG_vol_F$A1FREQ >= 0.5)] <- 1-PG_vol_F$A1FREQ[which(PG_vol_F$A1FREQ >= 0.5)]
colnames(PG_vol_F) <- paste(colnames(PG_vol_F), "brain", sep="_")
PG_vol_F$hg38pos_brain <- gsub("chrX", "chr23", PG_vol_F$hg38pos_brain)

HT_vol_F$hg38chr <- gsub(".*chr(.+):.*", "\\1", HT_vol_F$hg38pos)
HT_vol_F$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", HT_vol_F$hg38pos)))
HT_vol_F$pval <- 10^(-HT_vol_F$LOG10P)
HT_vol_F$varbeta <- (HT_vol_F$SE)^2
HT_vol_F$MAF <- NA
HT_vol_F$MAF[which(HT_vol_F$A1FREQ < 0.5)] <- HT_vol_F$A1FREQ[which(HT_vol_F$A1FREQ < 0.5)]
HT_vol_F$MAF[which(HT_vol_F$A1FREQ >= 0.5)] <- 1-HT_vol_F$A1FREQ[which(HT_vol_F$A1FREQ >= 0.5)]
colnames(HT_vol_F) <- paste(colnames(HT_vol_F), "brain", sep="_")
HT_vol_F$hg38pos_brain <- gsub("chrX", "chr23", HT_vol_F$hg38pos_brain)

LR_vol_F$hg38chr <- gsub(".*chr(.+):.*", "\\1", LR_vol_F$hg38pos)
LR_vol_F$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", LR_vol_F$hg38pos)))
LR_vol_F$pval <- 10^(-LR_vol_F$LOG10P)
LR_vol_F$varbeta <- (LR_vol_F$SE)^2
LR_vol_F$MAF <- NA
LR_vol_F$MAF[which(LR_vol_F$A1FREQ < 0.5)] <- LR_vol_F$A1FREQ[which(LR_vol_F$A1FREQ < 0.5)]
LR_vol_F$MAF[which(LR_vol_F$A1FREQ >= 0.5)] <- 1-LR_vol_F$A1FREQ[which(LR_vol_F$A1FREQ >= 0.5)]
colnames(LR_vol_F) <- paste(colnames(LR_vol_F), "brain", sep="_")
LR_vol_F$hg38pos_brain <- gsub("chrX", "chr23", LR_vol_F$hg38pos_brain)

HT_GM_vol_F$hg38chr <- gsub(".*chr(.+):.*", "\\1", HT_GM_vol_F$hg38pos)
HT_GM_vol_F$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", HT_GM_vol_F$hg38pos)))
HT_GM_vol_F$pval <- 10^(-HT_GM_vol_F$LOG10P)
HT_GM_vol_F$varbeta <- (HT_GM_vol_F$SE)^2
HT_GM_vol_F$MAF <- NA
HT_GM_vol_F$MAF[which(HT_GM_vol_F$A1FREQ < 0.5)] <- HT_GM_vol_F$A1FREQ[which(HT_GM_vol_F$A1FREQ < 0.5)]
HT_GM_vol_F$MAF[which(HT_GM_vol_F$A1FREQ >= 0.5)] <- 1-HT_GM_vol_F$A1FREQ[which(HT_GM_vol_F$A1FREQ >= 0.5)]
colnames(HT_GM_vol_F) <- paste(colnames(HT_GM_vol_F), "brain", sep="_")
HT_GM_vol_F$hg38pos_brain <- gsub("chrX", "chr23", HT_GM_vol_F$hg38pos_brain)

PG_vol_M$hg38chr <- gsub(".*chr(.+):.*", "\\1", PG_vol_M$hg38pos)
PG_vol_M$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", PG_vol_M$hg38pos)))
PG_vol_M$pval <- 10^(-PG_vol_M$LOG10P)
PG_vol_M$varbeta <- (PG_vol_M$SE)^2
PG_vol_M$MAF <- NA
PG_vol_M$MAF[which(PG_vol_M$A1FREQ < 0.5)] <- PG_vol_M$A1FREQ[which(PG_vol_M$A1FREQ < 0.5)]
PG_vol_M$MAF[which(PG_vol_M$A1FREQ >= 0.5)] <- 1-PG_vol_M$A1FREQ[which(PG_vol_M$A1FREQ >= 0.5)]
colnames(PG_vol_M) <- paste(colnames(PG_vol_M), "brain", sep="_")
PG_vol_M$hg38pos_brain <- gsub("chrX", "chr23", PG_vol_M$hg38pos_brain)

HT_vol_M$hg38chr <- gsub(".*chr(.+):.*", "\\1", HT_vol_M$hg38pos)
HT_vol_M$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", HT_vol_M$hg38pos)))
HT_vol_M$pval <- 10^(-HT_vol_M$LOG10P)
HT_vol_M$varbeta <- (HT_vol_M$SE)^2
HT_vol_M$MAF <- NA
HT_vol_M$MAF[which(HT_vol_M$A1FREQ < 0.5)] <- HT_vol_M$A1FREQ[which(HT_vol_M$A1FREQ < 0.5)]
HT_vol_M$MAF[which(HT_vol_M$A1FREQ >= 0.5)] <- 1-HT_vol_M$A1FREQ[which(HT_vol_M$A1FREQ >= 0.5)]
colnames(HT_vol_M) <- paste(colnames(HT_vol_M), "brain", sep="_")
HT_vol_M$hg38pos_brain <- gsub("chrX", "chr23", HT_vol_M$hg38pos_brain)

LR_vol_M$hg38chr <- gsub(".*chr(.+):.*", "\\1", LR_vol_M$hg38pos)
LR_vol_M$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", LR_vol_M$hg38pos)))
LR_vol_M$pval <- 10^(-LR_vol_M$LOG10P)
LR_vol_M$varbeta <- (LR_vol_M$SE)^2
LR_vol_M$MAF <- NA
LR_vol_M$MAF[which(LR_vol_M$A1FREQ < 0.5)] <- LR_vol_M$A1FREQ[which(LR_vol_M$A1FREQ < 0.5)]
LR_vol_M$MAF[which(LR_vol_M$A1FREQ >= 0.5)] <- 1-LR_vol_M$A1FREQ[which(LR_vol_M$A1FREQ >= 0.5)]
colnames(LR_vol_M) <- paste(colnames(LR_vol_M), "brain", sep="_")
LR_vol_M$hg38pos_brain <- gsub("chrX", "chr23", LR_vol_M$hg38pos_brain)


HT_GM_vol_M$hg38chr <- gsub(".*chr(.+):.*", "\\1", HT_GM_vol_M$hg38pos)
HT_GM_vol_M$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", HT_GM_vol_M$hg38pos)))
HT_GM_vol_M$pval <- 10^(-HT_GM_vol_M$LOG10P)
HT_GM_vol_M$varbeta <- (HT_GM_vol_M$SE)^2
HT_GM_vol_M$MAF <- NA
HT_GM_vol_M$MAF[which(HT_GM_vol_M$A1FREQ < 0.5)] <- HT_GM_vol_M$A1FREQ[which(HT_GM_vol_M$A1FREQ < 0.5)]
HT_GM_vol_M$MAF[which(HT_GM_vol_M$A1FREQ >= 0.5)] <- 1-HT_GM_vol_M$A1FREQ[which(HT_GM_vol_M$A1FREQ >= 0.5)]
colnames(HT_GM_vol_M) <- paste(colnames(HT_GM_vol_M), "brain", sep="_")
HT_GM_vol_M$hg38pos_brain <- gsub("chrX", "chr23", HT_GM_vol_M$hg38pos_brain)

test_all_eur$varbeta <- (test_all_eur$SE)^2
test_all_eur$N <- 69666
colnames(test_all_eur) <- paste(colnames(test_all_eur), "hormone", sep = "_")
fsh_all_eur$varbeta <- (fsh_all_eur$SE)^2
fsh_all_eur$N <- 32769
colnames(fsh_all_eur) <- paste(colnames(fsh_all_eur), "hormone", sep = "_")
lh_all_eur$varbeta <- (lh_all_eur$SE)^2
lh_all_eur$N <- 28015
colnames(lh_all_eur) <- paste(colnames(lh_all_eur), "hormone", sep = "_")
prog_all_eur$varbeta <- (prog_all_eur$SE)^2
prog_all_eur$N <- 16171
colnames(prog_all_eur) <- paste(colnames(prog_all_eur), "hormone", sep = "_")
oest_all_eur$varbeta <- (oest_all_eur$SE)^2
oest_all_eur$N <- 60249
colnames(oest_all_eur) <- paste(colnames(oest_all_eur), "hormone", sep = "_")

test_female_eur$varbeta <- (test_female_eur$SE)^2
test_female_eur$N <- 38541
colnames(test_female_eur) <- paste(colnames(test_female_eur), "hormone", sep = "_")
fsh_female_eur$varbeta <- (fsh_female_eur$SE)^2
fsh_female_eur$N <- 30596
colnames(fsh_female_eur) <- paste(colnames(fsh_female_eur), "hormone", sep = "_")
lh_female_eur$varbeta <- (lh_female_eur$SE)^2
lh_female_eur$N <- 25259
colnames(lh_female_eur) <- paste(colnames(lh_female_eur), "hormone", sep = "_")
oest_female_eur$varbeta <- (oest_female_eur$SE)^2
oest_female_eur$N <- 38360
colnames(oest_female_eur) <- paste(colnames(oest_female_eur), "hormone", sep = "_")
prog_female_eur$varbeta <- (prog_female_eur$SE)^2
prog_female_eur$N <- 14558
colnames(prog_female_eur) <- paste(colnames(prog_female_eur), "hormone", sep = "_")

test_male_eur$varbeta <- (test_male_eur$SE)^2
test_male_eur$N <- 37756
colnames(test_male_eur) <- paste(colnames(test_male_eur), "hormone", sep = "_")
fsh_male_eur$varbeta <- (fsh_male_eur$SE)^2
fsh_male_eur$N <- 3981
colnames(fsh_male_eur) <- paste(colnames(fsh_male_eur), "hormone", sep = "_")
lh_male_eur$varbeta <- (lh_male_eur$SE)^2
lh_male_eur$N <- 4532
colnames(lh_male_eur) <- paste(colnames(lh_male_eur), "hormone", sep = "_")
oest_male_eur$varbeta <- (oest_male_eur$SE)^2
oest_male_eur$N <- 21674
colnames(oest_male_eur) <- paste(colnames(oest_male_eur), "hormone", sep = "_")

female_infert1$varbeta <- (female_infert1$SE)^2
female_infert1$MAF <- NA
female_infert1$MAF[which(female_infert1$Freq1 < 0.5)] <- female_infert1$Freq1[which(female_infert1$Freq1 < 0.5)]
female_infert1$MAF[which(female_infert1$Freq1 >= 0.5)] <- 1-female_infert1$Freq1[which(female_infert1$Freq1 >= 0.5)]
colnames(female_infert1) <- paste(colnames(female_infert1), "infert", sep = "_")
female_infert1$ID_infert <- gsub("chrX", "chr23", female_infert1$ID_infert)
female_infert2$varbeta <- (female_infert2$SE)^2
female_infert2$MAF <- NA
female_infert2$MAF[which(female_infert2$Freq1 < 0.5)] <- female_infert2$Freq1[which(female_infert2$Freq1 < 0.5)]
female_infert2$MAF[which(female_infert2$Freq1 >= 0.5)] <- 1-female_infert2$Freq1[which(female_infert2$Freq1 >= 0.5)]
colnames(female_infert2) <- paste(colnames(female_infert2), "infert", sep = "_")
female_infert2$ID_infert <- gsub("chrX", "chr23", female_infert2$ID_infert)
female_infert3$varbeta <- (female_infert3$SE)^2
female_infert3$MAF <- NA
female_infert3$MAF[which(female_infert3$Freq1 < 0.5)] <- female_infert3$Freq1[which(female_infert3$Freq1 < 0.5)]
female_infert3$MAF[which(female_infert3$Freq1 >= 0.5)] <- 1-female_infert3$Freq1[which(female_infert3$Freq1 >= 0.5)]
colnames(female_infert3) <- paste(colnames(female_infert3), "infert", sep = "_")
female_infert3$ID_infert <- gsub("chrX", "chr23", female_infert3$ID_infert)
female_infert4$varbeta <- (female_infert4$SE)^2
female_infert4$MAF <- NA
female_infert4$MAF[which(female_infert4$Freq1 < 0.5)] <- female_infert4$Freq1[which(female_infert4$Freq1 < 0.5)]
female_infert4$MAF[which(female_infert4$Freq1 >= 0.5)] <- 1-female_infert4$Freq1[which(female_infert4$Freq1 >= 0.5)]
colnames(female_infert4) <- paste(colnames(female_infert4), "infert", sep = "_")
female_infert4$ID_infert <- gsub("chrX", "chr23", female_infert4$ID_infert)
female_infert5$varbeta <- (female_infert5$SE)^2
female_infert5$MAF <- NA
female_infert5$MAF[which(female_infert5$Freq1 < 0.5)] <- female_infert5$Freq1[which(female_infert5$Freq1 < 0.5)]
female_infert5$MAF[which(female_infert5$Freq1 >= 0.5)] <- 1-female_infert5$Freq1[which(female_infert5$Freq1 >= 0.5)]
colnames(female_infert5) <- paste(colnames(female_infert5), "infert", sep = "_")
female_infert5$ID_infert <- gsub("chrX", "chr23", female_infert5$ID_infert)
male_infert$varbeta <- (male_infert$SE)^2
male_infert$MAF <- NA
male_infert$MAF[which(male_infert$Freq1 < 0.5)] <- male_infert$Freq1[which(male_infert$Freq1 < 0.5)]
male_infert$MAF[which(male_infert$Freq1 >= 0.5)] <- 1-male_infert$Freq1[which(male_infert$Freq1 >= 0.5)]
colnames(male_infert) <- paste(colnames(male_infert), "infert", sep = "_")
male_infert$ID_infert <- gsub("chrX", "chr23", male_infert$ID_infert)


PG_vol_test <- merge(PG_vol, test_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_test <- merge(HT_vol, test_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_test <- merge(LR_vol, test_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_test <- merge(HT_GM_vol, test_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_F_test <- merge(PG_vol_F, test_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_F_test <- merge(HT_vol_F, test_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_F_test <- merge(LR_vol_F, test_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_F_test <- merge(HT_GM_vol_F, test_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_M_test <- merge(PG_vol_M, test_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_M_test <- merge(HT_vol_M, test_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_M_test <- merge(LR_vol_M, test_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_M_test <- merge(HT_GM_vol_M, test_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_fsh <- merge(PG_vol, fsh_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_fsh <- merge(HT_vol, fsh_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_fsh <- merge(LR_vol, fsh_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_fsh <- merge(HT_GM_vol, fsh_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_F_fsh <- merge(PG_vol_F, fsh_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_F_fsh <- merge(HT_vol_F, fsh_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_F_fsh <- merge(LR_vol_F, fsh_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_F_fsh <- merge(HT_GM_vol_F, fsh_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_M_fsh <- merge(PG_vol_M, fsh_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_M_fsh <- merge(HT_vol_M, fsh_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_M_fsh <- merge(LR_vol_M, fsh_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_M_fsh <- merge(HT_GM_vol_M, fsh_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_lh <- merge(PG_vol, lh_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_lh <- merge(HT_vol, lh_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_lh <- merge(LR_vol, lh_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_lh <- merge(HT_GM_vol, lh_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_F_lh <- merge(PG_vol_F, lh_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_F_lh <- merge(HT_vol_F, lh_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_F_lh <- merge(LR_vol_F, lh_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_F_lh <- merge(HT_GM_vol_F, lh_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_M_lh <- merge(PG_vol_M, lh_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_M_lh <- merge(HT_vol_M, lh_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_M_lh <- merge(LR_vol_M, lh_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_M_lh <- merge(HT_GM_vol_M, lh_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_oest <- merge(PG_vol, oest_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_oest <- merge(HT_vol, oest_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_oest <- merge(LR_vol, oest_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_oest <- merge(HT_GM_vol, oest_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_F_oest <- merge(PG_vol_F, oest_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_F_oest <- merge(HT_vol_F, oest_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_F_oest <- merge(LR_vol_F, oest_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_F_oest <- merge(HT_GM_vol_F, oest_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_M_oest <- merge(PG_vol_M, oest_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_M_oest <- merge(HT_vol_M, oest_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_M_oest <- merge(LR_vol_M, oest_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_M_oest <- merge(HT_GM_vol_M, oest_male_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_prog <- merge(PG_vol, prog_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_prog <- merge(HT_vol, prog_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_prog <- merge(LR_vol, prog_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_prog <- merge(HT_GM_vol, prog_all_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_F_prog <- merge(PG_vol_F, prog_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_vol_F_prog <- merge(HT_vol_F, prog_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
LR_vol_F_prog <- merge(LR_vol_F, prog_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")
HT_GM_vol_F_prog <- merge(HT_GM_vol_F, prog_female_eur, by.x="hg38pos_brain", by.y="ID_hormone")

PG_vol_F_infert1 <- merge(PG_vol_F, female_infert1, by.x="hg38pos_brain", by.y="ID_infert")
HT_vol_F_infert1 <- merge(HT_vol_F, female_infert1, by.x="hg38pos_brain", by.y="ID_infert")
LR_vol_F_infert1 <- merge(LR_vol_F, female_infert1, by.x="hg38pos_brain", by.y="ID_infert")
HT_GM_vol_F_infert1 <- merge(HT_GM_vol_F, female_infert1, by.x="hg38pos_brain", by.y="ID_infert")

PG_vol_F_infert2 <- merge(PG_vol_F, female_infert2, by.x="hg38pos_brain", by.y="ID_infert")
HT_vol_F_infert2 <- merge(HT_vol_F, female_infert2, by.x="hg38pos_brain", by.y="ID_infert")
LR_vol_F_infert2 <- merge(LR_vol_F, female_infert2, by.x="hg38pos_brain", by.y="ID_infert")
HT_GM_vol_F_infert2 <- merge(HT_GM_vol_F, female_infert2, by.x="hg38pos_brain", by.y="ID_infert")

PG_vol_F_infert3 <- merge(PG_vol_F, female_infert3, by.x="hg38pos_brain", by.y="ID_infert")
HT_vol_F_infert3 <- merge(HT_vol_F, female_infert3, by.x="hg38pos_brain", by.y="ID_infert")
LR_vol_F_infert3 <- merge(LR_vol_F, female_infert3, by.x="hg38pos_brain", by.y="ID_infert")
HT_GM_vol_F_infert3 <- merge(HT_GM_vol_F, female_infert3, by.x="hg38pos_brain", by.y="ID_infert")

PG_vol_F_infert4 <- merge(PG_vol_F, female_infert4, by.x="hg38pos_brain", by.y="ID_infert")
HT_vol_F_infert4 <- merge(HT_vol_F, female_infert4, by.x="hg38pos_brain", by.y="ID_infert")
LR_vol_F_infert4 <- merge(LR_vol_F, female_infert4, by.x="hg38pos_brain", by.y="ID_infert")
HT_GM_vol_F_infert4 <- merge(HT_GM_vol_F, female_infert4, by.x="hg38pos_brain", by.y="ID_infert")

PG_vol_F_infert5 <- merge(PG_vol_F, female_infert5, by.x="hg38pos_brain", by.y="ID_infert")
HT_vol_F_infert5 <- merge(HT_vol_F, female_infert5, by.x="hg38pos_brain", by.y="ID_infert")
LR_vol_F_infert5 <- merge(LR_vol_F, female_infert5, by.x="hg38pos_brain", by.y="ID_infert")
HT_GM_vol_F_infert5 <- merge(HT_GM_vol_F, female_infert5, by.x="hg38pos_brain", by.y="ID_infert")

PG_vol_M_infert <- merge(PG_vol_M, male_infert, by.x="hg38pos_brain", by.y="ID_infert")
HT_vol_M_infert <- merge(HT_vol_M, male_infert, by.x="hg38pos_brain", by.y="ID_infert")
LR_vol_M_infert <- merge(LR_vol_M, male_infert, by.x="hg38pos_brain", by.y="ID_infert")
HT_GM_vol_M_infert <- merge(HT_GM_vol_M, male_infert, by.x="hg38pos_brain", by.y="ID_infert")

PG_test_coloc <- coloc_run(PG_vol_cojo$SNP, PG_vol, PG_vol_test)
HT_test_coloc <- coloc_run(HT_vol_cojo$SNP, HT_vol, HT_vol_test)
LR_test_coloc <- coloc_run(LR_vol_cojo$SNP, LR_vol, LR_vol_test)
HT_GM_test_coloc <- coloc_run(HT_GM_vol_cojo$SNP, HT_GM_vol, HT_GM_vol_test)

PG_test_F_coloc <- coloc_run(PG_vol_cojo_F$SNP, PG_vol_F, PG_vol_F_test)
HT_test_F_coloc <- coloc_run(HT_vol_cojo_F$SNP, HT_vol_F, HT_vol_F_test)
LR_test_F_coloc <- coloc_run(LR_vol_cojo_F$SNP, LR_vol_F, LR_vol_F_test)
HT_GM_test_F_coloc <- coloc_run(HT_GM_vol_cojo_F$SNP, HT_GM_vol_F, HT_GM_vol_F_test)

PG_test_M_coloc <- coloc_run(PG_vol_cojo_M$SNP, PG_vol_M, PG_vol_M_test)
HT_test_M_coloc <- coloc_run(HT_vol_cojo_M$SNP, HT_vol_M, HT_vol_M_test)
LR_test_M_coloc <- coloc_run(LR_vol_cojo_M$SNP, LR_vol_M, LR_vol_M_test)
HT_GM_test_M_coloc <- coloc_run(HT_GM_vol_cojo_M$SNP, HT_GM_vol_M, HT_GM_vol_M_test)

PG_fsh_coloc <- coloc_run(PG_vol_cojo$SNP, PG_vol, PG_vol_fsh)
HT_fsh_coloc <- coloc_run(HT_vol_cojo$SNP, HT_vol, HT_vol_fsh)
LR_fsh_coloc <- coloc_run(LR_vol_cojo$SNP, LR_vol, LR_vol_fsh)
HT_GM_fsh_coloc <- coloc_run(HT_GM_vol_cojo$SNP, HT_GM_vol, HT_GM_vol_fsh)

PG_fsh_F_coloc <- coloc_run(PG_vol_cojo_F$SNP, PG_vol_F, PG_vol_F_fsh)
HT_fsh_F_coloc <- coloc_run(HT_vol_cojo_F$SNP, HT_vol_F, HT_vol_F_fsh)
LR_fsh_F_coloc <- coloc_run(LR_vol_cojo_F$SNP, LR_vol_F, LR_vol_F_fsh)
HT_GM_fsh_F_coloc <- coloc_run(HT_GM_vol_cojo_F$SNP, HT_GM_vol_F, HT_GM_vol_F_fsh)

PG_fsh_M_coloc <- coloc_run(PG_vol_cojo_M$SNP, PG_vol_M, PG_vol_M_fsh)
HT_fsh_M_coloc <- coloc_run(HT_vol_cojo_M$SNP, HT_vol_M, HT_vol_M_fsh)
LR_fsh_M_coloc <- coloc_run(LR_vol_cojo_M$SNP, LR_vol_M, LR_vol_M_fsh)
HT_GM_fsh_M_coloc <- coloc_run(HT_GM_vol_cojo_M$SNP, HT_GM_vol_M, HT_GM_vol_M_fsh)

PG_lh_coloc <- coloc_run(PG_vol_cojo$SNP, PG_vol, PG_vol_lh)
HT_lh_coloc <- coloc_run(HT_vol_cojo$SNP, HT_vol, HT_vol_lh)
LR_lh_coloc <- coloc_run(LR_vol_cojo$SNP, LR_vol, LR_vol_lh)
HT_GM_lh_coloc <- coloc_run(HT_GM_vol_cojo$SNP, HT_GM_vol, HT_GM_vol_lh)

PG_lh_F_coloc <- coloc_run(PG_vol_cojo_F$SNP, PG_vol_F, PG_vol_F_lh)
HT_lh_F_coloc <- coloc_run(HT_vol_cojo_F$SNP, HT_vol_F, HT_vol_F_lh)
LR_lh_F_coloc <- coloc_run(LR_vol_cojo_F$SNP, LR_vol_F, LR_vol_F_lh)
HT_GM_lh_F_coloc <- coloc_run(HT_GM_vol_cojo_F$SNP, HT_GM_vol_F, HT_GM_vol_F_lh)

PG_lh_M_coloc <- coloc_run(PG_vol_cojo_M$SNP, PG_vol_M, PG_vol_M_lh)
HT_lh_M_coloc <- coloc_run(HT_vol_cojo_M$SNP, HT_vol_M, HT_vol_M_lh)
LR_lh_M_coloc <- coloc_run(LR_vol_cojo_M$SNP, LR_vol_M, LR_vol_M_lh)
HT_GM_lh_M_coloc <- coloc_run(HT_GM_vol_cojo_M$SNP, HT_GM_vol_M, HT_GM_vol_M_lh)

PG_oest_coloc <- coloc_run(PG_vol_cojo$SNP, PG_vol, PG_vol_oest)
HT_oest_coloc <- coloc_run(HT_vol_cojo$SNP, HT_vol, HT_vol_oest)
LR_oest_coloc <- coloc_run(LR_vol_cojo$SNP, LR_vol, LR_vol_oest)
HT_GM_oest_coloc <- coloc_run(HT_GM_vol_cojo$SNP, HT_GM_vol, HT_GM_vol_oest)

PG_oest_F_coloc <- coloc_run(PG_vol_cojo_F$SNP, PG_vol_F, PG_vol_F_oest)
HT_oest_F_coloc <- coloc_run(HT_vol_cojo_F$SNP, HT_vol_F, HT_vol_F_oest)
LR_oest_F_coloc <- coloc_run(LR_vol_cojo_F$SNP, LR_vol_F, LR_vol_F_oest)
HT_GM_oest_F_coloc <- coloc_run(HT_GM_vol_cojo_F$SNP, HT_GM_vol_F, HT_GM_vol_F_oest)

PG_oest_M_coloc <- coloc_run(PG_vol_cojo_M$SNP, PG_vol_M, PG_vol_M_oest)
HT_oest_M_coloc <- coloc_run(HT_vol_cojo_M$SNP, HT_vol_M, HT_vol_M_oest)
LR_oest_M_coloc <- coloc_run(LR_vol_cojo_M$SNP, LR_vol_M, LR_vol_M_oest)
HT_GM_oest_M_coloc <- coloc_run(HT_GM_vol_cojo_M$SNP, HT_GM_vol_M, HT_GM_vol_M_oest)

PG_prog_coloc <- coloc_run(PG_vol_cojo$SNP, PG_vol, PG_vol_prog)
HT_prog_coloc <- coloc_run(HT_vol_cojo$SNP, HT_vol, HT_vol_prog)
LR_prog_coloc <- coloc_run(LR_vol_cojo$SNP, LR_vol, LR_vol_prog)
HT_GM_prog_coloc <- coloc_run(HT_GM_vol_cojo$SNP, HT_GM_vol, HT_GM_vol_prog)

PG_prog_F_coloc <- coloc_run(PG_vol_cojo_F$SNP, PG_vol_F, PG_vol_F_prog)
HT_prog_F_coloc <- coloc_run(HT_vol_cojo_F$SNP, HT_vol_F, HT_vol_F_prog)
LR_prog_F_coloc <- coloc_run(LR_vol_cojo_F$SNP, LR_vol_F, LR_vol_F_prog)
HT_GM_prog_F_coloc <- coloc_run(HT_GM_vol_cojo_F$SNP, HT_GM_vol_F, HT_GM_vol_F_prog)

PG_infert1_F_coloc <- coloc_infert_run(PG_vol_cojo_F$SNP, PG_vol_F, PG_vol_F_infert1)
HT_infert1_F_coloc <- coloc_infert_run(HT_vol_cojo_F$SNP, HT_vol_F, HT_vol_F_infert1)
LR_infert1_F_coloc <- coloc_infert_run(LR_vol_cojo_F$SNP, LR_vol_F, LR_vol_F_infert1)
HT_GM_infert1_F_coloc <- coloc_infert_run(HT_GM_vol_cojo_F$SNP, HT_GM_vol_F, HT_GM_vol_F_infert1)

PG_infert2_F_coloc <- coloc_infert_run(PG_vol_cojo_F$SNP, PG_vol_F, PG_vol_F_infert2)
HT_infert2_F_coloc <- coloc_infert_run(HT_vol_cojo_F$SNP, HT_vol_F, HT_vol_F_infert2)
LR_infert2_F_coloc <- coloc_infert_run(LR_vol_cojo_F$SNP, LR_vol_F, LR_vol_F_infert2)
HT_GM_infert2_F_coloc <- coloc_infert_run(HT_GM_vol_cojo_F$SNP, HT_GM_vol_F, HT_GM_vol_F_infert2)

PG_infert3_F_coloc <- coloc_infert_run(PG_vol_cojo_F$SNP, PG_vol_F, PG_vol_F_infert3)
HT_infert3_F_coloc <- coloc_infert_run(HT_vol_cojo_F$SNP, HT_vol_F, HT_vol_F_infert3)
LR_infert3_F_coloc <- coloc_infert_run(LR_vol_cojo_F$SNP, LR_vol_F, LR_vol_F_infert3)
HT_GM_infert3_F_coloc <- coloc_infert_run(HT_GM_vol_cojo_F$SNP, HT_GM_vol_F, HT_GM_vol_F_infert3)

PG_infert4_F_coloc <- coloc_infert_run(PG_vol_cojo_F$SNP, PG_vol_F, PG_vol_F_infert4)
HT_infert4_F_coloc <- coloc_infert_run(HT_vol_cojo_F$SNP, HT_vol_F, HT_vol_F_infert4)
LR_infert4_F_coloc <- coloc_infert_run(LR_vol_cojo_F$SNP, LR_vol_F, LR_vol_F_infert4)
HT_GM_infert4_F_coloc <- coloc_infert_run(HT_GM_vol_cojo_F$SNP, HT_GM_vol_F, HT_GM_vol_F_infert4)

PG_infert5_F_coloc <- coloc_infert_run(PG_vol_cojo_F$SNP, PG_vol_F, PG_vol_F_infert5)
HT_infert5_F_coloc <- coloc_infert_run(HT_vol_cojo_F$SNP, HT_vol_F, HT_vol_F_infert5)
LR_infert5_F_coloc <- coloc_infert_run(LR_vol_cojo_F$SNP, LR_vol_F, LR_vol_F_infert5)
HT_GM_infert5_F_coloc <- coloc_infert_run(HT_GM_vol_cojo_F$SNP, HT_GM_vol_F, HT_GM_vol_F_infert5)

PG_infert_M_coloc <- coloc_infert_run(PG_vol_cojo_M$SNP, PG_vol_M, PG_vol_M_infert)
HT_infert_M_coloc <- coloc_infert_run(HT_vol_cojo_M$SNP, HT_vol_M, HT_vol_M_infert)
LR_infert_M_coloc <- coloc_infert_run(LR_vol_cojo_M$SNP, LR_vol_M, LR_vol_M_infert)
HT_GM_infert_M_coloc <- coloc_infert_run(HT_GM_vol_cojo_M$SNP, HT_GM_vol_M, HT_GM_vol_M_infert)

write.csv(PG_test_coloc, "PG_test_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_test_coloc, "HT_test_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_test_coloc, "LR_test_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_test_coloc, "HT_GM_test_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_fsh_coloc, "PG_fsh_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_fsh_coloc, "HT_fsh_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_fsh_coloc, "LR_fsh_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_fsh_coloc, "HT_GM_fsh_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_lh_coloc, "PG_lh_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_lh_coloc, "HT_lh_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_lh_coloc, "LR_lh_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_lh_coloc, "HT_GM_lh_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_oest_coloc, "PG_oest_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_oest_coloc, "HT_oest_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_oest_coloc, "LR_oest_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_oest_coloc, "HT_GM_oest_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_prog_coloc, "PG_prog_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_prog_coloc, "HT_prog_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_prog_coloc, "LR_prog_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_prog_coloc, "HT_GM_prog_coloc.csv", quote=FALSE, row.names=FALSE)

system("dx upload *_coloc.csv --path data/coloc/")

write.csv(PG_test_F_coloc, "PG_test_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_test_F_coloc, "HT_test_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_test_F_coloc, "LR_test_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_test_F_coloc, "HT_GM_test_F_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_fsh_F_coloc, "PG_fsh_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_fsh_F_coloc, "HT_fsh_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_fsh_F_coloc, "LR_fsh_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_fsh_F_coloc, "HT_GM_fsh_F_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_lh_F_coloc, "PG_lh_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_lh_F_coloc, "HT_lh_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_lh_F_coloc, "LR_lh_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_lh_F_coloc, "HT_GM_lh_F_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_oest_F_coloc, "PG_oest_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_oest_F_coloc, "HT_oest_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_oest_F_coloc, "LR_oest_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_oest_F_coloc, "HT_GM_oest_F_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_prog_F_coloc, "PG_prog_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_prog_F_coloc, "HT_prog_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_prog_F_coloc, "LR_prog_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_prog_F_coloc, "HT_GM_prog_F_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_test_M_coloc, "PG_test_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_test_M_coloc, "HT_test_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_test_M_coloc, "LR_test_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_test_M_coloc, "HT_GM_test_M_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_fsh_M_coloc, "PG_fsh_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_fsh_M_coloc, "HT_fsh_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_fsh_M_coloc, "LR_fsh_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_fsh_M_coloc, "HT_GM_fsh_M_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_lh_M_coloc, "PG_lh_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_lh_M_coloc, "HT_lh_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_lh_M_coloc, "LR_lh_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_lh_M_coloc, "HT_GM_lh_M_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_oest_M_coloc, "PG_oest_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_oest_M_coloc, "HT_oest_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_oest_M_coloc, "LR_oest_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_oest_M_coloc, "HT_GM_oest_M_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_infert1_F_coloc, "PG_infert1_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_infert1_F_coloc, "HT_infert1_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_infert1_F_coloc, "LR_infert1_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_infert1_F_coloc, "HT_GM_infert1_F_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_infert2_F_coloc, "PG_infert2_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_infert2_F_coloc, "HT_infert2_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_infert2_F_coloc, "LR_infert2_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_infert2_F_coloc, "HT_GM_infert2_F_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_infert3_F_coloc, "PG_infert3_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_infert3_F_coloc, "HT_infert3_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_infert3_F_coloc, "LR_infert3_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_infert3_F_coloc, "HT_GM_infert3_F_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_infert4_F_coloc, "PG_infert4_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_infert4_F_coloc, "HT_infert4_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_infert4_F_coloc, "LR_infert4_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_infert4_F_coloc, "HT_GM_infert4_F_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_infert5_F_coloc, "PG_infert5_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_infert5_F_coloc, "HT_infert5_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_infert5_F_coloc, "LR_infert5_F_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_infert5_F_coloc, "HT_GM_infert5_F_coloc.csv", quote=FALSE, row.names=FALSE)

write.csv(PG_infert_M_coloc, "PG_infert_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_infert_M_coloc, "HT_infert_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_infert_M_coloc, "LR_infert_M_coloc.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_infert_M_coloc, "HT_GM_infert_M_coloc.csv", quote=FALSE, row.names=FALSE)


PG_test_coloc$ratio <- PG_test_coloc$PP.H4.abf/PG_test_coloc$PP.H3.abf
PG_test_coloc_sig <- subset(PG_test_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_test_coloc_sig, "PG_test_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_test_coloc$ratio <- HT_test_coloc$PP.H4.abf/HT_test_coloc$PP.H3.abf
HT_test_coloc_sig <- subset(HT_test_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_test_coloc_sig, "HT_test_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_test_coloc$ratio <- LR_test_coloc$PP.H4.abf/LR_test_coloc$PP.H3.abf
LR_test_coloc_sig <- subset(LR_test_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_test_coloc_sig, "LR_test_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_test_coloc$ratio <- HT_GM_test_coloc$PP.H4.abf/HT_GM_test_coloc$PP.H3.abf
HT_GM_test_coloc_sig <- subset(HT_GM_test_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_test_coloc_sig, "HT_GM_test_coloc_sig.csv", quote=FALSE, row.names=FALSE)

PG_fsh_coloc$ratio <- PG_fsh_coloc$PP.H4.abf/PG_fsh_coloc$PP.H3.abf
PG_fsh_coloc_sig <- subset(PG_fsh_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_fsh_coloc_sig, "PG_fsh_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_fsh_coloc$ratio <- HT_fsh_coloc$PP.H4.abf/HT_fsh_coloc$PP.H3.abf
HT_fsh_coloc_sig <- subset(HT_fsh_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_fsh_coloc_sig, "HT_fsh_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_fsh_coloc$ratio <- LR_fsh_coloc$PP.H4.abf/LR_fsh_coloc$PP.H3.abf
LR_fsh_coloc_sig <- subset(LR_fsh_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_fsh_coloc_sig, "LR_fsh_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_fsh_coloc$ratio <- HT_GM_fsh_coloc$PP.H4.abf/HT_GM_fsh_coloc$PP.H3.abf
HT_GM_fsh_coloc_sig <- subset(HT_GM_fsh_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_fsh_coloc_sig, "HT_GM_fsh_coloc_sig.csv", quote=FALSE, row.names=FALSE)

PG_lh_coloc$ratio <- PG_lh_coloc$PP.H4.abf/PG_lh_coloc$PP.H3.abf
PG_lh_coloc_sig <- subset(PG_lh_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_lh_coloc_sig, "PG_lh_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_lh_coloc$ratio <- HT_lh_coloc$PP.H4.abf/HT_lh_coloc$PP.H3.abf
HT_lh_coloc_sig <- subset(HT_lh_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_lh_coloc_sig, "HT_lh_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_lh_coloc$ratio <- LR_lh_coloc$PP.H4.abf/LR_lh_coloc$PP.H3.abf
LR_lh_coloc_sig <- subset(LR_lh_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_lh_coloc_sig, "LR_lh_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_lh_coloc$ratio <- HT_GM_lh_coloc$PP.H4.abf/HT_GM_lh_coloc$PP.H3.abf
HT_GM_lh_coloc_sig <- subset(HT_GM_lh_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_lh_coloc_sig, "HT_GM_lh_coloc_sig.csv", quote=FALSE, row.names=FALSE)

PG_oest_coloc$ratio <- PG_oest_coloc$PP.H4.abf/PG_oest_coloc$PP.H3.abf
PG_oest_coloc_sig <- subset(PG_oest_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_oest_coloc_sig, "PG_oest_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_oest_coloc$ratio <- HT_oest_coloc$PP.H4.abf/HT_oest_coloc$PP.H3.abf
HT_oest_coloc_sig <- subset(HT_oest_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_oest_coloc_sig, "HT_oest_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_oest_coloc$ratio <- LR_oest_coloc$PP.H4.abf/LR_oest_coloc$PP.H3.abf
LR_oest_coloc_sig <- subset(LR_oest_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_oest_coloc_sig, "LR_oest_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_oest_coloc$ratio <- HT_GM_oest_coloc$PP.H4.abf/HT_GM_oest_coloc$PP.H3.abf
HT_GM_oest_coloc_sig <- subset(HT_GM_oest_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_oest_coloc_sig, "HT_GM_oest_coloc_sig.csv", quote=FALSE, row.names=FALSE)

PG_prog_coloc$ratio <- PG_prog_coloc$PP.H4.abf/PG_prog_coloc$PP.H3.abf
PG_prog_coloc_sig <- subset(PG_prog_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_prog_coloc_sig, "PG_prog_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_prog_coloc$ratio <- HT_prog_coloc$PP.H4.abf/HT_prog_coloc$PP.H3.abf
HT_prog_coloc_sig <- subset(HT_prog_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_prog_coloc_sig, "HT_prog_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_prog_coloc$ratio <- LR_prog_coloc$PP.H4.abf/LR_prog_coloc$PP.H3.abf
LR_prog_coloc_sig <- subset(LR_prog_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_prog_coloc_sig, "LR_prog_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_prog_coloc$ratio <- HT_GM_prog_coloc$PP.H4.abf/HT_GM_prog_coloc$PP.H3.abf
HT_GM_prog_coloc_sig <- subset(HT_GM_prog_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_prog_coloc_sig, "HT_GM_prog_coloc_sig.csv", quote=FALSE, row.names=FALSE)

system("dx upload *_sig.csv --path data/coloc/sig/")

PG_test_F_coloc$ratio <- PG_test_F_coloc$PP.H4.abf/PG_test_F_coloc$PP.H3.abf
PG_test_F_coloc_sig <- subset(PG_test_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_test_F_coloc_sig, "PG_test_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_test_F_coloc$ratio <- HT_test_F_coloc$PP.H4.abf/HT_test_F_coloc$PP.H3.abf
HT_test_F_coloc_sig <- subset(HT_test_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_test_F_coloc_sig, "HT_test_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_test_F_coloc$ratio <- LR_test_F_coloc$PP.H4.abf/LR_test_F_coloc$PP.H3.abf
LR_test_F_coloc_sig <- subset(LR_test_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_test_F_coloc_sig, "LR_test_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_test_F_coloc$ratio <- HT_GM_test_F_coloc$PP.H4.abf/HT_GM_test_F_coloc$PP.H3.abf
HT_GM_test_F_coloc_sig <- subset(HT_GM_test_F_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_test_F_coloc_sig, "HT_GM_test_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)

PG_fsh_F_coloc$ratio <- PG_fsh_F_coloc$PP.H4.abf/PG_fsh_F_coloc$PP.H3.abf
PG_fsh_F_coloc_sig <- subset(PG_fsh_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_fsh_F_coloc_sig, "PG_fsh_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_fsh_F_coloc$ratio <- HT_fsh_F_coloc$PP.H4.abf/HT_fsh_F_coloc$PP.H3.abf
HT_fsh_F_coloc_sig <- subset(HT_fsh_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_fsh_F_coloc_sig, "HT_fsh_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_fsh_F_coloc$ratio <- LR_fsh_F_coloc$PP.H4.abf/LR_fsh_F_coloc$PP.H3.abf
LR_fsh_F_coloc_sig <- subset(LR_fsh_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_fsh_F_coloc_sig, "LR_fsh_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_fsh_F_coloc$ratio <- HT_GM_fsh_F_coloc$PP.H4.abf/HT_GM_fsh_F_coloc$PP.H3.abf
HT_GM_fsh_F_coloc_sig <- subset(HT_GM_fsh_F_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_fsh_F_coloc_sig, "HT_GM_fsh_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)

PG_lh_F_coloc$ratio <- PG_lh_F_coloc$PP.H4.abf/PG_lh_F_coloc$PP.H3.abf
PG_lh_F_coloc_sig <- subset(PG_lh_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_lh_F_coloc_sig, "PG_lh_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_lh_F_coloc$ratio <- HT_lh_F_coloc$PP.H4.abf/HT_lh_F_coloc$PP.H3.abf
HT_lh_F_coloc_sig <- subset(HT_lh_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_lh_F_coloc_sig, "HT_lh_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_lh_F_coloc$ratio <- LR_lh_F_coloc$PP.H4.abf/LR_lh_F_coloc$PP.H3.abf
LR_lh_F_coloc_sig <- subset(LR_lh_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_lh_F_coloc_sig, "LR_lh_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_lh_F_coloc$ratio <- HT_GM_lh_F_coloc$PP.H4.abf/HT_GM_lh_F_coloc$PP.H3.abf
HT_GM_lh_F_coloc_sig <- subset(HT_GM_lh_F_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_lh_F_coloc_sig, "HT_GM_lh_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)

PG_oest_F_coloc$ratio <- PG_oest_F_coloc$PP.H4.abf/PG_oest_F_coloc$PP.H3.abf
PG_oest_F_coloc_sig <- subset(PG_oest_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_oest_F_coloc_sig, "PG_oest_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_oest_F_coloc$ratio <- HT_oest_F_coloc$PP.H4.abf/HT_oest_F_coloc$PP.H3.abf
HT_oest_F_coloc_sig <- subset(HT_oest_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_oest_F_coloc_sig, "HT_oest_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_oest_F_coloc$ratio <- LR_oest_F_coloc$PP.H4.abf/LR_oest_F_coloc$PP.H3.abf
LR_oest_F_coloc_sig <- subset(LR_oest_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_oest_F_coloc_sig, "LR_oest_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_oest_F_coloc$ratio <- HT_GM_oest_F_coloc$PP.H4.abf/HT_GM_oest_F_coloc$PP.H3.abf
HT_GM_oest_F_coloc_sig <- subset(HT_GM_oest_F_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_oest_F_coloc_sig, "HT_GM_oest_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)

PG_prog_F_coloc$ratio <- PG_prog_F_coloc$PP.H4.abf/PG_prog_F_coloc$PP.H3.abf
PG_prog_F_coloc_sig <- subset(PG_prog_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_prog_F_coloc_sig, "PG_prog_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_prog_F_coloc$ratio <- HT_prog_F_coloc$PP.H4.abf/HT_prog_F_coloc$PP.H3.abf
HT_prog_F_coloc_sig <- subset(HT_prog_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_prog_F_coloc_sig, "HT_prog_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_prog_F_coloc$ratio <- LR_prog_F_coloc$PP.H4.abf/LR_prog_F_coloc$PP.H3.abf
LR_prog_F_coloc_sig <- subset(LR_prog_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_prog_F_coloc_sig, "LR_prog_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_prog_F_coloc$ratio <- HT_GM_prog_F_coloc$PP.H4.abf/HT_GM_prog_F_coloc$PP.H3.abf
HT_GM_prog_F_coloc_sig <- subset(HT_GM_prog_F_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_prog_F_coloc_sig, "HT_GM_prog_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)

PG_test_M_coloc$ratio <- PG_test_M_coloc$PP.H4.abf/PG_test_M_coloc$PP.H3.abf
PG_test_M_coloc_sig <- subset(PG_test_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_test_M_coloc_sig, "PG_test_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_test_M_coloc$ratio <- HT_test_M_coloc$PP.H4.abf/HT_test_M_coloc$PP.H3.abf
HT_test_M_coloc_sig <- subset(HT_test_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_test_M_coloc_sig, "HT_test_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_test_M_coloc$ratio <- LR_test_M_coloc$PP.H4.abf/LR_test_M_coloc$PP.H3.abf
LR_test_M_coloc_sig <- subset(LR_test_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_test_M_coloc_sig, "LR_test_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_test_M_coloc$ratio <- HT_GM_test_M_coloc$PP.H4.abf/HT_GM_test_M_coloc$PP.H3.abf
HT_GM_test_M_coloc_sig <- subset(HT_GM_test_M_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_test_M_coloc_sig, "HT_GM_test_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)

PG_fsh_M_coloc$ratio <- PG_fsh_M_coloc$PP.H4.abf/PG_fsh_M_coloc$PP.H3.abf
PG_fsh_M_coloc_sig <- subset(PG_fsh_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_fsh_M_coloc_sig, "PG_fsh_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_fsh_M_coloc$ratio <- HT_fsh_M_coloc$PP.H4.abf/HT_fsh_M_coloc$PP.H3.abf
HT_fsh_M_coloc_sig <- subset(HT_fsh_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_fsh_M_coloc_sig, "HT_fsh_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_fsh_M_coloc$ratio <- LR_fsh_M_coloc$PP.H4.abf/LR_fsh_M_coloc$PP.H3.abf
LR_fsh_M_coloc_sig <- subset(LR_fsh_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_fsh_M_coloc_sig, "LR_fsh_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_fsh_M_coloc$ratio <- HT_GM_fsh_M_coloc$PP.H4.abf/HT_GM_fsh_M_coloc$PP.H3.abf
HT_GM_fsh_M_coloc_sig <- subset(HT_GM_fsh_M_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_fsh_M_coloc_sig, "HT_GM_fsh_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)

PG_lh_M_coloc$ratio <- PG_lh_M_coloc$PP.H4.abf/PG_lh_M_coloc$PP.H3.abf
PG_lh_M_coloc_sig <- subset(PG_lh_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_lh_M_coloc_sig, "PG_lh_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_lh_M_coloc$ratio <- HT_lh_M_coloc$PP.H4.abf/HT_lh_M_coloc$PP.H3.abf
HT_lh_M_coloc_sig <- subset(HT_lh_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_lh_M_coloc_sig, "HT_lh_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_lh_M_coloc$ratio <- LR_lh_M_coloc$PP.H4.abf/LR_lh_M_coloc$PP.H3.abf
LR_lh_M_coloc_sig <- subset(LR_lh_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_lh_M_coloc_sig, "LR_lh_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_lh_M_coloc$ratio <- HT_GM_lh_M_coloc$PP.H4.abf/HT_GM_lh_M_coloc$PP.H3.abf
HT_GM_lh_M_coloc_sig <- subset(HT_GM_lh_M_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_lh_M_coloc_sig, "HT_GM_lh_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)

PG_oest_M_coloc$ratio <- PG_oest_M_coloc$PP.H4.abf/PG_oest_M_coloc$PP.H3.abf
PG_oest_M_coloc_sig <- subset(PG_oest_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_oest_M_coloc_sig, "PG_oest_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_oest_M_coloc$ratio <- HT_oest_M_coloc$PP.H4.abf/HT_oest_M_coloc$PP.H3.abf
HT_oest_M_coloc_sig <- subset(HT_oest_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_oest_M_coloc_sig, "HT_oest_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_oest_M_coloc$ratio <- LR_oest_M_coloc$PP.H4.abf/LR_oest_M_coloc$PP.H3.abf
LR_oest_M_coloc_sig <- subset(LR_oest_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_oest_M_coloc_sig, "LR_oest_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_oest_M_coloc$ratio <- HT_GM_oest_M_coloc$PP.H4.abf/HT_GM_oest_M_coloc$PP.H3.abf
HT_GM_oest_M_coloc_sig <- subset(HT_GM_oest_M_coloc, ratio >5 & PP.H4.abf)
write.csv(HT_GM_oest_M_coloc_sig, "HT_GM_oest_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)

PG_infert1_F_coloc$ratio <- PG_infert1_F_coloc$PP.H4.abf/PG_infert1_F_coloc$PP.H3.abf
PG_infert1_F_coloc_sig <- subset(PG_infert1_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_infert1_F_coloc_sig, "PG_infert1_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_infert1_F_coloc$ratio <- HT_infert1_F_coloc$PP.H4.abf/HT_infert1_F_coloc$PP.H3.abf
HT_infert1_F_coloc_sig <- subset(HT_infert1_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_infert1_F_coloc_sig, "HT_infert1_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_infert1_F_coloc$ratio <- LR_infert1_F_coloc$PP.H4.abf/LR_infert1_F_coloc$PP.H3.abf
LR_infert1_F_coloc_sig <- subset(LR_infert1_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_infert1_F_coloc_sig, "LR_infert1_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_infert1_F_coloc$ratio <- HT_GM_infert1_F_coloc$PP.H4.abf/HT_GM_infert1_F_coloc$PP.H3.abf
HT_GM_infert1_F_coloc_sig <- subset(HT_GM_infert1_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_GM_infert1_F_coloc_sig, "HT_GM_infert1_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
PG_infert2_F_coloc$ratio <- PG_infert2_F_coloc$PP.H4.abf/PG_infert2_F_coloc$PP.H3.abf
PG_infert2_F_coloc_sig <- subset(PG_infert2_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_infert2_F_coloc_sig, "PG_infert2_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_infert2_F_coloc$ratio <- HT_infert2_F_coloc$PP.H4.abf/HT_infert2_F_coloc$PP.H3.abf
HT_infert2_F_coloc_sig <- subset(HT_infert2_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_infert2_F_coloc_sig, "HT_infert2_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_infert2_F_coloc$ratio <- LR_infert2_F_coloc$PP.H4.abf/LR_infert2_F_coloc$PP.H3.abf
LR_infert2_F_coloc_sig <- subset(LR_infert2_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_infert2_F_coloc_sig, "LR_infert2_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_infert2_F_coloc$ratio <- HT_GM_infert2_F_coloc$PP.H4.abf/HT_GM_infert2_F_coloc$PP.H3.abf
HT_GM_infert2_F_coloc_sig <- subset(HT_GM_infert2_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_GM_infert2_F_coloc_sig, "HT_GM_infert2_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
PG_infert3_F_coloc$ratio <- PG_infert3_F_coloc$PP.H4.abf/PG_infert3_F_coloc$PP.H3.abf
PG_infert3_F_coloc_sig <- subset(PG_infert3_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_infert3_F_coloc_sig, "PG_infert3_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_infert3_F_coloc$ratio <- HT_infert3_F_coloc$PP.H4.abf/HT_infert3_F_coloc$PP.H3.abf
HT_infert3_F_coloc_sig <- subset(HT_infert3_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_infert3_F_coloc_sig, "HT_infert3_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_infert3_F_coloc$ratio <- LR_infert3_F_coloc$PP.H4.abf/LR_infert3_F_coloc$PP.H3.abf
LR_infert3_F_coloc_sig <- subset(LR_infert3_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_infert3_F_coloc_sig, "LR_infert3_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_infert3_F_coloc$ratio <- HT_GM_infert3_F_coloc$PP.H4.abf/HT_GM_infert3_F_coloc$PP.H3.abf
HT_GM_infert3_F_coloc_sig <- subset(HT_GM_infert3_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_GM_infert3_F_coloc_sig, "HT_GM_infert3_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
PG_infert4_F_coloc$ratio <- PG_infert4_F_coloc$PP.H4.abf/PG_infert4_F_coloc$PP.H3.abf
PG_infert4_F_coloc_sig <- subset(PG_infert4_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_infert4_F_coloc_sig, "PG_infert4_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_infert4_F_coloc$ratio <- HT_infert4_F_coloc$PP.H4.abf/HT_infert4_F_coloc$PP.H3.abf
HT_infert4_F_coloc_sig <- subset(HT_infert4_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_infert4_F_coloc_sig, "HT_infert4_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_infert4_F_coloc$ratio <- LR_infert4_F_coloc$PP.H4.abf/LR_infert4_F_coloc$PP.H3.abf
LR_infert4_F_coloc_sig <- subset(LR_infert4_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_infert4_F_coloc_sig, "LR_infert4_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_infert4_F_coloc$ratio <- HT_GM_infert4_F_coloc$PP.H4.abf/HT_GM_infert4_F_coloc$PP.H3.abf
HT_GM_infert4_F_coloc_sig <- subset(HT_GM_infert4_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_GM_infert4_F_coloc_sig, "HT_GM_infert4_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
PG_infert5_F_coloc$ratio <- PG_infert5_F_coloc$PP.H4.abf/PG_infert5_F_coloc$PP.H3.abf
PG_infert5_F_coloc_sig <- subset(PG_infert5_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_infert5_F_coloc_sig, "PG_infert5_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_infert5_F_coloc$ratio <- HT_infert5_F_coloc$PP.H4.abf/HT_infert5_F_coloc$PP.H3.abf
HT_infert5_F_coloc_sig <- subset(HT_infert5_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_infert5_F_coloc_sig, "HT_infert5_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_infert5_F_coloc$ratio <- LR_infert5_F_coloc$PP.H4.abf/LR_infert5_F_coloc$PP.H3.abf
LR_infert5_F_coloc_sig <- subset(LR_infert5_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_infert5_F_coloc_sig, "LR_infert5_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_infert5_F_coloc$ratio <- HT_GM_infert5_F_coloc$PP.H4.abf/HT_GM_infert5_F_coloc$PP.H3.abf
HT_GM_infert5_F_coloc_sig <- subset(HT_GM_infert5_F_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_GM_infert5_F_coloc_sig, "HT_GM_infert5_F_coloc_sig.csv", quote=FALSE, row.names=FALSE)


PG_infert_M_coloc$ratio <- PG_infert_M_coloc$PP.H4.abf/PG_infert_M_coloc$PP.H3.abf
PG_infert_M_coloc_sig <- subset(PG_infert_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(PG_infert_M_coloc_sig, "PG_infert_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_infert_M_coloc$ratio <- HT_infert_M_coloc$PP.H4.abf/HT_infert_M_coloc$PP.H3.abf
HT_infert_M_coloc_sig <- subset(HT_infert_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_infert_M_coloc_sig, "HT_infert_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
LR_infert_M_coloc$ratio <- LR_infert_M_coloc$PP.H4.abf/LR_infert_M_coloc$PP.H3.abf
LR_infert_M_coloc_sig <- subset(LR_infert_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(LR_infert_M_coloc_sig, "LR_infert_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)
HT_GM_infert_M_coloc$ratio <- HT_GM_infert_M_coloc$PP.H4.abf/HT_GM_infert_M_coloc$PP.H3.abf
HT_GM_infert_M_coloc_sig <- subset(HT_GM_infert_M_coloc, ratio >5 & PP.H4.abf >0.5)
write.csv(HT_GM_infert_M_coloc_sig, "HT_GM_infert_M_coloc_sig.csv", quote=FALSE, row.names=FALSE)

