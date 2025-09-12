###############
###LIBRARIES###
###############

install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages("R.utils")
install.packages("metafor")
library(TwoSampleMR)
library(data.table)
library(R.utils)
library(metafor)


##########
###DATA###
##########

dir.create("all")
dir.create("female_only")
dir.create("male_only")

system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o all/")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt -o all/")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o all/")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX.txt -o all/")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o all/")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt -o all/")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt -o all/")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt -o all/")

system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o female_only/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt -o female_only/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o female_only")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt -o female_only/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o female_only/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.txt -o female_only/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt -o female_only/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt -o female_only/")

system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o male_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt -o male_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o male_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt -o male_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o male_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.txt -o male_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt -o male_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt -o male_only/")

PG_vol <- fread("all/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
PG_vol_cojo <- fread("all/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
HT_vol <- fread("all/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_vol_cojo <- fread("all/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX.txt", data.table=FALSE)
LR_vol <- fread("all/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
LR_vol_cojo <- fread("all/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
HT_GM_vol <- fread("all/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM_vol_cojo <- fread("all/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)

PG_vol_F <- fread("female_only/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
PG_vol_cojo_F <- fread("female_only/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
HT_vol_F <- fread("female_only/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_vol_cojo_F <- fread("female_only/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
LR_vol_F <- fread("female_only/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
LR_vol_cojo_F <- fread("female_only/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)
HT_GM_vol_F <- fread("female_only/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM_vol_cojo_F <- fread("female_only/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)

PG_vol_M <- fread("male_only/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
PG_vol_cojo_M <- fread("male_only/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
HT_vol_M <- fread("male_only/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_vol_cojo_M <- fread("male_only/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
LR_vol_M <- fread("male_only/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
LR_vol_cojo_M <- fread("male_only/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)
HT_GM_vol_M <- fread("male_only/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM_vol_cojo_M <- fread("male_only/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)

system("dx download data/hormone_fertility_update_231125/0624_remeta/Testosterone_sex_comb_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/FSH_sex_comb_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/LH_sex_comb_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/Progesterone_sex_comb_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/Oestradiol_sex_comb_EUR_with_rsids.txt")

system("dx download data/hormone_fertility_update_231125/0624_remeta/lead_snps/all_lead_snp_sumstats_Testosterone_sex_comb_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/lead_snps/all_lead_snp_sumstats_FSH_sex_comb_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/lead_snps/all_lead_snp_sumstats_LH_sex_comb_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/lead_snps/all_lead_snp_sumstats_Oestradiol_sex_comb_EUR_with_rsids.txt")

system("dx download data/hormone_fertility_update_231125/0624_remeta/Testosterone_F_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/FSH_F_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/LH_F_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/Progesterone_F_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/Oestradiol_F_EUR_with_rsids.txt")

system("dx download data/hormone_fertility_update_231125/0624_remeta/lead_snps/all_lead_snp_sumstats_Testosterone_F_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/lead_snps/all_lead_snp_sumstats_FSH_F_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/lead_snps/all_lead_snp_sumstats_LH_F_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/lead_snps/all_lead_snp_sumstats_Oestradiol_F_EUR_with_rsids.txt")

system("dx download data/hormone_fertility_update_231125/0624_remeta/Testosterone_M_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/FSH_M_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/LH_M_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/Oestradiol_M_EUR_with_rsids.txt")

system("dx download data/hormone_fertility_update_231125/0624_remeta/lead_snps/all_lead_snp_sumstats_Testosterone_M_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/lead_snps/all_lead_snp_sumstats_FSH_M_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/lead_snps/all_lead_snp_sumstats_LH_M_EUR_with_rsids.txt")
system("dx download data/hormone_fertility_update_231125/0624_remeta/lead_snps/all_lead_snp_sumstats_Oestradiol_M_EUR_with_rsids.txt")

test_all_eur <- fread("Testosterone_sex_comb_EUR_with_rsids.txt", data.table=FALSE)
fsh_all_eur <- fread("FSH_sex_comb_EUR_with_rsids.txt", data.table=FALSE)
lh_all_eur <- fread("LH_sex_comb_EUR_with_rsids.txt", data.table=FALSE)
prog_all_eur <- fread("Progesterone_sex_comb_EUR_with_rsids.txt", data.table=FALSE)
oest_all_eur <- fread("Oestradiol_sex_comb_EUR_with_rsids.txt", data.table=FALSE)

test_all_eur_lead <- fread("all_lead_snp_sumstats_Testosterone_sex_comb_EUR_with_rsids.txt", data.table=FALSE)
fsh_all_eur_lead <- fread("all_lead_snp_sumstats_FSH_sex_comb_EUR_with_rsids.txt", data.table=FALSE)
lh_all_eur_lead <- fread("all_lead_snp_sumstats_LH_sex_comb_EUR_with_rsids.txt", data.table=FALSE)
#NO prog lead snps
oest_all_eur_lead <- fread("all_lead_snp_sumstats_Oestradiol_sex_comb_EUR_with_rsids.txt", data.table=FALSE)

test_female_eur <- fread("Testosterone_F_EUR_with_rsids.txt", data.table=FALSE)
fsh_female_eur <- fread("FSH_F_EUR_with_rsids.txt", data.table=FALSE)
lh_female_eur <- fread("LH_F_EUR_with_rsids.txt", data.table=FALSE)
prog_female_eur <- fread("Progesterone_F_EUR_with_rsids.txt", data.table=FALSE)
oest_female_eur <- fread("Oestradiol_F_EUR_with_rsids.txt", data.table=FALSE)

test_female_eur_lead <- fread("all_lead_snp_sumstats_Testosterone_F_EUR_with_rsids.txt", data.table=FALSE)
fsh_female_eur_lead <- fread("all_lead_snp_sumstats_FSH_F_EUR_with_rsids.txt", data.table=FALSE)
lh_female_eur_lead <- fread("all_lead_snp_sumstats_LH_F_EUR_with_rsids.txt", data.table=FALSE)
#NO prog lead snps
oest_female_eur_lead <- fread("all_lead_snp_sumstats_Oestradiol_F_EUR_with_rsids.txt", data.table=FALSE)

test_male_eur <- fread("Testosterone_M_EUR_with_rsids.txt", data.table=FALSE)
fsh_male_eur <- fread("FSH_M_EUR_with_rsids.txt", data.table=FALSE)
lh_male_eur <- fread("LH_M_EUR_with_rsids.txt", data.table=FALSE)
oest_male_eur <- fread("Oestradiol_M_EUR_with_rsids.txt", data.table=FALSE)

test_male_eur_lead <- fread("all_lead_snp_sumstats_Testosterone_M_EUR_with_rsids.txt", data.table=FALSE)
fsh_male_eur_lead <- fread("all_lead_snp_sumstats_FSH_M_EUR_with_rsids.txt", data.table=FALSE)
lh_male_eur_lead <- fread("all_lead_snp_sumstats_LH_M_EUR_with_rsids.txt", data.table=FALSE)
oest_male_eur_lead <- fread("all_lead_snp_sumstats_Oestradiol_M_EUR_with_rsids.txt", data.table=FALSE)

#female_infert1 <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/female_infertility_analysis1_eur_with_rsids.txt.gz"), header=TRUE)
#female_infert2 <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/female_infertility_analysis2_eur_with_rsids.txt.gz"), header=TRUE)
#female_infert3 <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/female_infertility_analysis3_eur_with_rsids.txt.gz"), header=TRUE)
#female_infert4 <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/female_infertility_analysis4_eur_with_rsids.txt.gz"), header=TRUE)
#female_infert5 <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/female_infertility_analysis5_eur_with_rsids.txt.gz"), header=TRUE)
#male_infert <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/full_sumstats_txt/male_infertility_eur_with_rsids.txt.gz"), header=TRUE)

#female_infert1_lead <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/lead_snps/all_lead_snps_female_infertility_analysis1_eur.txt.gz"), header=TRUE)
#female_infert2_lead <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/lead_snps/all_lead_snps_female_infertility_analysis2_eur.txt.gz"), header=TRUE)
#female_infert3_lead <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/lead_snps/all_lead_snps_female_infertility_analysis3_eur.txt.gz"), header=TRUE)
#female_infert4_lead <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/lead_snps/all_lead_snps_female_infertility_analysis4_eur.txt.gz"), header=TRUE)
#female_infert5_lead <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/lead_snps/all_lead_snps_female_infertility_analysis5_eur.txt.gz"), header=TRUE)
#male_infert_lead <- read.table(gzfile("/mnt/project/data/hormone_fertility_update_231125/lead_snps/all_lead_snps_male_infertility_eur.txt.gz"), header=TRUE)

##############
###ANALYSIS###
##############

#get data for lead snps pf brain structures
PG_sig <- subset(PG_vol, ID %in% PG_vol_cojo$SNP)
HT_sig <- subset(HT_vol, ID %in% HT_vol_cojo$SNP)
LR_sig <- subset(LR_vol, ID %in% LR_vol_cojo$SNP)
HT_GM_sig <- subset(HT_GM_vol, ID %in% HT_GM_vol_cojo$SNP)

PG_sig_F <- subset(PG_vol_F, ID %in% PG_vol_cojo_F$SNP)
HT_sig_F <- subset(HT_vol_F, ID %in% HT_vol_cojo_F$SNP)
LR_sig_F <- subset(LR_vol_F, ID %in% LR_vol_cojo_F$SNP)
HT_GM_sig_F <- subset(HT_GM_vol_F, ID %in% HT_GM_vol_cojo_F$SNP)

PG_sig_M <- subset(PG_vol_M, ID %in% PG_vol_cojo_M$SNP)
HT_sig_M <- subset(HT_vol_M, ID %in% HT_vol_cojo_M$SNP)
LR_sig_M <- subset(LR_vol_M, ID %in% LR_vol_cojo_M$SNP)
HT_GM_sig_M <- subset(HT_GM_vol_M, ID %in% HT_GM_vol_cojo_M$SNP)

#set up exposure data for brain structures
PG_exp <- format_data(
  PG_sig,
  type="exposure",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS",
)
PG_exp$exposure <- "PG"

HT_exp <- format_data(
  HT_sig,
  type="exposure",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_exp$exposure <- "HT"

LR_exp <- format_data(
  LR_sig,
  type="exposure",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_exp$exposure <- "LR"

HT_GM_exp <- format_data(
  HT_GM_sig,
  type="exposure",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_exp$exposure <- "HT_GM"

PG_F_exp <- format_data(
  PG_sig_F,
  type="exposure",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS",
)
PG_F_exp$exposure <- "PG_F"

HT_F_exp <- format_data(
  HT_sig_F,
  type="exposure",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_F_exp$exposure <- "HT_F"

LR_F_exp <- format_data(
  LR_sig_F,
  type="exposure",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_F_exp$exposure <- "LR_F"

HT_GM_F_exp <- format_data(
  HT_GM_sig_F,
  type="exposure",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_F_exp$exposure <- "HT_GM_F"

PG_M_exp <- format_data(
  PG_sig_M,
  type="exposure",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS",
)
PG_M_exp$exposure <- "PG_M"

HT_M_exp <- format_data(
  HT_sig_M,
  type="exposure",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_M_exp$exposure <- "HT_M"

LR_M_exp <- format_data(
  LR_sig_M,
  type="exposure",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_M_exp$exposure <- "LR_M"

HT_GM_M_exp <- format_data(
  HT_GM_sig_M,
  type="exposure",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_M_exp$exposure <- "HT_GM_M"

#Add N manually from Samvida table to all and lead hormone results
test_all_eur$N <- 69666
fsh_all_eur$N <- 32769
lh_all_eur$N <- 28015
prog_all_eur$N <- 16171
oest_all_eur$N <- 60249

test_all_eur_lead$N <- 69666
fsh_all_eur_lead$N <- 32769
lh_all_eur_lead$N <- 28015
#prog_all_eur_lead$N <- 16171
oest_all_eur_lead$N <- 60249

test_female_eur$N <- 38541
fsh_female_eur$N <- 30596
lh_female_eur$N <- 25259
oest_female_eur$N <- 38360
prog_female_eur$N <- 14558

test_female_eur_lead$N <- 38541
fsh_female_eur_lead$N <- 30596
lh_female_eur_lead$N <- 25259
oest_female_eur_lead$N <- 38360

test_male_eur$N <- 37756
fsh_male_eur$N <- 3981
lh_male_eur$N <- 4532
oest_male_eur$N <- 21674

test_male_eur_lead$N <- 37756
fsh_male_eur_lead$N <- 3981
lh_male_eur_lead$N <- 4532
oest_male_eur_lead$N <- 21674

#Format hormone exposire data
test_exp <- format_data(
  test_all_eur_lead,
  type="exposure",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_exp$exposure <- "test"

fsh_exp <- format_data(
  fsh_all_eur_lead,
  type="exposure",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_exp$exposure <- "fsh"

lh_exp <- format_data(
  lh_all_eur_lead,
  type="exposure",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
lh_exp$exposure <- "lh"

#PROG NO SIGNIFICANT SNPS

oest_exp <- format_data(
  oest_all_eur_lead,
  type="exposure",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
oest_exp$exposure <- "oest"

test_F_exp <- format_data(
  test_female_eur_lead,
  type="exposure",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_F_exp$exposure <- "test_F"

fsh_F_exp <- format_data(
  fsh_female_eur_lead,
  type="exposure",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_F_exp$exposure <- "fsh_F"

lh_F_exp <- format_data(
  lh_female_eur_lead,
  type="exposure",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
lh_F_exp$exposure <- "lh_F"

#PROG NO SIGNIFICANT SNPS

oest_F_exp <- format_data(
  oest_female_eur_lead,
  type="exposure",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
oest_F_exp$exposure <- "oest_F"

test_M_exp <- format_data(
  test_male_eur_lead,
  type="exposure",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_M_exp$exposure <- "test_M"

fsh_M_exp <- format_data(
  fsh_male_eur_lead,
  type="exposure",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_M_exp$exposure <- "fsh_M"

lh_M_exp <- format_data(
  lh_male_eur_lead,
  type="exposure",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
lh_M_exp$exposure <- "lh_M"

#PROG NO SIGNIFICANT SNPS

oest_M_exp <- format_data(
  oest_male_eur_lead,
  type="exposure",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
oest_M_exp$exposure <- "oest_M"


#female_infert1_exp <- format_data(
#  female_infert1_lead,
#  type="exposure",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col = "Allele1",
#  other_allele_col = "Allele2",
#  eaf_col = "Freq1",
#  pval_col = "PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#female_infert1_exp$exposure <- "female_infert1"

#female_infert2_exp <- format_data(
#  female_infert2_lead,
#  type="exposure",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col = "Allele1",
#  other_allele_col = "Allele2",
#  eaf_col = "Freq1",
#  pval_col = "PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#female_infert2_exp$exposure <- "female_infert2"

#female_infert3_exp <- format_data(
#  female_infert3_lead,
#  type="exposure",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col = "Allele1",
#  other_allele_col = "Allele2",
#  eaf_col = "Freq1",
#  pval_col = "PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#female_infert3_exp$exposure <- "female_infert3"

#female_infert4_exp <- format_data(
#  female_infert4_lead,
#  type="exposure",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col = "Allele1",
#  other_allele_col = "Allele2",
#  eaf_col = "Freq1",
#  pval_col = "PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#female_infert4_exp$exposure <- "female_infert4"

#female_infert5_exp <- format_data(
#  female_infert5_lead,
#  type="exposure",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col = "Allele1",
#  other_allele_col = "Allele2",
#  eaf_col = "Freq1",
#  pval_col = "PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#female_infert5_exp$exposure <- "female_infert5"

#male_infert_exp <- format_data(
#  male_infert_lead,
#  type="exposure",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col = "Allele1",
#  other_allele_col = "Allele2",
#  eaf_col = "Freq1",
#  pval_col = "PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#male_infert_exp$exposure <- "male_infert"

##################################
###RUN FORMATTING DATA TOGETHER###
##################################

#brain as exposure and hormone as outcome
PG_test_out <- format_data(
  snps=PG_exp$SNP,
  test_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
  )
PG_test_out$outcome <- "test"
PG_fsh_out <- format_data(
  snps=PG_exp$SNP,
  fsh_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
PG_fsh_out$outcome <- "fsh"
PG_lh_out <- format_data(
  snps=PG_exp$SNP,
  lh_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
PG_lh_out$outcome <- "lh"
PG_prog_out <- format_data(
  snps=PG_exp$SNP,
  prog_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
PG_prog_out$outcome <- "prog"
PG_oest_out <- format_data(
  snps=PG_exp$SNP,
  oest_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
PG_oest_out$outcome <- "oest"

HT_test_out <- format_data(
  snps=HT_exp$SNP,
  test_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_test_out$outcome <- "test"
HT_fsh_out <- format_data(
  snps=HT_exp$SNP,
  fsh_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_fsh_out$outcome <- "fsh"
HT_lh_out <- format_data(
  snps=HT_exp$SNP,
  lh_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_lh_out$outcome <- "lh"
HT_prog_out <- format_data(
  snps=HT_exp$SNP,
  prog_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_prog_out$outcome <- "prog"
HT_oest_out <- format_data(
  snps=HT_exp$SNP,
  oest_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_oest_out$outcome <- "oest"

LR_test_out <- format_data(
  snps=LR_exp$SNP,
  test_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_test_out$outcome <- "test"
LR_fsh_out <- format_data(
  snps=LR_exp$SNP,
  fsh_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_fsh_out$outcome <- "fsh"
LR_lh_out <- format_data(
  snps=LR_exp$SNP,
  lh_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_lh_out$outcome <- "lh"
LR_prog_out <- format_data(
  snps=LR_exp$SNP,
  prog_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_prog_out$outcome <- "prog"
LR_oest_out <- format_data(
  snps=LR_exp$SNP,
  oest_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_oest_out$outcome <- "oest"

HT_GM_test_out <- format_data(
  snps=HT_GM_exp$SNP,
  test_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_test_out$outcome <- "test"
HT_GM_fsh_out <- format_data(
  snps=HT_GM_exp$SNP,
  fsh_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_fsh_out$outcome <- "fsh"
HT_GM_lh_out <- format_data(
  snps=HT_GM_exp$SNP,
  lh_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_lh_out$outcome <- "lh"
HT_GM_prog_out <- format_data(
  snps=HT_GM_exp$SNP,
  prog_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_prog_out$outcome <- "prog"
HT_GM_oest_out <- format_data(
  snps=HT_GM_exp$SNP,
  oest_all_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_oest_out$outcome <- "oest"

#hormone as exposure, brain as outcome
test_HT_out <- format_data(
  snps=test_exp$SNP,
  HT_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_HT_out$outcome <- "HT"
test_PG_out <- format_data(
  snps=test_exp$SNP,
  PG_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_PG_out$outcome <- "PG"
test_LR_out <- format_data(
  snps=test_exp$SNP,
  LR_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_LR_out$outcome <- "LR"
test_HT_GM_out <- format_data(
  snps=test_exp$SNP,
  HT_GM_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_HT_GM_out$outcome <- "HT_GM"

fsh_HT_out <- format_data(
  snps=fsh_exp$SNP,
  HT_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_HT_out$outcome <- "HT"
fsh_PG_out <- format_data(
  snps=fsh_exp$SNP,
  PG_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_PG_out$outcome <- "PG"
fsh_LR_out <- format_data(
  snps=fsh_exp$SNP,
  LR_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_LR_out$outcome <- "LR"
fsh_HT_GM_out <- format_data(
  snps=fsh_exp$SNP,
  HT_GM_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_HT_GM_out$outcome <- "HT_GM"

lh_HT_out <- format_data(
  snps=lh_exp$SNP,
  HT_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
lh_HT_out$outcome <- "HT"
lh_PG_out <- format_data(
  snps=lh_exp$SNP,
  PG_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
lh_PG_out$outcome <- "PG"
lh_LR_out <- format_data(
  snps=lh_exp$SNP,
  LR_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
lh_LR_out$outcome <- "LR"
lh_HT_GM_out <- format_data(
  snps=lh_exp$SNP,
  HT_GM_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
lh_HT_GM_out$outcome <- "HT_GM"

oest_HT_out <- format_data(
  snps=oest_exp$SNP,
  HT_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
oest_HT_out$outcome <- "HT"
oest_PG_out <- format_data(
  snps=oest_exp$SNP,
  PG_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
oest_PG_out$outcome <- "PG"
oest_LR_out <- format_data(
  snps=oest_exp$SNP,
  LR_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
oest_LR_out$outcome <- "LR"
oest_HT_GM_out <- format_data(
  snps=oest_exp$SNP,
  HT_GM_vol,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
oest_HT_GM_out$outcome <- "HT_GM"

###FEMALE

#brain as exposure and hormone as outcome
PG_test_F_out <- format_data(
  snps=PG_F_exp$SNP,
  test_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
PG_test_F_out$outcome <- "test_F"
PG_fsh_F_out <- format_data(
  snps=PG_F_exp$SNP,
  fsh_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
PG_fsh_F_out$outcome <- "fsh_F"
PG_lh_F_out <- format_data(
  snps=PG_F_exp$SNP,
  lh_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
PG_lh_F_out$outcome <- "lh_F"
PG_prog_F_out <- format_data(
  snps=PG_F_exp$SNP,
  prog_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
PG_prog_F_out$outcome <- "prog_F"
PG_oest_F_out <- format_data(
  snps=PG_F_exp$SNP,
  oest_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
PG_oest_F_out$outcome <- "oest_F"

HT_test_F_out <- format_data(
  snps=HT_F_exp$SNP,
  test_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_test_F_out$outcome <- "test_F"
HT_fsh_F_out <- format_data(
  snps=HT_F_exp$SNP,
  fsh_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_fsh_F_out$outcome <- "fsh_F"
HT_lh_F_out <- format_data(
  snps=HT_F_exp$SNP,
  lh_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_lh_F_out$outcome <- "lh_F"
HT_prog_F_out <- format_data(
  snps=HT_F_exp$SNP,
  prog_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_prog_F_out$outcome <- "prog_F"
HT_oest_F_out <- format_data(
  snps=HT_F_exp$SNP,
  oest_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_oest_F_out$outcome <- "oest_F"

LR_test_F_out <- format_data(
  snps=LR_F_exp$SNP,
  test_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_test_F_out$outcome <- "test_F"
LR_fsh_F_out <- format_data(
  snps=LR_F_exp$SNP,
  fsh_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_fsh_F_out$outcome <- "fsh_F"
LR_lh_F_out <- format_data(
  snps=LR_F_exp$SNP,
  lh_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_lh_F_out$outcome <- "lh_F"
LR_prog_F_out <- format_data(
  snps=LR_F_exp$SNP,
  prog_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_prog_F_out$outcome <- "prog_F"
LR_oest_F_out <- format_data(
  snps=LR_F_exp$SNP,
  oest_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_oest_F_out$outcome <- "oest_F"

HT_GM_test_F_out <- format_data(
  snps=HT_GM_F_exp$SNP,
  test_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_test_F_out$outcome <- "test_F"
HT_GM_fsh_F_out <- format_data(
  snps=HT_GM_F_exp$SNP,
  fsh_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_fsh_F_out$outcome <- "fsh_F"
HT_GM_lh_F_out <- format_data(
  snps=HT_GM_F_exp$SNP,
  lh_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_lh_F_out$outcome <- "lh_F"
HT_GM_prog_F_out <- format_data(
  snps=HT_GM_F_exp$SNP,
  prog_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_prog_F_out$outcome <- "prog_F"
HT_GM_oest_F_out <- format_data(
  snps=HT_GM_F_exp$SNP,
  oest_female_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_oest_F_out$outcome <- "oest_F"

#hormone as exposure, brain as outcome
test_HT_F_out <- format_data(
  snps=test_F_exp$SNP,
  HT_vol_F,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_HT_F_out$outcome <- "HT_F"
test_PG_F_out <- format_data(
  snps=test_F_exp$SNP,
  PG_vol_F,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_PG_F_out$outcome <- "PG_F"
test_LR_F_out <- format_data(
  snps=test_F_exp$SNP,
  LR_vol_F,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_LR_F_out$outcome <- "LR_F"
test_HT_GM_F_out <- format_data(
  snps=test_F_exp$SNP,
  HT_GM_vol_F,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_HT_GM_F_out$outcome <- "HT_GM_F"

fsh_HT_F_out <- format_data(
  snps=fsh_F_exp$SNP,
  HT_vol_F,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_HT_F_out$outcome <- "HT_F"
fsh_PG_F_out <- format_data(
  snps=fsh_F_exp$SNP,
  PG_vol_F,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_PG_F_out$outcome <- "PG_F"
fsh_LR_F_out <- format_data(
  snps=fsh_F_exp$SNP,
  LR_vol_F,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_LR_F_out$outcome <- "LR_F"
fsh_HT_GM_F_out <- format_data(
  snps=fsh_F_exp$SNP,
  HT_GM_vol_F,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_HT_GM_F_out$outcome <- "HT_GM_F"

lh_HT_F_out <- format_data(
  snps=lh_F_exp$SNP,
  HT_vol_F,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
lh_HT_F_out$outcome <- "HT_F"
lh_PG_F_out <- format_data(
  snps=lh_F_exp$SNP,
  PG_vol_F,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
lh_PG_F_out$outcome <- "PG_F"
lh_LR_F_out <- format_data(
  snps=lh_F_exp$SNP,
  LR_vol_F,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
lh_LR_F_out$outcome <- "LR_F"
lh_HT_GM_F_out <- format_data(
  snps=lh_F_exp$SNP,
  HT_GM_vol_F,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
lh_HT_GM_F_out$outcome <- "HT_GM_F"

#DOES NOT WORK -  I think because of lack of crossover of genetics
#oest_HT_F_out <- format_data(
#  snps=oest_F_exp$SNP,
#  HT_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  chr_col="CHROM",
#  log_pval=TRUE,
#  pos_col="GENPOS"
#)
#oest_HT_F_out$outcome <- "HT_F"
#DOES NOT WORK -  I think because of lack of crossover of genetics
#oest_PG_F_out <- format_data(
#  snps=oest_F_exp$SNP,
#  PG_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
#)
#oest_PG_F_out$outcome <- "PG_F"
#DOES NOT WORK -  I think because of lack of crossover of genetics
#oest_LR_F_out <- format_data(
#  snps=oest_F_exp$SNP,
#  LR_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
#)
#oest_LR_F_out$outcome <- "LR_female"
#DOES NOT WORK -  I think because of lack of crossover of genetics
#oest_HT_GM_F_out <- format_data(
#  snps=oest_F_exp$SNP,
#  HT_GM_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
#)
#oest_HT_GM_F_out$outcome <- "HT_GM_female"

###MALE

#brain as exposure and hormone as outcome
PG_test_M_out <- format_data(
  snps=PG_M_exp$SNP,
  test_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
PG_test_M_out$outcome <- "test_M"
PG_fsh_M_out <- format_data(
  snps=PG_M_exp$SNP,
  fsh_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
PG_fsh_M_out$outcome <- "fsh_M"
PG_lh_M_out <- format_data(
  snps=PG_M_exp$SNP,
  lh_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
PG_lh_M_out$outcome <- "lh_M"

PG_oest_M_out <- format_data(
  snps=PG_M_exp$SNP,
  oest_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
PG_oest_M_out$outcome <- "oest_M"

HT_test_M_out <- format_data(
  snps=HT_M_exp$SNP,
  test_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_test_M_out$outcome <- "test_M"
HT_fsh_M_out <- format_data(
  snps=HT_M_exp$SNP,
  fsh_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_fsh_M_out$outcome <- "fsh_M"
HT_lh_M_out <- format_data(
  snps=HT_M_exp$SNP,
  lh_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_lh_M_out$outcome <- "lh_male"

HT_oest_M_out <- format_data(
  snps=HT_M_exp$SNP,
  oest_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_oest_M_out$outcome <- "oest_M"

LR_test_M_out <- format_data(
  snps=LR_M_exp$SNP,
  test_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_test_M_out$outcome <- "test_M"
LR_fsh_M_out <- format_data(
  snps=LR_M_exp$SNP,
  fsh_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_fsh_M_out$outcome <- "fsh_M"
LR_lh_M_out <- format_data(
  snps=LR_M_exp$SNP,
  lh_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_lh_M_out$outcome <- "lh_M"
LR_oest_M_out <- format_data(
  snps=LR_M_exp$SNP,
  oest_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
LR_oest_M_out$outcome <- "oest_M"

HT_GM_test_M_out <- format_data(
  snps=HT_GM_M_exp$SNP,
  test_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_test_M_out$outcome <- "test_M"
HT_GM_fsh_M_out <- format_data(
  snps=HT_GM_M_exp$SNP,
  fsh_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_fsh_M_out$outcome <- "fsh_M"
HT_GM_lh_M_out <- format_data(
  snps=HT_GM_M_exp$SNP,
  lh_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_lh_M_out$outcome <- "lh_M"
HT_GM_oest_M_out <- format_data(
  snps=HT_GM_M_exp$SNP,
  oest_male_eur,
  type="outcome",
  snp_col="RSID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="Allele1",
  other_allele_col="Allele2",
  eaf_col="Freq1",
  pval_col="PVALUE",
  samplesize_col="N",
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
HT_GM_oest_M_out$outcome <- "oest_M"

#hormone as exposure, brain as outcome
test_HT_M_out <- format_data(
  snps=test_M_exp$SNP,
  HT_vol_M,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_HT_M_out$outcome <- "HT_M"
test_PG_M_out <- format_data(
  snps=test_M_exp$SNP,
  PG_vol_M,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_PG_M_out$outcome <- "PG_M"
test_LR_M_out <- format_data(
  snps=test_M_exp$SNP,
  LR_vol_M,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_LR_M_out$outcome <- "LR_M"
test_HT_GM_M_out <- format_data(
  snps=test_M_exp$SNP,
  HT_GM_vol_M,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
test_HT_GM_M_out$outcome <- "HT_GM_M"

fsh_HT_M_out <- format_data(
  snps=fsh_M_exp$SNP,
  HT_vol_M,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_HT_M_out$outcome <- "HT_M"
fsh_PG_M_out <- format_data(
  snps=fsh_M_exp$SNP,
  PG_vol_M,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_PG_M_out$outcome <- "PG_M"
fsh_LR_M_out <- format_data(
  snps=fsh_M_exp$SNP,
  LR_vol_M,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_LR_M_out$outcome <- "LR_M"
fsh_HT_GM_M_out <- format_data(
  snps=fsh_M_exp$SNP,
  HT_GM_vol_M,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
fsh_HT_GM_M_out$outcome <- "HT_GM_M"

#DOES NOT WORK -  I think because of lack of crossover of genetics
#lh_HT_M_out <- format_data(
#  snps=lh_M_exp$SNP,
#  HT_vol_M,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
#)
#lh_HT_M_out$outcome <- "HT_M"
#DOES NOT WORK -  I think because of lack of crossover of genetics
#lh_PG_M_out <- format_data(
#  snps=lh_M_exp$SNP,
#  PG_vol_M,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
#)
#lh_PG_M_out$outcome <- "PG_male"
#lh_LR_M_out <- format_data(
#  snps=lh_M_exp$SNP,
#  LR_vol_M,
#  type="outcome",
#  snp_col="ID",
#  se_col="SE",
#  beta_col="BETA",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
#)
#lh_LR_M_out$outcome <- "LR_male"
#DOES NOT WORK -  I think because of lack of crossover of genetics
#lh_HT_GM_M_out <- format_data(
#  snps=lh_M_exp$SNP,
#  HT_GM_vol_M,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
#)
#lh_HT_GM_M_out$outcome <- "HT_GM_male"

oest_HT_M_out <- format_data(
  snps=oest_M_exp$SNP,
  HT_vol_M,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
oest_HT_M_out$outcome <- "HT_M"
oest_PG_M_out <- format_data(
  snps=oest_M_exp$SNP,
  PG_vol_M,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
oest_PG_M_out$outcome <- "PG_M"
oest_LR_M_out <- format_data(
  snps=oest_M_exp$SNP,
  LR_vol_M,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
oest_LR_M_out$outcome <- "LR_M"
oest_HT_GM_M_out <- format_data(
  snps=oest_M_exp$SNP,
  HT_GM_vol_M,
  type="outcome",
  snp_col="ID",
  beta_col="BETA",
  se_col="SE",
  effect_allele_col="ALLELE1",
  other_allele_col="ALLELE0",
  eaf_col="A1FREQ",
  pval_col="LOG10P",
  samplesize_col="N",
  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
)
oest_HT_GM_M_out$outcome <- "HT_GM_M"

###Brain as exposure, fertility as outcome
#brain as exposure and hormone as outcome
#PG_female_infert1_out <- format_data(
#  snps=PG_F_exp$SNP,
#  female_infert1,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#  )
#PG_female_infert1_out$outcome <- "female_infert1"
#PG_female_infert2_out <- format_data(
#  snps=PG_F_exp$SNP,
#  female_infert2,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#PG_female_infert2_out$outcome <- "female_infert2"
#PG_female_infert3_out <- format_data(
#  snps=PG_F_exp$SNP,
#  female_infert3,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
# ncontrol="N_CONTROLS"
#)
#PG_female_infert3_out$outcome <- "female_infert3"
#PG_female_infert4_out <- format_data(
#  snps=PG_F_exp$SNP,
#  female_infert4,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#PG_female_infert4_out$outcome <- "female_infert4"
#PG_female_infert5_out <- format_data(
#  snps=PG_F_exp$SNP,
#  female_infert5,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#PG_female_infert5_out$outcome <- "female_infert5"
#PG_male_infert_out <- format_data(
#  snps=PG_M_exp$SNP,
#  male_infert,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#PG_male_infert_out$outcome <- "male_infert"

#HT_female_infert1_out <- format_data(
#  snps=HT_F_exp$SNP,
#  female_infert1,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#HT_female_infert1_out$outcome <- "female_infert1"
#HT_female_infert2_out <- format_data(
#  snps=HT_F_exp$SNP,
#  female_infert2,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#HT_female_infert2_out$outcome <- "female_infert2"
#HT_female_infert3_out <- format_data(
#  snps=HT_F_exp$SNP,
#  female_infert3,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#HT_female_infert3_out$outcome <- "female_infert3"
#HT_female_infert4_out <- format_data(
#  snps=HT_F_exp$SNP,
#  female_infert4,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#HT_female_infert4_out$outcome <- "female_infert4"
#HT_female_infert5_out <- format_data(
#  snps=HT_F_exp$SNP,
#  female_infert5,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#HT_female_infert5_out$outcome <- "female_infert5"
#HT_male_infert_out <- format_data(
#  snps=HT_M_exp$SNP,
#  male_infert,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#HT_male_infert_out$outcome <- "male_infert"

#LR_female_infert1_out <- format_data(
#  snps=LR_F_exp$SNP,
#  female_infert1,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  effect_allele_col="Allele1",
#  se_col="SE",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#LR_female_infert1_out$outcome <- "female_infert1"
#LR_female_infert2_out <- format_data(
#  snps=LR_F_exp$SNP,
#  female_infert2,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#LR_female_infert2_out$outcome <- "female_infert2"
#LR_female_infert3_out <- format_data(
#  snps=LR_F_exp$SNP,
#  female_infert3,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#LR_female_infert3_out$outcome <- "female_infert3"
#LR_female_infert4_out <- format_data(
#  snps=LR_F_exp$SNP,
#  female_infert4,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#LR_female_infert4_out$outcome <- "female_infert4"
#LR_female_infert5_out <- format_data(
#  snps=LR_F_exp$SNP,
#  female_infert5,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#LR_female_infert5_out$outcome <- "female_infert5"
#LR_male_infert_out <- format_data(
#  snps=LR_M_exp$SNP,
#  male_infert,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#LR_male_infert_out$outcome <- "male_infert"

#HT_GM_female_infert1_out <- format_data(
#  snps=HT_GM_F_exp$SNP,
#  female_infert1,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#HT_GM_female_infert1_out$outcome <- "female_infert1"
#HT_GM_female_infert2_out <- format_data(
#  snps=HT_GM_F_exp$SNP,
#  female_infert2,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#HT_GM_female_infert2_out$outcome <- "female_infert2"
#HT_GM_female_infert3_out <- format_data(
#  snps=HT_GM_F_exp$SNP,
#  female_infert3,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#HT_GM_female_infert3_out$outcome <- "female_infert3"
#HT_GM_female_infert4_out <- format_data(
#  snps=HT_GM_F_exp$SNP,
#  female_infert4,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#HT_GM_female_infert4_out$outcome <- "female_infert4"
#HT_GM_female_infert5_out <- format_data(
#  snps=HT_GM_F_exp$SNP,
#  female_infert5,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#HT_GM_female_infert5_out$outcome <- "female_infert5"
#HT_GM_male_infert_out <- format_data(
#  snps=HT_GM_M_exp$SNP,
#  male_infert,
#  type="outcome",
#  snp_col="RSID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="Allele1",
#  other_allele_col="Allele2",
#  eaf_col="Freq1",
#  pval_col="PVALUE",
#  ncase="N_CASES",
#  ncontrol="N_CONTROLS"
#)
#HT_GM_male_infert_out$outcome <- "male_infert"

#female_infert1_HT_out <- format_data(
#  snps=female_infert1_exp$SNP,
#  HT_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert1_HT_out$outcome <- "HT_F"
#female_infert2_HT_out <- format_data(
#  HT_vol_F,
#  snps=female_infert2_exp$SNP,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert2_HT_out$outcome <- "HT_F"
#female_infert3_HT_out <- format_data(
#  snps=female_infert3_exp$SNP,
#  HT_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert3_HT_out$outcome <- "HT_F"
#female_infert4_HT_out <- format_data(
#  snps=female_infert4_exp$SNP,
#  HT_vol_F,
#  type="outcome",
#  beta_col="BETA",
#  snp_col="ID",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert4_HT_out$outcome <- "HT_F"
#female_infert5_HT_out <- format_data(
#  snps=female_infert5_exp$SNP,
#  HT_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert5_HT_out$outcome <- "HT_F"
#male_infert_HT_out <- format_data(
#  snps=male_infert_exp$SNP,
#  HT_vol_M,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#male_infert_HT_out$outcome <- "HT_M"

#female_infert1_PG_out <- format_data(
#  snps=female_infert1_exp$SNP,
#  PG_vol_F, 
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert1_PG_out$outcome <- "PG_F"
#female_infert2_PG_out <- format_data(
#  snps=female_infert2_exp$SNP,
#  PG_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert2_PG_out$outcome <- "PG_F"
#female_infert3_PG_out <- format_data(
#  snps=female_infert3_exp$SNP,
#  PG_vol_F,
#  snp_col="ID",
#  type="outcome",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  pval_col="LOG10P",
#  eaf_col="A1FREQ",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert3_PG_out$outcome <- "PG_F"
#female_infert4_PG_out <- format_data(
#  snps=female_infert4_exp$SNP,
#  PG_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
#  chr_col="CHROM",
#  pos_col="GENPOS"
#)
#female_infert4_PG_out$outcome <- "PG_F"
#female_infert5_PG_out <- format_data(
#  snps=female_infert5_exp$SNP,
#  PG_vol_F,
#  snp_col="ID",
#  beta_col="BETA",
#  type="outcome",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert5_PG_out$outcome <- "PG_F"
#male_infert_PG_out <- format_data(
#  snps=male_infert_exp$SNP,
#  PG_vol_M,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)

##male_infert_PG_out$outcome <- "PG_M"
#female_infert1_LR_out <- format_data(
#  snps=female_infert1_exp$SNP,
#  LR_vol_F,
#  type="outcome",
#  beta_col="BETA",
#  se_col="SE",
#  snp_col="ID",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert1_LR_out$outcome <- "LR_F"
#female_infert2_LR_out <- format_data(
#  snps=female_infert2_exp$SNP,
#  LR_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert2_LR_out$outcome <- "LR_F"
#female_infert3_LR_out <- format_data(
#  snps=female_infert3_exp$SNP,
#  LR_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert3_LR_out$outcome <- "LR_F"
#female_infert4_LR_out <- format_data(
#  snps=female_infert4_exp$SNP,
#  LR_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert4_LR_out$outcome <- "LR_F"
#female_infert5_LR_out <- format_data(
#  snps=female_infert5_exp$SNP,
#  LR_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert5_LR_out$outcome <- "LR_F"
#male_infert_LR_out <- format_data(
#  snps=male_infert_exp$SNP,
#  LR_vol_M,
#  type="outcome",
#  snp_col="ID",
#  se_col="SE",
#  beta_col="BETA",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#male_infert_LR_out$outcome <- "LR_M"

#female_infert1_HT_GM_out <- format_data(
#  snps=female_infert1_exp$SNP,
#  HT_GM_vol_F,
#  type="outcome",
#  beta_col="BETA",
#  snp_col="ID",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert1_HT_GM_out$outcome <- "HT_GM_F"
#female_infert2_HT_GM_out <- format_data(
#  snps=female_infert2_exp$SNP,
#  HT_GM_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert2_HT_GM_out$outcome <- "HT_GM_F"
#female_infert3_HT_GM_out <- format_data(
#  snps=female_infert3_exp$SNP,
#  HT_GM_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert3_HT_GM_out$outcome <- "HT_GM_F"
#female_infert4_HT_GM_out <- format_data(
#  snps=female_infert4_exp$SNP,
#  HT_GM_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert4_HT_GM_out$outcome <- "HT_GM_F"
#female_infert5_HT_GM_out <- format_data(
#  snps=female_infert5_exp$SNP,
#  HT_GM_vol_F,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#female_infert5_HT_GM_out$outcome <- "HT_GM_F"
#male_infert_HT_GM_out <- format_data(
#  snps=male_infert_exp$SNP,
#  HT_GM_vol_M,
#  type="outcome",
#  snp_col="ID",
#  beta_col="BETA",
#  se_col="SE",
#  effect_allele_col="ALLELE1",
#  other_allele_col="ALLELE0",
#  eaf_col="A1FREQ",
#  pval_col="LOG10P",
#  samplesize_col="N",
#  log_pval=TRUE,
##  chr_col="CHROM",
##  pos_col="GENPOS"
#)
#male_infert_HT_GM_out$outcome <- "HT_GM_M"

###HARMONISE DATA###

PG_test_harm <- harmonise_data(
  exposure_dat = PG_exp,
  outcome_dat = PG_test_out
)
PG_fsh_harm <- harmonise_data(
  exposure_dat = PG_exp,
  outcome_dat = PG_fsh_out
)
PG_lh_harm <- harmonise_data(
  exposure_dat = PG_exp,
  outcome_dat = PG_lh_out
)
PG_prog_harm <- harmonise_data(
  exposure_dat = PG_exp,
  outcome_dat = PG_prog_out
)
PG_oest_harm <- harmonise_data(
  exposure_dat = PG_exp,
  outcome_dat = PG_oest_out
)

HT_test_harm <- harmonise_data(
  exposure_dat = HT_exp,
  outcome_dat = HT_test_out
)
HT_fsh_harm <- harmonise_data(
  exposure_dat = HT_exp,
  outcome_dat = HT_fsh_out
)
HT_lh_harm <- harmonise_data(
  exposure_dat = HT_exp,
  outcome_dat = HT_lh_out
)
HT_prog_harm <- harmonise_data(
  exposure_dat = HT_exp,
  outcome_dat = HT_prog_out
)
HT_oest_harm <- harmonise_data(
  exposure_dat = HT_exp,
  outcome_dat = HT_oest_out
)

LR_test_harm <- harmonise_data(
  exposure_dat = LR_exp,
  outcome_dat = LR_test_out
)
LR_fsh_harm <- harmonise_data(
  exposure_dat = LR_exp,
  outcome_dat = LR_fsh_out
)
LR_lh_harm <- harmonise_data(
  exposure_dat = LR_exp,
  outcome_dat = LR_lh_out
)
LR_prog_harm <- harmonise_data(
  exposure_dat = LR_exp,
  outcome_dat = LR_prog_out
)
LR_oest_harm <- harmonise_data(
  exposure_dat = LR_exp,
  outcome_dat = LR_oest_out
)

HT_GM_test_harm <- harmonise_data(
  exposure_dat = HT_GM_exp,
  outcome_dat = HT_GM_test_out
)
HT_GM_fsh_harm <- harmonise_data(
  exposure_dat = HT_GM_exp,
  outcome_dat = HT_GM_fsh_out
)
HT_GM_lh_harm <- harmonise_data(
  exposure_dat = HT_GM_exp,
  outcome_dat = HT_GM_lh_out
)
HT_GM_prog_harm <- harmonise_data(
  exposure_dat = HT_GM_exp,
  outcome_dat = HT_GM_prog_out
)
HT_GM_oest_harm <- harmonise_data(
  exposure_dat = HT_GM_exp,
  outcome_dat = HT_GM_oest_out
)

test_HT_harm <- harmonise_data(
  exposure_dat = test_exp,
  outcome_dat = test_HT_out
)
test_PG_harm <- harmonise_data(
  exposure_dat = test_exp,
  outcome_dat = test_PG_out
)
test_LR_harm <- harmonise_data(
  exposure_dat = test_exp,
  outcome_dat = test_LR_out
)
test_HT_GM_harm <- harmonise_data(
  exposure_dat = test_exp,
  outcome_dat = test_HT_GM_out
)

fsh_HT_harm <- harmonise_data(
  exposure_dat = fsh_exp,
  outcome_dat = fsh_HT_out
)
fsh_PG_harm <- harmonise_data(
  exposure_dat = fsh_exp,
  outcome_dat = fsh_PG_out
)
fsh_LR_harm <- harmonise_data(
  exposure_dat = fsh_exp,
  outcome_dat = fsh_LR_out
)
fsh_HT_GM_harm <- harmonise_data(
  exposure_dat = fsh_exp,
  outcome_dat = fsh_HT_GM_out
)

lh_HT_harm <- harmonise_data(
  exposure_dat = lh_exp,
  outcome_dat = lh_HT_out
)
lh_PG_harm <- harmonise_data(
  exposure_dat = lh_exp,
  outcome_dat = lh_PG_out
)
lh_LR_harm <- harmonise_data(
  exposure_dat = lh_exp,
  outcome_dat = lh_LR_out
)
lh_HT_GM_harm <- harmonise_data(
  exposure_dat = lh_exp,
  outcome_dat = lh_HT_GM_out
)

oest_HT_harm <- harmonise_data(
  exposure_dat = oest_exp,
  outcome_dat = oest_HT_out
)
oest_PG_harm <- harmonise_data(
  exposure_dat = oest_exp,
  outcome_dat = oest_PG_out
)
oest_LR_harm <- harmonise_data(
  exposure_dat = oest_exp,
  outcome_dat = oest_LR_out
)
oest_HT_GM_harm <- harmonise_data(
  exposure_dat = oest_exp,
  outcome_dat = oest_HT_GM_out
)

###FEMALE
PG_test_F_harm <- harmonise_data(
  exposure_dat = PG_F_exp,
  outcome_dat = PG_test_F_out
)
PG_fsh_F_harm <- harmonise_data(
  exposure_dat = PG_F_exp,
  outcome_dat = PG_fsh_F_out
)
PG_lh_F_harm <- harmonise_data(
  exposure_dat = PG_F_exp,
  outcome_dat = PG_lh_F_out
)
PG_prog_F_harm <- harmonise_data(
  exposure_dat = PG_F_exp,
  outcome_dat = PG_prog_F_out
)
PG_oest_F_harm <- harmonise_data(
  exposure_dat = PG_F_exp,
  outcome_dat = PG_oest_F_out
)

HT_test_F_harm <- harmonise_data(
  exposure_dat = HT_F_exp,
  outcome_dat = HT_test_F_out
)
HT_fsh_F_harm <- harmonise_data(
  exposure_dat = HT_F_exp,
  outcome_dat = HT_fsh_F_out
)
HT_lh_F_harm <- harmonise_data(
  exposure_dat = HT_F_exp,
  outcome_dat = HT_lh_F_out
)
HT_prog_F_harm <- harmonise_data(
  exposure_dat = HT_F_exp,
  outcome_dat = HT_prog_F_out
)
HT_oest_F_harm <- harmonise_data(
  exposure_dat = HT_F_exp,
  outcome_dat = HT_oest_F_out
)

LR_test_F_harm <- harmonise_data(
  exposure_dat = LR_F_exp,
  outcome_dat = LR_test_F_out
)
LR_fsh_F_harm <- harmonise_data(
  exposure_dat = LR_F_exp,
  outcome_dat = LR_fsh_F_out
)
LR_lh_F_harm <- harmonise_data(
  exposure_dat = LR_F_exp,
  outcome_dat = LR_lh_F_out
)
LR_prog_F_harm <- harmonise_data(
  exposure_dat = LR_F_exp,
  outcome_dat = LR_prog_F_out
)
LR_oest_F_harm <- harmonise_data(
  exposure_dat = LR_F_exp,
  outcome_dat = LR_oest_F_out
)

HT_GM_test_F_harm <- harmonise_data(
  exposure_dat = HT_GM_F_exp,
  outcome_dat = HT_GM_test_F_out
)
HT_GM_fsh_F_harm <- harmonise_data(
  exposure_dat = HT_GM_F_exp,
  outcome_dat = HT_GM_fsh_F_out
)
HT_GM_lh_F_harm <- harmonise_data(
  exposure_dat = HT_GM_F_exp,
  outcome_dat = HT_GM_lh_F_out
)
HT_GM_prog_F_harm <- harmonise_data(
  exposure_dat = HT_GM_F_exp,
  outcome_dat = HT_GM_prog_F_out
)
HT_GM_oest_F_harm <- harmonise_data(
  exposure_dat = HT_GM_F_exp,
  outcome_dat = HT_GM_oest_F_out
)

test_HT_F_harm <- harmonise_data(
  exposure_dat = test_F_exp,
  outcome_dat = test_HT_F_out
)
test_PG_F_harm <- harmonise_data(
  exposure_dat = test_F_exp,
  outcome_dat = test_PG_F_out
)
test_LR_F_harm <- harmonise_data(
  exposure_dat = test_F_exp,
  outcome_dat = test_LR_F_out
)
test_HT_GM_F_harm <- harmonise_data(
  exposure_dat = test_F_exp,
  outcome_dat = test_HT_GM_F_out
)

fsh_HT_F_harm <- harmonise_data(
  exposure_dat = fsh_F_exp,
  outcome_dat = fsh_HT_F_out
)
fsh_PG_F_harm <- harmonise_data(
  exposure_dat = fsh_F_exp,
  outcome_dat = fsh_PG_F_out
)
fsh_LR_F_harm <- harmonise_data(
  exposure_dat = fsh_F_exp,
  outcome_dat = fsh_LR_F_out
)
fsh_HT_GM_F_harm <- harmonise_data(
  exposure_dat = fsh_F_exp,
  outcome_dat = fsh_HT_GM_F_out
)

lh_HT_F_harm <- harmonise_data(
  exposure_dat = lh_F_exp,
  outcome_dat = lh_HT_F_out
)
lh_PG_F_harm <- harmonise_data(
  exposure_dat = lh_F_exp,
  outcome_dat = lh_PG_F_out
)
lh_LR_F_harm <- harmonise_data(
  exposure_dat = lh_F_exp,
  outcome_dat = lh_LR_F_out
)
lh_HT_GM_F_harm <- harmonise_data(
  exposure_dat = lh_F_exp,
  outcome_dat = lh_HT_GM_F_out
)

#oest_HT_F_harm <- harmonise_data(
#  exposure_dat = oest_F_exp,
#  outcome_dat = oest_HT_F_out
#)
#oest_PG_F_harm <- harmonise_data(
#  exposure_dat = oest_F_exp,
#  outcome_dat = oest_PG_F_out
#)
#oest_LR_F_harm <- harmonise_data(
#  exposure_dat = oest_F_exp,
#  outcome_dat = oest_LR_F_out
#)
#oest_HT_GM_F_harm <- harmonise_data(
#  exposure_dat = oest_F_exp,
#  outcome_dat = oest_HT_GM_F_out
#)

###MALE
PG_test_M_harm <- harmonise_data(
  exposure_dat = PG_M_exp,
  outcome_dat = PG_test_M_out
)
PG_fsh_M_harm <- harmonise_data(
  exposure_dat = PG_M_exp,
  outcome_dat = PG_fsh_M_out
)
PG_lh_M_harm <- harmonise_data(
  exposure_dat = PG_M_exp,
  outcome_dat = PG_lh_M_out
)
PG_oest_M_harm <- harmonise_data(
  exposure_dat = PG_M_exp,
  outcome_dat = PG_oest_M_out
)

HT_test_M_harm <- harmonise_data(
  exposure_dat = HT_M_exp,
  outcome_dat = HT_test_M_out
)
HT_fsh_M_harm <- harmonise_data(
  exposure_dat = HT_M_exp,
  outcome_dat = HT_fsh_M_out
)
HT_lh_M_harm <- harmonise_data(
  exposure_dat = HT_M_exp,
  outcome_dat = HT_lh_M_out
)
HT_oest_M_harm <- harmonise_data(
  exposure_dat = HT_M_exp,
  outcome_dat = HT_oest_M_out
)

LR_test_M_harm <- harmonise_data(
  exposure_dat = LR_M_exp,
  outcome_dat = LR_test_M_out
)
LR_fsh_M_harm <- harmonise_data(
  exposure_dat = LR_M_exp,
  outcome_dat = LR_fsh_M_out
)
LR_lh_M_harm <- harmonise_data(
  exposure_dat = LR_M_exp,
  outcome_dat = LR_lh_M_out
)
LR_oest_M_harm <- harmonise_data(
  exposure_dat = LR_M_exp,
  outcome_dat = LR_oest_M_out
)

HT_GM_test_M_harm <- harmonise_data(
  exposure_dat = HT_GM_M_exp,
  outcome_dat = HT_GM_test_M_out
)
HT_GM_fsh_M_harm <- harmonise_data(
  exposure_dat = HT_GM_M_exp,
  outcome_dat = HT_GM_fsh_M_out
)
HT_GM_lh_M_harm <- harmonise_data(
  exposure_dat = HT_GM_M_exp,
  outcome_dat = HT_GM_lh_M_out
)
HT_GM_oest_M_harm <- harmonise_data(
  exposure_dat = HT_GM_M_exp,
  outcome_dat = HT_GM_oest_M_out
)

test_HT_M_harm <- harmonise_data(
  exposure_dat = test_M_exp,
  outcome_dat = test_HT_M_out
)
test_PG_M_harm <- harmonise_data(
  exposure_dat = test_M_exp,
  outcome_dat = test_PG_M_out
)
test_LR_M_harm <- harmonise_data(
  exposure_dat = test_M_exp,
  outcome_dat = test_LR_M_out
)
test_HT_GM_M_harm <- harmonise_data(
  exposure_dat = test_M_exp,
  outcome_dat = test_HT_GM_M_out
)

fsh_HT_M_harm <- harmonise_data(
  exposure_dat = fsh_M_exp,
  outcome_dat = fsh_HT_M_out
)
fsh_PG_M_harm <- harmonise_data(
  exposure_dat = fsh_M_exp,
  outcome_dat = fsh_PG_M_out
)
fsh_LR_M_harm <- harmonise_data(
  exposure_dat = fsh_M_exp,
  outcome_dat = fsh_LR_M_out
)
fsh_HT_GM_M_harm <- harmonise_data(
  exposure_dat = fsh_M_exp,
  outcome_dat = fsh_HT_GM_M_out
)

#lh_HT_M_harm <- harmonise_data(
#  exposure_dat = lh_M_exp,
#  outcome_dat = lh_HT_M_out
#)
#lh_PG_M_harm <- harmonise_data(
#  exposure_dat = lh_M_exp,
#  outcome_dat = lh_PG_M_out
#)
#lh_LR_M_harm <- harmonise_data(
#  exposure_dat = lh_M_exp,
#  outcome_dat = lh_LR_M_out
#)
#lh_HT_GM_M_harm <- harmonise_data(
#  exposure_dat = lh_M_exp,
#  outcome_dat = lh_HT_GM_M_out
#)

oest_HT_M_harm <- harmonise_data(
  exposure_dat = oest_M_exp,
  outcome_dat = oest_HT_M_out
)
oest_PG_M_harm <- harmonise_data(
  exposure_dat = oest_M_exp,
  outcome_dat = oest_PG_M_out
)
oest_LR_M_harm <- harmonise_data(
  exposure_dat = oest_M_exp,
  outcome_dat = oest_LR_M_out
)
oest_HT_GM_M_harm <- harmonise_data(
  exposure_dat = oest_M_exp,
  outcome_dat = oest_HT_GM_M_out
)

#HT_female_infert1_harm <- harmonise_data(
#  exposure_dat = HT_F_exp,
#  outcome_dat = HT_female_infert1_out
#)
#HT_female_infert2_harm <- harmonise_data(
#  exposure_dat = HT_F_exp,
#  outcome_dat = HT_female_infert2_out
#)
#HT_female_infert3_harm <- harmonise_data(
#  exposure_dat = HT_F_exp,
#  outcome_dat = HT_female_infert3_out
#)
#HT_female_infert4_harm <- harmonise_data(
#  exposure_dat = HT_F_exp,
#  outcome_dat = HT_female_infert4_out
#)
#HT_female_infert5_harm <- harmonise_data(
#  exposure_dat = HT_F_exp,
#  outcome_dat = HT_female_infert5_out
#)
#HT_male_infert_harm <- harmonise_data(
#  exposure_dat = HT_M_exp,
#  outcome_dat = HT_male_infert_out
#)

#PG_female_infert1_harm <- harmonise_data(
#  exposure_dat = PG_F_exp,
#  outcome_dat = PG_female_infert1_out
#)
#PG_female_infert2_harm <- harmonise_data(
#  exposure_dat = PG_F_exp,
#  outcome_dat = PG_female_infert2_out
#)
#PG_female_infert3_harm <- harmonise_data(
#  exposure_dat = PG_F_exp,
#  outcome_dat = PG_female_infert3_out
#)
#PG_female_infert4_harm <- harmonise_data(
#  exposure_dat = PG_F_exp,
#  outcome_dat = PG_female_infert4_out
#)
#PG_female_infert5_harm <- harmonise_data(
#  exposure_dat = PG_F_exp,
#  outcome_dat = PG_female_infert5_out
#)
#PG_male_infert_harm <- harmonise_data(
#  exposure_dat = PG_M_exp,
#  outcome_dat = PG_male_infert_out
#)

#LR_female_infert1_harm <- harmonise_data(
#  exposure_dat = LR_F_exp,
#  outcome_dat = LR_female_infert1_out
#)
#LR_female_infert2_harm <- harmonise_data(
#  exposure_dat = LR_F_exp,
#  outcome_dat = LR_female_infert2_out
#)
#LR_female_infert3_harm <- harmonise_data(
#  exposure_dat = LR_F_exp,
#  outcome_dat = LR_female_infert3_out
#)
#LR_female_infert4_harm <- harmonise_data(
#  exposure_dat = LR_F_exp,
#  outcome_dat = LR_female_infert4_out
#)
#LR_female_infert5_harm <- harmonise_data(
#  exposure_dat = LR_F_exp,
)
##  outcome_dat = LR_female_infert5_out
#LR_male_infert_harm <- harmonise_data(
#  exposure_dat = LR_M_exp,
#  outcome_dat = LR_male_infert_out
#)

#HT_GM_female_infert1_harm <- harmonise_data(
#  exposure_dat = HT_GM_F_exp,
#  outcome_dat = HT_GM_female_infert1_out
#)
#HT_GM_female_infert2_harm <- harmonise_data(
#  exposure_dat = HT_GM_F_exp,
#  outcome_dat = HT_GM_female_infert2_out
#)
#HT_GM_female_infert3_harm <- harmonise_data(
#  exposure_dat = HT_GM_F_exp,
#  outcome_dat = HT_GM_female_infert3_out
#)
#HT_GM_female_infert4_harm <- harmonise_data(
#  exposure_dat = HT_GM_F_exp,
#  outcome_dat = HT_GM_female_infert4_out
#)
#HT_GM_female_infert5_harm <- harmonise_data(
#  exposure_dat = HT_GM_F_exp,
#  outcome_dat = HT_GM_female_infert5_out
#)
#HT_GM_male_infert_harm <- harmonise_data(
#  exposure_dat = HT_GM_M_exp,
#  outcome_dat = HT_GM_male_infert_out

#)
#female_infert1_HT_harm <- harmonise_data(
#  exposure_dat = female_infert1_exp,
#  outcome_dat = female_infert1_HT_out
#)
#female_infert1_PG_harm <- harmonise_data(
#  exposure_dat = female_infert1_exp,
#  outcome_dat = female_infert1_PG_out
#)
#  exposure_dat = female_infert1_exp,
#female_infert1_LR_harm <- harmonise_data(
#  outcome_dat = female_infert1_LR_out
#female_infert1_HT_GM_harm <- harmonise_data(
##)
#  exposure_dat = female_infert1_exp,
#  outcome_dat = female_infert1_HT_GM_out
#)
#female_infert2_HT_harm <- harmonise_data(
 # exposure_dat = female_infert2_exp,
#  outcome_dat = female_infert2_HT_out
#female_infert2_PG_harm <- harmonise_data(
##)
#  exposure_dat = female_infert2_exp,
#  outcome_dat = female_infert2_PG_out
#)
#female_infert2_LR_harm <- harmonise_data(
#  exposure_dat = female_infert2_exp,
#  outcome_dat = female_infert2_LR_out
#)
#female_infert2_HT_GM_harm <- harmonise_data(
#  exposure_dat = female_infert2_exp,
#  outcome_dat = female_infert2_HT_GM_out
#)
#female_infert3_HT_harm <- harmonise_data(
#  exposure_dat = female_infert3_exp,
#  outcome_dat = female_infert3_HT_out
#)
#female_infert3_PG_harm <- harmonise_data(
#  exposure_dat = female_infert3_exp,
#  outcome_dat = female_infert3_PG_out
#)
#female_infert3_LR_harm <- harmonise_data(
#  exposure_dat = female_infert3_exp,
)
##  outcome_dat = female_infert3_LR_out
#female_infert3_HT_GM_harm <- harmonise_data(
#  exposure_dat = female_infert3_exp,
#  outcome_dat = female_infert3_HT_GM_out
#)
#female_infert4_HT_harm <- harmonise_data(
#  exposure_dat = female_infert4_exp,
#)
#female_infert4_PG_harm <- harmonise_data(
##  outcome_dat = female_infert4_HT_out
#  exposure_dat = female_infert4_exp,
#  outcome_dat = female_infert4_PG_out
#)
#female_infert4_LR_harm <- harmonise_data(
#  exposure_dat = female_infert4_exp,
#)
#female_infert4_HT_GM_harm <- harmonise_data(
##  outcome_dat = female_infert4_LR_out
#  exposure_dat = female_infert4_exp,
#  outcome_dat = female_infert4_HT_GM_out
#)
#female_infert5_HT_harm <- harmonise_data(
#  exposure_dat = female_infert5_exp,
#  outcome_dat = female_infert5_HT_out
#)
#female_infert5_PG_harm <- harmonise_data(
#  exposure_dat = female_infert5_exp,
#  outcome_dat = female_infert5_PG_out
#)
#female_infert5_LR_harm <- harmonise_data(
#  exposure_dat = female_infert5_exp,
#)
##  outcome_dat = female_infert5_LR_out
#  exposure_dat = female_infert5_exp,
##female_infert5_HT_GM_harm <- harmonise_data(
#  outcome_dat = female_infert5_HT_GM_out
#)
#male_infert_HT_harm <- harmonise_data(
#  exposure_dat = male_infert_exp,
#)
#male_infert_PG_harm <- harmonise_data(
##  outcome_dat = male_infert_HT_out
#  exposure_dat = male_infert_exp,
#  outcome_dat = male_infert_PG_out
#)
#male_infert_LR_harm <- harmonise_data(
#  exposure_dat = male_infert_exp,
#  outcome_dat = male_infert_LR_out
#)
#male_infert_HT_GM_harm <- harmonise_data(
#  exposure_dat = male_infert_exp,
#  outcome_dat = male_infert_HT_GM_out
#)

###RUN MR###

PG_test_mr_res <- mr(PG_test_harm)
PG_test_mr_het <- mr_heterogeneity(PG_test_harm)
PG_test_mr_plei <- mr_pleiotropy_test(PG_test_harm)
PG_test_mr_sing <- mr_singlesnp(PG_test_harm)

PG_fsh_mr_res <- mr(PG_fsh_harm)
PG_fsh_mr_het <- mr_heterogeneity(PG_fsh_harm)
PG_fsh_mr_plei <- mr_pleiotropy_test(PG_fsh_harm)
PG_fsh_mr_sing <- mr_singlesnp(PG_fsh_harm)

PG_lh_mr_res <- mr(PG_lh_harm)
PG_lh_mr_het <- mr_heterogeneity(PG_lh_harm)
PG_lh_mr_plei <- mr_pleiotropy_test(PG_lh_harm)
PG_lh_mr_sing <- mr_singlesnp(PG_lh_harm)

PG_prog_mr_res <- mr(PG_prog_harm)
PG_prog_mr_het <- mr_heterogeneity(PG_prog_harm)
PG_prog_mr_plei <- mr_pleiotropy_test(PG_prog_harm)
PG_prog_mr_sing <- mr_singlesnp(PG_prog_harm)

PG_oest_mr_res <- mr(PG_oest_harm)
PG_oest_mr_het <- mr_heterogeneity(PG_oest_harm)
PG_oest_mr_plei <- mr_pleiotropy_test(PG_oest_harm)
PG_oest_mr_sing <- mr_singlesnp(PG_oest_harm)

HT_test_mr_res <- mr(HT_test_harm)
HT_test_mr_het <- mr_heterogeneity(HT_test_harm)
HT_test_mr_plei <- mr_pleiotropy_test(HT_test_harm)
HT_test_mr_sing <- mr_singlesnp(HT_test_harm)
HT_test_direction <- directionality_test(HT_test_harm)

HT_fsh_mr_res <- mr(HT_fsh_harm)
HT_fsh_mr_het <- mr_heterogeneity(HT_fsh_harm)
HT_fsh_mr_plei <- mr_pleiotropy_test(HT_fsh_harm)
HT_fsh_mr_sing <- mr_singlesnp(HT_fsh_harm)

HT_lh_mr_res <- mr(HT_lh_harm)
HT_lh_mr_het <- mr_heterogeneity(HT_lh_harm)
HT_lh_mr_plei <- mr_pleiotropy_test(HT_lh_harm)
HT_lh_mr_sing <- mr_singlesnp(HT_lh_harm)

HT_prog_mr_res <- mr(HT_prog_harm)
HT_prog_mr_het <- mr_heterogeneity(HT_prog_harm)
HT_prog_mr_plei <- mr_pleiotropy_test(HT_prog_harm)
HT_prog_mr_sing <- mr_singlesnp(HT_prog_harm)

HT_oest_mr_res <- mr(HT_oest_harm)
HT_oest_mr_het <- mr_heterogeneity(HT_oest_harm)
HT_oest_mr_plei <- mr_pleiotropy_test(HT_oest_harm)
HT_oest_mr_sing <- mr_singlesnp(HT_oest_harm)

LR_test_mr_res <- mr(LR_test_harm)
LR_test_mr_het <- mr_heterogeneity(LR_test_harm)
LR_test_mr_plei <- mr_pleiotropy_test(LR_test_harm)
LR_test_mr_sing <- mr_singlesnp(LR_test_harm)

LR_fsh_mr_res <- mr(LR_fsh_harm)
LR_fsh_mr_het <- mr_heterogeneity(LR_fsh_harm)
LR_fsh_mr_plei <- mr_pleiotropy_test(LR_fsh_harm)
LR_fsh_mr_sing <- mr_singlesnp(LR_fsh_harm)

LR_lh_mr_res <- mr(LR_lh_harm)
LR_lh_mr_het <- mr_heterogeneity(LR_lh_harm)
LR_lh_mr_plei <- mr_pleiotropy_test(LR_lh_harm)
LR_lh_mr_sing <- mr_singlesnp(LR_lh_harm)

LR_prog_mr_res <- mr(LR_prog_harm)
LR_prog_mr_het <- mr_heterogeneity(LR_prog_harm)
LR_prog_mr_plei <- mr_pleiotropy_test(LR_prog_harm)
LR_prog_mr_sing <- mr_singlesnp(LR_prog_harm)

LR_oest_mr_res <- mr(LR_oest_harm)
LR_oest_mr_het <- mr_heterogeneity(LR_oest_harm)
LR_oest_mr_plei <- mr_pleiotropy_test(LR_oest_harm)
LR_oest_mr_sing <- mr_singlesnp(LR_oest_harm)

HT_GM_test_mr_res <- mr(HT_GM_test_harm)
HT_GM_test_mr_het <- mr_heterogeneity(HT_GM_test_harm)
HT_GM_test_mr_plei <- mr_pleiotropy_test(HT_GM_test_harm)
HT_GM_test_mr_sing <- mr_singlesnp(HT_GM_test_harm)
HT_GM_test_direction <- directionality_test(HT_GM_test_harm)

HT_GM_fsh_mr_res <- mr(HT_GM_fsh_harm)
HT_GM_fsh_mr_het <- mr_heterogeneity(HT_GM_fsh_harm)
HT_GM_fsh_mr_plei <- mr_pleiotropy_test(HT_GM_fsh_harm)
HT_GM_fsh_mr_sing <- mr_singlesnp(HT_GM_fsh_harm)

HT_GM_lh_mr_res <- mr(HT_GM_lh_harm)
HT_GM_lh_mr_het <- mr_heterogeneity(HT_GM_lh_harm)
HT_GM_lh_mr_plei <- mr_pleiotropy_test(HT_GM_lh_harm)
HT_GM_lh_mr_sing <- mr_singlesnp(HT_GM_lh_harm)

HT_GM_prog_mr_res <- mr(HT_GM_prog_harm)
HT_GM_prog_mr_het <- mr_heterogeneity(HT_GM_prog_harm)
HT_GM_prog_mr_plei <- mr_pleiotropy_test(HT_GM_prog_harm)
HT_GM_prog_mr_sing <- mr_singlesnp(HT_GM_prog_harm)

HT_GM_oest_mr_res <- mr(HT_GM_oest_harm)
HT_GM_oest_mr_het <- mr_heterogeneity(HT_GM_oest_harm)
HT_GM_oest_mr_plei <- mr_pleiotropy_test(HT_GM_oest_harm)
HT_GM_oest_mr_sing <- mr_singlesnp(HT_GM_oest_harm)
HT_GM_oest_mr_loo <- mr_leaveoneout(HT_GM_oest_harm)
HT_GM_oest_mr_dir <- directionality_test(HT_GM_oest_harm)


test_HT_mr_res <- mr(test_HT_harm)
test_HT_mr_het <- mr_heterogeneity(test_HT_harm)
test_HT_mr_plei <- mr_pleiotropy_test(test_HT_harm)
test_HT_mr_sing <- mr_singlesnp(test_HT_harm)
#test_HT_direction <- directionality_test(test_HT_harm)

test_PG_mr_res <- mr(test_PG_harm)
test_PG_mr_het <- mr_heterogeneity(test_PG_harm)
test_PG_mr_plei <- mr_pleiotropy_test(test_PG_harm)
test_PG_mr_sing <- mr_singlesnp(test_PG_harm)

test_LR_mr_res <- mr(test_LR_harm)
test_LR_mr_het <- mr_heterogeneity(test_LR_harm)
test_LR_mr_plei <- mr_pleiotropy_test(test_LR_harm)
test_LR_mr_sing <- mr_singlesnp(test_LR_harm)

test_HT_GM_mr_res <- mr(test_HT_GM_harm)
test_HT_GM_mr_het <- mr_heterogeneity(test_HT_GM_harm)
test_HT_GM_mr_plei <- mr_pleiotropy_test(test_HT_GM_harm)
test_HT_GM_mr_sing <- mr_singlesnp(test_HT_GM_harm)

fsh_HT_mr_res <- mr(fsh_HT_harm)
fsh_HT_mr_het <- mr_heterogeneity(fsh_HT_harm)
fsh_HT_mr_plei <- mr_pleiotropy_test(fsh_HT_harm)
fsh_HT_mr_sing <- mr_singlesnp(fsh_HT_harm)
fsh_HT_direction <- directionality_test(fsh_HT_harm)

fsh_PG_mr_res <- mr(fsh_PG_harm)
fsh_PG_mr_het <- mr_heterogeneity(fsh_PG_harm)
fsh_PG_mr_plei <- mr_pleiotropy_test(fsh_PG_harm)
fsh_PG_mr_sing <- mr_singlesnp(fsh_PG_harm)

fsh_LR_mr_res <- mr(fsh_LR_harm)
fsh_LR_mr_het <- mr_heterogeneity(fsh_LR_harm)
fsh_LR_mr_plei <- mr_pleiotropy_test(fsh_LR_harm)
fsh_LR_mr_sing <- mr_singlesnp(fsh_LR_harm)

fsh_HT_GM_mr_res <- mr(fsh_HT_GM_harm)
fsh_HT_GM_mr_het <- mr_heterogeneity(fsh_HT_GM_harm)
fsh_HT_GM_mr_plei <- mr_pleiotropy_test(fsh_HT_GM_harm)
fsh_HT_GM_mr_sing <- mr_singlesnp(fsh_HT_GM_harm)

lh_HT_mr_res <- mr(lh_HT_harm)
lh_HT_mr_het <- mr_heterogeneity(lh_HT_harm)
lh_HT_mr_plei <- mr_pleiotropy_test(lh_HT_harm)
lh_HT_mr_sing <- mr_singlesnp(lh_HT_harm)
lh_HT_direction <- directionality_test(lh_HT_harm)

lh_PG_mr_res <- mr(lh_PG_harm)
lh_PG_mr_het <- mr_heterogeneity(lh_PG_harm)
lh_PG_mr_plei <- mr_pleiotropy_test(lh_PG_harm)
lh_PG_mr_sing <- mr_singlesnp(lh_PG_harm)

lh_LR_mr_res <- mr(lh_LR_harm)
lh_LR_mr_het <- mr_heterogeneity(lh_LR_harm)
lh_LR_mr_plei <- mr_pleiotropy_test(lh_LR_harm)
lh_LR_mr_sing <- mr_singlesnp(lh_LR_harm)

lh_HT_GM_mr_res <- mr(lh_HT_GM_harm)
lh_HT_GM_mr_het <- mr_heterogeneity(lh_HT_GM_harm)
lh_HT_GM_mr_plei <- mr_pleiotropy_test(lh_HT_GM_harm)
lh_HT_GM_mr_sing <- mr_singlesnp(lh_HT_GM_harm)

oest_HT_mr_res <- mr(oest_HT_harm)
oest_HT_mr_het <- mr_heterogeneity(oest_HT_harm)
oest_HT_mr_plei <- mr_pleiotropy_test(oest_HT_harm)
oest_HT_mr_sing <- mr_singlesnp(oest_HT_harm)
oest_HT_direction <- directionality_test(oest_HT_harm)

oest_PG_mr_res <- mr(oest_PG_harm)
oest_PG_mr_het <- mr_heterogeneity(oest_PG_harm)
oest_PG_mr_plei <- mr_pleiotropy_test(oest_PG_harm)
oest_PG_mr_sing <- mr_singlesnp(oest_PG_harm)

oest_LR_mr_res <- mr(oest_LR_harm)
oest_LR_mr_het <- mr_heterogeneity(oest_LR_harm)
oest_LR_mr_plei <- mr_pleiotropy_test(oest_LR_harm)
oest_LR_mr_sing <- mr_singlesnp(oest_LR_harm)

oest_HT_GM_mr_res <- mr(oest_HT_GM_harm)
oest_HT_GM_mr_het <- mr_heterogeneity(oest_HT_GM_harm)
oest_HT_GM_mr_plei <- mr_pleiotropy_test(oest_HT_GM_harm)
oest_HT_GM_mr_sing <- mr_singlesnp(oest_HT_GM_harm)

###FEMALE 

PG_test_F_mr_res <- mr(PG_test_F_harm)
PG_test_F_mr_het <- mr_heterogeneity(PG_test_F_harm)
PG_test_F_mr_plei <- mr_pleiotropy_test(PG_test_F_harm)
PG_test_F_mr_sing <- mr_singlesnp(PG_test_F_harm)

PG_fsh_F_mr_res <- mr(PG_fsh_F_harm)
PG_fsh_F_mr_het <- mr_heterogeneity(PG_fsh_F_harm)
PG_fsh_F_mr_plei <- mr_pleiotropy_test(PG_fsh_F_harm)
PG_fsh_F_mr_sing <- mr_singlesnp(PG_fsh_F_harm)
PG_fsh_f_mr_loo <- mr_leaveoneout(PG_fsh_F_harm)
PG_fsh_f_mr_dir <- directionality_test(PG_fsh_F_harm)

PG_lh_F_mr_res <- mr(PG_lh_F_harm)
PG_lh_F_mr_het <- mr_heterogeneity(PG_lh_F_harm)
PG_lh_F_mr_plei <- mr_pleiotropy_test(PG_lh_F_harm)
PG_lh_F_mr_sing <- mr_singlesnp(PG_lh_F_harm)
PG_lh_F_mr_loo <- mr_leaveoneout(PG_lh_F_harm)
PG_lh_F_mr_dir <- directionality_test(PG_lh_F_harm)

PG_prog_F_mr_res <- mr(PG_prog_F_harm)
PG_prog_F_mr_het <- mr_heterogeneity(PG_prog_F_harm)
PG_prog_F_mr_plei <- mr_pleiotropy_test(PG_prog_F_harm)
PG_prog_F_mr_sing <- mr_singlesnp(PG_prog_F_harm)

PG_oest_F_mr_res <- mr(PG_oest_F_harm)
PG_oest_F_mr_het <- mr_heterogeneity(PG_oest_F_harm)
PG_oest_F_mr_plei <- mr_pleiotropy_test(PG_oest_F_harm)
PG_oest_F_mr_sing <- mr_singlesnp(PG_oest_F_harm)

HT_test_F_mr_res <- mr(HT_test_F_harm)
HT_test_F_mr_het <- mr_heterogeneity(HT_test_F_harm)
HT_test_F_mr_plei <- mr_pleiotropy_test(HT_test_F_harm)
HT_test_F_mr_sing <- mr_singlesnp(HT_test_F_harm)
HT_test_F_direction <- directionality_test(HT_test_F_harm)

HT_fsh_F_mr_res <- mr(HT_fsh_F_harm)
HT_fsh_F_mr_het <- mr_heterogeneity(HT_fsh_F_harm)
HT_fsh_F_mr_plei <- mr_pleiotropy_test(HT_fsh_F_harm)
HT_fsh_F_mr_sing <- mr_singlesnp(HT_fsh_F_harm)

HT_lh_F_mr_res <- mr(HT_lh_F_harm)
HT_lh_F_mr_het <- mr_heterogeneity(HT_lh_F_harm)
HT_lh_F_mr_plei <- mr_pleiotropy_test(HT_lh_F_harm)
HT_lh_F_mr_sing <- mr_singlesnp(HT_lh_F_harm)

HT_prog_F_mr_res <- mr(HT_prog_F_harm)
HT_prog_F_mr_het <- mr_heterogeneity(HT_prog_F_harm)
HT_prog_F_mr_plei <- mr_pleiotropy_test(HT_prog_F_harm)
HT_prog_F_mr_sing <- mr_singlesnp(HT_prog_F_harm)

HT_oest_F_mr_res <- mr(HT_oest_F_harm)
HT_oest_F_mr_het <- mr_heterogeneity(HT_oest_F_harm)
HT_oest_F_mr_plei <- mr_pleiotropy_test(HT_oest_F_harm)
HT_oest_F_mr_sing <- mr_singlesnp(HT_oest_F_harm)

LR_test_F_mr_res <- mr(LR_test_F_harm)
LR_test_F_mr_het <- mr_heterogeneity(LR_test_F_harm)
LR_test_F_mr_plei <- mr_pleiotropy_test(LR_test_F_harm)
LR_test_F_mr_sing <- mr_singlesnp(LR_test_F_harm)

LR_fsh_F_mr_res <- mr(LR_fsh_F_harm)
LR_fsh_F_mr_het <- mr_heterogeneity(LR_fsh_F_harm)
LR_fsh_F_mr_plei <- mr_pleiotropy_test(LR_fsh_F_harm)
LR_fsh_F_mr_sing <- mr_singlesnp(LR_fsh_F_harm)

LR_lh_F_mr_res <- mr(LR_lh_F_harm)
LR_lh_F_mr_het <- mr_heterogeneity(LR_lh_F_harm)
LR_lh_F_mr_plei <- mr_pleiotropy_test(LR_lh_F_harm)
LR_lh_F_mr_sing <- mr_singlesnp(LR_lh_F_harm)

LR_prog_F_mr_res <- mr(LR_prog_F_harm)
LR_prog_F_mr_het <- mr_heterogeneity(LR_prog_F_harm)
LR_prog_F_mr_plei <- mr_pleiotropy_test(LR_prog_F_harm)
LR_prog_F_mr_sing <- mr_singlesnp(LR_prog_F_harm)

LR_oest_F_mr_res <- mr(LR_oest_F_harm)
LR_oest_F_mr_het <- mr_heterogeneity(LR_oest_F_harm)
LR_oest_F_mr_plei <- mr_pleiotropy_test(LR_oest_F_harm)
LR_oest_F_mr_sing <- mr_singlesnp(LR_oest_F_harm)

HT_GM_test_F_mr_res <- mr(HT_GM_test_F_harm)
HT_GM_test_F_mr_het <- mr_heterogeneity(HT_GM_test_F_harm)
HT_GM_test_F_mr_plei <- mr_pleiotropy_test(HT_GM_test_F_harm)
HT_GM_test_F_mr_sing <- mr_singlesnp(HT_GM_test_F_harm)
HT_GM_test_F_direction <- directionality_test(HT_GM_test_F_harm)

HT_GM_fsh_F_mr_res <- mr(HT_GM_fsh_F_harm)
HT_GM_fsh_F_mr_het <- mr_heterogeneity(HT_GM_fsh_F_harm)
HT_GM_fsh_F_mr_plei <- mr_pleiotropy_test(HT_GM_fsh_F_harm)
HT_GM_fsh_F_mr_sing <- mr_singlesnp(HT_GM_fsh_F_harm)

HT_GM_lh_F_mr_res <- mr(HT_GM_lh_F_harm)
HT_GM_lh_F_mr_het <- mr_heterogeneity(HT_GM_lh_F_harm)
HT_GM_lh_F_mr_plei <- mr_pleiotropy_test(HT_GM_lh_F_harm)
HT_GM_lh_F_mr_sing <- mr_singlesnp(HT_GM_lh_F_harm)

HT_GM_prog_F_mr_res <- mr(HT_GM_prog_F_harm)
HT_GM_prog_F_mr_het <- mr_heterogeneity(HT_GM_prog_F_harm)
HT_GM_prog_F_mr_plei <- mr_pleiotropy_test(HT_GM_prog_F_harm)
HT_GM_prog_F_mr_sing <- mr_singlesnp(HT_GM_prog_F_harm)

HT_GM_oest_F_mr_res <- mr(HT_GM_oest_F_harm)
HT_GM_oest_F_mr_het <- mr_heterogeneity(HT_GM_oest_F_harm)
HT_GM_oest_F_mr_plei <- mr_pleiotropy_test(HT_GM_oest_F_harm)
HT_GM_oest_F_mr_sing <- mr_singlesnp(HT_GM_oest_F_harm)

test_HT_F_mr_res <- mr(test_HT_F_harm)
test_HT_F_mr_het <- mr_heterogeneity(test_HT_F_harm)
test_HT_F_mr_plei <- mr_pleiotropy_test(test_HT_F_harm)
test_HT_F_mr_sing <- mr_singlesnp(test_HT_F_harm)
#test_HT_F_direction <- directionality_test(test_HT_F_harm)

test_PG_F_mr_res <- mr(test_PG_F_harm)
test_PG_F_mr_het <- mr_heterogeneity(test_PG_F_harm)
test_PG_F_mr_plei <- mr_pleiotropy_test(test_PG_F_harm)
test_PG_F_mr_sing <- mr_singlesnp(test_PG_F_harm)

test_LR_F_mr_res <- mr(test_LR_F_harm)
test_LR_F_mr_het <- mr_heterogeneity(test_LR_F_harm)
test_LR_F_mr_plei <- mr_pleiotropy_test(test_LR_F_harm)
test_LR_F_mr_sing <- mr_singlesnp(test_LR_F_harm)

test_HT_GM_F_mr_res <- mr(test_HT_GM_F_harm)
test_HT_GM_F_mr_het <- mr_heterogeneity(test_HT_GM_F_harm)
test_HT_GM_F_mr_plei <- mr_pleiotropy_test(test_HT_GM_F_harm)
test_HT_GM_F_mr_sing <- mr_singlesnp(test_HT_GM_F_harm)

fsh_HT_F_mr_res <- mr(fsh_HT_F_harm)
fsh_HT_F_mr_het <- mr_heterogeneity(fsh_HT_F_harm)
fsh_HT_F_mr_plei <- mr_pleiotropy_test(fsh_HT_F_harm)
fsh_HT_F_mr_sing <- mr_singlesnp(fsh_HT_F_harm)
fsh_HT_F_direction <- directionality_test(fsh_HT_F_harm)

fsh_PG_F_mr_res <- mr(fsh_PG_F_harm)
fsh_PG_F_mr_het <- mr_heterogeneity(fsh_PG_F_harm)
fsh_PG_F_mr_plei <- mr_pleiotropy_test(fsh_PG_F_harm)
fsh_PG_F_mr_sing <- mr_singlesnp(fsh_PG_F_harm)

fsh_LR_F_mr_res <- mr(fsh_LR_F_harm)
fsh_LR_F_mr_het <- mr_heterogeneity(fsh_LR_F_harm)
fsh_LR_F_mr_plei <- mr_pleiotropy_test(fsh_LR_F_harm)
fsh_LR_F_mr_sing <- mr_singlesnp(fsh_LR_F_harm)

fsh_HT_GM_F_mr_res <- mr(fsh_HT_GM_F_harm)
fsh_HT_GM_F_mr_het <- mr_heterogeneity(fsh_HT_GM_F_harm)
fsh_HT_GM_F_mr_plei <- mr_pleiotropy_test(fsh_HT_GM_F_harm)
fsh_HT_GM_F_mr_sing <- mr_singlesnp(fsh_HT_GM_F_harm)

lh_HT_F_mr_res <- mr(lh_HT_F_harm)
lh_HT_F_mr_het <- mr_heterogeneity(lh_HT_F_harm)
lh_HT_F_mr_plei <- mr_pleiotropy_test(lh_HT_F_harm)
lh_HT_F_mr_sing <- mr_singlesnp(lh_HT_F_harm)
lh_HT_F_direction <- directionality_test(lh_HT_F_harm)

lh_PG_F_mr_res <- mr(lh_PG_F_harm)
lh_PG_F_mr_het <- mr_heterogeneity(lh_PG_F_harm)
lh_PG_F_mr_plei <- mr_pleiotropy_test(lh_PG_F_harm)
lh_PG_F_mr_sing <- mr_singlesnp(lh_PG_F_harm)

lh_LR_F_mr_res <- mr(lh_LR_F_harm)
lh_LR_F_mr_het <- mr_heterogeneity(lh_LR_F_harm)
lh_LR_F_mr_plei <- mr_pleiotropy_test(lh_LR_F_harm)
lh_LR_F_mr_sing <- mr_singlesnp(lh_LR_F_harm)

lh_HT_GM_F_mr_res <- mr(lh_HT_GM_F_harm)
lh_HT_GM_F_mr_het <- mr_heterogeneity(lh_HT_GM_F_harm)
lh_HT_GM_F_mr_plei <- mr_pleiotropy_test(lh_HT_GM_F_harm)
lh_HT_GM_F_mr_sing <- mr_singlesnp(lh_HT_GM_F_harm)

#oest_HT_F_mr_res <- mr(oest_HT_F_harm)
#oest_HT_F_mr_het <- mr_heterogeneity(oest_HT_F_harm)
#oest_HT_F_mr_plei <- mr_pleiotropy_test(oest_HT_F_harm)
#oest_HT_F_mr_sing <- mr_singlesnp(oest_HT_F_harm)
#oest_HT_F_direction <- directionality_test(oest_HT_F_harm)

#oest_PG_F_mr_res <- mr(oest_PG_F_harm)
#oest_PG_F_mr_het <- mr_heterogeneity(oest_PG_F_harm)
#oest_PG_F_mr_plei <- mr_pleiotropy_test(oest_PG_F_harm)
#oest_PG_F_mr_sing <- mr_singlesnp(oest_PG_F_harm)

#oest_LR_F_mr_res <- mr(oest_LR_F_harm)
#oest_LR_F_mr_het <- mr_heterogeneity(oest_LR_F_harm)
#oest_LR_F_mr_plei <- mr_pleiotropy_test(oest_LR_F_harm)
#oest_LR_F_mr_sing <- mr_singlesnp(oest_LR_F_harm)

#oest_HT_GM_F_mr_res <- mr(oest_HT_GM_F_harm)
#oest_HT_GM_F_mr_het <- mr_heterogeneity(oest_HT_GM_F_harm)
#oest_HT_GM_F_mr_plei <- mr_pleiotropy_test(oest_HT_GM_F_harm)
#oest_HT_GM_F_mr_sing <- mr_singlesnp(oest_HT_GM_F_harm)

###MALE 

PG_test_M_mr_res <- mr(PG_test_M_harm)
PG_test_M_mr_het <- mr_heterogeneity(PG_test_M_harm)
PG_test_M_mr_plei <- mr_pleiotropy_test(PG_test_M_harm)
PG_test_M_mr_sing <- mr_singlesnp(PG_test_M_harm)

PG_fsh_M_mr_res <- mr(PG_fsh_M_harm)
PG_fsh_M_mr_het <- mr_heterogeneity(PG_fsh_M_harm)
PG_fsh_M_mr_plei <- mr_pleiotropy_test(PG_fsh_M_harm)
PG_fsh_M_mr_sing <- mr_singlesnp(PG_fsh_M_harm)

PG_lh_M_mr_res <- mr(PG_lh_M_harm)
PG_lh_M_mr_het <- mr_heterogeneity(PG_lh_M_harm)
PG_lh_M_mr_plei <- mr_pleiotropy_test(PG_lh_M_harm)
PG_lh_M_mr_sing <- mr_singlesnp(PG_lh_M_harm)

PG_oest_M_mr_res <- mr(PG_oest_M_harm)
PG_oest_M_mr_het <- mr_heterogeneity(PG_oest_M_harm)
PG_oest_M_mr_plei <- mr_pleiotropy_test(PG_oest_M_harm)
PG_oest_M_mr_sing <- mr_singlesnp(PG_oest_M_harm)

HT_test_M_mr_res <- mr(HT_test_M_harm)
HT_test_M_mr_het <- mr_heterogeneity(HT_test_M_harm)
HT_test_M_mr_plei <- mr_pleiotropy_test(HT_test_M_harm)
HT_test_M_mr_sing <- mr_singlesnp(HT_test_M_harm)
HT_test_M_direction <- directionality_test(HT_test_M_harm)

HT_fsh_M_mr_res <- mr(HT_fsh_M_harm)
HT_fsh_M_mr_het <- mr_heterogeneity(HT_fsh_M_harm)
HT_fsh_M_mr_plei <- mr_pleiotropy_test(HT_fsh_M_harm)
HT_fsh_M_mr_sing <- mr_singlesnp(HT_fsh_M_harm)

HT_lh_M_mr_res <- mr(HT_lh_M_harm)
HT_lh_M_mr_het <- mr_heterogeneity(HT_lh_M_harm)
HT_lh_M_mr_plei <- mr_pleiotropy_test(HT_lh_M_harm)
HT_lh_M_mr_sing <- mr_singlesnp(HT_lh_M_harm)

HT_oest_M_mr_res <- mr(HT_oest_M_harm)
HT_oest_M_mr_het <- mr_heterogeneity(HT_oest_M_harm)
HT_oest_M_mr_plei <- mr_pleiotropy_test(HT_oest_M_harm)
HT_oest_M_mr_sing <- mr_singlesnp(HT_oest_M_harm)

LR_test_M_mr_res <- mr(LR_test_M_harm)
LR_test_M_mr_het <- mr_heterogeneity(LR_test_M_harm)
LR_test_M_mr_plei <- mr_pleiotropy_test(LR_test_M_harm)
LR_test_M_mr_sing <- mr_singlesnp(LR_test_M_harm)

LR_fsh_M_mr_res <- mr(LR_fsh_M_harm)
LR_fsh_M_mr_het <- mr_heterogeneity(LR_fsh_M_harm)
LR_fsh_M_mr_plei <- mr_pleiotropy_test(LR_fsh_M_harm)
LR_fsh_M_mr_sing <- mr_singlesnp(LR_fsh_M_harm)

LR_lh_M_mr_res <- mr(LR_lh_M_harm)
LR_lh_M_mr_het <- mr_heterogeneity(LR_lh_M_harm)
LR_lh_M_mr_plei <- mr_pleiotropy_test(LR_lh_M_harm)
LR_lh_M_mr_sing <- mr_singlesnp(LR_lh_M_harm)

LR_oest_M_mr_res <- mr(LR_oest_M_harm)
LR_oest_M_mr_het <- mr_heterogeneity(LR_oest_M_harm)
LR_oest_M_mr_plei <- mr_pleiotropy_test(LR_oest_M_harm)
LR_oest_M_mr_sing <- mr_singlesnp(LR_oest_M_harm)

HT_GM_test_M_mr_res <- mr(HT_GM_test_M_harm)
HT_GM_test_M_mr_het <- mr_heterogeneity(HT_GM_test_M_harm)
HT_GM_test_M_mr_plei <- mr_pleiotropy_test(HT_GM_test_M_harm)
HT_GM_test_M_mr_sing <- mr_singlesnp(HT_GM_test_M_harm)
HT_GM_test_M_direction <- directionality_test(HT_GM_test_M_harm)

HT_GM_fsh_M_mr_res <- mr(HT_GM_fsh_M_harm)
HT_GM_fsh_M_mr_het <- mr_heterogeneity(HT_GM_fsh_M_harm)
HT_GM_fsh_M_mr_plei <- mr_pleiotropy_test(HT_GM_fsh_M_harm)
HT_GM_fsh_M_mr_sing <- mr_singlesnp(HT_GM_fsh_M_harm)

HT_GM_lh_M_mr_res <- mr(HT_GM_lh_M_harm)
HT_GM_lh_M_mr_het <- mr_heterogeneity(HT_GM_lh_M_harm)
HT_GM_lh_M_mr_plei <- mr_pleiotropy_test(HT_GM_lh_M_harm)
HT_GM_lh_M_mr_sing <- mr_singlesnp(HT_GM_lh_M_harm)

HT_GM_oest_M_mr_res <- mr(HT_GM_oest_M_harm)
HT_GM_oest_M_mr_het <- mr_heterogeneity(HT_GM_oest_M_harm)
HT_GM_oest_M_mr_plei <- mr_pleiotropy_test(HT_GM_oest_M_harm)
HT_GM_oest_M_mr_sing <- mr_singlesnp(HT_GM_oest_M_harm)

test_HT_M_mr_res <- mr(test_HT_M_harm)
test_HT_M_mr_het <- mr_heterogeneity(test_HT_M_harm)
test_HT_M_mr_plei <- mr_pleiotropy_test(test_HT_M_harm)
test_HT_M_mr_sing <- mr_singlesnp(test_HT_M_harm)
test_HT_M_direction <- directionality_test(test_HT_M_harm)

test_PG_M_mr_res <- mr(test_PG_M_harm)
test_PG_M_mr_het <- mr_heterogeneity(test_PG_M_harm)
test_PG_M_mr_plei <- mr_pleiotropy_test(test_PG_M_harm)
test_PG_M_mr_sing <- mr_singlesnp(test_PG_M_harm)

test_LR_M_mr_res <- mr(test_LR_M_harm)
test_LR_M_mr_het <- mr_heterogeneity(test_LR_M_harm)
test_LR_M_mr_plei <- mr_pleiotropy_test(test_LR_M_harm)
test_LR_M_mr_sing <- mr_singlesnp(test_LR_M_harm)

test_HT_GM_M_mr_res <- mr(test_HT_GM_M_harm)
test_HT_GM_M_mr_het <- mr_heterogeneity(test_HT_GM_M_harm)
test_HT_GM_M_mr_plei <- mr_pleiotropy_test(test_HT_GM_M_harm)
test_HT_GM_M_mr_sing <- mr_singlesnp(test_HT_GM_M_harm)

fsh_HT_M_mr_res <- mr(fsh_HT_M_harm)
fsh_HT_M_mr_het <- mr_heterogeneity(fsh_HT_M_harm)
fsh_HT_M_mr_plei <- mr_pleiotropy_test(fsh_HT_M_harm)
fsh_HT_M_mr_sing <- mr_singlesnp(fsh_HT_M_harm)
fsh_HT_M_direction <- directionality_test(fsh_HT_M_harm)

fsh_PG_M_mr_res <- mr(fsh_PG_M_harm)
fsh_PG_M_mr_het <- mr_heterogeneity(fsh_PG_M_harm)
fsh_PG_M_mr_plei <- mr_pleiotropy_test(fsh_PG_M_harm)
fsh_PG_M_mr_sing <- mr_singlesnp(fsh_PG_M_harm)

fsh_LR_M_mr_res <- mr(fsh_LR_M_harm)
fsh_LR_M_mr_het <- mr_heterogeneity(fsh_LR_M_harm)
fsh_LR_M_mr_plei <- mr_pleiotropy_test(fsh_LR_M_harm)
fsh_LR_M_mr_sing <- mr_singlesnp(fsh_LR_M_harm)

fsh_HT_GM_M_mr_res <- mr(fsh_HT_GM_M_harm)
fsh_HT_GM_M_mr_het <- mr_heterogeneity(fsh_HT_GM_M_harm)
fsh_HT_GM_M_mr_plei <- mr_pleiotropy_test(fsh_HT_GM_M_harm)
fsh_HT_GM_M_mr_sing <- mr_singlesnp(fsh_HT_GM_M_harm)

#lh_HT_M_mr_res <- mr(lh_HT_M_harm)
#lh_HT_M_mr_het <- mr_heterogeneity(lh_HT_M_harm)
#lh_HT_M_mr_plei <- mr_pleiotropy_test(lh_HT_M_harm)
#lh_HT_M_mr_sing <- mr_singlesnp(lh_HT_M_harm)
#lh_HT_M_direction <- directionality_test(lh_HT_M_harm)

#lh_PG_M_mr_res <- mr(lh_PG_M_harm)
#lh_PG_M_mr_het <- mr_heterogeneity(lh_PG_M_harm)
#lh_PG_M_mr_plei <- mr_pleiotropy_test(lh_PG_M_harm)
#lh_PG_M_mr_sing <- mr_singlesnp(lh_PG_M_harm)

#lh_LR_M_mr_res <- mr(lh_LR_M_harm)
#lh_LR_M_mr_het <- mr_heterogeneity(lh_LR_M_harm)
#lh_LR_M_mr_plei <- mr_pleiotropy_test(lh_LR_M_harm)
#lh_LR_M_mr_sing <- mr_singlesnp(lh_LR_M_harm)

#lh_HT_GM_M_mr_res <- mr(lh_HT_GM_M_harm)
#lh_HT_GM_M_mr_het <- mr_heterogeneity(lh_HT_GM_M_harm)
#lh_HT_GM_M_mr_plei <- mr_pleiotropy_test(lh_HT_GM_M_harm)
#lh_HT_GM_M_mr_sing <- mr_singlesnp(lh_HT_GM_M_harm)

oest_HT_M_mr_res <- mr(oest_HT_M_harm)
oest_HT_M_mr_het <- mr_heterogeneity(oest_HT_M_harm)
oest_HT_M_mr_plei <- mr_pleiotropy_test(oest_HT_M_harm)
oest_HT_M_mr_sing <- mr_singlesnp(oest_HT_M_harm)
oest_HT_M_direction <- directionality_test(oest_HT_M_harm)

oest_PG_M_mr_res <- mr(oest_PG_M_harm)
oest_PG_M_mr_het <- mr_heterogeneity(oest_PG_M_harm)
oest_PG_M_mr_plei <- mr_pleiotropy_test(oest_PG_M_harm)
oest_PG_M_mr_sing <- mr_singlesnp(oest_PG_M_harm)

oest_LR_M_mr_res <- mr(oest_LR_M_harm)
oest_LR_M_mr_het <- mr_heterogeneity(oest_LR_M_harm)
oest_LR_M_mr_plei <- mr_pleiotropy_test(oest_LR_M_harm)
oest_LR_M_mr_sing <- mr_singlesnp(oest_LR_M_harm)

oest_HT_GM_M_mr_res <- mr(oest_HT_GM_M_harm)
oest_HT_GM_M_mr_het <- mr_heterogeneity(oest_HT_GM_M_harm)
oest_HT_GM_M_mr_plei <- mr_pleiotropy_test(oest_HT_GM_M_harm)
oest_HT_GM_M_mr_sing <- mr_singlesnp(oest_HT_GM_M_harm)

HT_female_infert1_mr_res <- mr(HT_female_infert1_harm)
HT_female_infert1_mr_het <- mr_heterogeneity(HT_female_infert1_harm)
HT_female_infert1_mr_plei <- mr_pleiotropy_test(HT_female_infert1_harm)
HT_female_infert1_mr_sing <- mr_singlesnp(HT_female_infert1_harm)
HT_female_infert2_mr_res <- mr(HT_female_infert2_harm)
HT_female_infert2_mr_het <- mr_heterogeneity(HT_female_infert2_harm)
HT_female_infert2_mr_plei <- mr_pleiotropy_test(HT_female_infert2_harm)
HT_female_infert2_mr_sing <- mr_singlesnp(HT_female_infert2_harm)
HT_female_infert3_mr_res <- mr(HT_female_infert3_harm)
HT_female_infert3_mr_het <- mr_heterogeneity(HT_female_infert3_harm)
HT_female_infert3_mr_plei <- mr_pleiotropy_test(HT_female_infert3_harm)
HT_female_infert3_mr_sing <- mr_singlesnp(HT_female_infert3_harm)
HT_female_infert4_mr_res <- mr(HT_female_infert4_harm)
HT_female_infert4_mr_het <- mr_heterogeneity(HT_female_infert4_harm)
HT_female_infert4_mr_plei <- mr_pleiotropy_test(HT_female_infert4_harm)
HT_female_infert4_mr_sing <- mr_singlesnp(HT_female_infert4_harm)
HT_female_infert5_mr_res <- mr(HT_female_infert5_harm)
HT_female_infert5_mr_het <- mr_heterogeneity(HT_female_infert5_harm)
HT_female_infert5_mr_plei <- mr_pleiotropy_test(HT_female_infert5_harm)
HT_female_infert5_mr_sing <- mr_singlesnp(HT_female_infert5_harm)
HT_male_infert_mr_res <- mr(HT_male_infert_harm)
HT_male_infert_mr_het <- mr_heterogeneity(HT_male_infert_harm)
HT_male_infert_mr_plei <- mr_pleiotropy_test(HT_male_infert_harm)
HT_male_infert_mr_sing <- mr_singlesnp(HT_male_infert_harm)

PG_female_infert1_mr_res <- mr(PG_female_infert1_harm)
PG_female_infert1_mr_het <- mr_heterogeneity(PG_female_infert1_harm)
PG_female_infert1_mr_plei <- mr_pleiotropy_test(PG_female_infert1_harm)
PG_female_infert1_mr_sing <- mr_singlesnp(PG_female_infert1_harm)
PG_female_infert2_mr_res <- mr(PG_female_infert2_harm)
PG_female_infert2_mr_het <- mr_heterogeneity(PG_female_infert2_harm)
PG_female_infert2_mr_plei <- mr_pleiotropy_test(PG_female_infert2_harm)
PG_female_infert2_mr_sing <- mr_singlesnp(PG_female_infert2_harm)
PG_female_infert3_mr_res <- mr(PG_female_infert3_harm)
PG_female_infert3_mr_het <- mr_heterogeneity(PG_female_infert3_harm)
PG_female_infert3_mr_plei <- mr_pleiotropy_test(PG_female_infert3_harm)
PG_female_infert3_mr_sing <- mr_singlesnp(PG_female_infert3_harm)
PG_female_infert4_mr_res <- mr(PG_female_infert4_harm)
PG_female_infert4_mr_het <- mr_heterogeneity(PG_female_infert4_harm)
PG_female_infert4_mr_plei <- mr_pleiotropy_test(PG_female_infert4_harm)
PG_female_infert4_mr_sing <- mr_singlesnp(PG_female_infert4_harm)
PG_female_infert5_mr_res <- mr(PG_female_infert5_harm)
PG_female_infert5_mr_het <- mr_heterogeneity(PG_female_infert5_harm)
PG_female_infert5_mr_plei <- mr_pleiotropy_test(PG_female_infert5_harm)
PG_female_infert5_mr_sing <- mr_singlesnp(PG_female_infert5_harm)
PG_male_infert_mr_res <- mr(PG_male_infert_harm)
PG_male_infert_mr_het <- mr_heterogeneity(PG_male_infert_harm)
PG_male_infert_mr_plei <- mr_pleiotropy_test(PG_male_infert_harm)
PG_male_infert_mr_sing <- mr_singlesnp(PG_male_infert_harm)

LR_female_infert1_mr_res <- mr(LR_female_infert1_harm)
LR_female_infert1_mr_het <- mr_heterogeneity(LR_female_infert1_harm)
LR_female_infert1_mr_plei <- mr_pleiotropy_test(LR_female_infert1_harm)
LR_female_infert1_mr_sing <- mr_singlesnp(LR_female_infert1_harm)
LR_female_infert2_mr_res <- mr(LR_female_infert2_harm)
LR_female_infert2_mr_het <- mr_heterogeneity(LR_female_infert2_harm)
LR_female_infert2_mr_plei <- mr_pleiotropy_test(LR_female_infert2_harm)
LR_female_infert2_mr_sing <- mr_singlesnp(LR_female_infert2_harm)
LR_female_infert3_mr_res <- mr(LR_female_infert3_harm)
LR_female_infert3_mr_het <- mr_heterogeneity(LR_female_infert3_harm)
LR_female_infert3_mr_plei <- mr_pleiotropy_test(LR_female_infert3_harm)
LR_female_infert3_mr_sing <- mr_singlesnp(LR_female_infert3_harm)
LR_female_infert4_mr_res <- mr(LR_female_infert4_harm)
LR_female_infert4_mr_het <- mr_heterogeneity(LR_female_infert4_harm)
LR_female_infert4_mr_plei <- mr_pleiotropy_test(LR_female_infert4_harm)
LR_female_infert4_mr_sing <- mr_singlesnp(LR_female_infert4_harm)
LR_female_infert5_mr_res <- mr(LR_female_infert5_harm)
LR_female_infert5_mr_het <- mr_heterogeneity(LR_female_infert5_harm)
LR_female_infert5_mr_plei <- mr_pleiotropy_test(LR_female_infert5_harm)
LR_female_infert5_mr_sing <- mr_singlesnp(LR_female_infert5_harm)
LR_male_infert_mr_res <- mr(LR_male_infert_harm)
LR_male_infert_mr_het <- mr_heterogeneity(LR_male_infert_harm)
LR_male_infert_mr_plei <- mr_pleiotropy_test(LR_male_infert_harm)
LR_male_infert_mr_sing <- mr_singlesnp(LR_male_infert_harm)

HT_GM_female_infert1_mr_res <- mr(HT_GM_female_infert1_harm)
HT_GM_female_infert1_mr_het <- mr_heterogeneity(HT_GM_female_infert1_harm)
HT_GM_female_infert1_mr_plei <- mr_pleiotropy_test(HT_GM_female_infert1_harm)
HT_GM_female_infert1_mr_sing <- mr_singlesnp(HT_GM_female_infert1_harm)
HT_GM_female_infert2_mr_res <- mr(HT_GM_female_infert2_harm)
HT_GM_female_infert2_mr_het <- mr_heterogeneity(HT_GM_female_infert2_harm)
HT_GM_female_infert2_mr_plei <- mr_pleiotropy_test(HT_GM_female_infert2_harm)
HT_GM_female_infert2_mr_sing <- mr_singlesnp(HT_GM_female_infert2_harm)
HT_GM_female_infert3_mr_res <- mr(HT_GM_female_infert3_harm)
HT_GM_female_infert3_mr_het <- mr_heterogeneity(HT_GM_female_infert3_harm)
HT_GM_female_infert3_mr_plei <- mr_pleiotropy_test(HT_GM_female_infert3_harm)
HT_GM_female_infert3_mr_sing <- mr_singlesnp(HT_GM_female_infert3_harm)
HT_GM_female_infert4_mr_res <- mr(HT_GM_female_infert4_harm)
HT_GM_female_infert4_mr_het <- mr_heterogeneity(HT_GM_female_infert4_harm)
HT_GM_female_infert4_mr_plei <- mr_pleiotropy_test(HT_GM_female_infert4_harm)
HT_GM_female_infert4_mr_sing <- mr_singlesnp(HT_GM_female_infert4_harm)
HT_GM_female_infert5_mr_res <- mr(HT_GM_female_infert5_harm)
HT_GM_female_infert5_mr_het <- mr_heterogeneity(HT_GM_female_infert5_harm)
HT_GM_female_infert5_mr_plei <- mr_pleiotropy_test(HT_GM_female_infert5_harm)
HT_GM_female_infert5_mr_sing <- mr_singlesnp(HT_GM_female_infert5_harm)
HT_GM_male_infert_mr_res <- mr(HT_GM_male_infert_harm)
HT_GM_male_infert_mr_het <- mr_heterogeneity(HT_GM_male_infert_harm)
HT_GM_male_infert_mr_plei <- mr_pleiotropy_test(HT_GM_male_infert_harm)
HT_GM_male_infert_mr_sing <- mr_singlesnp(HT_GM_male_infert_harm)

female_infert1_HT_mr_res <- mr(female_infert1_HT_harm)
female_infert1_HT_mr_het <- mr_heterogeneity(female_infert1_HT_harm)
female_infert1_HT_mr_plei <- mr_pleiotropy_test(female_infert1_HT_harm)
female_infert1_HT_mr_sing <- mr_singlesnp(female_infert1_HT_harm)
female_infert2_HT_mr_res <- mr(female_infert2_HT_harm)
female_infert2_HT_mr_het <- mr_heterogeneity(female_infert2_HT_harm)
female_infert2_HT_mr_plei <- mr_pleiotropy_test(female_infert2_HT_harm)
female_infert2_HT_mr_sing <- mr_singlesnp(female_infert2_HT_harm)
female_infert3_HT_mr_res <- mr(female_infert3_HT_harm)
female_infert3_HT_mr_het <- mr_heterogeneity(female_infert3_HT_harm)
female_infert3_HT_mr_plei <- mr_pleiotropy_test(female_infert3_HT_harm)
female_infert3_HT_mr_sing <- mr_singlesnp(female_infert3_HT_harm)
female_infert4_HT_mr_res <- mr(female_infert4_HT_harm)
female_infert4_HT_mr_het <- mr_heterogeneity(female_infert4_HT_harm)
female_infert4_HT_mr_plei <- mr_pleiotropy_test(female_infert4_HT_harm)
female_infert4_HT_mr_sing <- mr_singlesnp(female_infert4_HT_harm)
female_infert5_HT_mr_res <- mr(female_infert5_HT_harm)
female_infert5_HT_mr_het <- mr_heterogeneity(female_infert5_HT_harm)
female_infert5_HT_mr_plei <- mr_pleiotropy_test(female_infert5_HT_harm)
female_infert5_HT_mr_sing <- mr_singlesnp(female_infert5_HT_harm)
male_infert_HT_mr_res <- mr(male_infert_HT_harm)
male_infert_HT_mr_het <- mr_heterogeneity(male_infert_HT_harm)
male_infert_HT_mr_plei <- mr_pleiotropy_test(male_infert_HT_harm)
male_infert_HT_mr_sing <- mr_singlesnp(male_infert_HT_harm)

female_infert1_PG_mr_res <- mr(female_infert1_PG_harm)
female_infert1_PG_mr_het <- mr_heterogeneity(female_infert1_PG_harm)
female_infert1_PG_mr_plei <- mr_pleiotropy_test(female_infert1_PG_harm)
female_infert1_PG_mr_sing <- mr_singlesnp(female_infert1_PG_harm)
female_infert2_PG_mr_res <- mr(female_infert2_PG_harm)
female_infert2_PG_mr_het <- mr_heterogeneity(female_infert2_PG_harm)
female_infert2_PG_mr_plei <- mr_pleiotropy_test(female_infert2_PG_harm)
female_infert2_PG_mr_sing <- mr_singlesnp(female_infert2_PG_harm)
female_infert3_PG_mr_res <- mr(female_infert3_PG_harm)
female_infert3_PG_mr_het <- mr_heterogeneity(female_infert3_PG_harm)
female_infert3_PG_mr_plei <- mr_pleiotropy_test(female_infert3_PG_harm)
female_infert3_PG_mr_sing <- mr_singlesnp(female_infert3_PG_harm)
female_infert4_PG_mr_res <- mr(female_infert4_PG_harm)
female_infert4_PG_mr_het <- mr_heterogeneity(female_infert4_PG_harm)
female_infert4_PG_mr_plei <- mr_pleiotropy_test(female_infert4_PG_harm)
female_infert4_PG_mr_sing <- mr_singlesnp(female_infert4_PG_harm)
female_infert5_PG_mr_res <- mr(female_infert5_PG_harm)
female_infert5_PG_mr_het <- mr_heterogeneity(female_infert5_PG_harm)
female_infert5_PG_mr_plei <- mr_pleiotropy_test(female_infert5_PG_harm)
female_infert5_PG_mr_sing <- mr_singlesnp(female_infert5_PG_harm)
male_infert_PG_mr_res <- mr(male_infert_PG_harm)
male_infert_PG_mr_het <- mr_heterogeneity(male_infert_PG_harm)
male_infert_PG_mr_plei <- mr_pleiotropy_test(male_infert_PG_harm)
male_infert_PG_mr_sing <- mr_singlesnp(male_infert_PG_harm)

female_infert1_LR_mr_res <- mr(female_infert1_LR_harm)
female_infert1_LR_mr_het <- mr_heterogeneity(female_infert1_LR_harm)
female_infert1_LR_mr_plei <- mr_pleiotropy_test(female_infert1_LR_harm)
female_infert1_LR_mr_sing <- mr_singlesnp(female_infert1_LR_harm)
female_infert2_LR_mr_res <- mr(female_infert2_LR_harm)
female_infert2_LR_mr_het <- mr_heterogeneity(female_infert2_LR_harm)
female_infert2_LR_mr_plei <- mr_pleiotropy_test(female_infert2_LR_harm)
female_infert2_LR_mr_sing <- mr_singlesnp(female_infert2_LR_harm)
female_infert3_LR_mr_res <- mr(female_infert3_LR_harm)
female_infert3_LR_mr_het <- mr_heterogeneity(female_infert3_LR_harm)
female_infert3_LR_mr_plei <- mr_pleiotropy_test(female_infert3_LR_harm)
female_infert3_LR_mr_sing <- mr_singlesnp(female_infert3_LR_harm)
female_infert4_LR_mr_res <- mr(female_infert4_LR_harm)
female_infert4_LR_mr_het <- mr_heterogeneity(female_infert4_LR_harm)
female_infert4_LR_mr_plei <- mr_pleiotropy_test(female_infert4_LR_harm)
female_infert4_LR_mr_sing <- mr_singlesnp(female_infert4_LR_harm)
female_infert5_LR_mr_res <- mr(female_infert5_LR_harm)
female_infert5_LR_mr_het <- mr_heterogeneity(female_infert5_LR_harm)
female_infert5_LR_mr_plei <- mr_pleiotropy_test(female_infert5_LR_harm)
female_infert5_LR_mr_sing <- mr_singlesnp(female_infert5_LR_harm)
male_infert_LR_mr_res <- mr(male_infert_LR_harm)
male_infert_LR_mr_het <- mr_heterogeneity(male_infert_LR_harm)
male_infert_LR_mr_plei <- mr_pleiotropy_test(male_infert_LR_harm)
male_infert_LR_mr_sing <- mr_singlesnp(male_infert_LR_harm)

female_infert1_HT_GM_mr_res <- mr(female_infert1_HT_GM_harm)
female_infert1_HT_GM_mr_het <- mr_heterogeneity(female_infert1_HT_GM_harm)
female_infert1_HT_GM_mr_plei <- mr_pleiotropy_test(female_infert1_HT_GM_harm)
female_infert1_HT_GM_mr_sing <- mr_singlesnp(female_infert1_HT_GM_harm)
female_infert2_HT_GM_mr_res <- mr(female_infert2_HT_GM_harm)
female_infert2_HT_GM_mr_het <- mr_heterogeneity(female_infert2_HT_GM_harm)
female_infert2_HT_GM_mr_plei <- mr_pleiotropy_test(female_infert2_HT_GM_harm)
female_infert2_HT_GM_mr_sing <- mr_singlesnp(female_infert2_HT_GM_harm)
female_infert3_HT_GM_mr_res <- mr(female_infert3_HT_GM_harm)
female_infert3_HT_GM_mr_het <- mr_heterogeneity(female_infert3_HT_GM_harm)
female_infert3_HT_GM_mr_plei <- mr_pleiotropy_test(female_infert3_HT_GM_harm)
female_infert3_HT_GM_mr_sing <- mr_singlesnp(female_infert3_HT_GM_harm)
female_infert4_HT_GM_mr_res <- mr(female_infert4_HT_GM_harm)
female_infert4_HT_GM_mr_het <- mr_heterogeneity(female_infert4_HT_GM_harm)
female_infert4_HT_GM_mr_plei <- mr_pleiotropy_test(female_infert4_HT_GM_harm)
female_infert4_HT_GM_mr_sing <- mr_singlesnp(female_infert4_HT_GM_harm)
female_infert5_HT_GM_mr_res <- mr(female_infert5_HT_GM_harm)
female_infert5_HT_GM_mr_het <- mr_heterogeneity(female_infert5_HT_GM_harm)
female_infert5_HT_GM_mr_plei <- mr_pleiotropy_test(female_infert5_HT_GM_harm)
female_infert5_HT_GM_mr_sing <- mr_singlesnp(female_infert5_HT_GM_harm)
male_infert_HT_GM_mr_res <- mr(male_infert_HT_GM_harm)
male_infert_HT_GM_mr_het <- mr_heterogeneity(male_infert_HT_GM_harm)
male_infert_HT_GM_mr_plei <- mr_pleiotropy_test(male_infert_HT_GM_harm)
male_infert_HT_GM_mr_sing <- mr_singlesnp(male_infert_HT_GM_harm)



###PLOTS###
PG_test_scat <- mr_scatter_plot(PG_test_mr_res, PG_test_harm)
PG_test_for <- mr_forest_plot(PG_test_mr_sing)
PG_fsh_scat <- mr_scatter_plot(PG_fsh_mr_res, PG_fsh_harm)
PG_fsh_for <- mr_forest_plot(PG_fsh_mr_sing)
PG_lh_scat <- mr_scatter_plot(PG_lh_mr_res, PG_lh_harm)
PG_lh_for <- mr_forest_plot(PG_lh_mr_sing)
PG_prog_scat <- mr_scatter_plot(PG_prog_mr_res, PG_prog_harm)
PG_prog_for <- mr_forest_plot(PG_prog_mr_sing)
PG_oest_scat <- mr_scatter_plot(PG_oest_mr_res, PG_oest_harm)
PG_oest_for <- mr_forest_plot(PG_oest_mr_sing)

HT_test_scat <- mr_scatter_plot(HT_test_mr_res, HT_test_harm)
HT_test_for <- mr_forest_plot(HT_test_mr_sing)
HT_fsh_scat <- mr_scatter_plot(HT_fsh_mr_res, HT_fsh_harm)
HT_fsh_for <- mr_forest_plot(HT_fsh_mr_sing)
HT_lh_scat <- mr_scatter_plot(HT_lh_mr_res, HT_lh_harm)
HT_lh_for <- mr_forest_plot(HT_lh_mr_sing)
HT_prog_scat <- mr_scatter_plot(HT_prog_mr_res, HT_prog_harm)
HT_prog_for <- mr_forest_plot(HT_prog_mr_sing)
HT_oest_scat <- mr_scatter_plot(HT_oest_mr_res, HT_oest_harm)
HT_oest_for <- mr_forest_plot(HT_oest_mr_sing)

pdf("HT_test_scatter.pdf")
HT_test_scat
dev.off()

LR_test_scat <- mr_scatter_plot(LR_test_mr_res, LR_test_harm)
LR_test_for <- mr_forest_plot(LR_test_mr_sing)
LR_fsh_scat <- mr_scatter_plot(LR_fsh_mr_res, LR_fsh_harm)
LR_fsh_for <- mr_forest_plot(LR_fsh_mr_sing)
LR_lh_scat <- mr_scatter_plot(LR_lh_mr_res, LR_lh_harm)
LR_lh_for <- mr_forest_plot(LR_lh_mr_sing)
LR_prog_scat <- mr_scatter_plot(LR_prog_mr_res, LR_prog_harm)
LR_prog_for <- mr_forest_plot(LR_prog_mr_sing)
LR_oest_scat <- mr_scatter_plot(LR_oest_mr_res, LR_oest_harm)
LR_oest_for <- mr_forest_plot(LR_oest_mr_sing)

HT_GM_test_scat <- mr_scatter_plot(HT_GM_test_mr_res, HT_GM_test_harm)
HT_GM_test_for <- mr_forest_plot(HT_GM_test_mr_sing)
HT_GM_fsh_scat <- mr_scatter_plot(HT_GM_fsh_mr_res, HT_GM_fsh_harm)
HT_GM_fsh_for <- mr_forest_plot(HT_GM_fsh_mr_sing)
HT_GM_lh_scat <- mr_scatter_plot(HT_GM_lh_mr_res, HT_GM_lh_harm)
HT_GM_lh_for <- mr_forest_plot(HT_GM_lh_mr_sing)
HT_GM_prog_scat <- mr_scatter_plot(HT_GM_prog_mr_res, HT_GM_prog_harm)
HT_GM_prog_for <- mr_forest_plot(HT_GM_prog_mr_sing)
HT_GM_oest_scat <- mr_scatter_plot(HT_GM_oest_mr_res, HT_GM_oest_harm)
HT_GM_oest_for <- mr_forest_plot(HT_GM_oest_mr_sing)
pdf("HT_GM_oest_scatter.pdf")
HT_GM_oest_scat
dev.off()
pdf("HT_GM_oest_forest.pdf")
HT_GM_oest_for
dev.off()

test_HT_scat <- mr_scatter_plot(test_HT_mr_res, test_HT_harm)
test_HT_for <- mr_forest_plot(test_HT_mr_sing)
test_PG_scat <- mr_scatter_plot(test_PG_mr_res, test_PG_harm)
test_PG_for <- mr_forest_plot(test_PG_mr_sing)
test_LR_scat <- mr_scatter_plot(test_LR_mr_res, test_LR_harm)
test_LR_for <- mr_forest_plot(test_LR_mr_sing)
test_HT_GM_scat <- mr_scatter_plot(test_HT_GM_mr_res, test_HT_GM_harm)
test_HT_GM_for <- mr_forest_plot(test_HT_GM_mr_sing)

fsh_HT_scat <- mr_scatter_plot(fsh_HT_mr_res, fsh_HT_harm)
fsh_HT_for <- mr_forest_plot(fsh_HT_mr_sing)
fsh_PG_scat <- mr_scatter_plot(fsh_PG_mr_res, fsh_PG_harm)
fsh_PG_for <- mr_forest_plot(fsh_PG_mr_sing)
fsh_LR_scat <- mr_scatter_plot(fsh_LR_mr_res, fsh_LR_harm)
fsh_LR_for <- mr_forest_plot(fsh_LR_mr_sing)
fsh_HT_GM_scat <- mr_scatter_plot(fsh_HT_GM_mr_res, fsh_HT_GM_harm)
fsh_HT_GM_for <- mr_forest_plot(fsh_HT_GM_mr_sing)

lh_HT_scat <- mr_scatter_plot(lh_HT_mr_res, lh_HT_harm)
lh_HT_for <- mr_forest_plot(lh_HT_mr_sing)
lh_PG_scat <- mr_scatter_plot(lh_PG_mr_res, lh_PG_harm)
lh_PG_for <- mr_forest_plot(lh_PG_mr_sing)
lh_LR_scat <- mr_scatter_plot(lh_LR_mr_res, lh_LR_harm)
lh_LR_for <- mr_forest_plot(lh_LR_mr_sing)
lh_HT_GM_scat <- mr_scatter_plot(lh_HT_GM_mr_res, lh_HT_GM_harm)
lh_HT_GM_for <- mr_forest_plot(lh_HT_GM_mr_sing)

oest_HT_scat <- mr_scatter_plot(oest_HT_mr_res, oest_HT_harm)
oest_HT_for <- mr_forest_plot(oest_HT_mr_sing)
oest_PG_scat <- mr_scatter_plot(oest_PG_mr_res, oest_PG_harm)
oest_PG_for <- mr_forest_plot(oest_PG_mr_sing)
oest_LR_scat <- mr_scatter_plot(oest_LR_mr_res, oest_LR_harm)
oest_LR_for <- mr_forest_plot(oest_LR_mr_sing)
oest_HT_GM_scat <- mr_scatter_plot(oest_HT_GM_mr_res, oest_HT_GM_harm)
oest_HT_GM_for <- mr_forest_plot(oest_HT_GM_mr_sing)

###FEMALES

PG_test_F_scat <- mr_scatter_plot(PG_test_F_mr_res, PG_test_F_harm)
PG_test_F_for <- mr_forest_plot(PG_test_F_mr_sing)
PG_fsh_F_scat <- mr_scatter_plot(PG_fsh_F_mr_res, PG_fsh_F_harm)
PG_fsh_F_for <- mr_forest_plot(PG_fsh_F_mr_sing)
pdf("PG_fsh_F_scatter.pdf")
PG_fsh_F_scat
dev.off()

PG_fsh_F_for2 <- forest_plot(PG_fsh_F_mr_res)

pdf("PG_fsh_F_forest.pdf")
PG_fsh_F_for
dev.off()
PG_lh_F_scat <- mr_scatter_plot(PG_lh_F_mr_res, PG_lh_F_harm)
PG_lh_F_for <- mr_forest_plot(PG_lh_F_mr_sing)
pdf("PG_lh_F_scatter.pdf")
PG_lh_F_scat
dev.off()
pdf("PG_lh_F_forest.pdf")
PG_lh_F_for
dev.off()
PG_prog_F_scat <- mr_scatter_plot(PG_prog_F_mr_res, PG_prog_F_harm)
PG_prog_F_for <- mr_forest_plot(PG_prog_F_mr_sing)
PG_oest_F_scat <- mr_scatter_plot(PG_oest_F_mr_res, PG_oest_F_harm)
PG_oest_F_for <- mr_forest_plot(PG_oest_F_mr_sing)

HT_test_F_scat <- mr_scatter_plot(HT_test_F_mr_res, HT_test_F_harm)
HT_test_F_for <- mr_forest_plot(HT_test_F_mr_sing)
HT_fsh_F_scat <- mr_scatter_plot(HT_fsh_F_mr_res, HT_fsh_F_harm)
HT_fsh_F_for <- mr_forest_plot(HT_fsh_F_mr_sing)
HT_lh_F_scat <- mr_scatter_plot(HT_lh_F_mr_res, HT_lh_F_harm)
HT_lh_F_for <- mr_forest_plot(HT_lh_F_mr_sing)
HT_prog_F_scat <- mr_scatter_plot(HT_prog_F_mr_res, HT_prog_F_harm)
HT_prog_F_for <- mr_forest_plot(HT_prog_F_mr_sing)
HT_oest_F_scat <- mr_scatter_plot(HT_oest_F_mr_res, HT_oest_F_harm)
HT_oest_F_for <- mr_forest_plot(HT_oest_F_mr_sing)

pdf("HT_test_F_scatter.pdf")
HT_test_F_scat
dev.off()

LR_test_F_scat <- mr_scatter_plot(LR_test_F_mr_res, LR_test_F_harm)
LR_test_F_for <- mr_forest_plot(LR_test_F_mr_sing)
LR_fsh_F_scat <- mr_scatter_plot(LR_fsh_F_mr_res, LR_fsh_F_harm)
LR_fsh_F_for <- mr_forest_plot(LR_fsh_F_mr_sing)
LR_lh_F_scat <- mr_scatter_plot(LR_lh_F_mr_res, LR_lh_F_harm)
LR_lh_F_for <- mr_forest_plot(LR_lh_F_mr_sing)
LR_prog_F_scat <- mr_scatter_plot(LR_prog_F_mr_res, LR_prog_F_harm)
LR_prog_F_for <- mr_forest_plot(LR_prog_F_mr_sing)
LR_oest_F_scat <- mr_scatter_plot(LR_oest_F_mr_res, LR_oest_F_harm)
LR_oest_F_for <- mr_forest_plot(LR_oest_F_mr_sing)

HT_GM_test_F_scat <- mr_scatter_plot(HT_GM_test_F_mr_res, HT_GM_test_F_harm)
HT_GM_test_F_for <- mr_forest_plot(HT_GM_test_F_mr_sing)
HT_GM_fsh_F_scat <- mr_scatter_plot(HT_GM_fsh_F_mr_res, HT_GM_fsh_F_harm)
HT_GM_fsh_F_for <- mr_forest_plot(HT_GM_fsh_F_mr_sing)
HT_GM_lh_F_scat <- mr_scatter_plot(HT_GM_lh_F_mr_res, HT_GM_lh_F_harm)
HT_GM_lh_F_for <- mr_forest_plot(HT_GM_lh_F_mr_sing)
HT_GM_prog_F_scat <- mr_scatter_plot(HT_GM_prog_F_mr_res, HT_GM_prog_F_harm)
HT_GM_prog_F_for <- mr_forest_plot(HT_GM_prog_F_mr_sing)
HT_GM_oest_F_scat <- mr_scatter_plot(HT_GM_oest_F_mr_res, HT_GM_oest_F_harm)
HT_GM_oest_F_for <- mr_forest_plot(HT_GM_oest_F_mr_sing)

test_HT_F_scat <- mr_scatter_plot(test_HT_F_mr_res, test_HT_F_harm)
test_HT_F_for <- mr_forest_plot(test_HT_F_mr_sing)
test_PG_F_scat <- mr_scatter_plot(test_PG_F_mr_res, test_PG_F_harm)
test_PG_F_for <- mr_forest_plot(test_PG_F_mr_sing)
test_LR_F_scat <- mr_scatter_plot(test_LR_F_mr_res, test_LR_F_harm)
test_LR_F_for <- mr_forest_plot(test_LR_F_mr_sing)
test_HT_GM_F_scat <- mr_scatter_plot(test_HT_GM_F_mr_res, test_HT_GM_F_harm)
test_HT_GM_F_for <- mr_forest_plot(test_HT_GM_F_mr_sing)

fsh_HT_F_scat <- mr_scatter_plot(fsh_HT_F_mr_res, fsh_HT_F_harm)
fsh_HT_F_for <- mr_forest_plot(fsh_HT_F_mr_sing)
fsh_PG_F_scat <- mr_scatter_plot(fsh_PG_F_mr_res, fsh_PG_F_harm)
fsh_PG_F_for <- mr_forest_plot(fsh_PG_F_mr_sing)
fsh_LR_F_scat <- mr_scatter_plot(fsh_LR_F_mr_res, fsh_LR_F_harm)
fsh_LR_F_for <- mr_forest_plot(fsh_LR_F_mr_sing)
fsh_HT_GM_F_scat <- mr_scatter_plot(fsh_HT_GM_F_mr_res, fsh_HT_GM_F_harm)
fsh_HT_GM_F_for <- mr_forest_plot(fsh_HT_GM_F_mr_sing)

lh_HT_F_scat <- mr_scatter_plot(lh_HT_F_mr_res, lh_HT_F_harm)
lh_HT_F_for <- mr_forest_plot(lh_HT_F_mr_sing)
lh_PG_F_scat <- mr_scatter_plot(lh_PG_F_mr_res, lh_PG_F_harm)
lh_PG_F_for <- mr_forest_plot(lh_PG_F_mr_sing)
lh_LR_F_scat <- mr_scatter_plot(lh_LR_F_mr_res, lh_LR_F_harm)
lh_LR_F_for <- mr_forest_plot(lh_LR_F_mr_sing)
lh_HT_GM_F_scat <- mr_scatter_plot(lh_HT_GM_F_mr_res, lh_HT_GM_F_harm)
lh_HT_GM_F_for <- mr_forest_plot(lh_HT_GM_F_mr_sing)

oest_HT_F_scat <- mr_scatter_plot(oest_HT_F_mr_res, oest_HT_F_harm)
oest_HT_F_for <- mr_forest_plot(oest_HT_F_mr_sing)
oest_PG_F_scat <- mr_scatter_plot(oest_PG_F_mr_res, oest_PG_F_harm)
oest_PG_F_for <- mr_forest_plot(oest_PG_F_mr_sing)
oest_LR_F_scat <- mr_scatter_plot(oest_LR_F_mr_res, oest_LR_F_harm)
oest_LR_F_for <- mr_forest_plot(oest_LR_F_mr_sing)
oest_HT_GM_F_scat <- mr_scatter_plot(oest_HT_GM_F_mr_res, oest_HT_GM_F_harm)
oest_HT_GM_F_for <- mr_forest_plot(oest_HT_GM_F_mr_sing)

###MALES
PG_test_M_scat <- mr_scatter_plot(PG_test_M_mr_res, PG_test_M_harm)
PG_test_M_for <- mr_forest_plot(PG_test_M_mr_sing)
PG_fsh_M_scat <- mr_scatter_plot(PG_fsh_M_mr_res, PG_fsh_M_harm)
PG_fsh_M_for <- mr_forest_plot(PG_fsh_M_mr_sing)
PG_lh_M_scat <- mr_scatter_plot(PG_lh_M_mr_res, PG_lh_M_harm)
PG_lh_M_for <- mr_forest_plot(PG_lh_M_mr_sing)
PG_prog_M_scat <- mr_scatter_plot(PG_prog_M_mr_res, PG_prog_M_harm)
PG_prog_M_for <- mr_forest_plot(PG_prog_M_mr_sing)
PG_oest_M_scat <- mr_scatter_plot(PG_oest_M_mr_res, PG_oest_M_harm)
PG_oest_M_for <- mr_forest_plot(PG_oest_M_mr_sing)

HT_test_M_scat <- mr_scatter_plot(HT_test_M_mr_res, HT_test_M_harm)
HT_test_M_for <- mr_forest_plot(HT_test_M_mr_sing)
HT_fsh_M_scat <- mr_scatter_plot(HT_fsh_M_mr_res, HT_fsh_M_harm)
HT_fsh_M_for <- mr_forest_plot(HT_fsh_M_mr_sing)
HT_lh_M_scat <- mr_scatter_plot(HT_lh_M_mr_res, HT_lh_M_harm)
HT_lh_M_for <- mr_forest_plot(HT_lh_M_mr_sing)
HT_prog_M_scat <- mr_scatter_plot(HT_prog_M_mr_res, HT_prog_M_harm)
HT_prog_M_for <- mr_forest_plot(HT_prog_M_mr_sing)
HT_oest_M_scat <- mr_scatter_plot(HT_oest_M_mr_res, HT_oest_M_harm)
HT_oest_M_for <- mr_forest_plot(HT_oest_M_mr_sing)

pdf("HT_test_M_scatter.pdf")
HT_test_M_scat
dev.off()

LR_test_M_scat <- mr_scatter_plot(LR_test_M_mr_res, LR_test_M_harm)
LR_test_M_for <- mr_forest_plot(LR_test_M_mr_sing)
LR_fsh_M_scat <- mr_scatter_plot(LR_fsh_M_mr_res, LR_fsh_M_harm)
LR_fsh_M_for <- mr_forest_plot(LR_fsh_M_mr_sing)
LR_lh_M_scat <- mr_scatter_plot(LR_lh_M_mr_res, LR_lh_M_harm)
LR_lh_M_for <- mr_forest_plot(LR_lh_M_mr_sing)
LR_prog_M_scat <- mr_scatter_plot(LR_prog_M_mr_res, LR_prog_M_harm)
LR_prog_M_for <- mr_forest_plot(LR_prog_M_mr_sing)
LR_oest_M_scat <- mr_scatter_plot(LR_oest_M_mr_res, LR_oest_M_harm)
LR_oest_M_for <- mr_forest_plot(LR_oest_M_mr_sing)

HT_GM_test_M_scat <- mr_scatter_plot(HT_GM_test_M_mr_res, HT_GM_test_M_harm)
HT_GM_test_M_for <- mr_forest_plot(HT_GM_test_M_mr_sing)
HT_GM_fsh_M_scat <- mr_scatter_plot(HT_GM_fsh_M_mr_res, HT_GM_fsh_M_harm)
HT_GM_fsh_M_for <- mr_forest_plot(HT_GM_fsh_M_mr_sing)
HT_GM_lh_M_scat <- mr_scatter_plot(HT_GM_lh_M_mr_res, HT_GM_lh_M_harm)
HT_GM_lh_M_for <- mr_forest_plot(HT_GM_lh_M_mr_sing)
HT_GM_prog_M_scat <- mr_scatter_plot(HT_GM_prog_M_mr_res, HT_GM_prog_M_harm)
HT_GM_prog_M_for <- mr_forest_plot(HT_GM_prog_M_mr_sing)
HT_GM_oest_M_scat <- mr_scatter_plot(HT_GM_oest_M_mr_res, HT_GM_oest_M_harm)
HT_GM_oest_M_for <- mr_forest_plot(HT_GM_oest_M_mr_sing)

test_HT_M_scat <- mr_scatter_plot(test_HT_M_mr_res, test_HT_M_harm)
test_HT_M_for <- mr_forest_plot(test_HT_M_mr_sing)
test_PG_M_scat <- mr_scatter_plot(test_PG_M_mr_res, test_PG_M_harm)
test_PG_M_for <- mr_forest_plot(test_PG_M_mr_sing)
test_LR_M_scat <- mr_scatter_plot(test_LR_M_mr_res, test_LR_M_harm)
test_LR_M_for <- mr_forest_plot(test_LR_M_mr_sing)
test_HT_GM_M_scat <- mr_scatter_plot(test_HT_GM_M_mr_res, test_HT_GM_M_harm)
test_HT_GM_M_for <- mr_forest_plot(test_HT_GM_M_mr_sing)

fsh_HT_M_scat <- mr_scatter_plot(fsh_HT_M_mr_res, fsh_HT_M_harm)
fsh_HT_M_for <- mr_forest_plot(fsh_HT_M_mr_sing)
fsh_PG_M_scat <- mr_scatter_plot(fsh_PG_M_mr_res, fsh_PG_M_harm)
fsh_PG_M_for <- mr_forest_plot(fsh_PG_M_mr_sing)
fsh_LR_M_scat <- mr_scatter_plot(fsh_LR_M_mr_res, fsh_LR_M_harm)
fsh_LR_M_for <- mr_forest_plot(fsh_LR_M_mr_sing)
fsh_HT_GM_M_scat <- mr_scatter_plot(fsh_HT_GM_M_mr_res, fsh_HT_GM_M_harm)
fsh_HT_GM_M_for <- mr_forest_plot(fsh_HT_GM_M_mr_sing)

lh_HT_M_scat <- mr_scatter_plot(lh_HT_M_mr_res, lh_HT_M_harm)
lh_HT_M_for <- mr_forest_plot(lh_HT_M_mr_sing)
lh_PG_M_scat <- mr_scatter_plot(lh_PG_M_mr_res, lh_PG_M_harm)
lh_PG_M_for <- mr_forest_plot(lh_PG_M_mr_sing)
lh_LR_M_scat <- mr_scatter_plot(lh_LR_M_mr_res, lh_LR_M_harm)
lh_LR_M_for <- mr_forest_plot(lh_LR_M_mr_sing)
lh_HT_GM_M_scat <- mr_scatter_plot(lh_HT_GM_M_mr_res, lh_HT_GM_M_harm)
lh_HT_GM_M_for <- mr_forest_plot(lh_HT_GM_M_mr_sing)

oest_HT_M_scat <- mr_scatter_plot(oest_HT_M_mr_res, oest_HT_M_harm)
oest_HT_M_for <- mr_forest_plot(oest_HT_M_mr_sing)
oest_PG_M_scat <- mr_scatter_plot(oest_PG_M_mr_res, oest_PG_M_harm)
oest_PG_M_for <- mr_forest_plot(oest_PG_M_mr_sing)
oest_LR_M_scat <- mr_scatter_plot(oest_LR_M_mr_res, oest_LR_M_harm)
oest_LR_M_for <- mr_forest_plot(oest_LR_M_mr_sing)
oest_HT_GM_M_scat <- mr_scatter_plot(oest_HT_GM_M_mr_res, oest_HT_GM_M_harm)
oest_HT_GM_M_for <- mr_forest_plot(oest_HT_GM_M_mr_sing)


################
###COMPILE TABLES

all_res1 <- do.call(rbind, list(test_HT_mr_res, test_PG_mr_res, test_LR_mr_res, test_HT_GM_mr_res,
                               test_HT_F_mr_res, test_PG_F_mr_res, test_LR_F_mr_res, test_HT_GM_F_mr_res,
                               test_HT_M_mr_res, test_PG_M_mr_res, test_LR_M_mr_res, test_HT_GM_M_mr_res,
                               fsh_HT_mr_res, fsh_PG_mr_res, fsh_LR_mr_res, fsh_HT_GM_mr_res,
                               fsh_HT_F_mr_res, fsh_PG_F_mr_res, fsh_LR_F_mr_res, fsh_HT_GM_F_mr_res,
                               fsh_HT_M_mr_res, fsh_PG_M_mr_res, fsh_LR_M_mr_res, fsh_HT_GM_M_mr_res,
                               lh_HT_mr_res, lh_PG_mr_res, lh_LR_mr_res, lh_HT_GM_mr_res,
                               lh_HT_F_mr_res, lh_PG_F_mr_res, lh_LR_F_mr_res, lh_HT_GM_F_mr_res,
                               oest_HT_mr_res, oest_PG_mr_res, oest_LR_mr_res, oest_HT_GM_mr_res,
                               oest_HT_M_mr_res, oest_PG_M_mr_res, oest_LR_M_mr_res, oest_HT_GM_M_mr_res))

all_res2 <- do.call(rbind, list(HT_test_mr_res, HT_fsh_mr_res, HT_lh_mr_res, HT_prog_mr_res, HT_oest_mr_res,
                               PG_test_mr_res, PG_fsh_mr_res, PG_lh_mr_res, PG_prog_mr_res, PG_oest_mr_res,
                               LR_test_mr_res, LR_fsh_mr_res, LR_lh_mr_res, LR_prog_mr_res, LR_oest_mr_res,
                               HT_GM_test_mr_res, HT_GM_fsh_mr_res, HT_GM_lh_mr_res, HT_GM_prog_mr_res, HT_GM_oest_mr_res,
                               HT_test_F_mr_res, HT_fsh_F_mr_res, HT_lh_F_mr_res, HT_prog_F_mr_res, HT_oest_F_mr_res,
                               PG_test_F_mr_res, PG_fsh_F_mr_res, PG_lh_F_mr_res, PG_prog_F_mr_res, PG_oest_F_mr_res,
                               LR_test_F_mr_res, LR_fsh_F_mr_res, LR_lh_F_mr_res, LR_prog_F_mr_res, LR_oest_F_mr_res,
                               HT_GM_test_F_mr_res, HT_GM_fsh_F_mr_res, HT_GM_lh_F_mr_res, HT_GM_prog_F_mr_res, HT_GM_oest_F_mr_res,
                               HT_test_M_mr_res, HT_fsh_M_mr_res, HT_lh_M_mr_res, HT_oest_M_mr_res,
                               PG_test_M_mr_res, PG_fsh_M_mr_res, PG_lh_M_mr_res, PG_oest_M_mr_res,
                               LR_test_M_mr_res, LR_fsh_M_mr_res, LR_lh_M_mr_res, LR_oest_M_mr_res,
                               HT_GM_test_M_mr_res, HT_GM_fsh_M_mr_res, HT_GM_lh_M_mr_res,  HT_GM_oest_M_mr_res))

all_res_fertility <- do.call(rbind, list(HT_female_infert1_mr_res, PG_female_infert1_mr_res, LR_female_infert1_mr_res, HT_GM_female_infert1_mr_res,
                                         HT_female_infert2_mr_res, PG_female_infert2_mr_res, LR_female_infert2_mr_res, HT_GM_female_infert2_mr_res,
                                         HT_female_infert3_mr_res, PG_female_infert3_mr_res, LR_female_infert3_mr_res, HT_GM_female_infert3_mr_res,
                                         HT_female_infert4_mr_res, PG_female_infert4_mr_res, LR_female_infert4_mr_res, HT_GM_female_infert4_mr_res,
                                         HT_female_infert5_mr_res, PG_female_infert5_mr_res, LR_female_infert5_mr_res, HT_GM_female_infert5_mr_res,
                                         HT_male_infert_mr_res, PG_male_infert_mr_res, LR_male_infert_mr_res, HT_GM_male_infert_mr_res))

all_res1_het <- do.call(rbind, list(test_HT_mr_het, test_PG_mr_het, test_LR_mr_het, test_HT_GM_mr_het,
                                    test_HT_F_mr_het, test_PG_F_mr_het, test_LR_F_mr_het, test_HT_GM_F_mr_het,
                                    test_HT_M_mr_het, test_PG_M_mr_het, test_LR_M_mr_het, test_HT_GM_M_mr_het,
                                    fsh_HT_mr_het, fsh_PG_mr_het, fsh_LR_mr_het, fsh_HT_GM_mr_het,
                                    fsh_HT_F_mr_het, fsh_PG_F_mr_het, fsh_LR_F_mr_het, fsh_HT_GM_F_mr_het,
                                    fsh_HT_M_mr_het, fsh_PG_M_mr_het, fsh_LR_M_mr_het, fsh_HT_GM_M_mr_het,
                                    lh_HT_mr_het, lh_PG_mr_het, lh_LR_mr_het, lh_HT_GM_mr_het,
                                    lh_HT_F_mr_het, lh_PG_F_mr_het, lh_LR_F_mr_het, lh_HT_GM_F_mr_het,
                                    oest_HT_mr_het, oest_PG_mr_het, oest_LR_mr_het, oest_HT_GM_mr_het,
                                    oest_HT_M_mr_het, oest_PG_M_mr_het, oest_LR_M_mr_het, oest_HT_GM_M_mr_het))

all_res1_plei <- do.call(rbind, list(test_HT_mr_plei, test_PG_mr_plei, test_LR_mr_plei, test_HT_GM_mr_plei,
                                  test_HT_F_mr_plei, test_PG_F_mr_plei, test_LR_F_mr_plei, test_HT_GM_F_mr_plei,
                                  test_HT_M_mr_plei, test_PG_M_mr_plei, test_LR_M_mr_plei, test_HT_GM_M_mr_plei,
                                  fsh_HT_mr_plei, fsh_PG_mr_plei, fsh_LR_mr_plei, fsh_HT_GM_mr_plei,
                                  fsh_HT_F_mr_plei, fsh_PG_F_mr_plei, fsh_LR_F_mr_plei, fsh_HT_GM_F_mr_plei,
                                  fsh_HT_M_mr_plei, fsh_PG_M_mr_plei, fsh_LR_M_mr_plei, fsh_HT_GM_M_mr_plei,
                                  lh_HT_mr_plei, lh_PG_mr_plei, lh_LR_mr_plei, lh_HT_GM_mr_plei,
                                  lh_HT_F_mr_plei, lh_PG_F_mr_plei, lh_LR_F_mr_plei, lh_HT_GM_F_mr_plei,
                                  oest_HT_mr_plei, oest_PG_mr_plei, oest_LR_mr_plei, oest_HT_GM_mr_plei,
                                  oest_HT_M_mr_plei, oest_PG_M_mr_plei, oest_LR_M_mr_plei, oest_HT_GM_M_mr_plei))

all_res2_het <- do.call(rbind, list(HT_test_mr_het, HT_fsh_mr_het, HT_lh_mr_het, HT_prog_mr_het, HT_oest_mr_het,
                                 PG_test_mr_het, PG_fsh_mr_het, PG_lh_mr_het, PG_prog_mr_het, PG_oest_mr_het,
                                 LR_test_mr_het, LR_fsh_mr_het, LR_lh_mr_het, LR_prog_mr_het, LR_oest_mr_het,
                                 HT_GM_test_mr_het, HT_GM_fsh_mr_het, HT_GM_lh_mr_het, HT_GM_prog_mr_het, HT_GM_oest_mr_het,
                                 HT_test_F_mr_het, HT_fsh_F_mr_het, HT_lh_F_mr_het, HT_prog_F_mr_het, HT_oest_F_mr_het,
                                 PG_test_F_mr_het, PG_fsh_F_mr_het, PG_lh_F_mr_het, PG_prog_F_mr_het, PG_oest_F_mr_het,
                                 LR_test_F_mr_het, LR_fsh_F_mr_het, LR_lh_F_mr_het, LR_prog_F_mr_het, LR_oest_F_mr_het,
                                 HT_GM_test_F_mr_het, HT_GM_fsh_F_mr_het, HT_GM_lh_F_mr_het, HT_GM_prog_F_mr_het, HT_GM_oest_F_mr_het,
                                 HT_test_M_mr_het, HT_fsh_M_mr_het, HT_lh_M_mr_het, HT_oest_M_mr_het,
                                 PG_test_M_mr_het, PG_fsh_M_mr_het, PG_lh_M_mr_het, PG_oest_M_mr_het,
                                 LR_test_M_mr_het, LR_fsh_M_mr_het, LR_lh_M_mr_het, LR_oest_M_mr_het,
                                 HT_GM_test_M_mr_het, HT_GM_fsh_M_mr_het, HT_GM_lh_M_mr_het,  HT_GM_oest_M_mr_het))

all_res2_plei <- do.call(rbind, list(HT_test_mr_plei, HT_fsh_mr_plei, HT_lh_mr_plei, HT_prog_mr_plei, HT_oest_mr_plei,
                                     PG_test_mr_plei, PG_fsh_mr_plei, PG_lh_mr_plei, PG_prog_mr_plei, PG_oest_mr_plei,
                                     LR_test_mr_plei, LR_fsh_mr_plei, LR_lh_mr_plei, LR_prog_mr_plei, LR_oest_mr_plei,
                                     HT_GM_test_mr_plei, HT_GM_fsh_mr_plei, HT_GM_lh_mr_plei, HT_GM_prog_mr_plei, HT_GM_oest_mr_plei,
                                     HT_test_F_mr_plei, HT_fsh_F_mr_plei, HT_lh_F_mr_plei, HT_prog_F_mr_plei, HT_oest_F_mr_plei,
                                     PG_test_F_mr_plei, PG_fsh_F_mr_plei, PG_lh_F_mr_plei, PG_prog_F_mr_plei, PG_oest_F_mr_plei,
                                     LR_test_F_mr_plei, LR_fsh_F_mr_plei, LR_lh_F_mr_plei, LR_prog_F_mr_plei, LR_oest_F_mr_plei,
                                     HT_GM_test_F_mr_plei, HT_GM_fsh_F_mr_plei, HT_GM_lh_F_mr_plei, HT_GM_prog_F_mr_plei, HT_GM_oest_F_mr_plei,
                                     HT_test_M_mr_plei, HT_fsh_M_mr_plei, HT_lh_M_mr_plei, HT_oest_M_mr_plei,
                                     PG_test_M_mr_plei, PG_fsh_M_mr_plei, PG_lh_M_mr_plei, PG_oest_M_mr_plei,
                                     LR_test_M_mr_plei, LR_fsh_M_mr_plei, LR_lh_M_mr_plei, LR_oest_M_mr_plei,
                                     HT_GM_test_M_mr_plei, HT_GM_fsh_M_mr_plei, HT_GM_lh_M_mr_plei,  HT_GM_oest_M_mr_plei))

write.csv(all_res1, "MR_brain_out_0624_table.csv", quote=FALSE, row.names=FALSE)
write.csv(all_res1_het, "MR_het_brain_out_0624_table.csv", quote=FALSE, row.names=FALSE)
write.csv(all_res1_plei, "MR_plei_brain_out_0624_table.csv", quote=FALSE, row.names=FALSE)

write.csv(all_res2, "MR_brain_exp_0624_table.csv", quote=FALSE, row.names=FALSE)
write.csv(all_res2_het, "MR_het_brain_exp_0624_table.csv", quote=FALSE, row.names=FALSE)
write.csv(all_res2_plei, "MR_plei_brain_exp_0624_table.csv", quote=FALSE, row.names=FALSE)

write.csv(all_res_fertility, "MR_brain_fertility.csv", quote=FALSE, row.names=FALSE)

system("dx upload MR_brain_out_0624_table.csv --path data/")
system("dx upload MR_het_brain_out_0624_table.csv --path data/")
system("dx upload MR_plei_brain_out_0624_table.csv --path data/")

system("dx upload MR_brain_exp_0624_table.csv --path data/")
system("dx upload MR_het_brain_exp_0624_table.csv --path data/")
system("dx upload MR_plei_brain_exp_0624_table.csv --path data/")

system("dx upload MR_brain_fertility.csv --path data/")

#####################
###POST PROCESSING###
#####################

system("dx download data/MR_brain_out_0624_table.csv")
system("dx download data/MR_het_brain_out_0624_table.csv")
system("dx download data/MR_plei_brain_out_0624_table.csv")

system("dx download data/MR_brain_exp_0624_table.csv")
system("dx download data/MR_het_brain_exp_0624_table.csv")
system("dx download data/MR_plei_brain_exp_0624_table.csv")

system("dx download data/MR_brain_fertility.csv")


all_res1 <- fread("MR_brain_out_0624_table.csv", data.table=FALSE)
all_res2 <- fread("MR_brain_exp_0624_table.csv", data.table=FALSE)
all_res_fertility <- fread("MR_brain_fertility.csv", data.table=FALSE)

all_res <- do.call(rbind, list(all_res1, all_res2, all_res_fertility))
all_res$test <- paste(all_res$outcome, all_res$exposure, sep="_")
p <- 0.05/(length(unique(all_res$test)))

all_res$id.exposure <- NULL
all_res$id.outcome <-NULL
write.csv(all_res, "MR_all_results_0624.csv", quote=FALSE, row.names=FALSE)
system("dx upload MR_all_results_0624.csv --path data/")

all_res_sig <- subset(all_res, pval < p)
all_res_sig_nom <- subset(all_res, pval < 0.05 )
all_res_sig_nom <- subset(all_res_sig_nom, nsnp >4)
all_res_sig_001 <- subset(all_res_sig_nom, pval < 0.01)

all_res_sig_nom_order <- all_res_sig_nom[order(all_res_sig_nom$pval),]

n_p=length(unique(all_res$test))
all_res$p.fdr <- p.adjust(all_res$pval, method="fdr")
all_res_sig_fdr <- subset(all_res, p.fdr < 0.05)

#write.csv(all_res_sig, "MR_all_sig.csv", quote=FALSE, row.names=FALSE)
#system("dx upload MR_all_sig.csv --path data/")

pdf("test_LR_all_scat.pdf")
test_LR_scat
dev.off()
pdf("test_HT_F_scat.pdf")
test_HT_F_scat
dev.off()
pdf("test_PG_F_scat.pdf")
test_PG_F_scat
dev.off()
pdf("test_LR_F_scat.pdf")
test_LR_F_scat
dev.off()


sig_relations = subset(all_res_sig, nsnp > 4)
sig_relations$test <- paste(sig_relations$outcome, sig_relations$exposure, sep="_")
sig_all_relations <- subset(all_res, test %in% sig_relations$test)
df = sig_all_relations
df$sig <- df$pval <p
cols <- c("gray40", "hotpink")[match(df$sig, c("FALSE", "TRUE"))]
z <- qnorm((1 + 0.95) / 2, 0, 1)

pdf("all_sig_MR_forest.pdf", width = 12, height = 10)
forest(
  df$b,
  sei = df$se,
  slab = sprintf("    %s", df$method),
  xlab = "MR Estimate of causal effect",
  annotate = FALSE,
  ilab = data.frame(
    sprintf("%.2f", df$b),
    sprintf("(±%.2f)", df$se),
    sprintf("%.2e", df$pval)),
  ilab.xpos = c(1, 1.15, 1.5),
  pch = 16,
#  atransf = exp,
#  at = log(c(0.5, 1, 2, 4, 8, 16)),
  rows = c(1,2,3,4,5,8,9,10,11,12,15,16,17,18,19),
  xlim = c(-1, 2),
  ylim = c(0, 23),
  shade = TRUE,
#  header=TRUE,
  col=cols)

# plot disease labels
par(font=2)
text(-1, c(6,13,20), pos=4, c("HT-GM and Oestradiol in sex-combined", "PG and FSH in females", "PG and LH in females"))
text(0.85, 20, pos=4, "Effect Estimate ± SE")
text(1.39, 20, pos=4, "P-value")

dev.off()

system("dx upload *.pdf --path plots/")
