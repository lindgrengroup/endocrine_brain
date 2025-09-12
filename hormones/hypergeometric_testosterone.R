###############
###LIBRARIES###
###############

library(data.table)
install.packages("Rgb")
library(Rgb)

##########
###DATA###
##########

system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX.txt")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt")

HT_vol <- fread("assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX.txt", data.table=FALSE)
PG_vol <- fread("assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
LR_vol <- fread("assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
HT_GM <- fread("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)

system("dx download data/hormone_fertility_update_231125/0624_remeta/lead_snps/withukb/all_lead_snp_sumstats_Testosterone_sex_comb_EUR_with_rsids.txt")
test_lead <- fread("all_lead_snp_sumstats_Testosterone_sex_comb_EUR_with_rsids.txt", data.table=FALSE)

##############
###ANALYSIS###
##############

all_results <- rbind(HT_vol, PG_vol, LR_vol, HT_GM)

all_results$testosterone <-FALSE
all_results$testosterone[which(all_results$SNP %in% c("rs11941568", "rs2184968", "rs76895963", "rs7808966", "rs10220706"))] <- TRUE

dat <- data.frame(
  "testosterone" = c(sum(all_results$testosterone==TRUE), length(test_lead$RSID)), 
  "non_test" = c(sum(all_results$testosterone==FALSE), 1000000-length(test_lead$RSID)),
  row.names=c("expected", "observed"),
  stringsAsFactors = FALSE
)

x <- mosaicplot(dat, color = TRUE)
png("test_enrichment_fisher.png")
x
dev.off()

fisher_exact <- fisher.test(dat, alternative="greater")


