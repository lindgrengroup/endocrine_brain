###############
###LIBRARIES###
###############
#install.packages("devtools")
#library(devtools)
#install_version("rlang", version = "1.0.5", repos = "http://cran.us.r-project.org")
library(data.table)
library(ggplot2)
#devtools::install_github("NightingaleHealth/ggforestplot")
#library(ggforestplot)
#library(dplyr)
install.packages("remotes")
remotes::install_github("const-ae/ggsignif")
library(ggsignif)
install.packages("ggpubr")
library(ggpubr)
install.packages("readxl")
library("readxl")

##########
###DATA###
##########

dir.create("all")
dir.create("female_only")
dir.create("male_only")

system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX.txt -o all/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt -o female_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt -o male_only/")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o all/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o female_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o male_only/")

system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt -o all/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt -o female_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt -o male_only/")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt -o all/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt -o female_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt -o male_only/")

system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt -o all/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt -o female_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt -o male_only/")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o all/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o female_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o male_only/")

system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt -o all/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.txt -o female_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/cojo_results/merged/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.txt -o male_only/")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o all/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o female_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o male_only/")

HT_all_cojo <- fread("all/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX.txt", data.table=FALSE)
HT_female_cojo <- fread("female_only/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
HT_male_cojo <- fread("male_only/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
HT_all <- fread("all/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_female <- fread("female_only/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_male <- fread("male_only/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

HT_GM_all_cojo <- fread("all/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)
HT_GM_female_cojo <- fread("female_only/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)
HT_GM_male_cojo <- fread("male_only/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)
HT_GM_all <- fread("all/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM_female <- fread("female_only/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM_male <- fread("male_only/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

PG_all_cojo <- fread("all/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
PG_female_cojo <- fread("female_only/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
PG_male_cojo <- fread("male_only/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
PG_all <- fread("all/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
PG_female <- fread("female_only/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
PG_male <- fread("male_only/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

OB_all_cojo <- fread("all/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
OB_female_cojo <- fread("female_only/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)
OB_male_cojo <- fread("male_only/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)
OB_all <- fread("all/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table = FALSE)
OB_female <- fread("female_only/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
OB_male <- fread("male_only/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

system("dx download data/sex_diff/update0724/Sex_diff_july24.xlsx")
sex_diff_results <- read_excel("Sex_diff_july24.xlsx")
write.csv(sex_diff_results, "sex_diff_july24.csv", quote=FALSE, row.names=FALSE)
system("dx upload sex_diff_july24.csv --path data/sex_diff/update0724/")
##############
###ANALYSIS###
##############

HT_combined_snps <- unique(c(HT_all_cojo$SNP, HT_female_cojo$SNP, HT_male_cojo$SNP))
HT_all_select_snps <- subset(HT_all, ID %in% HT_combined_snps)
HT_female_select_snps <- subset(HT_female, ID %in% HT_combined_snps)
HT_male_select_snps <- subset(HT_male, ID %in% HT_combined_snps)

HT_all_select_snps$group <- "all"
HT_female_select_snps$group <- "female"
HT_male_select_snps$group <- "male"
HT_select_snps <- rbind(HT_all_select_snps, HT_female_select_snps, HT_male_select_snps)

PG_combined_snps <- unique(c(PG_all_cojo$SNP, PG_female_cojo$SNP, PG_male_cojo$SNP))
PG_all_select_snps <- subset(PG_all, ID %in% PG_combined_snps)
PG_female_select_snps <- subset(PG_female, ID %in% PG_combined_snps)
PG_male_select_snps <- subset(PG_male, ID %in% PG_combined_snps)

PG_all_select_snps$group <- "all"
PG_female_select_snps$group <- "female"
PG_male_select_snps$group <- "male"
PG_select_snps <- rbind(PG_all_select_snps, PG_female_select_snps, PG_male_select_snps)

OB_combined_snps <- unique(c(OB_all_cojo$SNP, OB_female_cojo$SNP, OB_male_cojo$SNP))
OB_all_select_snps <- subset(OB_all, ID %in% OB_combined_snps)
OB_female_select_snps <- subset(OB_female, ID %in% OB_combined_snps)
OB_male_select_snps <- subset(OB_male, ID %in% OB_combined_snps)

OB_all_select_snps$group <- "all"
OB_female_select_snps$group <- "female"
OB_male_select_snps$group <- "male"
OB_select_snps <- rbind(OB_all_select_snps, OB_female_select_snps, OB_male_select_snps)

HT_GM_combined_snps <- unique(c(HT_GM_all_cojo$SNP, HT_GM_female_cojo$SNP, HT_GM_male_cojo$SNP))
HT_GM_all_select_snps <- subset(HT_GM_all, ID %in% HT_GM_combined_snps)
HT_GM_female_select_snps <- subset(HT_GM_female, ID %in% HT_GM_combined_snps)
HT_GM_male_select_snps <- subset(HT_GM_male, ID %in% HT_GM_combined_snps)

HT_GM_all_select_snps$group <- "all"
HT_GM_female_select_snps$group <- "female"
HT_GM_male_select_snps$group <- "male"
HT_GM_select_snps <- rbind(HT_GM_all_select_snps, HT_GM_female_select_snps, HT_GM_male_select_snps)


HT_select_snps$group <- as.factor(HT_select_snps$group)
PG_select_snps$group <- as.factor(PG_select_snps$group)
OB_select_snps$group <- as.factor(OB_select_snps$group)
HT_GM_select_snps$group <- as.factor(HT_GM_select_snps$group)

sex_diff_results$pheno <- paste(sex_diff_results$GLAND, sex_diff_results$PHENOTYPE, sep="_")

HT_select_snps_sig <- subset(HT_select_snps, ID %in% sex_diff_results$SNP)
PG_select_snps_sig <- subset(PG_select_snps, ID %in% sex_diff_results$SNP)
OB_select_snps_sig <- subset(OB_select_snps, ID %in% sex_diff_results$SNP)
HT_GM_select_snps_sig <- subset(HT_GM_select_snps, ID %in% sex_diff_results$SNP)

HT_select_snps_sig$pheno <- "HT"
PG_select_snps_sig$pheno <- "PG"
OB_select_snps_sig$pheno <- "OB"
HT_GM_select_snps_sig$pheno <- "HT_GM"

select_snps_sig <- rbind(HT_select_snps_sig, PG_select_snps_sig, OB_select_snps_sig, HT_GM_select_snps_sig)
select_snps_MF <- subset(select_snps_sig, group !="all")
select_snps_MF$pvalue <- 10^(-select_snps_MF$LOG10P)

write.csv(select_snps_MF, "sex_diff_sig_snps.csv", quote=FALSE, row.names=FALSE)
system("dx upload sex_diff_sig_snps.csv --path data/sex_diff/")

##############
###PLOTTING###
##############

select_snps_MF$ID <- factor(select_snps_MF$ID, levels=c("rs144968764", "rs7749444", "rs186990314", "rs6544040"))
p <- ggplot(data=select_snps_MF, aes(y=ID, x=BETA, xmin=(BETA-SE), xmax=(BETA+SE), col=as.factor(pheno), fill=pheno, pch=group)) +
  geom_point(position = position_dodge(width=0.5), size=4) +
  geom_errorbarh(height=.1, position=position_dodge(width=0.5)) +
  labs(x='Effect Size', y = 'Genetic Variant') +
  scale_colour_manual(values=c("#5bbda3","#9d83d6","#EEA243", "#b51a0e"), limits=c("HT", "PG", "OB", "HT_GM"), labels = c("Hypothalamus", "Pituitary\nGland", "Olfactory\nBulb", "Hypothalamus\nGrey Matter"))+
  theme_minimal(base_size=16)
p <- p + geom_vline(xintercept=0, colour="grey40")

############
###OUTPUT###
############

pdf("sex_diff_forest_plot_paper_update0724.pdf", width=8, height=6)
p
dev.off()





