###############
###LIBRARIES###
###############

library(data.table)

##########
###DATA###
##########

PG_vol <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_vol <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
LR_vol <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM_vol <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

PG_vol_F <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_vol_F <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
LR_vol_F <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM_vol_F <- fread("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

PG_vol_M <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_vol_M <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
LR_vol_M <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM_vol_M <- fread("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)


##############
###ANALYSIS###
##############

PG_vol$chr_format <- paste("chr", PG_vol$CHROM, sep="")
PG_vol$start <- PG_vol$GENPOS -1
PG_vol$end <- PG_vol$GENPOS
PG_vol_bed <- PG_vol[,c("chr_format", "ID", "start", "end")]
write.table(PG_vol_bed, "PG_vol_all_hg19.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

HT_vol$chr_format <- paste("chr", HT_vol$CHROM, sep="")
HT_vol$start <- HT_vol$GENPOS -1
HT_vol$end <- HT_vol$GENPOS
HT_vol_bed <- HT_vol[,c("chr_format", "ID", "start", "end")]
write.table(HT_vol_bed, "HT_vol_all_hg19.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

LR_vol$chr_format <- paste("chr", LR_vol$CHROM, sep="")
LR_vol$start <- LR_vol$GENPOS -1
LR_vol$end <- LR_vol$GENPOS
LR_vol_bed <- LR_vol[,c("chr_format", "ID", "start", "end")]
write.table(LR_vol_bed, "LR_vol_all_hg19.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

HT_GM_vol$chr_format <- paste("chr", HT_GM_vol$CHROM, sep="")
HT_GM_vol$start <- HT_GM_vol$GENPOS -1
HT_GM_vol$end <- HT_GM_vol$GENPOS
HT_GM_vol_bed <- HT_GM_vol[,c("chr_format", "ID", "start", "end")]
write.table(HT_GM_vol_bed, "HT_GM_vol_all_hg19.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

PG_vol_F$chr_format <- paste("chr", PG_vol_F$CHROM, sep="")
PG_vol_F$start <- PG_vol_F$GENPOS -1
PG_vol_F$end <- PG_vol_F$GENPOS
PG_vol_F_bed <- PG_vol_F[,c("chr_format", "ID", "start", "end")]
write.table(PG_vol_F_bed, "PG_vol_female_hg19.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

HT_vol_F$chr_format <- paste("chr", HT_vol_F$CHROM, sep="")
HT_vol_F$start <- HT_vol_F$GENPOS -1
HT_vol_F$end <- HT_vol_F$GENPOS
HT_vol_F_bed <- HT_vol_F[,c("chr_format", "ID", "start", "end")]
write.table(HT_vol_F_bed, "HT_vol_female_hg19.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

LR_vol_F$chr_format <- paste("chr", LR_vol_F$CHROM, sep="")
LR_vol_F$start <- LR_vol_F$GENPOS -1
LR_vol_F$end <- LR_vol_F$GENPOS
LR_vol_F_bed <- LR_vol_F[,c("chr_format", "ID", "start", "end")]
write.table(LR_vol_F_bed, "LR_vol_female_hg19.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

HT_GM_vol_F$chr_format <- paste("chr", HT_GM_vol_F$CHROM, sep="")
HT_GM_vol_F$start <- HT_GM_vol_F$GENPOS -1
HT_GM_vol_F$end <- HT_GM_vol_F$GENPOS
HT_GM_vol_F_bed <- HT_GM_vol_F[,c("chr_format", "ID", "start", "end")]
write.table(HT_GM_vol_F_bed, "HT_GM_vol_female_hg19.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

PG_vol_M$chr_format <- paste("chr", PG_vol_M$CHROM, sep="")
PG_vol_M$start <- PG_vol_M$GENPOS -1
PG_vol_M$end <- PG_vol_M$GENPOS
PG_vol_M_bed <- PG_vol_M[,c("chr_format", "ID", "start", "end")]
write.table(PG_vol_M_bed, "PG_vol_male_hg19.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

HT_vol_M$chr_format <- paste("chr", HT_vol_M$CHROM, sep="")
HT_vol_M$start <- HT_vol_M$GENPOS -1
HT_vol_M$end <- HT_vol_M$GENPOS
HT_vol_M_bed <- HT_vol_M[,c("chr_format", "ID", "start", "end")]
write.table(HT_vol_M_bed, "HT_vol_male_hg19.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

LR_vol_M$chr_format <- paste("chr", LR_vol_M$CHROM, sep="")
LR_vol_M$start <- LR_vol_M$GENPOS -1
LR_vol_M$end <- LR_vol_M$GENPOS
LR_vol_M_bed <- LR_vol_M[,c("chr_format", "ID", "start", "end")]
write.table(LR_vol_M_bed, "LR_vol_male_hg19.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

HT_GM_vol_M$chr_format <- paste("chr", HT_GM_vol_M$CHROM, sep="")
HT_GM_vol_M$start <- HT_GM_vol_M$GENPOS -1
HT_GM_vol_M$end <- HT_GM_vol_M$GENPOS
HT_GM_vol_M_bed <- HT_GM_vol_M[,c("chr_format", "ID", "start", "end")]
write.table(HT_GM_vol_M_bed, "HT_GM_vol_male_hg19.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

system("dx upload *_hg19.bed --path data/brain_all/genotype_process/liftover/input/")
