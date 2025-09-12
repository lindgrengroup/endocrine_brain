#############
###LIBRARY###
#############

library(data.table)

###############
###FUNCTIONS###
###############

get_genes <- function(snp_list, sumstats){
  gene_results <- lapply(snp_list, function(i){
    lead_chr <- sumstats$hg38chr[which(sumstats$ID==i)]
    lead_pos <- sumstats$hg38_genpos[which(sumstats$ID==i)]
    gene_info_select <- subset(gene_info_subset, `Chromosome/scaffold name`==lead_chr & `Transcription start site (TSS)`>lead_pos-1000000 & `Transcription start site (TSS)`<lead_pos+1000000)
    if(nrow(gene_info_select)>0){
      gene_info_select$SNP <- i
    }
    return(gene_info_select)
  })
  gene_results_table <- do.call(rbind, gene_results)
  return(gene_results_table)
}

##########
###DATA###
##########

gene_info <- fread("/mnt/project/data/resources/ensembl_build38_all_gene_transcript_names_positions.txt", data.table=FALSE)

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


##############
###ANALYSIS###
##############

PG_vol_select <- subset(PG_vol, ID %in% PG_vol_cojo$SNP)
PG_vol$hg38chr <- gsub(".*chr(.+):.*", "\\1", PG_vol$hg38pos)
PG_vol$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", PG_vol$hg38pos)))

HT_vol_select <- subset(HT_vol, ID %in% HT_vol_cojo$SNP)
HT_vol$hg38chr <- gsub(".*chr(.+):.*", "\\1", HT_vol$hg38pos)
HT_vol$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", HT_vol$hg38pos)))

LR_vol_select <- subset(LR_vol, ID %in% LR_vol_cojo$SNP)
LR_vol$hg38chr <- gsub(".*chr(.+):.*", "\\1", LR_vol$hg38pos)
LR_vol$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", LR_vol$hg38pos)))

HT_GM_vol_select <- subset(HT_GM_vol, ID %in% HT_GM_vol_cojo$SNP)
HT_GM_vol$hg38chr <- gsub(".*chr(.+):.*", "\\1", HT_GM_vol$hg38pos)
HT_GM_vol$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", HT_GM_vol$hg38pos)))

PG_vol_select_F <- subset(PG_vol_F, ID %in% PG_vol_cojo_F$SNP)
PG_vol_F$hg38chr <- gsub(".*chr(.+):.*", "\\1", PG_vol_F$hg38pos)
PG_vol_F$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", PG_vol_F$hg38pos)))

HT_vol_select_F <- subset(HT_vol_F, ID %in% HT_vol_cojo_F$SNP)
HT_vol_F$hg38chr <- gsub(".*chr(.+):.*", "\\1", HT_vol_F$hg38pos)
HT_vol_F$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", HT_vol_F$hg38pos)))

LR_vol_select_F <- subset(LR_vol_F, ID %in% LR_vol_cojo_F$SNP)
LR_vol_F$hg38chr <- gsub(".*chr(.+):.*", "\\1", LR_vol_F$hg38pos)
LR_vol_F$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", LR_vol_F$hg38pos)))

HT_GM_vol_select_F <- subset(HT_GM_vol_F, ID %in% HT_GM_vol_cojo_F$SNP)
HT_GM_vol_F$hg38chr <- gsub(".*chr(.+):.*", "\\1", HT_GM_vol_F$hg38pos)
HT_GM_vol_F$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", HT_GM_vol_F$hg38pos)))

PG_vol_select_M <- subset(PG_vol_M, ID %in% PG_vol_cojo_M$SNP)
PG_vol_M$hg38chr <- gsub(".*chr(.+):.*", "\\1", PG_vol_M$hg38pos)
PG_vol_M$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", PG_vol_M$hg38pos)))

HT_vol_select_M <- subset(HT_vol_M, ID %in% HT_vol_cojo_M$SNP)
HT_vol_M$hg38chr <- gsub(".*chr(.+):.*", "\\1", HT_vol_M$hg38pos)
HT_vol_M$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", HT_vol_M$hg38pos)))

LR_vol_select_M <- subset(LR_vol_M, ID %in% LR_vol_cojo_M$SNP)
LR_vol_M$hg38chr <- gsub(".*chr(.+):.*", "\\1", LR_vol_M$hg38pos)
LR_vol_M$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", LR_vol_M$hg38pos)))

HT_GM_vol_select_M <- subset(HT_GM_vol_M, ID %in% HT_GM_vol_cojo_M$SNP)
HT_GM_vol_M$hg38chr <- gsub(".*chr(.+):.*", "\\1", HT_GM_vol_M$hg38pos)
HT_GM_vol_M$hg38_genpos <- as.numeric(gsub("\\D", "", sub(".*:", "", HT_GM_vol_M$hg38pos)))


gene_info_subset <- subset(gene_info, `Chromosome/scaffold name` %in% c("1","2","3","4","5","6","7","8","9","10","11", "12","13","14","15","16","17","18","19","20","21","22","X"))
gene_info_subset <- subset(gene_info_subset, `Ensembl Canonical`==1 & is.na(`Ensembl Canonical`) ==FALSE)

PG_genes <- get_genes(PG_vol_cojo$SNP, PG_vol)
HT_genes <- get_genes(HT_vol_cojo$SNP, HT_vol)
LR_genes <- get_genes(LR_vol_cojo$SNP, LR_vol)
HT_GM_genes <- get_genes(HT_GM_vol_cojo$SNP, HT_GM_vol)

PG_genes_F <- get_genes(PG_vol_cojo_F$SNP, PG_vol_F)
HT_genes_F <- get_genes(HT_vol_cojo_F$SNP, HT_vol_F)
LR_genes_F <- get_genes(LR_vol_cojo_F$SNP, LR_vol_F)
HT_GM_genes_F <- get_genes(HT_GM_vol_cojo_F$SNP, HT_GM_vol_F)

PG_genes_M <- get_genes(PG_vol_cojo_M$SNP, PG_vol_M)
HT_genes_M <- get_genes(HT_vol_cojo_M$SNP, HT_vol_M)
LR_genes_M <- get_genes(LR_vol_cojo_M$SNP, LR_vol_M)
HT_GM_genes_M <- get_genes(HT_GM_vol_cojo_M$SNP, HT_GM_vol_M)

write.csv(PG_genes, "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_genes, "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_genes, "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_genes, "assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm_TSS_1Mb_genes.csv", quote=FALSE, row.names=FALSE)

system("dx upload *genes.csv --path data/brain_all/genotype_process/all/regenie_step2/filtered/coloc/genes/")

write.csv(PG_genes_F, "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_genes_F, "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_genes_F, "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_genes_F, "assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm_TSS_1Mb_genes.csv", quote=FALSE, row.names=FALSE)

system("dx upload *genes.csv --path data/brain_all/genotype_process/female_only/regenie_step2/filtered/coloc/genes/")


write.csv(PG_genes_M, "assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_genes_M, "assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", quote=FALSE, row.names=FALSE)
write.csv(LR_genes_M, "assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", quote=FALSE, row.names=FALSE)
write.csv(HT_GM_genes_M, "assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm_TSS_1Mb_genes.csv", quote=FALSE, row.names=FALSE)

system("dx upload *genes.csv --path data/brain_all/genotype_process/male_only/regenie_step2/filtered/coloc/genes/")

################
###SPARE CODE###
################

genes$keep<-0
for (i in 1:nrow(genes)){
  CHR<-genes$Chromosome.scaffold.name[i]
  TSS<-genes$Transcription.start.site..TSS.[i]
  subset<-sentinels[sentinels$chr==CHR & sentinels$pos>TSS-1000000 & sentinels$pos<TSS+1000000,]
  if (nrow(subset)>0){
    genes$keep[i]<-1
  }
}

genes<-genes[genes$keep==1,]
write.table(genes,"genes_with_TSS_within_1Mb_of_poi_sentinel_SNP.txt",sep="\t",quote=F,row.names=F)