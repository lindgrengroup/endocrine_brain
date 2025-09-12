###############
###LIBRARIES###
###############

library(data.table)


##########
###DATA###
##########

HT_sig_pattern = "^assoc.resid_HT.volume.*\\.sig.txt$"
HT_sig_files <- list.files("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/coloc/gtex_results", pattern=HT_sig_pattern)
HT_sig <- lapply(HT_sig_files, function(x){
  y <- paste("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/coloc/gtex_results/", x, sep="")
  z <- fread(y, data.table=FALSE)
  cell_type <- sub(".*\\.", "", sub("\\.coloc.sig.txt.*", "", x))
  if(dim(z)[1] > 0) 
  {z$tissue <- cell_type}
  return(z)
})
HT_F_sig_files <- list.files("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/coloc/gtex_results", pattern=HT_sig_pattern)
HT_F_sig <- lapply(HT_F_sig_files, function(x){
  y <- paste("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/coloc/gtex_results/", x, sep="")
  z <- fread(y, data.table=FALSE)
  cell_type <- sub(".*\\.", "", sub("\\.coloc.sig.txt.*", "", x))
  if(dim(z)[1] > 0) 
  {z$tissue <- cell_type}
  return(z)
})
HT_M_sig_files <- list.files("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/coloc/gtex_results", pattern=HT_sig_pattern)
HT_M_sig <- lapply(HT_M_sig_files, function(x){
  y <- paste("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/coloc/gtex_results/", x, sep="")
  z <- fread(y, data.table=FALSE)
  cell_type <- sub(".*\\.", "", sub("\\.coloc.sig.txt.*", "", x))
  if(dim(z)[1] > 0) 
  {z$tissue <- cell_type}
  return(z)
})

PG_sig_pattern = "^assoc.resid_PG.volume.*\\.sig.txt$"
PG_sig_files <- list.files("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/coloc/gtex_results", pattern=PG_sig_pattern)
PG_sig <- lapply(PG_sig_files, function(x){
  y <- paste("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/coloc/gtex_results/", x, sep="")
  z <- fread(y, data.table=FALSE)
  cell_type <- sub(".*\\.", "", sub("\\.coloc.sig.txt.*", "", x))
  if(dim(z)[1] > 0) 
  {z$tissue <- cell_type}
  return(z)
})
PG_F_sig_files <- list.files("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/coloc/gtex_results", pattern=PG_sig_pattern)
PG_F_sig <- lapply(PG_F_sig_files, function(x){
  y <- paste("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/coloc/gtex_results/", x, sep="")
  z <- fread(y, data.table=FALSE)
  cell_type <- sub(".*\\.", "", sub("\\.coloc.sig.txt.*", "", x))
  if(dim(z)[1] > 0) 
  {z$tissue <- cell_type}
  return(z)
})
PG_M_sig_files <- list.files("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/coloc/gtex_results", pattern=PG_sig_pattern)
PG_M_sig <- lapply(PG_M_sig_files, function(x){
  y <- paste("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/coloc/gtex_results/", x, sep="")
  z <- fread(y, data.table=FALSE)
  cell_type <- sub(".*\\.", "", sub("\\.coloc.sig.txt.*", "", x))
  if(dim(z)[1] > 0) 
  {z$tissue <- cell_type}
  return(z)
})

LR_sig_pattern = "^assoc.resid_LR.volume.*\\.sig.txt$"
LR_sig_files <- list.files("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/coloc/gtex_results", pattern=LR_sig_pattern)
LR_sig <- lapply(LR_sig_files, function(x){
  y <- paste("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/coloc/gtex_results/", x, sep="")
  z <- fread(y, data.table=FALSE)
  cell_type <- sub(".*\\.", "", sub("\\.coloc.sig.txt.*", "", x))
  if(dim(z)[1] > 0) 
  {z$tissue <- cell_type}
  return(z)
})
LR_F_sig_files <- list.files("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/coloc/gtex_results", pattern=LR_sig_pattern)
LR_F_sig <- lapply(LR_F_sig_files, function(x){
  y <- paste("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/coloc/gtex_results/", x, sep="")
  z <- fread(y, data.table=FALSE)
  cell_type <- sub(".*\\.", "", sub("\\.coloc.sig.txt.*", "", x))
  if(dim(z)[1] > 0) 
  {z$tissue <- cell_type}
  return(z)
})
LR_M_sig_files <- list.files("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/coloc/gtex_results", pattern=LR_sig_pattern)
LR_M_sig <- lapply(LR_M_sig_files, function(x){
  y <- paste("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/coloc/gtex_results/", x, sep="")
  z <- fread(y, data.table=FALSE)
  cell_type <- sub(".*\\.", "", sub("\\.coloc.sig.txt.*", "", x))
  if(dim(z)[1] > 0) 
  {z$tissue <- cell_type}
  return(z)
})

HT_GM_sig_pattern = "^assoc.resid_HT-SumGM-threshold=0.*\\.sig.txt$"
HT_GM_sig_files <- list.files("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/coloc/gtex_results", pattern=HT_GM_sig_pattern)
HT_GM_sig <- lapply(HT_GM_sig_files, function(x){
  y <- paste("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/coloc/gtex_results/", x, sep="")
  z <- fread(y, data.table=FALSE)
  cell_type <- sub(".*\\.", "", sub("\\.coloc.sig.txt.*", "", x))
  if(dim(z)[1] > 0) 
  {z$tissue <- cell_type}
  return(z)
})
HT_GM_F_sig_files <- list.files("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/coloc/gtex_results", pattern=HT_GM_sig_pattern)
HT_GM_F_sig <- lapply(HT_GM_F_sig_files, function(x){
  y <- paste("/mnt/project/data/brain_all/genotype_process/female_only/regenie_step2/filtered/coloc/gtex_results/", x, sep="")
  z <- fread(y, data.table=FALSE)
  cell_type <- sub(".*\\.", "", sub("\\.coloc.sig.txt.*", "", x))
  if(dim(z)[1] > 0) 
  {z$tissue <- cell_type}
  return(z)
})
HT_GM_M_sig_files <- list.files("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/coloc/gtex_results", pattern=HT_GM_sig_pattern)
HT_GM_M_sig <- lapply(HT_GM_M_sig_files, function(x){
  y <- paste("/mnt/project/data/brain_all/genotype_process/male_only/regenie_step2/filtered/coloc/gtex_results/", x, sep="")
  z <- fread(y, data.table=FALSE)
  cell_type <- sub(".*\\.", "", sub("\\.coloc.sig.txt.*", "", x))
  if(dim(z)[1] > 0) 
  {z$tissue <- cell_type}
  return(z)
})

dir.create("all")
dir.create("female_only")
dir.create("male_only")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/coloc/genes/assoc.* -o all/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/coloc/genes/assoc.* -o female_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/coloc/genes/assoc.* -o male_only/")

HT_all_genes <- fread("all/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", data.table=FALSE)
HT_female_genes <- fread("female_only/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", data.table=FALSE)
HT_male_genes <- fread("male_only/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", data.table=FALSE)

PG_all_genes <- fread("all/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", data.table=FALSE)
PG_female_genes <- fread("female_only/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", data.table=FALSE)
PG_male_genes <- fread("male_only/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", data.table=FALSE)

LR_all_genes <- fread("all/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", data.table=FALSE)
LR_female_genes <- fread("female_only/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", data.table=FALSE)
LR_male_genes <- fread("male_only/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm_TSS_1Mb_genes.csv", data.table = FALSE)

HT_GM_all_genes <- fread("all/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm_TSS_1Mb_genes.csv", data.table=FALSE)
HT_GM_female_genes <- fread("female_only/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm_TSS_1Mb_genes.csv", data.table=FALSE)
HT_GM_male_genes <- fread("male_only/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm_TSS_1Mb_genes.csv", data.table=FALSE)


##############
###ANALYSIS###
##############

HT_sig_df <- data.frame(do.call(rbind, HT_sig))
HT_F_sig_df <- data.frame(do.call(rbind, HT_F_sig))
HT_M_sig_df <- data.frame(do.call(rbind, HT_M_sig))

PG_sig_df <- data.frame(do.call(rbind, PG_sig))
PG_F_sig_df <- data.frame(do.call(rbind, PG_F_sig))
PG_M_sig_df <- data.frame(do.call(rbind, PG_M_sig))

LR_sig_df <- data.frame(do.call(rbind, LR_sig))
LR_F_sig_df <- data.frame(do.call(rbind, LR_F_sig))
LR_M_sig_df <- data.frame(do.call(rbind, LR_M_sig))

HT_GM_sig_df <- data.frame(do.call(rbind, HT_GM_sig))
HT_GM_F_sig_df <- data.frame(do.call(rbind, HT_GM_F_sig))
HT_GM_M_sig_df <- data.frame(do.call(rbind, HT_GM_M_sig))

HT_sig_df$pheno <- "HT"
HT_F_sig_df$pheno <- "HT_F"
HT_M_sig_df$pheno <- "HT_M"

PG_sig_df$pheno <- "PG"
PG_F_sig_df$pheno <- "PG_F"
PG_M_sig_df$pheno <- "PG_M"

LR_sig_df$pheno <- "LR"
LR_F_sig_df$pheno <- "LR_F"
LR_M_sig_df$pheno <- "LR_M"

HT_GM_sig_df$pheno <- "HT_GM"
HT_GM_F_sig_df$pheno <- "HT_GM_F"
HT_GM_M_sig_df$pheno <- "HT_GM_M"

HT_sig_df_merge <- merge(HT_sig_df, HT_all_genes, by.x="i", by.y="Gene stable ID")
HT_F_sig_df_merge <- merge(HT_F_sig_df, HT_female_genes, by.x="i", by.y="Gene stable ID")
HT_M_sig_df_merge <- merge(HT_M_sig_df, HT_male_genes, by.x="i", by.y = "Gene stable ID")

PG_sig_df_merge <- merge(PG_sig_df, PG_all_genes, by.x="i", by.y="Gene stable ID")
PG_F_sig_df_merge <- merge(PG_F_sig_df, PG_female_genes, by.x="i", by.y="Gene stable ID")
PG_M_sig_df_merge <- merge(PG_M_sig_df, PG_male_genes, by.x="i", by.y="Gene stable ID")

LR_sig_df_merge <- merge(LR_sig_df, LR_all_genes, by.x="i", by.y="Gene stable ID")
LR_F_sig_df_merge <- merge(LR_F_sig_df, LR_female_genes, by.x="i", by.y="Gene stable ID")
LR_M_sig_df_merge <- merge(LR_M_sig_df, LR_male_genes, by.x="i", by.y="Gene stable ID")

HT_GM_sig_df_merge <- merge(HT_GM_sig_df, HT_GM_all_genes, by.x="i", by.y="Gene stable ID")
HT_GM_F_sig_df_merge <- merge(HT_GM_F_sig_df, HT_GM_female_genes, by.x="i", by.y="Gene stable ID")
HT_GM_M_sig_df_merge <- merge(HT_GM_M_sig_df, HT_GM_male_genes, by.x="i", by.y="Gene stable ID")

all_sig <- data.frame(do.call(rbind, list(HT_sig_df, HT_F_sig_df, HT_M_sig_df,
                                          PG_sig_df, PG_F_sig_df, PG_M_sig_df,
                                          LR_sig_df, LR_F_sig_df, LR_M_sig_df,
                                          HT_GM_sig_df, HT_GM_F_sig_df, HT_GM_M_sig_df)))

all_sig_merge <- data.frame(do.call(rbind, list(HT_sig_df_merge, HT_F_sig_df_merge, HT_M_sig_df_merge,
                                          PG_sig_df_merge, PG_F_sig_df_merge, PG_M_sig_df_merge,
                                          LR_sig_df_merge, LR_F_sig_df_merge, LR_M_sig_df_merge,
                                          HT_GM_sig_df_merge, HT_GM_F_sig_df_merge, HT_GM_M_sig_df_merge)))
all_sig_order <- all_sig[order(all_sig$tissue),]
all_sig_order_merge <- all_sig_merge[order(all_sig_merge$tissue),]
write.csv(all_sig_order, "all_sig_gtex_coloc_results.csv", quote=FALSE, row.names=FALSE)
write.csv(all_sig_order_merge, "all_sig_gtex_coloc_results_merge.csv", quote=FALSE, row.names=FALSE)
system("dx upload all_sig_gtex_coloc_results.csv --path data/brain_all/genotype_process/all/regenie_step2/filtered/coloc/")
system("dx upload all_sig_gtex_coloc_results_merge.csv --path data/brain_all/genotype_process/all/regenie_step2/filtered/coloc/")
