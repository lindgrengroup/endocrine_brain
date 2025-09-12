###############
###LIBRARIES###
###############

library(data.table)

##########
###DATA###
##########

HT_vol <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/liftover/output/formatted/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.hg38.txt", data.table = FALSE)

##############
###ANALYSIS###
##############

HT_vol$hg38chr <- sub("\\:.*", "", HT_vol$hg38pos)
HT_vol$hg38genpos <- sub('.*:', '', HT_vol$hg38pos)

# ENSG00000010327 location in hg38

HT_vol_exome_res <- subset(HT_vol, hg38chr=="chr3" & hg38genpos>52495338 & hg38genpos<52524495)

write.csv(HT_vol_exome_res, "exome_crossover_HT_ENSG00000010327.csv", quote=FALSE, row.names=FALSE)

system("dx upload exome_crossover_HT_ENSG00000010327.csv --path data/")

