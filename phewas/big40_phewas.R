#############
###LIBRARY###
#############

library(data.table)
library(ggplot2)
install.packages('pals')
library(pals)

##########
###DATA###
##########

big40 <- fread("/mnt/project/data/2024-02-19_BIG40_filtered.txt", data.table=FALSE) #a1 is ref
HT_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.withX.txt", data.table=FALSE)
PG_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
OB_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.filtered.merged.txt", data.table=FALSE)
HT_GM_cojo <- fread("/mnt/project/data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.merged.txt", data.table=FALSE)

idps <- fread("/mnt/project/data/BIG40-IDPs_v4_IDPs.csv", data.table = FALSE)

##############
###ANLAYSIS###
##############

all_cojo <- rbind(HT_cojo, PG_cojo, OB_cojo, HT_GM_cojo)

colnames(big40) <- c("IDP", "chr", "rsid", "pos", "a1", "a2", "beta", "se", "pval(-log10)")

#Check allele alignment
big_40_check <- merge(big40, all_cojo[,c("SNP", "Chr", "bp", "refA")], by.x="rsid", by.y="SNP", suffixes=c(".40", "all"))
big_40_check$matchallele <- ifelse(big_40_check$a1 == big_40_check$refA, "yes", "no")
big_40_check$matchbp <- ifelse(big_40_check$pos==big_40_check$bp, "yes", "no")

big40 <- merge(big40, idps[,c("Pheno", "IDP short name", "Category name")], by.x="IDP", by.y="Pheno")

big40_HT <- subset(big40, rsid %in% HT_cojo$SNP)
big40_HT_order <- big40_HT[order(big40_HT$`pval(-log10)`, decreasing = TRUE),]
big40_HT_order_sig <- subset(big40_HT_order, `pval(-log10)` >-log10(cut_off))
big40_HT_order_sig$pvalue <- 10^(-big40_HT_order_sig$`pval(-log10)`)
big40_HT_order_sig_paper <- big40_HT_order_sig
big40_HT_order_sig_paper$`pval(-log10)` <- NULL

big40_PG <- subset(big40, rsid %in% PG_cojo$SNP)
big40_PG_order <- big40_PG[order(big40_PG$`pval(-log10)`, decreasing = TRUE),]
big40_PG_order_sig <- subset(big40_PG_order, `pval(-log10)` >-log10(cut_off))
big40_PG_order_sig$pvalue <- 10^(-big40_PG_order_sig$`pval(-log10)`)
big40_PG_order_sig_paper <- big40_PG_order_sig
big40_PG_order_sig_paper$`pval(-log10)` <- NULL

big40_OB <- subset(big40, rsid %in% OB_cojo$SNP)
big40_OB_order <- big40_OB[order(big40_OB$`pval(-log10)`, decreasing = TRUE),]
big40_OB_order_sig <- subset(big40_OB_order, `pval(-log10)` >-log10(cut_off))
big40_OB_order_sig$pvalue <- 10^(-big40_OB_order_sig$`pval(-log10)`)
big40_OB_order_sig_paper <- big40_OB_order_sig
big40_OB_order_sig_paper$`pval(-log10)` <- NULL

big40_HT_GM <- subset(big40, rsid %in% HT_GM_cojo$SNP)
big40_HT_GM_order <- big40_HT_GM[order(big40_HT_GM$`pval(-log10)`, decreasing = TRUE),]
big40_HT_GM_order_sig <- subset(big40_HT_GM_order, `pval(-log10)` >-log10(cut_off))
big40_HT_GM_order_sig$pvalue <- 10^(-big40_HT_GM_order_sig$`pval(-log10)`)
big40_HT_GM_order_sig_paper <- big40_HT_GM_order_sig
big40_HT_GM_order_sig_paper$`pval(-log10)` <- NULL

cut_off <- 0.05/(length(unique(big40$IDP))*length(all_cojo$SNP))

sig40 <- subset(big40, `pval(-log10)` > -log10(cut_off))
unique_sig_idps <- unique(sig40$`IDP short name`)

sig40_vol_intense <- subset(sig40, `Category name` =="regional and tissue volume" | `Category name`=="regional and tissue intensity")
unique_sig_tiss_vol <- unique(sig40_vol_intense$`IDP short name`)

##############
###PLOTTING###
##############

big40_temp <- big40
big40_temp$IDP <- factor(big40_temp$IDP, levels=unique(big40_temp$IDP[order(big40_temp$`Category name`)]))

p <- ggplot(big40_temp, aes(x=IDP, y=`pval(-log10)`, color=`Category name`)) +
  geom_point() +
  scale_color_manual(values= c("gold", "purple", "orange", "lightblue", "red", "seagreen", 
          "lightpink1", "blue", "lightsalmon", "violet", "gold", "maroon", "palegreen",
          "orangered4", "yellowgreen", "sienna1", "olivedrab")) +
  scale_x_discrete(limits=levels(big40_temp$IDP), labels=function(x) big40_temp$`Category name`[match(x, big40_temp$IDP)], expand=c(.04,0)) +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) +
  geom_hline(yintercept=-log10(cut_off), color="grey57")



############
###OUTPUT###
############

write.csv(big40_HT_order, "big40_HT_update_200224.csv", quote=FALSE, row.names=FALSE)
write.csv(big40_PG_order, "big40_PG_update_200224.csv", quote=FALSE, row.names=FALSE)
write.csv(big40_OB_order, "big40_OB_update_200224.csv", quote=FALSE, row.names=FALSE)
write.csv(big40_HT_GM_order, "big40_HT_GM_updaye_200224.csv", quote=FALSE, row.names=FALSE)

write.csv(big40_HT_order_sig_paper, "big40_HT_update_200224_sig.csv", quote=FALSE, row.names=FALSE)
write.csv(big40_PG_order_sig_paper, "big40_PG_update_200224_sig.csv", quote=FALSE, row.names=FALSE)
write.csv(big40_OB_order_sig_paper, "big40_OB_update_200224_sig.csv", quote=FALSE, row.names=FALSE)
write.csv(big40_HT_GM_order_sig_paper, "big40_HT_GM_update_200224_sig.csv", quote=FALSE, row.names=FALSE)

pdf("big40_update_200224_manhattan.pdf",width=8, height=4.8)
p
dev.off()
png("big40_update_200224_manhattan.png",width=1600, height=960)
p
dev.off()

ggsave(file="big40_update_200224_manhattan.svg", plot=p, width=8, height=4.8)


system("dx upload *.csv --path data/")
system("dx upload *.pdf --path plots/")


