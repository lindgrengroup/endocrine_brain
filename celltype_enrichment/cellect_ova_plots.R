###############
###LIBRARIES###
###############

library(ggplot2)
library(stringr)
library(dplyr)
library(data.table)

###############
###FUNCTIONS###
###############

##########
###DATA###
##########

system("dx download data/prioritization*")
cellect <- fread("prioritization_ova.csv", data.table=FALSE)
cellect_magma <- fread("prioritization_magma_ova.csv", data.table=FALSE)

  
##############
###ANALYSIS###
##############

# cellect_fdr <- cellect %>% 
#    group_by(gwas) %>%
#    mutate(n=length(pvalue)) %>% 
#    mutate(fdr_p=p.adjust(pvalue, method="fdr")) %>% 
#    ungroup()
#cellect_fdr <- as.data.frame(cellect_fdr)
#cellect_fdr$fdr_stat <-  cellect_fdr$fdr_p <0.05

cellect_magma$group <- unlist(lapply(cellect_magma$gwas, function (x) str_split(x, "\\_|\\.")[[1]][3]))
cellect_magma$pheno <- unlist(lapply(cellect_magma$gwas, function (x) {
  tmp <- str_split(x, "\\_|\\.")[[1]][c(1,2)]
  return(paste(tmp, collapse = "_"))}))
cellect_magma <- subset(cellect_magma, group %in% c("female", "comb"))

cellect_magma_fdr <- cellect_magma %>%
  group_by(gwas) %>%
  mutate(n=length(pvalue)) %>% 
  mutate(fdr_p=p.adjust(pvalue, method="fdr")) %>% 
  ungroup()
cellect_magma_fdr <- as.data.frame(cellect_magma_fdr)
cellect_magma_fdr$fdr_stat <-  cellect_magma_fdr$fdr_p <0.05

cellect_magma_fdr$log10p <- -log10(cellect_magma_fdr$pvalue)
cellect_magma_fdr$log10fdr <- -log10(cellect_magma_fdr$fdr_p)

cellect_magma_fdr_short <- cellect_magma_fdr[,c("pheno", "group","annotation", "beta", "beta_se", "n", "pvalue", "fdr_p", "fdr_stat")]
write.csv(cellect_magma_fdr_short, "ova_cellect_magma_fdr_short.csv", quote=FALSE, row.names=FALSE)
system("dx upload ova_cellect_magma_fdr_short.csv --path data/")

cellect_magma_fdr_sig <- subset(cellect_magma_fdr_short, fdr_stat==TRUE)
write.csv(cellect_magma_fdr_sig, "cellect_magma_fdr_sig.csv", quote=FALSE, row.names=FALSE)
system("dx upload cellect_magma_fdr_sig.csv --path data/")

p <- ggplot(cellect_magma_fdr, aes(y=annotation, x=log10fdr, color=group)) + 
  geom_point(position = position_dodge(width=0.5), size=2) + facet_wrap(~pheno) + 
  labs(x='-log10(FDR corrected p-value)', y = 'Ovarian cell type') +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", color="gray") +
  theme_minimal(base_size=8)

p_1 <- ggplot(cellect_magma_fdr, aes(y=annotation, x=log10fdr, color=pheno, shape=group)) + 
  geom_point(position = position_dodge(width=.5), size=3) + 
  labs(x='-log10(FDR corrected p-value)', y = 'Ovarian cell type') +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", color="gray") +
  theme_minimal(base_size=12) +
  scale_color_manual(values=c("#5bbda3","#9d83d6","#EEA243", "#b51a0e"), limits=c("HT_Vol", "PG_Vol", "LR_Vol", "HT_SumGM"), labels = c("Hypothalamus", "Pituitary\nGland", "Olfactory\nBulb", "Hypothalamus\nGrey Matter")) +
  scale_shape(limits=c("female", "male", "comb"), labels = c("Female", "Male", "Sex\ncombined"))  +
  guides(color = guide_legend(override.aes = list(size = 6)), shape= guide_legend(override.aes = list(size = 6))) 

cellect_magma_fdr$group <- factor(cellect_magma_fdr$group, levels = c("female", "male", "comb"), 
                  labels = c("Female", "Male", "Sex-combined"))

p_2 <- ggplot(cellect_magma_fdr, aes(y=annotation, x=log10fdr, color=pheno, shape=group)) + 
  geom_point(position = position_dodge(width=.5), size=3) + 
  labs(x='-log10(FDR corrected p-value)', y = 'Ovarian cell type') +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", color="gray45") +
  theme_minimal(base_size=12) +
  scale_color_manual(values=c("#5bbda3","#9d83d6","#EEA243", "#b51a0e"), limits=c("HT_Vol", "PG_Vol", "LR_Vol", "HT_SumGM"), labels = c("Hypothalamus", "Pituitary\nGland", "Olfactory\nBulb", "Hypothalamus\nGrey Matter")) +
  scale_shape(limits=c("Female", "Male", "Sex-combined"), labels = c("Female", "Male", "Sex\ncombined"))  +
  guides(color = guide_legend(override.aes = list(size = 6)), shape= guide_legend(override.aes = list(size = 6))) +
  facet_wrap(~group, nrow=1) +
  theme(strip.text.x=element_text(face="bold", size=12))

pdf("cellect_ova.pdf", width=8, height=6)
p
dev.off()
pdf("cellect_ovarian_color_combine_update.pdf", width=8, height=6)
p_1
dev.off()
pdf("cellect_ovarian_color_combine_update_facet.pdf", width=10, height=6)
p_2
dev.off()

system("dx upload cellect_ovarian_color_combine_update.pdf --path plots/")
system("dx upload cellect_ovarian_color_combine_update_facet.pdf --path plots/")
