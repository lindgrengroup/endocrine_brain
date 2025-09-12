###############
###LIBRARIES###
###############

install.packages("readxl")
library(readxl)
library(ggplot2)
library(stringr)
library(dplyr)
library(data.table)

###############
###FUNCTIONS###
###############

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
##########
###DATA###
##########

cellect <- read_excel_allsheets("/mnt/project/cellect/CELLECT_HCL_GTEx.xlsx")
cellect_extra <- fread("/mnt/project/cellect/prioritization_ldsc.csv")
#cellect_extra_magma <- fread("/mnt/project/cellect/prioritization_magma.csv")
  
##############
###ANALYSIS###
##############

cellect_res <- cellect[which(names(cellect) != "Columns")]
cellect_res_tab <- do.call(rbind, cellect_res)
cellect_res_tab$annotation_short <- str_replace_all(cellect_res_tab$annotation, "[[:punct:]]", "")
rownames(cellect_res_tab) <- NULL
cellect_res_tab_small <- cellect_res_tab[,c("gwas", "specificity_id", "annotation", "beta", "beta_se", "pvalue", "annotation_short")]

annotation_dict <- distinct(cellect_res_tab[,c("annotation", "annotation_short")])

cellect_extra$gwas <- c("infertility", "HT.volume.female", "HT.volume.male")[match(cellect_extra$gwas, c("Inf_Eur_1", "female", "male"))]
cellect_extra$specificity_id <- "Hypothalamus.siletti.hypomap"
cellect_extra$annotation_short <- cellect_extra$annotation
cellect_extra$annotation <- NULL
cellect_extra <- left_join(cellect_extra, annotation_dict, by="annotation_short")
cellect_extra_HT_male <- subset(cellect_extra, gwas=="HT.volume.male")

cellect_all <- rbind(cellect_res_tab_small, cellect_extra_HT_male)
cellect_all_fdr <- cellect_all %>% 
    group_by(gwas) %>%
    mutate(n=length(pvalue)) %>% 
    mutate(fdr_p=p.adjust(pvalue, method="fdr")) %>% 
    ungroup()
cellect_all_fdr <- as.data.frame(cellect_all_fdr)
cellect_all_fdr$fdr_stat <-  cellect_all_fdr$fdr_p <0.05

cellect_all_fdr$log10p <- -log10(cellect_all_fdr$pvalue)
cellect_all_fdr$group <- unlist(lapply(cellect_all_fdr$gwas, function (x) str_split(x, "\\-|\\.")[[1]][3]))
cellect_all_fdr$pheno <- unlist(lapply(cellect_all_fdr$gwas, function (x) {
  tmp <- str_split(x, "\\-|\\.")[[1]][c(1,2)]
  return(paste(tmp, collapse = "_"))}))
cellect_all_fdr$log10fdr <- -log10(cellect_all_fdr$fdr_p)
cellect_all_fdr$annotation_brief <- gsub("^.*\\.\\.", "", cellect_all_fdr$annotation)

cellect_all_fdr_short <- cellect_all_fdr[,c("pheno", "group","annotation_brief", "beta", "beta_se", "n", "pvalue", "fdr_p", "fdr_stat")]
write.csv(cellect_all_fdr_short, "cellect_all_fdr_short.csv", quote=FALSE, row.names=FALSE)
system("dx upload cellect_all_fdr_short.csv --path data/")

p <- ggplot(cellect_all_fdr, aes(y=annotation, x=log10fdr, color=group)) + 
  geom_point(position = position_dodge(width=0.5), size=2) + facet_wrap(~pheno) + 
  labs(x='-log10(FDR corrected p-value)', y = 'Hypothalmic cell type') +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", color="gray") +
  theme_minimal(base_size=8)

p_1 <- ggplot(cellect_all_fdr, aes(y=annotation_brief, x=log10fdr, color=pheno, shape=group)) + 
  geom_point(position = position_dodge(width=2), size=3) + 
  labs(x='-log10(FDR corrected p-value)', y = 'Hypothalmic cell type') +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", color="gray") +
  theme_minimal(base_size=12) +
  scale_color_manual(values=c("#5bbda3","#9d83d6","#EEA243", "#b51a0e"), limits=c("HT_volume", "PG_volume", "LR_volume", "HT_SumGM"), labels = c("Hypothalamus", "Pituitary\nGland", "Olfactory\nBulb", "Hypothalamus\nGrey Matter")) +
  scale_shape(limits=c("female", "male", "sex_combined"), labels = c("Female", "Male", "Sex\ncombined"))  +
  guides(color = guide_legend(override.aes = list(size = 6)), shape= guide_legend(override.aes = list(size = 6))) 

#cellect_all_fdr$group <- factor(cellect_all_fdr$group, levels = c("male", "female", "sex_combined"), 
#                                  labels = c("Male", "Female", "Sex-combined"))

group.labs <- c("Male", "Female", "Sex-combined")
names(group.labs) <- c("male", "female", "sex_combined")

p_3 <- ggplot(cellect_all_fdr, aes(y=annotation_brief, x=log10fdr, color=pheno, shape=group)) + 
  geom_point(position = position_dodge(width=.5), size=3) + 
  labs(x='-log10(FDR corrected p-value)', y = 'Hypothalmic cell type') +
  geom_vline(xintercept = -log10(0.05), linetype="dashed", color="gray45") +
  theme_minimal(base_size=12) +
  scale_color_manual(values=c("#5bbda3","#9d83d6","#EEA243", "#b51a0e"), limits=c("HT_volume", "PG_volume", "LR_volume", "HT_SumGM"), labels = c("Hypothalamus", "Pituitary\nGland", "Olfactory\nBulb", "Hypothalamus\nGrey Matter")) +
  scale_shape(limits=c("female", "male", "sex_combined"), labels = c("Female", "Male", "Sex\ncombined"))  +
  guides(color = guide_legend(override.aes = list(size = 6)), shape= guide_legend(override.aes = list(size = 6))) +
  facet_wrap(~factor(group, levels=c("male", "female", "sex_combined", labels=c("Male", "Female", "Sex-combined"))), nrow=1) +
  theme(strip.text.x=element_text(face="bold", size=12))




pdf("cellect_hypothalmic.pdf", width=8, height=6)
p
dev.off()
pdf("cellect_hypothalmic_color_combine_update.pdf", width=8, height=6)
p_1
dev.off()
pdf("cellect_hypothalmic_color_combine_update_facet.pdf", width=15, height=6)
p_3
dev.off()

system("dx upload cellect_hypothalmic_color_combine_update.pdf --path plots/")
system("dx upload cellect_hypothalmic_color_combine_update_facet.pdf --path plots/")

