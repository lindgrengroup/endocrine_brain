###############
###LIBRARIES###
###############

install.packages("docopt")
if(!require("remotes"))
  install.packages("remotes") # if necessary
install.packages("R.utils")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("dplyr")
install.packages("tidyr")
install.packages("data.table")
#devtools::install_github("boxiangliu/locuscomparer")
BiocManager::install("snpStats")
#install_github("chr1swallace/coloc@main",build_vignettes=TRUE)


library(docopt)
library(remotes)
library(R.utils)
library(dplyr)
library(tidyr)
library(data.table)
library(devtools)
library(locuscomparer)
library(snpStats)
library(coloc)

###############
###FUNCTIONS###
###############




################
###DATA INPUT###
################


doc <- "Usage:
    colocalisation_parse.R --i <input> --g <genes> --t <tissue> --o <output> --s <significant>
Arguments:
    --i Path to input full GWAS sumstats of results to test [default: NULL]
    --g Path to file of all genes within distance of top GWAS hits [default: NULL]
    --t Path to tissue gtex files [default: NULL]
    --o Path to output table of results [default: NULL]
    --s Path to output significant results table [default: NULL]
"

opts <- docopt(doc)

##########
###DATA###
##########

gwas_full <- fread(opts$i, data.table=FALSE)
tissue <- fread(opts$t, data.table=FALSE)
genes <- fread(opts$g, data.table=FALSE)

##############
###ANALYSIS###
##############

tissue$gene_id<-gsub("\\..*","",tissue$gene_id)

tissue_subset <-subset(tissue, gene_id %in% genes$`Gene stable ID`)

tissue_subset_ext <- tissue_subset %>% 
  separate(variant_id, c("CHROM", "GENPOS", "REF", "ALT", "BUILD"))
tissue_subset_ext$GENPOS <- as.numeric(tissue_subset_ext$GENPOS)
#tissue_subset_ext$CHROM <- replace(tissue_subset_ext$CHROM, tissue_subset_ext$CHROM=="chrX", "chr23")

tissue_subset_ext$alleles <- apply(cbind(tissue_subset_ext$REF, tissue_subset_ext$ALT), 1, function(x) paste(sort(x), collapse="_"))
tissue_subset_ext$chrpos <- paste(tissue_subset_ext$CHROM, tissue_subset_ext$GENPOS, sep=":")
tissue_subset_ext$MarkerName <- paste(tissue_subset_ext$chrpos, tissue_subset_ext$alleles, sep=":")
tissue_subset_ext$N <- tissue_subset_ext$ma_count/(2*tissue_subset_ext$maf)
tissue_subset_ext$varbeta <- (tissue_subset_ext$slope_se)^2

#PG_vol$chrpos <- paste(PG_vol$CHROM, PG_vol$GENPOS, sep=":")
gwas_full$alleles <- apply(cbind(gwas_full$ALLELE0, gwas_full$ALLELE1), 1, function(x) paste(sort(x), collapse="_"))
gwas_full$MarkerName <- paste(gwas_full$hg38pos, gwas_full$alleles, sep=":")
gwas_full$varbeta <- (gwas_full$SE)^2
gwas_full$MAF <- NA
gwas_full$MAF[which(gwas_full$A1FREQ < 0.5)] <- gwas_full$A1FREQ[which(gwas_full$A1FREQ < 0.5)]
gwas_full$MAF[which(gwas_full$A1FREQ >= 0.5)] <- 1-gwas_full$A1FREQ[which(gwas_full$A1FREQ >= 0.5)]
gwas_full$pval <- 10^(-gwas_full$LOG10P)
gwas_full$EXTRA <- NULL

gwas_full_subset <- subset(gwas_full, MarkerName %in% tissue_subset_ext$MarkerName)
tissue_subset_ext_subset <- subset(tissue_subset_ext, MarkerName %in% gwas_full_subset$MarkerName)

colnames(gwas_full_subset) <- paste(colnames(gwas_full_subset), "gwas", sep="_")
colnames(tissue_subset_ext_subset) <- paste(colnames(tissue_subset_ext_subset), "gtex", sep="_")

gwas_gtex <- merge(tissue_subset_ext_subset, gwas_full_subset, by.x="MarkerName_gtex", by.y="MarkerName_gwas", all=TRUE)

gwas_gtex <- gwas_gtex[complete.cases(gwas_gtex),]

#Align alleles

gwas_gtex$beta_align_gwas <- NA
gwas_gtex$beta_align_gwas[gwas_gtex$ALLELE1_gwas==gwas_gtex$ALT_gtex] <- gwas_gtex$BETA_gwas[gwas_gtex$ALLELE1_gwas==gwas_gtex$ALT_gtex]

gwas_gtex$beta_align_gwas[gwas_gtex$ALLELE1_gwas==gwas_gtex$REF_gtex] <- gwas_gtex$BETA_gwas[gwas_gtex$ALLELE1_gwas==gwas_gtex$REF_gtex]*(-1)

gwas_gtex_genes <- unique(gwas_gtex$gene_id_gtex)

gwas_gtex_results <- lapply(gwas_gtex_genes, function(i){
  data_subset <- subset(gwas_gtex, gene_id_gtex==i)
  coloc_df <- coloc.abf(dataset1 = list(snp=data_subset$MarkerName_gtex, pvalues=data_subset$pval_gwas, beta=data_subset$beta_align_gwas, varbeta=data_subset$varbeta_gwas, MAF=data_subset$MAF_gwas, N=data_subset$N_gwas, type="quant"),
                        dataset2 = list(snp=data_subset$MarkerName_gtex, pvalues=data_subset$pval_nominal_gtex, beta=data_subset$slope_gtex, varbeta=data_subset$varbeta_gtex, MAF=data_subset$maf_gtex, N=data_subset$N_gtex, sdY=1, type="quant"), 
                        p1 = 1e-04, p2 = 1e-04, p12 = 1e-06)
  coloc_df_res <- as.list(coloc_df$summary)
  coloc_df_res$i <- i
  return(coloc_df_res)
})
gwas_gtex_results_table <- do.call(rbind, gwas_gtex_results)
gwas_gtex_results_df <- as.data.frame(gwas_gtex_results_table)
gwas_gtex_results_df <- unnest(gwas_gtex_results_df, cols = c(nsnps, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf, i))

write.table(gwas_gtex_results_df, opts$o, quote=FALSE, row.names = FALSE)

gwas_gtex_results_df$ratio <- gwas_gtex_results_df$PP.H4.abf/gwas_gtex_results_df$PP.H3.abf
gwas_gtex_coloc_sig <- subset(gwas_gtex_results_df, ratio >5 & PP.H4.abf >0.5)

write.table(gwas_gtex_coloc_sig, opts$s, quote=FALSE, row.names = FALSE)


