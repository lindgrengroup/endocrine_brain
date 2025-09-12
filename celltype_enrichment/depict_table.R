###############
###LIBRARIES###
###############

library(data.table)

##########
###DATA###
##########

system("dx download DEPICT_analysis/MRI_depict_0724/*")

HT_all_gene_prior <- fread("HT_volume_comb_geneprioritization.txt", data.table=FALSE)
HT_all_gene_enrich <- fread("HT_volume_comb_genesetenrichment.txt", data.table=FALSE)
HT_all_tissue <- fread("HT_volume_comb_tissueenrichment.txt", data.table=FALSE)

HT_female_gene_prior <- fread("HT_volume_females_geneprioritization.txt", data.table=FALSE)
HT_female_gene_enrich <- fread("HT_volume_females_genesetenrichment.txt", data.table=FALSE)
HT_female_tissue <- fread("HT_volume_females_tissueenrichment.txt", data.table=FALSE)

HT_male_gene_prior <- fread("HT_volume_males_geneprioritization.txt", data.table=FALSE)
HT_male_gene_enrich <- fread("HT_volume_males_genesetenrichment.txt", data.table=FALSE)
HT_male_tissue <- fread("HT_volume_males_tissueenrichment.txt", data.table=FALSE)

PG_all_gene_prior <- fread("PG_volume_comb_geneprioritization.txt", data.table=FALSE)
PG_all_gene_enrich <- fread("PG_volume_comb_genesetenrichment.txt", data.table=FALSE)
PG_all_tissue <- fread("PG_volume_comb_tissueenrichment.txt", data.table=FALSE)

PG_female_gene_prior <- fread("PG_volume_females_geneprioritization.txt", data.table=FALSE)
PG_female_gene_enrich <- fread("PG_volume_females_genesetenrichment.txt", data.table=FALSE)
PG_female_tissue <- fread("PG_volume_females_tissueenrichment.txt", data.table=FALSE)

PG_male_gene_prior <- fread("PG_volume_males_geneprioritization.txt", data.table=FALSE)
PG_male_gene_enrich <- fread("PG_volume_males_genesetenrichment.txt", data.table=FALSE)
PG_male_tissue <- fread("PG_volume_males_tissueenrichment.txt", data.table=FALSE)

LR_all_gene_prior <- fread("LR_volume_comb_geneprioritization.txt", data.table=FALSE)
LR_all_gene_enrich <- fread("LR_volume_comb_genesetenrichment.txt", data.table=FALSE)
LR_all_tissue <- fread("LR_volume_comb_tissueenrichment.txt", data.table=FALSE)

LR_female_gene_prior <- fread("LR_volume_females_geneprioritization.txt", data.table=FALSE)
LR_female_gene_enrich <- fread("LR_volume_females_genesetenrichment.txt", data.table=FALSE)
LR_female_tissue <- fread("LR_volume_females_tissueenrichment.txt", data.table=FALSE)

LR_male_gene_prior <- fread("LR_volume_males_geneprioritization.txt", data.table=FALSE)
LR_male_gene_enrich <- fread("LR_volume_males_genesetenrichment.txt", data.table=FALSE)
LR_male_tissue <- fread("LR_volume_males_tissueenrichment.txt", data.table=FALSE)

HT_GM_all_gene_prior <- fread("HT_SumGM_comb_geneprioritization.txt", data.table=FALSE)
HT_GM_all_gene_enrich <- fread("HT_SumGM_comb_genesetenrichment.txt", data.table=FALSE)
HT_GM_all_tissue <- fread("HT_SumGM_comb_tissueenrichment.txt", data.table=FALSE)

HT_GM_female_gene_prior <- fread("HT_SumGM_females_geneprioritization.txt", data.table=FALSE)
HT_GM_female_gene_enrich <- fread("HT_SumGM_females_genesetenrichment.txt", data.table=FALSE)
HT_GM_female_tissue <- fread("HT_SumGM_females_tissueenrichment.txt", data.table=FALSE)

HT_GM_male_gene_prior <- fread("HT_SumGM_males_geneprioritization.txt", data.table=FALSE)
HT_GM_male_gene_enrich <- fread("HT_SumGM_males_genesetenrichment.txt", data.table=FALSE)
HT_GM_male_tissue <- fread("HT_SumGM_males_tissueenrichment.txt", data.table=FALSE)

##############
###ANALYSIS###
##############

#GENE ENRICHMENT
HT_all_gene_enrich$pheno <- "HT"
HT_all_gene_enrich$group <- "all"
HT_female_gene_enrich$pheno <- "HT"
HT_female_gene_enrich$group <- "female"
HT_male_gene_enrich$pheno <- "HT"
HT_male_gene_enrich$group <- "male"

PG_all_gene_enrich$pheno <- "PG"
PG_all_gene_enrich$group <- "all"
PG_female_gene_enrich$pheno <- "PG"
PG_female_gene_enrich$group <- "female"
PG_male_gene_enrich$pheno <- "PG"
PG_male_gene_enrich$group <- "male"

LR_all_gene_enrich$pheno <- "OB"
LR_all_gene_enrich$group <- "all"
LR_female_gene_enrich$pheno <- "OB"
LR_female_gene_enrich$group <- "female"
LR_male_gene_enrich$pheno <- "OB"
LR_male_gene_enrich$group <- "male"

HT_GM_all_gene_enrich$pheno <- "HT_GM"
HT_GM_all_gene_enrich$group <- "all"
HT_GM_female_gene_enrich$pheno <- "HT_GM"
HT_GM_female_gene_enrich$group <- "female"
HT_GM_male_gene_enrich$pheno <- "HT_GM"
HT_GM_male_gene_enrich$group <- "male"

all_gene_enrich <- rbind(HT_all_gene_enrich, HT_female_gene_enrich, HT_male_gene_enrich,
                         PG_all_gene_enrich, PG_female_gene_enrich, PG_male_gene_enrich,
                         LR_all_gene_enrich, LR_female_gene_enrich, LR_male_gene_enrich,
                         HT_GM_all_gene_enrich, HT_GM_female_gene_enrich, HT_GM_male_gene_enrich)

#TISSUE ENRICHMENT
HT_all_tissue$pheno <- "HT"
HT_all_tissue$group <- "all"
HT_female_tissue$pheno <- "HT"
HT_female_tissue$group <- "female"
HT_male_tissue$pheno <- "HT"
HT_male_tissue$group <- "male"

PG_all_tissue$pheno <- "PG"
PG_all_tissue$group <- "all"
PG_female_tissue$pheno <- "PG"
PG_female_tissue$group <- "female"
PG_male_tissue$pheno <- "PG"
PG_male_tissue$group <- "male"

LR_all_tissue$pheno <- "OB"
LR_all_tissue$group <- "all"
LR_female_tissue$pheno <- "OB"
LR_female_tissue$group <- "female"
LR_male_tissue$pheno <- "OB"
LR_male_tissue$group <- "male"

HT_GM_all_tissue$pheno <- "HT_GM"
HT_GM_all_tissue$group <- "all"
HT_GM_female_tissue$pheno <- "HT_GM"
HT_GM_female_tissue$group <- "female"
HT_GM_male_tissue$pheno <- "HT_GM"
HT_GM_male_tissue$group <- "male"

all_tissue <- rbind(HT_all_tissue, HT_female_tissue, HT_male_tissue,
                    PG_all_tissue, PG_female_tissue, PG_male_tissue,
                    LR_all_tissue, LR_female_tissue, LR_male_tissue,
                    HT_GM_all_tissue, HT_GM_female_tissue, HT_GM_male_tissue)

#GENE PRIORITISATION
HT_all_gene_prior$pheno <- "HT"
HT_all_gene_prior$group <- "all"
HT_female_gene_prior$pheno <- "HT"
HT_female_gene_prior$group <- "female"
HT_male_gene_prior$pheno <- "HT"
HT_male_gene_prior$group <- "male"

PG_all_gene_prior$pheno <- "PG"
PG_all_gene_prior$group <- "all"
PG_female_gene_prior$pheno <- "PG"
PG_female_gene_prior$group <- "female"
PG_male_gene_prior$pheno <- "PG"
PG_male_gene_prior$group <- "male"

LR_all_gene_prior$pheno <- "OB"
LR_all_gene_prior$group <- "all"
LR_female_gene_prior$pheno <- "OB"
LR_female_gene_prior$group <- "female"
LR_male_gene_prior$pheno <- "OB"
LR_male_gene_prior$group <- "male"

HT_GM_all_gene_prior$pheno <- "HT_GM"
HT_GM_all_gene_prior$group <- "all"
HT_GM_female_gene_prior$pheno <- "HT_GM"
HT_GM_female_gene_prior$group <- "female"
HT_GM_male_gene_prior$pheno <- "HT_GM"
HT_GM_male_gene_prior$group <- "male"

all_gene_prior <- rbind(HT_all_gene_prior, HT_female_gene_prior, HT_male_gene_prior,
                        PG_all_gene_prior, PG_female_gene_prior, PG_male_gene_prior,
                        LR_all_gene_prior, LR_female_gene_prior, LR_male_gene_prior,
                        HT_GM_all_gene_prior, HT_GM_female_gene_prior, HT_GM_male_gene_prior)

write.csv(all_gene_enrich, "all_gene_enrichment.csv", row.names=FALSE, quote=FALSE)
write.csv(all_gene_prior, "all_gene_prioritisation.csv", row.names=FALSE, quote=FALSE)
write.csv(all_tissue, "all_tissue_prioritisation.csv", row.names=FALSE, quote=FALSE)

system("dx upload *.csv --path DEPICT_analysis/MRI_depict_0724/")
