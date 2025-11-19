######################
###COMMON LOCATIONS###
######################

genotype_file_dir="Bulk/Imputation/UKB imputation from genotype"
script_loc="scripts"
software_loc="software"

############################
###ALL BRAIN PHENOS - M&F###
############################

pheno_input="Phenotypes_MRI/regenie_input/phenotype_files/replication_resid_norm_euro_all.phe"
project_out="data/brain_all"
pheno_file="replication_resid_norm_euro_all.phe"
pheno_group="all"
pheno_name="replication_resid_norm_euro_all"
covar_file="replication_resid_norm_euro_all.cov"
covar_input="Phenotypes_MRI/regenie_input/covar_files/replication_resid_norm_euro_all.cov"
snplist_qc="${project_out}/genotype_process/${pheno_group}/genotype_array_snps_qc_pass.snplist"

#############################
###ALL BRAIN PHENOS FEMALE###
#############################
pheno_input="Phenotypes_MRI/regenie_input/phenotype_files/replication_resid_norm_euro_female.phe"
project_out="data/brain_all"
pheno_file="replication_resid_norm_euro_female.phe"
pheno_group="female_only"
pheno_name="replication_resid_norm_euro_female"
covar_file="replication_resid_norm_euro_female.cov"
covar_input="Phenotypes_MRI/regenie_input/covar_files/replication_resid_norm_euro_female.cov"
snplist_qc="${project_out}/genotype_process/${pheno_group}/genotype_array_snps_qc_pass.snplist"

###########################
###ALL BRAIN PHENOS MALE###
###########################
pheno_input="Phenotypes_MRI/regenie_input/phenotype_files/replication_resid_norm_euro_male.phe"
project_out="data/brain_all"
pheno_file="replication_resid_norm_euro_male.phe"
pheno_group="male_only"
pheno_name="replication_resid_norm_euro_male"
covar_file="replication_resid_norm_euro_male.cov"
covar_input="Phenotypes_MRI/regenie_input/covar_files/replication_resid_norm_euro_male.cov"
snplist_qc="${project_out}/genotype_process/${pheno_group}/genotype_array_snps_qc_pass.snplist"

