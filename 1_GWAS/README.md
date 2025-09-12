# Genome wide association studies and associated processing

Scripts in this folder to peform GWAS and downstream processing:

1. regenie1.sh - run regenie step 1
2. regenie2.sh - run regenie step 2
3. xchrom_info_merge.R - parsable file to add info score from original genotyping files to Xchrom REGENIE output
4. merge_regenie.sh merge by chromosome regenie output and add X chrom info
5. hw_imputed.sh - run plink to calculayte hardy-weinberg equilibrium for imputed data for usie in gwas_filter.R
6. gwas_filter.R - parsable file to be used in bash scripting for filtering output of REGENIE GWAS for MAF, info etc.
7. gwas_plot.R - parasble file to make manhattan and qq plot
8. gwas_ldsc_format.R - parsable script for formatting output of REGENIE GWAS for input into LDSC
9. run_ldsc.sh - bash script to run LDSC
10. ldsc_calc_pz.R - calculate the p-value and z-score for the LDSC intercept
11. gwas_cojo_format.R - parsable script for formatting output of REGENIE GWAS for input into conditional joint analysis using GCTA-COJO for selection of independent loci
12. gwas_miami.R - generate miami plots from sex-specific GWAS
13. sex_diff_forest_plot.R - make plots showing sexual dimorphism in genetic effect
14. heritability_plot.R - test for sex differences in heritability and plot
15. cojo_cond_pre_format.R - extract SNP list from GWAS for use in conditional join analysis using GCTA-COJO acorss snps discovered in other phenotypes
16. cojo_cond_filter.R - filter the output from cojo conditional joint analysis using GCTA-COJO to genome-wide signficiant and those not in LD.
17. run_cojo_across_gwas.sh - run gwas_cojo_format.R, cojo_cond_pre_format, GCTA-COJO and cojo_cond_filter.R to perform cojo conditional formatting agasint other phenotypes
18. run_gwas_filter_plot_Xchr.sh - run gwas_filter.R file and gwas_plot.R
19. meta_manhattan_qq_upset.R - make the combined manhattanupset plot from figure 1. Also make associated qq plots.
20. results_all_mega_table_create.R - compile GWAS results into a table for manuscript
21. format_X_chr_plink.R - change X to 23 in ChrX plink file
22. convert_to_dosage_withX.sh - run conversion of bgen to dosage data
23. dosage_prep.R - parsable file to prepare input for qctools to convert bgen to dosage from genotype and phenotype files.
24. dosage_format_withX.R - parsable file to format the output from qctool conversion of bgen to dosage of genotype.
25. hard_dosage_format.R - parsable file to process output from hard dosage conversions. 
26. gwas_sumstats_liftover_input_prep.R - prepare files ready to perfom liftover between genome builds
27. liftover_prep_parse.R - parsable file to prepare files ready to perform liftover between genome builds
28. run_liftover_prep.sh - script to run liftover_prep_parse.R
29. liftover_post_format.R - parsable file to format the output of perfomring genome build liftover

