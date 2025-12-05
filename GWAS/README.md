# Genome-wide Association Studies (GWAS) and associated analyses

Scripts in this folder to perform GWAS and downstream associated analyses:

1.  regenie_step1.sh - run step 1 of GWAS in REGENIE.
2.  regenie_step2.sh - run step 2 of GWAS in REGENIE.
3.  xchrom_info_merge.R - merge the info score information per SNP available in imputed files with the output of X chromosome REGENIE run.
4.  merge_regenie.sh - merge the per chromomse output of REGENIE.
5.  gwas_plot.R - parsable script to plot manhattan plot and qqplot from filtered GWAS output.
6.  gwas_miami.R - script to plot miami plots to display sex-specific results for each of the four neuroendocrine phenotypes.
7.  hw_imputed.sh - calculate the hardy-weinberg equilibrium across the imputed genotype files
8.  gwas_filter.R - parsable script to filter REGENIE GWAS results for info score, Hardy-Weinberg equilibrium, minor allele frequency and multiallelic alleles.
9.  cut_off_adjust_reviews.R -  calculate different lists of significant hits based on different mutliple testing method specifications
10.  meta_manhattan_qq_upset_paper_review.R - plot to generate figure one of paper including manhattan plot and aligned upset plot and colour coded qq plots
11.  run_gwas_filter_plot_Xchr.sh - filter and plot the results of REGENIE
12.  gwas_cojo_format.R - format the output of REGENIE for input to cojo
13.  format_X_chr_plink.R - change plink file X encoding from "X" to "23"
14.  run_cojo.sh - run cojo across each chromosome and merge
15.  run_cojo_Xchr.sh - convert imputed bgen to plink and run COJO conditional SNP selection across each of the GWASn for files once converged with X chromosome
16.  cojo_cond_filter.R -  filter based on the output of the cojo conditional analysis
17.  results_all_mega_table_create.R - create a meta table that compiles the significantly associated genetic variants found across the four neuroendocrine phenotypes
18.  cojo_cond_pre_format.R - format files for cojo conditional analysis
19.  run_cojo_across_gwas.sh - condition each all GWAS on the SNPS discovered in the GWAS of each of the neuroendocrine phenotypes
20.  mtag_format.R -  parsable script to prepare filtered summary statistics for MTAG meta-analysis
21.  mtag_post_process.R - parsable script to format the output of MTAG meta-analysis
22.  run_mtag_review.sh - run MTAG to meta analyse the four different the summary statistics of the four different neuroendocrine phenotypes. Parses mtag_format.R and mtag_post_process.R
23.  gwas_cojo_format_mtag.R - preapre the output of MTAG meta-analysis for input to cojo
24.  run_cojo_mtag_review.sh -  apply cojo conditional independent SNP selection on the output of MTAG meta-analysed summary statistics
25.  sex_diff_forest_plot.R - generate a forest plot of sex-specific betas of the SNPS that have been identified as being sexually dimorphic.
26.  heritability_plot.R - plot the sex-specific heritability.
27.  gwas_sumstats_liftover_input_prep.R - interactive script to prepare summary statistics for liftover.
28.  liftover_prep_parse.R - prepare summary statistics files for liftover.
29.  liftover_post_format.R - format the output of liftover.
30.  run_liftover_prep.sh - file to prepare GWAS output files and run liftover from hg19 to hg38. Parses liftover prep_parse.R and liftover_post_format.R
31.  run_ldsc.sh - run LDSC to estimate heritability. Then run genetic correlation against hormone and infertility summary statitsics from Venkatesh et al.
32.  gwas_ldsc_format.R - parsable script to prepare summary statistics for LDSC
33.  ldsc_calc_pz.R - calculate the z-score for heritability estimates
34.  dosage_prep.R - prepare files for input to qctools to extract hard calls of dosage
35.  convert_to_dosage_withX.sh -  convert genotyping BGEN files to hard dosage files
36.  hard_dosage_format.R - parsable script that formats the output of qctools creation of hard genotype calls
37.  dosage_format_withX.R -  parsable script that formats the output of qctools creation of hard genotype calls including X chromosome
