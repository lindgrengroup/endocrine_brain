# Exome-wide association analysis

Scripts in this folder to oricess results of exome-wide association analyses:

1. exome_config.sh - config file to run exome related scripts
2. prep_exome_input.R - prepare phenotype and covariate files for use in exome-wide association studies.
3. exome_merge.sh - merge the exome output files
4. run_exome_sig.sh - run exome_sig_gene_skat.R and exome_sig_sing_var.R
5. run_exome_plot.sh - run exome_plot_SKAT.R and exome_plot_sing_var.R
6. exome_crossover.R - Check whether there was any common geentic associations discoevred during GWAS within the _STAB1_ locus, the gene discovered during the exom-wide association analysis.
