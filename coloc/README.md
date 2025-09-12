# Colocalisation analysis

Scripts in this folder to peform colocalisation analysis:

1. get_genes.R - file to extract genes within range of genetic variants for use in colocalisation anaylsis
2. coloc_gtex_parse.R - parsable file to be set using bash script to run colocalisation for all GWAS sumstats across all tissues specific eQTLs within the GTEx consortium.
3. run_eqtl_coloc.sh - run coloc_gtex_parse.R across tissues from GTEx data
4. coloc_gtex_process.R - process the output from running colocalisation across GTEx eQTLs and to make one table of results per phenotype.
5. colocalisation_repro.R - R script to perform colocalisation analysis between neuroendocrine phenotypes and reproductive traits taken from Venkatesh SS et al., : Genome-wide analyses identify 21 infertility loci and over 400 reproductive hormone loci across the allele frequency spectrum (2024)
