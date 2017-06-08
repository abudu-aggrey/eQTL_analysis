# eQTL_analysis
eQTL analysis performed with TwinsUK data

run_eQTL_analysis.sh: contains main script to run the analysis


get_expression_sample_id.R: pull out ID's of genotyped individuals that have expression data avialable


relatedness.sh: script to create your relatedness matrix


merge_exp_data.R: transforms expression data and merges with .sample file 


GEMMA_snp_file.R: creates SNP annotation file for GEMMA


sort_eQTL_results.R: sorts the results by p-value


run_coloc.R: perform colocalisation analysis with results from eQTL analysis and GWAS summary stats
