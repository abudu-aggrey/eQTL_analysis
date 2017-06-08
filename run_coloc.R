#generate input files for coloc analysis.

#read in results from eQTL analysis and pull out rsID, beta and se

eQTL_data <- read.table("./TwinsUK/CD207_eQTL_analysis_all.assoc.txt", header=T)
eQTL_data <- eQTL_data[ ,c(2,8,9)]
colnames(eQTL_data)[1] <- "rsID"

#read in GWAS data and subset for SNPs from eQTL analysis
ECZ_GWAS <- read.table("./EAGLE_AD_GWAS/EAGLE_AD_GWAS_results_2015.txt", header=T)
GWAS_data <- merge(ECZ_GWAS, eQTL_data, by="rsID")

#extract snp, beta, se and eaf
GWAS_data <- GWAS_data[ ,c(1, 8,9,6)]
colnames(GWAS_data)[2] <- "beta"
colnames(GWAS_data)[3] <- "se"

#exclude SNPs in expression data with no genotype data
colnames(eQTL_data)[2] <- "exp_beta"
colnames(eQTL_data)[3] <- "exp_se"
eQTL_GWAS_merge <- merge(eQTL_data, GWAS_data, by = "rsID")
eQTL_data <- eQTL_GWAS_merge[ ,1:3]


#Perform coloc
#install.packages("coloc")
library(coloc)

my.res <- coloc.abf(dataset1=list(beta=GWAS_data$beta, varbeta=GWAS_data$se, N=nrow(GWAS_data), snp=GWAS_data$rsID, type="quant"),
            dataset2=list(beta=eQTL_data$exp_beta, varbeta=eQTL_data$exp_se, N=nrow(eQTL_data), snp=eQTL_data$rsID, type="quant"),
            MAF=GWAS_data$eaf)

coloc_results <- my.res$results
coloc_results <- coloc_results[order(coloc_results$SNP.PP.H4,decreasing=T),]
