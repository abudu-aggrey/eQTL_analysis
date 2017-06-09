#generate input files for coloc analysis.

#read in results from eQTL analysis 

eQTL_data <- read.table("./TwinsUK/CD207_eQTL_analysis_all.assoc.txt", header=T)

#define region around leand SNP
SNP_pos <- eQTL_data$ps[which(eQTL_data$rs == "rs112111458")]
upp_pos <- SNP_pos + 100000
low_pos <- SNP_pos - 100000


#read in GWAS data and subset for SNPs from eQTL analysis
ECZ_GWAS <- read.table("./EAGLE_AD_GWAS/EAGLE_AD_GWAS_results_2015.txt", header=T)
GWAS_data <- ECZ_GWAS[which(ECZ_GWAS$rsID %in% eQTL_data$rs), ]
locus_snps <- GWAS_data[which(GWAS_data$position > low_pos & GWAS_data$position < upp_pos), ]
locus_eQTL <- eQTL_data[which(eQTL_data$rs %in% locus_snps$rsID), ]


#Perform coloc
#install.packages("coloc")
library(coloc)

#N = 102762
#N minus 23andme = 102762 - 62231 = 40531
#ratio of cases to controls = 0.150
#using p-values

my.res <- coloc.abf(dataset1=list(pvalues=locus_snps$p.value, N=40531, snp=locus_snps$rsID, type="cc", s=0.150),
            dataset2=list(pvalues=locus_eQTL$p_score, N=672, snp=locus_eQTL$rs, type="quant"),
            MAF=locus_eQTL$af)


coloc_results <- my.res$results
coloc_results <- coloc_results[order(coloc_results$SNP.PP.H4,decreasing=T),]

write.table(coloc_results, file="./TwinsUK/coloc/results_CD207_minus_23andme.txt", quote=F, sep=" ", row.names=F)
