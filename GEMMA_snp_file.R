#create SNP annotation file 


#read in original impuatation ".gen" file
gen <- read.table("./TwinsUK/Eurobats_genotypes/chr2_Eurobats_Public.gen")

#read in text file of SNP names
CD207_snps <- read.table("./TwinsUK/SNP_list.txt")

#name the SNP column of both files
colnames(gen)[2] <- "RSID"
colnames(CD207_snps)[1] <- "RSID"

#merge CD207_snps with gen file
snp_list <- merge(CD207_snps, gen, by="RSID")

#subset "RSID" and "BP"
snp_file <- snp_list[ ,c(1,3)]
colnames(snp_file)[2] <- "BP"

#create a CHR column
snp_file$CHR <- rep(2,length(CD207_snps$RSID))

#remove headers and write to table
names(snp_file) <- NULL

write.table(snp_file, "./TwinsUK/Eurobats_genotypes/SNP_annotation_file.txt", sep=" ", quote=F, row.names=F)
