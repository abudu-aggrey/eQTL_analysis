#read in text file with expression data, transpose, and create ID column.

CD207 <- read.table("./TwinsUK/CD207_transcript_data.txt")
CD207_trans <- as.data.frame(t(CD207))

#row names = sample id's
exp_id <- rownames(CD207_trans)
CD207_trans$ID_1 <- exp_id
colnames(CD207_trans)[1] <- "CD207_exp"

#transform expression data with rntransform

#install.packages("GenABEL", "./R/x86_64-pc-linux-gnu-library/3.3", repos='http://cran.us.r-project.org')
library(GenABEL)

CD207_trans$CD207_exp <- as.numeric(CD207_trans$CD207_exp)
CD207_exp_rntransform <- rntransform(CD207_exp, CD207_trans)
CD207_trans$CD207_exp <- CD207_exp_rntransform

#make same format as .sample file

CD207_trans <- as.matrix(CD207_trans)
header <- c("P", 0)
CD207_trans <- rbind(header, CD207_trans)
CD207_trans <- as.data.frame(CD207_trans)


#read in .sample file and merge

sample_file <- read.table("./TwinsUK/Eurobats_genotypes/chr2_Eurobats_Public_subset_samples.sample", header=T)

#install "plyr" package to use "join" to merge
#install.packages("plyr")
library(plyr)

exp_sample <- join(sample_file, CD207_trans, by="ID_1")


#read in text file with sample age
age <- read.table("./TwinsUK/twin_age.txt", header=T)

header_1 <- c(0, "P")
age <- as.matrix(age)
age <- rbind(header_1, age)
age <- as.data.frame(age)

#rename first column:
colnames(age)[1] <- "ID_1"

#merge with sample data:

exp_sample <- join(exp_sample, age, by="ID_1")

#create intercept column (use when creating covariate file)
exp_sample$intercept <- c("B", rep(1,672))

write.table(exp_sample, "./TwinsUK/Eurobats_genotypes/chr2_Eurobats_Public_exp_pheno.sample", row.names=F, quote=F, sep=" ")


