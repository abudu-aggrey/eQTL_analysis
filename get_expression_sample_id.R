#read in text file with expression data, transpose, and create ID column.

CD207 <- read.table("./TwinsUK/CD207_transcript_data.txt")
CD207_trans <- as.data.frame(t(CD207))
exp_id <- rownames(CD207_trans)

#write list of ID's with expression data to text file
exp_id <- as.data.frame(exp_id)
colnames(exp_id) <- NULL

write.table(exp_id, "./TwinsUK/expression_sample_ids.txt", row.names=F, quote=F)