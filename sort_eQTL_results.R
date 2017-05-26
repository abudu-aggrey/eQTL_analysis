#read in results file

res <- read.table("./TwinsUK/CD207_eQTL_analysis_test.assoc.txt", header=T)

#order by p-score

res_sorted <- res[order(res$p_score), ]

#write to table

write.table(res_sorted, "./TwinsUK/CD207_eQTL_analysis_test.assoc.txt", sep=" ", quote=F, row.names=F)
