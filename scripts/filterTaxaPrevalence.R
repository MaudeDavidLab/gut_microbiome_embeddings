seqtab <- readRDS("../dada2_output/filtered_150/seqtab_final.rds")

#Only keep samples with more than 5000 reads
seqtab <- seqtab[rowSums(seqtab) > 5000, ]

#filter: keep taxa present in more than .07% of samples, or 13 samples
keep.07 <- apply(seqtab, 2, function(vec) return(sum(vec > 0) > 13))


print("Number of taxa at each level")
print(paste(".07 %", sum(keep.07)))

seqtab.07 <- seqtab[ , keep.07]
write.table(seqtab, "../dada2_output/filtered_150/seqtab_final_filter.07.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

