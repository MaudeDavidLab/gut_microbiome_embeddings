library(dada2); packageVersion("dada2")

seqtabs = list()
for(i in seq(40)){
	if(i < 10){
		i = paste("0", i, sep = "")
	}
	file = paste("../dada2_output/filtered_150/seqtab_0", i, ".rds",  sep = "")
	print(file)
	seqtabs[[i]] <- readRDS(file)

}

all <- mergeSequenceTables(tables = seqtabs)

seqtab <- removeBimeraDenovo(all, method = "consensus", multithread = TRUE)

saveRDS(seqtab, "../dada2_output/filtered_150/seqtab_final.rds")
write.table("../dada2_output/filtered_150/seqtab_final.txt")



############################################################
##### Filter while we have the object loaded  ##############
############################################################
seqtab <- readRDS("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/seqtab_final_filter.07.txt")
#Only keep samples with more than 5000 reads
seqtab <- seqtab[rowSums(seqtab) > 5000, ]

#filter: keep taxa present in more than .07% of samples, or 13 samples
nsamples <- nrow(seqtab)

prevalenceFilter <- function(rate = 0.07){
  threshold <- nsamples * rate
  keep <- apply(seqtab, 2, function(vec) return(sum(vec > 0) > threshold))
  
  print("Number of taxa at each level")
  print(paste(rate,  "%", sum(keep)))
  
  seqtab_filt <- seqtab[ , keep]
  saveRDS(seqtab, paste("../data/AG_new/seqtab_filter", rate, ".rds", sep = "" ))
  write.table(seqtab, paste("../data/AG_new/seqtab_filter", rate, ".txt", sep = "" ), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
}

rates <- c(0.01, 0.03, 0.05, 0.07, 0.09)
for(rate in rates){
  prevalenceFilter(rate)
}
