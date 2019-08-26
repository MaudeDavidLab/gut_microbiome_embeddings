library(dada2)
#path <- "../fastqs" # CHANGE ME to the directory containing the fastq files after unzipping.
data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/halfvarson/"
path = paste(data_dir, "fastqs/", sep = "")
fns <- list.files(path, pattern="fastq.gz")
filtpath <- paste(data_dir , "fastqs/filtered", sep = "") # CHANGE ME to the directory containing your filtered fastq files
out <- filterAndTrim(file.path(path,fns), file.path(filtpath,fns),
                     truncLen=99, maxEE=5, truncQ=11, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=F)


filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE) # CHANGE if different file extensions
sample.names <- sapply(strsplit(basename(filts), "\\."), `[`, 1) # Assumes filename = sample_XXX.fastq.gz
names(filts) <- sample.names

keep <- out[,"reads.out"] > 10000
filts <- filts[keep]
sample.names <- sample.names[keep]

##########################
## Learn error rates #####
##########################

set.seed(100)
err <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE)
# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
}
# Construct sequence table and write to disk
dds <- readRDS(paste(data_dir, "dada2_output/dds.rds", sep = ""))
keep <- sapply(dds, function(x) return(!is.null(x)))
seqtab <- makeSequenceTable(dds[keep])

##############################################
## Remove samples with too many or too #######
##  few reads ################################
##############################################

seqtab <- seqtab[rowSums(seqtab) < 1000000, ]
seqtab <- seqtab[rowSums(seqtab) > 10000, ]


#### These few lines label samples by their qid (matching with mapping file) instead of their 
#### err id. Specific to the way I downloaded this dataset

err_qid <- read.delim(paste(data_dir, "/err-to-qid.txt", sep = ""))
rownames(seqtab) <- as.character(err_qid$sample_title[match(rownames(seqtab), err_qid$run_accession)])
#saveRDS(seqtab, "dada2_output/seqtab.rds") # CHANGE ME to where you want sequence table saved

#################################################
#### write sequence table #######################
#################################################

write.table(seqtab, paste(data_dir, "seqtab.txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
saveRDS(seqtab, paste(data_dir, "seqtab.rds", sep = ""))

##################################################
###  Write representative sequences  #############
##################################################

fasta_file_name <- paste(data_dir, "queryseqs.fasta", sep = "")
taxa_seqs <- colnames(seqtab)
headers <- paste(">seq", seq(1, length(taxa_seqs)), sep = "")
fasta <- paste(headers, taxa_seqs, sep = "\n", collapse = "\n")
write(fasta, fasta_file_name)


