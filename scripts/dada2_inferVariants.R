library(dada2); packageVersion("dada2")
args = commandArgs(trailingOnly=TRUE)
filtpath = args[1]
errorFile = args[2]
print(paste("Filt path: ", filtpath))
print(paste("Error File: ", errorFile))
outFile = args[3]
print(paste("Outfile: ", outFile))

# File parsing
#filtpath <- "/nfs3/PHARM/David_Lab/christine/public_data/AG_new/fastqs/runs_use/filtered_150/subset2" # CHANGE ME to the directory containing your filtered fastq files
print(filtpath)
filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE) # CHANGE if different file extensions
sample.names <- unlist(lapply(strsplit(basename(filts), "\\."), function(x) return(x[1])))
names(filts) <- sample.names
# Learn error rates
set.seed(100)


#err <- learnErrors(filts, nbases = 1e8, multithread=TRUE, randomize=TRUE)
#saveRDS(err, "/nfs3/PHARM/David_Lab/christine/public_data/AG_new/dada2_output/filtered_100/err.rds")
#errorFile = "/nfs3/PHARM/David_Lab/christine/public_data/AG_new/dada2_output/filtered_150/err.rds"
err <- readRDS(errorFile)
print(paste("Loaded error rates from ", errorFile))

# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  if(length(derep$uniques) < 100000){
      dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
  }
  else{
    print(paste("sample has ", length(derep$uniques), " unique reads"))
  }
  #saveRDS(dds, "/nfs3/PHARM/David_Lab/christine/public_data/AG_new/dada2_output/filtered_100/dds.rds")
}
# Construct sequence table and write to disk
seqtab <- makeSequenceTable(dds)
saveRDS(seqtab, outFile) # CHANGE ME to where you want sequence table saved
