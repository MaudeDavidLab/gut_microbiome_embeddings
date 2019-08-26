data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/" ##CHANGE TO YOUR DATA DIRECTORY

#########################
### Read ASV table ######
#########################
asv_file = paste(data_dir, "/halfvarson/seqtab.rds", sep = "") ## CHANGE TO YOUR DATASET OF INTEREST
seqtab <- readRDS(asv_file)
#seqtab <- read.table(asv_file, quote = F, sep = "\t", row.names = T, col.names = T)

#### change column names to id
query_seqs_file = paste(data_dir, "/halfvarson/queryseqs.fasta", sep ="") ## CHANGE TO YOUR FASTA SEQUENCE FILE
query_fasta <- read.table(query_seqs_file)
query_headers <- query_fasta[seq(1, nrow(query_fasta), by = 2), 1]
query_headers <- gsub(">", "", query_headers)
query_seqs <- as.character(query_fasta[seq(2, nrow(query_fasta), by = 2), 1])

query_headers <- query_headers[query_seqs %in% colnames(seqtab)]
query_seqs <- query_seqs[query_seqs %in% colnames(seqtab)]

colnames(seqtab)[match(colnames(seqtab), query_seqs)] == query_seqs #check
colnames(seqtab)[match(colnames(seqtab), query_seqs)] <- query_headers

##############################
## Read quality vector  ######
## transformation table ######
##############################
transform_mat_file <- paste(data_dir, "/embed/embed_.07_100dim.txt", sep = "") #CHANGE TO EMBED IN A DIFFERENT DIMENSIONAL SPACE
qual_vecs <- read.table(transform_mat_file)

#############################
## Read best hits table #####
## from blast ###############
#############################

blast_output_file = paste(data_dir, "halfvarson/embed/best_hits.tsv", sep = "")  ## CHANGE TO BLAST OUTPUT - SEE README
best_hits <- read.delim(blast_output_file, header=FALSE, row.names = 1)
colnames(best_hits) <- c("hit_id", "query_seq", "hit_seq", "evalue", "bitscore")

#Filter best_hits table to only include hits that pass the e-value threshold
best_hits <- best_hits[best_hits$evalue < 1*10^(-29), ]

#################################
### Assign nearest neighbor id ##
#################################

#Drop any ASVs from the table that don't have near enough hits in the transformation matrix
seqtab <- seqtab[ , colnames(seqtab) %in% rownames(best_hits)] #17784 taxa left

#Assign the id of each ASV's nearest hit in the embedding transformation table.
colnames(seqtab) <- best_hits[colnames(seqtab), "hit_id"]


###############################
### Match orders ##############
###############################
qual_vecs <- qual_vecs[colnames(seqtab), ]


#############################################
### Take dot product and save file ##########
#############################################

embedded <- as.matrix(seqtab) %*% as.matrix(qual_vecs)
embedded_file = paste(data_dir, "halfvarson/embed/seqtab_embedded_.07_100dim", sep = "")
saveRDS(embedded, paste(embedded_file, ".rds", sep = ""))
write.table(embedded, paste(embedded_file, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)








