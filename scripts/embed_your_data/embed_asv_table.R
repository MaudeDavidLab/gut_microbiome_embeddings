library(readr)
#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
#if (length(args)==0) {
#  stop("At least one argument must be supplied (input file).n", call.=FALSE)
#} 

####### Pass appropriate arguments through embedDataset.sh ###############
##########################################################################

#asv_file = args[1]
#repseqs_file = args[2]
#best_hits = args[3]
#qual_vec_file = args[4]
#out_dir = args[5]
setwd("C:/Users/ctata/Documents/Lab/quality_vectors_final/data/halfvarson/")
asv_file = "seqtab_t.txt"
repseqs_file = "queryseqs.fasta"
best_hits_file = "embed_2/best_hits.tsv"
qual_vec_file = "../embed/embed_.07_100dim.txt"
out_dir = "embed_2/"

id_thresh = 99
#########################
### Read ASV table ######
#########################
#taxa are rows when read in, but we transform the matrix so that the dot product makes sense

seqtab <- read.table(asv_file, header = T, row.names = 1)
seqtab <- t(seqtab)
#seqtab <- readRDS(asv_file)
colnames(seqtab) <- gsub("X","", colnames(seqtab))


#############################
## Read best hits table #####
## from blast ###############
#############################

best_hits <- read.delim(best_hits_file, header=FALSE, row.names = 1)  ## CHANGE TO BLAST OUTPUT - SEE README
colnames(best_hits) <- c("hit_id", "query_seq", "hit_seq", "evalue", "bitscore", "length", "piden")

#Filter best_hits table to only include hits that pass the e-value threshold
best_hits <- best_hits[best_hits$evalue < 1*10^(-29), ]
best_hits <- best_hits[best_hits$piden >= id_thresh, ]


##M3 specific
#asvid_map <- read.csv("ASVid_map.csv")
# write fasta
#fasta <- paste(">", asvid_map$ASVid, sep = "")
#fasta <- paste(fasta, asvid_map$Sequence, collapse = "\n", sep = "\n")
#write(fasta, "M3_fasta.fasta")
#fasta_df <- asvid_map
#rownames(fasta_df) <- fasta_df$ASVid
#asvids <- intersect(colnames(seqtab), rownames(best_hits))
#best_hits <- best_hits[asvids, ]
#seqtab <- seqtab[ , asvids]
#colnames(seqtab) <- best_hits$hit_id

#Match up the query id in best_hits to the full ASV in colnames(seqtab) using the repseqs.fasta file
fasta <- read.table(repseqs_file)
headers <- gsub(">", "", fasta[seq(1, nrow(fasta), by = 2), 1])
query_seqs <- as.character(fasta[seq(2, nrow(fasta), by = 2), 1])
fasta_df <- data.frame(query_seqs, row.names = headers)
best_hits$full_query_seq <- as.character(fasta_df[rownames(best_hits), 1])
best_hits <- best_hits[as.character(best_hits$full_query_seq) %in% colnames(seqtab), ]

print(dim(best_hits))

#################################
### Assign nearest neighbor id ##
#################################

#Drop any ASVs from the table that don't have near enough hits in the transformation matrix
seqtab_hits <- seqtab[ , as.character(best_hits$full_query_seq) ]
print(dim(seqtab_hits))

#Assign the id of each ASV's nearest hit in the embedding transformation table.
seqtab_out <- seqtab_hits
colnames(seqtab_out) <- best_hits$query_seq
colnames(seqtab_hits) <- best_hits$hit_id


##############################
## Read quality vector  ######
## transformation table ######
##############################
#transform_mat_file <- paste(data_dir, "/embed/embed_.07_100dim.txt", sep = "") #CHANGE TO EMBED IN A DIFFERENT DIMENSIONAL SPACE
qual_vecs <- read.table(qual_vec_file)
qual_vecs_ord <- qual_vecs[colnames(seqtab_hits), ]

#############################################
### Take dot product and save file ##########
#############################################

embedded <- as.matrix(asinh(seqtab_hits)) %*% as.matrix(qual_vecs_ord)
#
embedded_out_file = paste(out_dir, "embedded_.07_100dim_",id_thresh, sep = "")
saveRDS(embedded, paste(embedded_out_file, ".rds", sep = ""))
write.table(embedded, paste(embedded_out_file, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)


seqtab_file = paste(out_dir, "seqtab_.07_100dim_", id_thresh, sep = "")
write.table(seqtab_out, paste(seqtab_file, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)






