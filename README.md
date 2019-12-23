# Data files: http://files.cgrb.oregonstate.edu/David_Lab/microbiome_embeddings/

# How to embed your data
## All code in this tutorial is meant as a guide/example. It's not a complete wrapper, and will likely not be plug-and-play for you and your system. A complete package is in the works, however, it is not available at this time.

## 1. Begin with a sample by ASV table as processed by Dada2, and a fasta file of the representative sequences.
Use this tutorial: https://benjjneb.github.io/dada2/tutorial.html to process your raw fastq reads. At the end of the tutorial, you are left with a variable called seqtab.nochim which is a sample by ASV table. Write this table to a file. Also write a fasta file with each full length sequence. See the end of scripts/embed_your_data/dada2.R or below:
To write a fasta file containing the sequences: 
 ```
 Rscript
    #################################################
    #### write sequence table #######################
    #################################################
  ```
    write.table(seqtab, paste(data_dir, "seqtab.txt", sep = "\t"),
                quote = F, sep = "\t", row.names = T, col.names = T)

    ##################################################
    ###  Write representative sequences  #############
    ##################################################

    fasta_file_name <- paste(data_dir, "repseqs.fasta", sep = "")
    taxa_seqs <- colnames(seqtab)
    headers <- paste(">seq", seq(1, length(taxa_seqs)), sep = "")
    fasta <- paste(headers, taxa_seqs, sep = "\n", collapse = "\n")
    write(fasta, fasta_file_name)
  ```
 
  ```
## 2. Blastn against the embedding database, and filter for best results
#### See steps 2 and 3 of scripts/embed_your_data/embedDataset.sh
```
"$blast_software_dir/blastn" -db "$blast_db/embedding_db_.07" -query "$data_dir/repseqs.fasta" -out "$out_dir/blast_hits.tsv"  -outfmt "6 qseqid sseqid qseq sseq evalue bitscore length pident"

cat "$out_dir/blast_hits.tsv" | sort -k1,1 -k5,5g -k6,6nr | sort -u -k1,1 --merge > "$out_dir/best_hits.tsv"

```

## 3. Assign the appropriate embedding sequence id to each of your sequences. If there is no match, toss out that sequence. See scripts/embed_your_data/embed_asv_table.R 
#### See step 4 of scripts/embed_your_data/embedDataset.sh

```
#########################
### Read ASV table ######
#########################
#taxa are rows when read in, but we transform the matrix so that the dot product makes sense

#seqtab <- read.table(asv_file, header = T, row.names = 1)
seqtab <- readRDS(asv_file)
colnames(seqtab) <- gsub("X","", colnames(seqtab))


#############################
## Read best hits table #####
## from blast ###############
#############################

id_thresh = 99
best_hits <- read.delim(best_hits_file, header=FALSE, row.names = 1)  ## CHANGE TO BLAST OUTPUT - SEE README
colnames(best_hits) <- c("hit_id", "query_seq", "hit_seq", "evalue", "bitscore", "length", "piden")

#Filter best_hits table to only include hits that pass the e-value threshold
best_hits <- best_hits[best_hits$evalue < 1*10^(-29), ]
best_hits <- best_hits[best_hits$piden >= id_thresh, ]



fasta <- read.table(repseqs_file)
headers <- gsub(">", "", fasta[seq(1, nrow(fasta), by = 2), 1])
query_seqs <- as.character(fasta[seq(2, nrow(fasta), by = 2), 1])
fasta_df <- data.frame(query_seqs, row.names = headers)
best_hits$full_query_seq <- as.character(fasta_df[rownames(best_hits), 1])
best_hits <- best_hits[as.character(best_hits$full_query_seq) %in% colnames(seqtab), ]

print(dim(best_hits))

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
embedded_out_file = paste(out_dir, "embedded_.07_besthit_100dim", sep = "")
saveRDS(embedded, paste(embedded_out_file, ".rds", sep = ""))
write.table(embedded, paste(embedded_out_file, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)


seqtab_file = paste(out_dir, "seqtab_.07_97id_100dim", sep = "")
write.table(seqtab_out, paste(seqtab_file, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)

```

#### You have successfully embedded your data! Each row represents that sample's microbiome signature.
