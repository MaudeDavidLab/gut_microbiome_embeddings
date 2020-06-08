
library(Biostrings)
library(rBLAST)
library(seqinr)
library(dplyr)
library(data.table)
library(pbapply)
source("scripts/embed_your_data/convert_ASV_embedIDs.R")



#' @export
getExampleSeqtab <- function(){
  seqtab <- utils::read.csv("data/test_files/seqtab_test.csv")
  rownames(seqtab) <- seqtab$X
  rownames(seqtab) <- gsub("X", "", rownames(seqtab))
  seqtab <- seqtab[, 2:ncol(seqtab)]
  return(seqtab)
}

getExampleSeqtab_AG <- function(){
  seqtab <- data.table::fread("data/seqtab_final_filter.07.txt", header = F, sep = '\t')
  seqs <- scan("data/sequences_.07.txt", what = character())
  seqtab <- as.data.frame(seqtab)
  rownames(seqtab) <- seqtab$V1
  rownames(seqtab) <- gsub("X", "", rownames(seqtab))
  seqtab <- seqtab[, 2:ncol(seqtab)]
  colnames(seqtab) <- seqs
}

getExampleSeqtab_Halfvarson <- function(){
  seqtab <- data.table::fread("data/halfvarson/seqtab.txt", header = F, sep = '\t')
  fasta <- seqinr::read.fasta("data/halfvarson/sequences.fasta")
  seqs <- toupper(unlist(seqinr::getSequence(fasta, as.string = T)))
  seqtab <- as.data.frame(seqtab)
  rownames(seqtab) <- seqtab$V1
  rownames(seqtab) <- gsub("X", "", rownames(seqtab))
  seqtab <- seqtab[, 2:ncol(seqtab)]
  colnames(seqtab) <- seqs
  return(seqtab)
}


#' @export
getExampleFasta <- function(){
  fasta_file <- "data/test_files/fasta_test.fasta"
  return(fasta_file)
}

# This function takes in the best_hits datatable from the blast output, and creates an asv by embedding database asv
# matrix where each element is 1 over the number of hits for that ASV in the entire database. Instead of breaking ties arbitrarily
# we allow each ASV to hit multiple embedding ASV
getHitTransformation <- function(best_hits){
  asv2embed_split <- split(best_hits, best_hits$qseqid)
  asv2embed_list <- lapply(asv2embed_split, function(x) {
    num <- data.frame(matrix(rep(1, nrow(x)), nrow = 1))
    colnames(num) <- x$sseqid
    return(num)
  })

  asv2embed_short <- rbind.fill(asv2embed_list)
  asv2embed_short[is.na(asv2embed_short)] <- 0
  asv2embed_short <- t(apply(asv2embed_short, 1, function(x) return(x / sum(x))))
  rownames(asv2embed_short) <- names(asv2embed_list)
  return(asv2embed_short)
}

#This function runs blast to find the closest hits to each query sequence in the embedding database
#If blast has been run externally for increased speed, the function may take in the blast output table
#If this function is taking more than a few minutes to run, consider running BLAST externally
#"$blast_software_dir/blastn" -db "$blast_db/embedding_db_.07" -query "$data_dir/repseqs.fasta" -out "$out_dir/blast_hits.tsv"  -outfmt "6 qseqid sseqid qseq sseq evalue bitscore length pident"
#cat "$out_dir/blast_hits.tsv" | sort -k1,1 -k5,5g -k6,6nr | sort -u -k1,1 --merge > "$out_dir/best_hits.tsv"

getBestHits <- function(blast_hits_file, id_thresh = 99){
  #save best hits per query ASV
  print("Filtering best hits from blast output")
  blast_hits <- read.delim(blast_hits_file, header= TRUE, sep = "\t")
  colnames(blast_hits) <- c('qseqid', 'sseqid', 'qseq', 'sseq', 'evalue', 'bitscore', 'length', 'pident')
  best_hits_list <- pblapply(unique(blast_hits$qseqid), function(asv_name){
    tmp <- blast_hits[blast_hits$qseqid == asv_name, ]
    tmp <- tmp[tmp$evalue == min(tmp$evalue), ]
    tmp <- tmp[tmp$pident >= id_thresh, ]
    return(tmp)
  })
  best_hits <- do.call(rbind, best_hits_list)
  return(best_hits)
}

getEmbeddingHits <- function(fasta_file, blast_hits_file = NA, id_thresh = 99, out_dir = "data/blast_output/"){
  if(is.na(blast_hits_file)){
    runBlast(fasta_file = fasta_file, out_dir = out_dir)
    blast_hits_file = file.path(out_dir, "blast_hits.tsv")
  }
  return(getBestHits(blast_hits_file, id_thresh))
}


getFastaDF <- function(fasta_file){
  fasta <- read.fasta(fasta_file)
  asv_ids <- names(fasta)
  asv_seqs <- toupper(unlist(seqinr::getSequence(fasta, as.string = T)))
  fasta_df <- data.frame(asv_ids, row.names = asv_seqs)
  return(fasta_df)
}

#For each query sequence, assigns the name of the closest embedding sequence. If there are multiple closest hits, splits the abundance
#of the query sequence evenly among all closest hits
transformSeqtab <- function(seqtab, best_hits, fasta_file_name ){

  #sample by ASV
  print('converting ids')
  colnames(seqtab) <- convertIDs(ids = colnames(seqtab), from_id = "ASV", to_id = "embedIDs", fasta_file_name) #fasta should contain all sequences in seqtab
  seqtab <- seqtab[, colnames(seqtab) %in% as.character(best_hits$qseqid)]

  # Split counts among all best hits
  print('renaming ASVs by embedding database')
  asv2embed_short <- getHitTransformation(best_hits)
  asv2embed_short <- asv2embed_short[colnames(seqtab), ]
  colnames(seqtab) == rownames(asv2embed_short)

  seqtab_transformed <- as.matrix(seqtab) %*% as.matrix(asv2embed_short)
  return(seqtab_transformed)
}

#' @export
embedSeqtab <- function(seqtab, best_hits, embedding_file_name, fasta_file_name){
  seqtab_transformed <- transformSeqtab(seqtab = seqtab, best_hits = best_hits, fasta_file_name = fasta_file_name)
  qual_vecs <- read.csv(embedding_file_name, row.names=1, sep="")
  qual_vecs <- qual_vecs[colnames(seqtab_transformed), ]
  embedded <- as.matrix(seqtab_transformed) %*% as.matrix(qual_vecs)
  return(embedded)
}



#Instructions:
#1. Collect your sequence table and fasta file. Seqtab will come from ASV processing like Dada2, and by in a sample by ASV orientation
seqtab = getExampleSeqtab()
fasta_file_name = getExampleFasta()


#2. Download all files in the folder http://files.cgrb.oregonstate.edu/David_Lab/microbiome_embeddings/data/blastdb/. Set variable below equal to path name of the download on your machine
blast_db = "data/blast_db"


#3. Run blast to rename your sequences with the nearest sequence available in the embedding matrix
output_file_name = "data/blast_output/blast_hits.tsv"
#Use the following command template to run blast using the command line
#blast_software_dir/blastn -db path_to_blast_db -query fasta.fasta -out output_file_name  -outfmt "6 qseqid sseqid qseq sseq evalue bitscore length pident"

best_hits = getBestHits(blast_hits_file = output_file_name, id_thresh = 99)


#4. Download the embedding transformation matrix from here: http://files.cgrb.oregonstate.edu/David_Lab/microbiome_embeddings/data/embed/embed_.07_100dim.txt
#set the variable below to the file name on your local machine
embedding_file_name = "data/embed/embed_.07_100dim.txt"


#5. Embed your data
embedded = embedSeqtab(seqtab, best_hits, embedding_file_name, fasta_file_name)







