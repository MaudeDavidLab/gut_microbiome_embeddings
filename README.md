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

## 3. Assign the appropriate embedding sequence id to each of your sequences. If there is no match, toss out that sequence. See scripts/embed_your_data/embed_asv_table.R 
#### See step 4 of scripts/embed_your_data/embedDataset.sh

#### You have successfully embedded your data! Each row represents that sample's microbiome signature.
