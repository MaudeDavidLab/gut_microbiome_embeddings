

setwd("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/")

combineKOlists <- function(l){
  combo <- do.call(rbind, lapply(lapply(l, unlist), "[",
                                 unique(unlist(c(sapply(l,names))))))
  return(combo)
}


getKOLists <- function(genomes){
  ko_lists <- list()
  
  
  #ko_lists <- ko_lists[1:365]
  #genomes <- genomes[!(genomes %in% names(ko_lists))]
  
  
  for(i in seq(along = genomes)){
    cat(paste(i, " ") )
    genome <- genomes[i]
    print(genome)
    
    errorFlag <- FALSE
    genome_entry <- tryCatch({
      keggGet(paste("genome:", genome, sep = ""))
    }, error = function(e){
      cat("error finding ")
      print(as.character(genome))
      errorFlag <<- TRUE
    })
    if(errorFlag){
      next
    }
    
    taxa <- genome_entry[[1]]$TAXONOMY
    taxa_name <- taxa$LINEAGE
    taxa_id <- taxa$TAXONOMY
    taxa_subspecies <- genome_entry[[1]]$DEFINITION
    taxa_keywords <- genome_entry[[1]]$KEYWORDS
    ko_list <- table(names(keggLink(genome_entry[[1]]$ENTRY, "ko")))
    ko_lists[[as.character(genome)]] <- ko_list
  }
  return(ko_lists)
}

#we have to do pathway presence/absence, because abundance doesn't make sense
getPathwayTable <- function(otu_genome_hit_table){
  pathway_lists <- list()
  org_names <- c()
  genomes <- unique(otu_genome_hit_table$NNgenome)
  
  for(i in seq(along= genomes)){
    print(paste(i, genomes[i]))
    
    errorFlag <- FALSE
    pathways <- tryCatch({
      unique(names(keggLink(genomes[i], "path")))
      
    }, error = function(e){
      cat("error finding ")
      print(as.character(genomes[i]))
      errorFlag <<- TRUE
    })
    if(errorFlag){
      next
    }
    pathways <- gsub(paste("path:", as.character(genomes[i]), sep = ""), "", pathways)
    df = data.frame(matrix(rep(1, length(pathways)), nrow = 1))
    
    colnames(df) <- pathways
    pathway_lists[[i]] <- df
    org_names <- c(org_names, as.character(genomes[i]))
  }
  
  pathway_table <- rbind.fill(pathway_lists)
  pathway_table[is.na(pathway_table)] <- 0
  rownames(pathway_table) <- org_names
  return(pathway_table)
}

############################################
#####  Read in Piphillan output table ######
############################################

otu_genome_hit_table <-read_delim("otu_genome_hit_table_2.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
otu_genome_hit_table <- otu_genome_hit_table[!duplicated(otu_genome_hit_table$OTU), ]


####################################################################
##### Read in sequence fasta to rename with full length sequence ###
####################################################################
seqs_fasta <- read.table("pathways/seqs_filter.07_piph.fasta", quote="\"", comment.char="")
headers <- gsub(">", "", as.character(seqs_fasta[seq(1, nrow(seqs_fasta), 2), ]))
seqs <- as.character(seqs_fasta[seq(2, nrow(seqs_fasta), 2), ])
seq_id_df <- data.frame(seqs)
rownames(seq_id_df) <- headers

otu_genome_hit_table$OTU <- seq_id_df[otu_genome_hit_table$OTU, ]

##################################################
###### Creates genomecode to pathway table #######
###### Table provided, no need to rerun    #######
##################################################

pathway_table <- getPathwayTable(otu_genome_hit_table)
saveRDS(pathway_table, "genomecode_pathway_table_2.rds")


#######################################
#### Construct ASV x pathway table ####
#######################################

pathway_table <- readRDS("pathways/genomecode_pathway_table.rds")


otu_genome_pruned <- otu_genome_hit_table[otu_genome_hit_table$NNgenome %in% rownames(pathway_table), ]
otu_pathway_table <- pathway_table[as.character(otu_genome_pruned$NNgenome), ]
rownames(otu_pathway_table) <- otu_genome_pruned$OTU

saveRDS(otu_pathway_table, "pathways/otu_pathway_table.RDS")
write.table(otu_pathway_table, "pathways/otu_pathway_table.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE, quote = FALSE)







