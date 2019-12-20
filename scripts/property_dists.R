library(vegan)
library(KEGGREST)
library(gtools)
library(RColorBrewer)
library(cowplot)
library(pheatmap)
library(lsa)

####This script should answer the question:
## "Do similar properties have similar metabolic pathways". 
## We find distance between each property based on the taxa vectors
## THEN find distance between each property using the metabolic pathway
## correlation vector. We then apply a Mantel test to ascertain the similarity
## of the two distance matrices calculated using the two different metrics.

setwd("C:/Users/ctata/Documents/Lab/quality_vectors_final/data")
data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_final/data/"
pathway_dir = paste(data_dir, "/pathways/", sep = "")


####################################
###  Read qual vecs  ###############
####################################

glove_emb <- read.table(paste(data_dir, "/embed/glove_emb_AG_newfilter.07_100.txt", sep = ""),
                        quote="\"", comment.char="", row.names = 1, sep = " ", header = F)

glove_emb = glove_emb[-which(rownames(glove_emb) == '<unk>'), ]
colnames(glove_emb) <- paste("property_100_", seq(1, 100), sep = "")
  
###################################
###   Get pathway table   #########
###################################

pathway_table <- readRDS(paste(data_dir, "/pathways/otu_pathway_table.RDS", sep = ""))
keep <- colSums(pathway_table) > 0
keep2 <- colSums(pathway_table) < nrow(pathway_table)
pathway_table <- pathway_table[, keep & keep2]


#####################################
### Match, clean, and center   ######
#####################################

embed_table_glove <- glove_emb
taxa_names <- intersect(rownames(pathway_table), rownames(embed_table_glove)) #should be the same regardless of the embedding table
pathway_table <- pathway_table[taxa_names, ]
embed_table_glove <- embed_table_glove[taxa_names, ]
embed_table_glove <- apply(embed_table_glove, 2, function(x) return((x - mean(x) ) / sd(x)))


#####################################
## Get correlation matrics between ##
##  properties and pathways  ########
#####################################

getCorMat <- function(embedding_table, pathway_table){
  cor_list <- list()
  for(i in seq(1, ncol(embedding_table))){
    cor = apply(pathway_table, 2, function(pathway_vec) return(cor(pathway_vec, embedding_table[ ,i])))
    cor_list[[i]] <- cor
  }
  cor_mat <- data.frame(matrix(unlist(cor_list), byrow = T, nrow = length(cor_list)))
  colnames(cor_mat) <- colnames(pathway_table)
  rownames(cor_mat) <- colnames(embedding_table)
  return(cor_mat)
}

corrmat <- getCorMat(embed_table_glove, pathway_table)


########################################################
### Find dist between properties using taxa vectors  ###
########################################################
prop_dists_taxa <- cosine(embed_table_glove)


########################################################
### Find dist between properties using taxa vectors  ###
########################################################
prop_dists_corrs <- cosine(t(corrmat))


###############################
### Mantel Test ###############
###############################
man <- mantel(prop_dists_taxa, prop_dists_corrs, permutations = 10000)


#####################
### procrustes#######
#####################
prop_pca_taxa = prcomp(embed_table_glove)
prop_pca_corr = prcomp(t(corrmat))

pro <- procrustes(prop_pca_taxa$rotation, prop_pca_corr$rotation, scores = "sites")
plot(pro, kind = 1)

protest(prop_pca_taxa$rotation, prop_pca_corr$rotation, scores = "sites", permutations = 10000) 




### Find properties that are closest together
inx = which(prop_dists_taxa == min(abs(prop_dists_taxa[prop_dists_taxa > 0])), arr.ind = T)

### Find the pathways that these properties most closesly relate to
corr1 <- corrmat[ , inx[1,1]]
corr2 <- corrmat[ , inx[1, 2]]


plot(seq(1, length(corr1)), corr1, type = "l", col = "red")
lines(seq(1, length(corr2)), corr2)


plot(prop_dists_taxa, prop_dists_corrs)
#similar properties largely correlate with the same pathways



