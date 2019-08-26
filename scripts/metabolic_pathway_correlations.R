library(vegan)
library(KEGGREST)
library(gtools)
library(RColorBrewer)
library(cowplot)
library(pheatmap)

setwd("C:/Users/ctata/Documents/Lab/quality_vectors_git/data")
data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/"
pathway_dir = paste(data_dir, "/pathways/", sep = "")


v <- readRDS(paste(pathway_dir, "pca_embedded_taxa.rds", sep = ""))



####################################
###  Read qual vecs  ###############
####################################

glove_emb <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/embed/glove_emb_AG_newfilter.07_100.txt",
                        quote="\"", comment.char="", row.names = 1, sep = " ", header = F)

glove_emb = glove_emb[-which(rownames(glove_emb) == '<unk>'), ]
colnames(glove_emb) <- paste("property_100_", seq(1, 100), sep = "")

###################################
###   Get pathway table   #########
###################################

pathway_table <- readRDS("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/pathways/otu_pathway_table.RDS")
keep <- colSums(pathway_table) > 0
keep2 <- colSums(pathway_table) < nrow(pathway_table)
pathway_table <- pathway_table[, keep & keep2]


#####################################
### Match, clean, and center   ######
#####################################

embed_table_glove <- glove_emb
embed_table_pca <- v[,1:100]
colnames(embed_table_pca) <- paste("pca", seq(1, ncol(embed_table_pca)), sep = "")


taxa_names <- intersect(rownames(pathway_table), rownames(embed_table_glove)) #should be the same regardless of the embedding table
pathway_table <- pathway_table[taxa_names, ]
embed_table_glove <- embed_table_glove[taxa_names, ]
embed_table_pca <- embed_table_pca[taxa_names, ]
embed_table_pca <- apply(embed_table_pca, 2, function(x) return((x - mean(x) ) / sd(x)))
embed_table_glove <- apply(embed_table_glove, 2, function(x) return((x - mean(x) ) / sd(x)))



####################################
###  Functions #####################
####################################

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


getMaxCorr <- function(embed_table, pathway_table){
  cor_mat <- getCorMat(embed_table, pathway_table)
  max_corr_inx <- apply(cor_mat, 1, function(cor_vec){
    return(which(cor_vec == max(cor_vec))[1])
  })
  return(colnames(pathway_table)[max_corr_inx])
}


makeNullPathwayTable_reorder <- function(pathway_table){
  new_order <- sample(seq(1, nrow(pathway_table)), size = nrow(pathway_table))
  return(pathway_table[new_order, ])
}

getTopNPathways <- function(cor_mat, n = 20){
  topN <- apply(cor_mat, 1, function(prop_path_cor){
    sorted <- sort(prop_path_cor, decreasing = T)
    return(names(sorted)[1:n])
  })
  return(topN)
}


permutationTest <- function(embed_vec, pathway_table){
  set.seed(123)
  #Find pathway with the highest correlation
  corrs <- apply(pathway_table, 2, cor, embed_vec)
  max_corr <- max(corrs)
  max_inx <- which(corrs == max_corr)
  
  #If we do the exact same process with randomly generated data multiple times, what are the chances we see a correlation has 
  #high as we did?
  maxPerm = 1000
  null_max_corrs <- c()
  for(iter in seq(1, maxPerm)){
    nullTable <- makeNullPathwayTable_reorder(pathway_table)
    corrs_null <- apply(nullTable, 2, cor, embed_vec)
    max_corr_null <- max(abs(corrs_null), na.rm = T)
    max_inx <- which(abs(corrs_null) == abs(max_corr_null))
    null_max_corrs <- c(null_max_corrs, max_corr_null)
  }
  pval <- sum(abs(null_max_corrs) >= abs(max_corr), na.rm=T) / maxPerm
  pval
  
  return(list(max_corr = max_corr,
              max_inx = max_inx,
              null_dist = null_max_corrs,
              pval = pval))
}


plotHeatmap <- function(cor_mat_list){
  breaksList = seq(-0.35, 0.35, by = .01)
  colors<-colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(length(breaksList))
  labs<- c("pca", "embed", "pca_null", "embed_null")
  i = 1
  heatmaps <- list()
  for(mat in cor_mat_list){
    if(i == 2){
      legend = T
    }else{
      legend = F
    }
    if(i ==1  | i == 2 | i ==3 | i == 4){
      col_labels = rep(" ", ncol(mat))
    }
    row_labels <- rep(" ", nrow(mat))
    
    heatmaps[[i]] <- pheatmap(mat, color = colors, breaks = breaksList,
                              treeheight_col = 0, treeheight_row = 0, 
                              labels_row = row_labels, labels_col = col_labels, 
                              border_color = NA, legend = F)
    #dev.off()
    i = i + 1
  }
  
  p <- plot_grid(heatmaps[[1]][[4]], heatmaps[[2]][[4]],
                 heatmaps[[3]][[4]], heatmaps[[4]][[4]])
  p
  
  pdf("../figures/cor_metabolic_pathways_grid.pdf", width = 5, height = 5)
  p
  dev.off()
}


#############################################################
######################## Plot heatmaps  ######################
##############################################################

null_pathway_table <- makeNullPathwayTable_reorder(pathway_table)
cor_mat_pca <- getCorMat(embed_table_pca, pathway_table)
cor_mat_glove <- getCorMat(embed_table_glove, pathway_table)
cor_mat_pca_null <- getCorMat(embed_table_pca, null_pathway_table)
cor_mat_glove_null <- getCorMat(embed_table_glove, null_pathway_table)

plotHeatmap(list(cor_mat_pca, cor_mat_glove, cor_mat_pca_null, cor_mat_glove_null))




###################################################################
############# get top pathway for each property ###################
###################################################################

corr_matches <- getMaxCorr(embed_table_glove, pathway_table)
path_names <- sapply(corr_matches, function(path_id){
  entry <- keggGet(paste("map", path_id, sep = ""))
  path_name <- entry[[1]]$NAME
  return(path_name)
})
df <- data.frame(colnames(embed_table_glove), corr_matches, path_names)
colnames(df) <- c("dim",	"pathway_id", 	"pathway_names")
write.table(df, "pathways/property_pathway_dict.txt", sep = "\t", quote = F, row.names = F)


#which pathway got picked a lot?
pathway_hits <- table(names(unlist(lapply(corr_matches, function(x) return(x$max_inx)))))
#out of 148 possible biological pathways, 54 had high correspondence with an embedding dimension or more. 
#There are 78 Desert pathways that are almost always present or almost always absent


####################################################################
#############  Get top 20 pathways for each property  ##############
####################################################################

top20Paths <- getTopNPathways(cor_mat_glove, n = 20)
top20Paths_named <- apply(top20Paths, 2, function(path_ids){
  lentries <- lapply(paste("map", path_ids, sep = ""), function(path) return(keggGet(path)))
  pathway_names <- unlist(lapply(lentries, function(entry) return(entry[[1]]$NAME)))
  return(pathway_names)
})
#Get pathway id to name dictionary

top20Paths_df <- top20Paths_named
colnames(top20Paths_df) <- colnames(embed_table_glove)
colnames(top20Paths) <- colnames(embed_table_glove)
top20Paths_ids <- apply(top20Paths, 2, function(x) return(paste("ko", x, sep = "")))

write.table(top20Paths_df, "pathways/top20Paths_per_property_names.csv", sep = ",", row.names = F)
write.table(top20Paths_ids, "pathways/top20Paths_per_property_ids.csv", sep = ",", row.names = F)

#####################################################################
###########  Calculate and save p values from permuation test #######
#####################################################################
increment <- 10
numGroups <- 100 / increment
#for(i in seq(1, numGroups)){
#  start <- (i-1) * increment + 1
#  end <- (i * increment)
#  corr_matches <- lapply(seq(start,end), function(i){
#    print(i)
#    return(permutationTest(embedding_table[ , i], pathway_table))
#  })
#  file <- paste("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/feces/piph/corr_matches_pca_", start, "_", end, ".RDS", sep = "")
#  print(file)
#  saveRDS(corr_matches, file)
#}



####################################################################################
#################   load corr_matches, which were calculated in increments, ########
#################    merge, and resave as one object                       #########
####################################################################################

corr_matches <- list()
for(i in seq(1, numGroups)){
  start <- (i-1) * increment + 1
  end <- (i * increment)
  file <- paste("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/AG_new/pathways/corr_matches_", start, "_", end, ".RDS", sep = "")
  print(file)
  corr_matches_tmp <- readRDS(file)
  corr_matches <- c(corr_matches, corr_matches_tmp)
}

saveRDS(corr_matches, paste(pathway_dir, "corr_matches.rds", sep = ""))
saveRDS(pathway_table, paste(pathway_dir, "pathway_table.RDS", sep = ""))


##########################################################################
#########   Load objects to get stats without recalculating  #############
##########################################################################
corr_matches_glove <- readRDS(paste(pathway_dir, "corr_matches_glove.rds", sep = ""))
corr_matches_pca <- readRDS(paste(pathway_dir, "corr_matches_pca.rds", sep = ""))
unlist(lapply(corr_matches_glove, function(x) return(x$pval)))
unlist(lapply(corr_matches_pca, function(x) return(x$pval))) 


#Check out the desert
tmp <- pheatmap(cor_mat_glove)
splits <- cutree(tmp$tree_col, h =  sort(tmp$tree_col$height, decreasing = T)[7])
desert_pathways <- names(splits[splits == 2])

hist(colSums(pathway_table[ , desert_pathways]), col = "red", breaks = 50)
hist(colSums(pathway_table[ , -which(colnames(pathway_table) %in% desert_pathways)]), col = "blue", add = T, breaks = 50)





###########################################################################
## For each property, find the lowest correlation that's still sig. #######
###########################################################################
library(pbapply)

findAllSig <- function(embed_vec, pathway_table){
  set.seed(123)
  #Find pathway with the highest correlation
  corrs <- apply(pathway_table, 2, cor, embed_vec)
  corrs <- sort(corrs, decreasing = T)
  #If we do the exact same process with randomly generated data multiple times, what are the chances we see a correlation has 
  #high as we did?
  maxPerm = 1000
  null_max_corrs <- c()
  for(iter in seq(1, maxPerm)){
    nullTable <- makeNullPathwayTable_reorder(pathway_table)
    corrs_null <- apply(nullTable, 2, cor, embed_vec)
    max_corr_null <- max(abs(corrs_null), na.rm = T)
    max_inx <- which(abs(corrs_null) == abs(max_corr_null))
    null_max_corrs <- c(null_max_corrs, max_corr_null)
  }
  return(corrs[corrs > max(null_max_corrs)])
}

getPathwayNames <- function(pathway_ids){
  lentries <- lapply(paste("map", pathway_ids, sep = ""), function(path) return(keggGet(path)))
  pathway_names <- unlist(lapply(lentries, function(entry) return(entry[[1]]$NAME)))
  return(pathway_names)
}

allSigCorrs <- pbapply(embed_table_glove, 2, findAllSig, pathway_table)
sigCorrs_ids <- lapply(allSigCorrs, function(i) return(names(i)))
sigCorrs_vals <- lapply(allSigCorrs, function(i) return(as.numeric(i)))
allSigPathNames <- lapply(allSigCorrs_ids, function(i) return(getPathwayNames(i)))

prop_label <- c()
for(i in seq(1, length(allSigCorrs))){
  prop_label <- c(prop_label, rep(paste("property_100_", i, sep = ""), length(allSigCorrs[[i]])))
}


df <- data.frame(property = prop_label, path_id = as.character(unlist(sigCorrs_ids)), path_name = as.character(unlist(allSigPathNames)), corr_val = unlist(sigCorrs_vals))

write.table(df, "pathways/property_pathway_dict_allsig.txt", row.names = F, col.names = T, quote = F, sep = "\t")
