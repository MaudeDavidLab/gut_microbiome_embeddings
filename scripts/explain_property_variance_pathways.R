setwd("C:/Users/ctata/Documents/Lab/quality_vectors_git/data")
data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/"
pathway_dir = paste(data_dir, "/pathways/", sep = "")


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
taxa_names <- intersect(rownames(pathway_table), rownames(embed_table_glove)) #should be the same regardless of the embedding table
pathway_table <- pathway_table[taxa_names, ]
embed_table_glove <- embed_table_glove[taxa_names, ]
embed_table_glove <- apply(embed_table_glove, 2, function(x) return((x - mean(x) ) / sd(x)))



i = 1
r2 <- c()
for(i in seq(1, ncol(embed_table_glove))){
  vec_y <- embed_table_glove[ , i]
  vecs_x <- pathway_table
  df <- cbind(vec_y, vecs_x)
  linemod <- lm("vec_y ~ .", data = df)
  r2 <- c(r2, summary(linemod)$r.squared)
  #res <- predict.lm(linemod, vecs_x)
  #var(vec_y) - var(vec_y - res) #variance explained
  i = i + 1
}



#apparently I just re-derived an R-squared statistic, so that's cool

  