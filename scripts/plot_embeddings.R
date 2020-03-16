glove_emb <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors_final/data/embed/glove_emb_AG_newfilter.07_100.txt",
                        quote="\"", comment.char="", row.names = 1, sep = " ", header = F)

glove_emb = glove_emb[-which(rownames(glove_emb) == '<unk>'), ]
colnames(glove_emb) <- paste("property_100_", seq(1, 100), sep = "")

library(Rtsne)

Rtsne(glove_emb)
