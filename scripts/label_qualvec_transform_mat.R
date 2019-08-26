# 

qual_vecs_50 <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/embed/glove_emb_AG_newfilter.07_50.txt", row.names = 1)
qual_vecs_100 <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/embed/glove_emb_AG_newfilter.07_100.txt", row.names = 1)
qual_vecs_250 <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/embed/glove_emb_AG_newfilter.07_250.txt", row.names = 1)
qual_vecs_500 <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/embed/glove_emb_AG_newfilter.07_500.txt", row.names = 1)
qual_vecs_750 <- read.table("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/embed/glove_emb_AG_newfilter.07_750.txt", row.names = 1)

headers <- paste("embed", seq(1, nrow(qual_vecs_50)), sep = "")
seqs <- rownames(qual_vecs_50)

rownames(qual_vecs_50) <- headers
rownames(qual_vecs_100) <- headers
rownames(qual_vecs_250) <- headers
rownames(qual_vecs_500) <- headers
rownames(qual_vecs_750) <- headers
write.table(qual_vecs_50, "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/embed/embed_.07_50dim.txt", quote = F, col.names = F)
write.table(qual_vecs_100, "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/embed/embed_.07_100dim.txt", quote = F, col.names = F)
write.table(qual_vecs_250, "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/embed/embed_.07_250dim.txt", quote = F, col.names = F)
write.table(qual_vecs_500, "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/embed/embed_.07_500dim.txt", quote = F, col.names = F)
write.table(qual_vecs_750, "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/embed/embed_.07_750dim.txt", quote = F, col.names = F)

headers_new <- paste(">", headers, sep= "")
fasta <- paste(headers_new, seqs, sep = "\n", collapse = "\n")
write(fasta, "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/embed/seqs_.07.fasta")
