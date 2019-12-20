library(msa)
library(phangorn)
library(ape)
library(phyloseq)

data_file = "/nfs3/PHARM/David_Lab/christine/quality_vectors_git/data/embed/embed_.07_100dim.txt"
qual_vecs <- read.table(data_file,
                        quote="\"", comment.char="", row.names = 1)
qual_vecs <- qual_vecs[rownames(qual_vecs) != "<unk>", ]

#2. read phylogenetic tree build using fasttree and clustalo msa
phy_tree = read.tree("../data/embed/tree.nwk")
print("read phy tree")

qual_vecs = qual_vecs[rownames(qual_vecs) %in% phy_tree$tip.label, ]

#1. write function to build tree (hierarchical clustering) from embeddings

hierarchical_cluster <- function(mat){
  dist_mat <- dist(mat, method = 'euclidean')
  hclust_avg <- hclust(dist_mat, method = 'average')
  return(hclust_avg)
}

#hclust_true <- hierarchical_cluster(qual_vecs)
#hclust_tree <- as.phylo(hclust_true)
#print("Tree built from embeddings")
#saveRDS(hclust_tree, "../data/hclust_embed_true.rds")
hclust_tree = readRDS("../data/hclust_embed_true.rds")
print("read hclust tree")


#3. write function to calculate distance between two trees

dist_true <- phangorn::treedist(hclust_tree, phy_tree)
#print("Calculated true distance")
saveRDS(dist_true, "../data/dist_true.rds")
#print("Calculated distance")


#4. write function to permute rownames on matrix

null_dists <- list()
num_iter = 200
phy_tree_permute = phy_tree
for(i in seq(1, num_iter)){
  #if(i %% 100 ==0){
  print(i)
  #}
  permute <- sample(phy_tree$tip.label, length(phy_tree$tip.label), replace = F)
  phy_tree_permute$tip.label <- permute
  dist_perm <- phangorn::treedist(hclust_tree, phy_tree_permute)
  print(dist_perm)
  null_dists[[i]] <- dist_perm
}
#plot(density(null_dists), xlim = c(min(null_dists) - 1, max(null_dists) + 1))
#abline(v = dist[2])
#p = sum(null_dists < dist_true[1]) / length(null_dists)
#p
print("calculated null distances")
saveRDS(null_dists, "../data/null_dists3.rds")

#Checks the distance between two trees in comparison to random trees


#####################################################
########### load objects and check results ##########
#####################################################
library(ggplot2)
data_dir = "C:/Users/ctata/Documents/Lab/quality_vectors_git/data/"
null_dists1 <- readRDS(paste(data_dir, "null_dists.rds", sep = ""))
null_dists2 <- readRDS(paste(data_dir, "null_dists2.rds", sep = ""))
null_dists3 <- readRDS(paste(data_dir, "null_dists3.rds", sep = ""))
null_dists4 <- readRDS(paste(data_dir, "null_dists4.rds", sep = ""))
true_dist <- readRDS(paste(data_dir, "dist_true.rds", sep = ""))

i = 1
sym_diffs_total <- c()
for(null_dist_list in list(null_dists1, null_dists2, null_dists3, null_dists4)){
  sym_diffs_total <- c(sym_diffs_total, unlist(lapply(null_dist_list, function(x) return(x[[i]]))))
}
ggplot(data.frame(sym_diffs_total))+
  geom_histogram(aes(x = sym_diffs_total, colour = "blue"), fill = "lightblue", color = "black", bins = 25)+
  geom_vline(aes(xintercept = true_dist[[i]], colour = "red"), linetype = 4, size = 2)+ theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size=25))


p = sum(sym_diffs_total < true_dist[[i]]) / length(sym_diffs_total)
p



i = 2
sym_diffs_total <- c()
for(null_dist_list in list(null_dists1, null_dists2, null_dists3, null_dists4)){
  sym_diffs_total <- c(sym_diffs_total, unlist(lapply(null_dist_list, function(x) return(x[[i]]))))
}

ggplot(data.frame(sym_diffs_total))+
  geom_histogram(aes(x = sym_diffs_total, colour = "blue"), fill = "lightblue", color = "black")+
  geom_vline(aes(xintercept = true_dist[[i]], colour = "red"), linetype = 4, size = 2)+ theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), text = element_text(size=25))
  


hist(sym_diffs_total, xlim = c(true_dist[[i]], max(sym_diffs_total) ))
abline(v = true_dist[[i]], col = "red")
p = sum(sym_diffs_total < true_dist[[i]]) / length(sym_diffs_total)
p