library(KEGGREST)
setwd("C:/Users/ctata/Documents/Lab/quality_vectors_git/scripts")
allSigCorrs <- readRDS("../data/pathways/allSigCorrs.rds")
imp_half <- read.csv("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/halfvarson/metabolic_pathways_importance_half.csv")
imp_ag <- read.csv("C:/Users/ctata/Documents/Lab/quality_vectors_git/data/AG_new/metabolic_pathways_importance_ag.csv")
corrs <- read.csv("../data/pathways/top20Paths_per_property_ids.csv")


imp_half <- imp_half[order(imp_half$X), ]
imp_ag <- imp_ag[order(imp_ag$X), ]

n = 0
health_props <- as.character(imp_half$X[imp_half$cum_timepoint_numtrees < -n & imp_ag$cum_timepoint_numtrees < -n])
ibd_props <- as.character(imp_half$X[imp_half$cum_timepoint_numtrees > n & imp_ag$cum_timepoint_numtrees > n])
health_prop_id <- as.numeric(gsub("property_100_", "", health_props))
ibd_prop_id <- as.numeric(gsub("property_100_", "", ibd_props))


pathway_names <- unique(unlist(lapply(allSigCorrs, function(x) return(names(x)))))
health_pathways <- sapply(pathway_names,function(x) 0)
ibd_pathways <- sapply(pathway_names,function(x) 0)

for(id in health_prop_id){
  for(path in names(allSigCorrs[[id]])){
    health_pathways[[path]] = health_pathways[[path]] + 1
  }
}

for(id in ibd_prop_id){
  for(path in names(allSigCorrs[[id]])){
    ibd_pathways[[path]] = ibd_pathways[[path]] + 1
  }
}


health_pathways_imp <- health_pathways[health_pathways > 0]
ibd_pathways_imp <- ibd_pathways[ibd_pathways > 0]



##Mantel test
cosineSim <- function(x, y){
  apply()
  dot = x %*% y
  magx = sqrt(sum(x^2))
  magy = sqrt(sum(y^2))
  return(dot / (magx * magy))
}
health_sim <- cosine(as.matrix(glove_emb[ , health_prop_id]))
ibd_sim <- cosine(as.matrix(glove_emb[ , ibd_prop_id]))
mantel.test(health_prop_id, ibd_prop_id)


names(health_pathways_imp) %in% names(ibd_pathways_imp)

health_pathway_names <- names(health_pathways_imp[! names(health_pathways_imp) %in% names(ibd_pathways_imp)])
ibd_pathway_names <- names(ibd_pathways_imp[!names(ibd_pathways_imp) %in% names(health_pathways_imp)])


health_fullpath_names <- sapply(health_pathway_names, function(path_id){
  print(path_id)
  entry <- keggGet(paste("map", path_id, sep = ""))
  path_name <- entry[[1]]$NAME
  class_name <- entry[[1]]$CLASS
  return(path_name)
})


ibd_fullpath_names <- sapply(ibd_pathway_names, function(path_id){
  print(path_id)
  entry <- keggGet(paste("map", path_id, sep = ""))
  path_name <- entry[[1]]$NAME
  class_name <- entry[[1]]$CLASS
  return(list(path_name, class_name))
})


getPathwayName <- function(path_id){
  entry <- keggGet(paste("map", path_id, sep = ""))
  path_name <- entry[[1]]$NAME
  class_name <- entry[[1]]$CLASS
  return(path_name)
}



allSigCorrs_pathNames <- lapply(allSigCorrs, function(x){
  path_names <- sapply(names(x), function(path_id) return(getPathwayName(path_id)))
  return(paste(names(path_names), path_names, sep = ": "))
})


lapply(names(allSigCorrs_pathNames), function(prop_name){
  x <- c(prop_name, allSigCorrs_pathNames[[prop_name]])
  write(x, "../figures/final_figures/supp_table_1.txt", append = TRUE)
} )
