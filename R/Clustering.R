#' Clustering of the solutions of the genetic algorithm
#' 
#' @param genetic_results The return value of the genetic algorithm
#' @param ATCtree ATC tree with upper bound of the DFS
#' @param method (from hclust function) the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param dist.normalize Do we normalize the distance (so it bellongs to [0;1])
#' @return the hierarchical clustering of the results of the genetic algorithm
#' 
#' @export
hclust_genetic_solution <- function(genetic_results,ATCtree, dist.normalize = T,
                                    method = "complete"){
  
  if(dist.normalize){
    divergence <- get_dissimilarity(genetic_results, ATCtree)
  }else{
    divergence <- get_dissimilarity(genetic_results, ATCtree, F)
  }
  divergence <- do.call(rbind,divergence)
  divergence <- as.dist(divergence)
  hc <- hclust(divergence, method = method)
  
  return (hc)
}

#' @export
tsne_genetic <- function(genetic_results,ATCtree, true_solutions,dim=2, dist.normalize = T){
  library(tsne)
  if(dist.normalize){
    divergence <- get_dissimilarity(genetic_results, ATCtree)
  }else{
    divergence <- get_dissimilarity(genetic_results, ATCtree, F)
  }
  divergence <- do.call(rbind,divergence)
  #divergence <- as.dist(divergence)
  
  tsn <- tsne(divergence, k=dim)
  tsne_df <- data.frame(tsn)
  
  info_df <- get_answer_class(genetic_results,
                              true_solutions)
  colnames(tsne_df) <- c("x","y")
  merged_df <- cbind(tsne_df, info_df)
  
  return (merged_df)
}

clustering_from_df <- function(df_genetic_results, ATCtree, 
                               dim=2, dist.normalize = T){
  library(tsne)
  
}