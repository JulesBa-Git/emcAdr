#' Calculate p-value of sampled value
#'
#' @param empirical_distribution A numeric vector of values representing the empirical distribution (return value of DistributionAproximation function)
#' @param sampled_values A scalar or a vector of real valued number representing the sampled value (score to be tested)
#' @param isFiltered A boolean representing if we want to use the filtered distribution or the distribution as is (False by default)
#' @param includeZeroValue A boolean that indicate if you want to take into account the null score (False by default)
#' @return A numeric value representing the empirical p-value
#'
#' @export
p_value_on_sampled <- function(empirical_distribution, sampled_values, isFiltered = F, includeZeroValue = F) {
  # Sort empirical distribution in ascending order (if the distribution comes. from 
  # the histogramToDitribution function it should already be sorted)
  if(isFiltered){
    if(includeZeroValue){
      empirical_distribution_array <- histogramToDitribution(empirical_distribution$FilteredDistribution)
      }else{
      empirical_distribution_array <- histogramToDitribution(empirical_distribution$FilteredDistribution[2:length(
        empirical_distribution$FilteredDistribution)])
      }
    empirical_distribution_array <- append(empirical_distribution_array, empirical_distribution$OutstandingRR)
  }
  else{
    if(includeZeroValue){
      empirical_distribution_array <- histogramToDitribution(empirical_distribution$Distribution)
    }else{
      empirical_distribution_array <- histogramToDitribution(empirical_distribution$Distribution[2:length(
        empirical_distribution$Distribution)])
    }
    empirical_distribution_array <- append(empirical_distribution_array, empirical_distribution$OutstandingRR)
  }
  
  empirical_distribution_array <- sort(empirical_distribution_array)
  
  p_values <- numeric(length(sampled_values))
  
  # Iterate over each sampled value
  for (i in seq_along(sampled_values)) {
    # Calculate ECDF value for the current sampled value
    ecdf_value <- sum(empirical_distribution_array <= sampled_values[i]) / length(empirical_distribution_array)
    
    # Calculate p-value
    p_values[i] <- 1 - ecdf_value
  }
  
  return(p_values) 
}

#' Calculate the divergence between 2 distributions (the true Distribution and the learned one)
#' @param empirical_distribution A numeric vector of values representing the empirical distribution (return value of DistributionAproximation function)
#' @param true_distribution A numeric vector of values representing the true distribution computed by the trueDistributionSizeTwoCocktail function
#' @param method A string, either "TV" or "KL" to respectively use the total variation distance or the Kullback-Leibler divergence. (default = "TV")
#' @return A numeric value representing the divergence of the 2 distributions
#' 
#' @export
calculate_divergence <- function(empirical_distribution, true_distribution, method = "TV", Filtered = F){
  RRmax <- (length(empirical_distribution$Distribution) - 1) / 10
  
  if(RRmax > 30){
    if(Filtered){
      dist_ouststandingRR <- incorporateOustandingRRToDistribution(true_distribution$OutstandingRR, RRmax)
      length(true_distribution$FilteredDistribution) <- length(dist_ouststandingRR)
      true_distribution$FilteredDistribution[is.na(true_distribution$FilteredDistribution)] <- 0
      true_distribution$FilteredDistribution <- true_distribution$FilteredDistribution + dist_ouststandingRR
    }
    else{
      dist_ouststandingRR <- incorporateOustandingRRToDistribution(true_distribution$OutstandingRR, RRmax)
      length(true_distribution$Distribution) <- length(dist_ouststandingRR)
      true_distribution$Distribution[is.na(true_distribution$Distribution)] <- 0
      true_distribution$Distribution <- true_distribution$Distribution + dist_ouststandingRR
    }
  }
  else if(RRmax < 30){
    if(Filtered){
      true_distribution$FilteredDistribution[((RRmax*10)+1)] <- sum(true_distribution$FilteredDistribution[((RRmax*10)+1):length(true_distribution$FilteredDistribution)]) + length(true_distribution$OutstandingRR)
      length(true_distribution$FilteredDistribution) <- ((RRmax*10)+1)
    }else{
      true_distribution$Distribution[((RRmax*10)+1)] <- sum(true_distribution$Distribution[((RRmax*10)+1):length(true_distribution$Distribution)]) + length(true_distribution$OutstandingRR)
      length(true_distribution$Distribution) <- ((RRmax*10)+1)
    }
  }
  else{
    if (Filtered) {
      true_distribution$FilteredDistribution[((RRmax*10)+1)] <- length(true_distribution$OutstandingRR)
    }else{
      true_distribution$Distribution[((RRmax*10)+1)] <- length(true_distribution$OutstandingRR)
    }
  }
  
  #empirical_distribution$Distribution <- empirical_distribution$Distribution / sum(empirical_distribution$Distribution)
  #true_distribution$Distribution <- true_distribution$Distribution / sum(true_distribution$Distribution)
  adjusted_empirical <- empirical_distribution$Distribution[2:length(empirical_distribution$Distribution)] / sum(empirical_distribution$Distribution[2:length(empirical_distribution$Distribution)])
  adjusted_true <- true_distribution$Distribution[2:length(true_distribution$Distribution)] / sum(true_distribution$Distribution[2:length(true_distribution$Distribution)])
  
  if(method == "TV"){
    #return(sum(abs(empirical_distribution$Distribution - true_distribution$Distribution)))
    return (sum(abs(adjusted_empirical - adjusted_true)))
  }
  else if(method == "KL"){
    #return(sum(true_distribution$Distribution * log(true_distribution$Distribution / empirical_distribution$Distribution)))
    return(sum(adjusted_true * log(adjusted_true / adjusted_empirical)))
  }
  else{
    stop("The method should be either \"TV\" for Total variation distance or \"KL\" for Kullback-Leibler divergence")
  }
}