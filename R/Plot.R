#' Plot the evolution of the mean and the best value of the population used by the GeneticAlgorithm
#'
#' @param list A list with 2 elements returned by the GeneticAlgorithm: "mean" and "best", containing the numeric vectors representing the mean and best fitness of the population
#' @param mean_color A string specifying the color of the mean values
#' @param best_color A string specifying the color of the best values
#' @param xlab A string specifying the label for the x-axis
#' @param ylab A string specifying the label for the y-axis
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom grid units
#' @export
#'
plot_evolution <- function(list, mean_color = "#F2A900", best_color = "#008080", xlab = "Epochs", ylab = "Score") {
  requireNamespace(dplyr)
  requireNamespace(ggplot2)
  
  epochs <- seq_along(list$meanFitnesses)
  
  data <- data.frame(epochs = epochs, mean = list$meanFitnesses, best = list$BestFitnesses)  
  
  ggplot2::ggplot(data, aes(x = epochs)) +
    geom_point(aes(y = best, color = "Best"), size = 2, shape = 16, fill = "white") +
    geom_point(aes(y = mean, color = "Mean"), size = 2, shape = 16, fill = "white") +
    geom_segment(aes(x = epochs, xend = dplyr::lead(epochs, default = last(epochs)), y = mean, yend = lead(mean, default = last(mean))), color = mean_color, size = 0.5, linetype = "dashed") +
    geom_segment(aes(x = epochs, xend = dplyr::lead(epochs, default = last(epochs)), y = best, yend = lead(best, default = last(best))), color = best_color, size = 0.5, linetype = "dashed") +
    scale_color_manual(values = c(best_color, mean_color)) +
    labs(title = "Evolution of the population", x = xlab, y = ylab) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = unit(16, "pt"), margin = margin(b = 10)),
          axis.title = element_text(face = "bold", size = unit(14, "pt")),
          axis.text = element_text(size = unit(12, "pt")),
          legend.title = element_blank(),
          legend.text = element_text(face = "bold", size = unit(12, "pt")),
          legend.position = "top")
}

#' Make a Quantile-Quantile diagram from the output of the MCMC algorithm (DistributionAproximation)
#' and the algorithm that exhaustively calculates the distribution
#' 
#' @param estimated Outputed object of DistributionApproximation function
#' @param true Outputed object of either DistributionApproximation function or True distribution
#' computation function
#' @param filtered Make use of the classic distributuion estimation or of the filtred one
#' (number of patient taking the cocktail > beta)
#' @param color The color of the dashed line of the qq-plot
#' 
#' @import ggplot2
#' @export
qq_plot_output <- function(estimated, true, filtered = F, color = "steelblue"){
  requireNamespace(ggplot2)
  if(filtered){
    estimated_distribution <- histogramToDitribution(estimated$FilteredDistribution[2:length(estimated$FilteredDistribution)])
    true_distribution <- histogramToDitribution(true$FilteredDistribution[2:length(true$FilteredDistribution)])
  }
  else{
    estimated_distribution <- histogramToDitribution(estimated$Distribution[2:length(estimated$Distribution)])
    true_distribution <- histogramToDitribution(true$Distribution[2:length(true$Distribution)])
  }
  
  num_quantiles <- min(length(estimated_distribution), length(true_distribution))
  probs <- seq(0,1, length.out = num_quantiles)
  
  quantiles_estim <- quantile(estimated_distribution, probs)
  quantiles_true <- quantile(true_distribution, probs)
  
  qq_df <- data.frame(estimated_quantiles = quantiles_estim, true_quantiles = quantiles_true)
  View(qq_df)
  ggplot(qq_df, aes(x = estimated_quantiles, y = true_quantiles)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = color) + # Adds a reference line y = x
    theme_minimal() +
    labs(x = "Estimated Distribution", y = "True Distribution", title = "QQ Plot of Estimated vs True distribution")
}

#' Plot the histogram of the approximation of the RR distribution 
#' @param freq_array The returned value of DistributionApproximation function
#' @param binwidth The width of the histogram bins
#' @param hist_color The fill color for the histogram bars
#' @param density_color The color for the density curve
#' @param sqrt A Boolean to specify whether we normalize the freq_array or not, it is recommended on large random walk.
#' @param xlab Label of X axis
#'
#' @import ggplot2
#' @importFrom dplyr data_frame
#' @export
plot_frequency <- function(freq_array, sqrt = F, binwidth = .1, hist_color = "#69b3a2", density_color = "#FF5733",
                           xlab = "Score") {
  requireNamespace(dplyr)
  requireNamespace(ggplot2)
  
  # Create a data frame from the returned value array
  if(sqrt){
    df <- data.frame(x = histogramToDitribution(sqrt(freq_array)))
    y_lab = "sqrt(frequency)"
  }else{
    df <- data.frame(x = histogramToDitribution(freq_array))
    y_lab = "frequency"
  }
  
  # Create histogram plot
  ggplot(df, aes(x = x, y = after_stat(density))) +
    geom_histogram(binwidth = binwidth, fill = hist_color, color = "#e9ecef") +
    labs(title = "The approximation of the distribution histogram", x = xlab, y = y_lab) +
    theme_minimal()# +
    #geom_density(color = density_color, size = 1, alpha = .2, fill = density_color)
}


