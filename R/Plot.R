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
plot_evolution <- function(list, mean_color = "#F2A900", best_color = "#008080", xlab = "Epochs", ylab = "RR") {
  
  epochs <- seq_along(list$meanFitnesses)
  
  data <- data.frame(epochs = epochs, mean = list$meanFitnesses, best = list$BestFitnesses)  
  
  ggplot(data, aes(x = epochs)) +
    geom_point(aes(y = best, color = "Best"), size = 2, shape = 16, fill = "white") +
    geom_point(aes(y = mean, color = "Mean"), size = 2, shape = 16, fill = "white") +
    geom_segment(aes(x = epochs, xend = lead(epochs, default = last(epochs)), y = mean, yend = lead(mean, default = last(mean))), color = mean_color, size = 0.5, linetype = "dashed") +
    geom_segment(aes(x = epochs, xend = lead(epochs, default = last(epochs)), y = best, yend = lead(best, default = last(best))), color = best_color, size = 0.5, linetype = "dashed") +
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

#' Plot the histogram of the approximation of the RR distribution 
#' @param freq_array The returned value of DistributionApproximation function
#' @param binwidth The width of the histogram bins
#' @param hist_color The fill color for the histogram bars
#' @param density_color The color for the density curve
#' @param sqrt A Boolean to specify whether we normalize the freq_array or not, it is recommended on large random walk.
#'
#' @import ggplot2
#' @importFrom dplyr data_frame
#' @export
plot_frequency <- function(freq_array, sqrt = T,binwidth = .3, hist_color = "#69b3a2", density_color = "#FF5733") {
  
  # Create a data frame from the returned value array
  if(sqrt){
    df <- data.frame(x = histogramToDitribution(sqrt(freq_array)))
    y_lab = "sqrt(frequency)"
  }else{
    df <- data.frame(x = histogramToDitribution(freq_array))
    y_lab = "frequency"
  }
  
  # Create histogram plot
  ggplot(df, aes(x = x, y = ..density..)) +
    geom_histogram(binwidth = binwidth, fill = hist_color, color = "#e9ecef") +
    labs(title = "The approximation of the distribution histogram", x = "RR", y = y_lab) +
    theme_minimal() +
    geom_density(color = density_color, size = 1, alpha = .2, fill = density_color)
}

