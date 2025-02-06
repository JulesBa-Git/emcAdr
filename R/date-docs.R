#' ATC Tree Upper Bound 2014
#'
#' The expected form of the ATC tree, this version of the ATC tree is from the WHO
#' website (2024-02-23).
#'
#' @format A data frame 4 variables:
#' \describe{
#'   \item{ATCCode}{The code of ATC nodes}
#'   \item{Name}{The name of ATC nodes}
#'   \item{ATC_length}{The number of character in the ATCCode}
#'   \item{upperBound}{The index of the last son of the current node in the tree}
#' }
#' @source World Health Organization, ATC classification register
"ATC_Tree_UpperBound_2024"

#' FAERS myopathy
#'
#' The dataframe corresponding to patients drug intake from the FAERS dataset. It
#' shows expected input for the genetic algorithm as well as the MCMC algorithm.
#'
#' @format A data frame with 2 columns:
#' \describe{
#'   \item{patientATC}{Drug intake for each patient in the form of vector of index of the ATC tree}
#'   \item{patientADR}{Does the patient experience the considered Adverse Event ? (here : myopathy)}
#' }
#' @source Food & Drug Administration Event Reporting System data
"FAERS_myopathy"