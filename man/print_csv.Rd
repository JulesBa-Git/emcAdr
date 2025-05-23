% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{print_csv}
\alias{print_csv}
\title{Print every cocktails found during the genetic algorithm when used with the 
hyperparam_test_genetic_algorithm function. This enables to condense the solutions 
found in each files by collapsing similar cocktail in a single row by cocktail.}
\usage{
print_csv(
  input_filenames,
  observations,
  repetition,
  ATCtree,
  csv_filename = "solutions.csv"
)
}
\arguments{
\item{input_filenames}{: A List containing filename of hyperparam_test_genetic_algorithm output file}

\item{observations}{: observation of the AE based on the medications of each patients
(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)
on which we want to compute the risk distribution}

\item{repetition}{: The parameter nb_test_desired used in the hyperparam test function}

\item{ATCtree}{: ATC tree with upper bound of the DFS (without the root)}

\item{csv_filename}{: Name of the output file, "solutions.csv" by default}
}
\value{
No return value, should process the output of the genetic algorithm in 
files produced by hyperparam_test_genetic_algorithm and output a summary csv file.
The csv file is outputed in current directory and named after the csv_filename
variable (solutions.csv by default).
}
\description{
Print every cocktails found during the genetic algorithm when used with the 
hyperparam_test_genetic_algorithm function. This enables to condense the solutions 
found in each files by collapsing similar cocktail in a single row by cocktail.
}
\examples{
\donttest{
 data("ATC_Tree_UpperBound_2024")
 data("FAERS_myopathy")
 files = c('250e_700ind_0.2mr_0ne_2alpha.txt') # results of hyperparam_test_genetic_algorithm

 print_csv(input_filenames = files, observations = FAERS_myopathy,
          repetition = 5, ATCtree = ATC_Tree_UpperBound_2024)
}
}
