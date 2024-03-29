% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{DistributionApproximation}
\alias{DistributionApproximation}
\title{The Evolutionary MCMC method that runs the random walk on a single cocktail}
\usage{
DistributionApproximation(
  epochs,
  ATCtree,
  observations,
  temperature = 1L,
  nbResults = 5L,
  Smax = 4L,
  p_type1 = 0.01,
  beta = 4L,
  RRmax = 100L
)
}
\arguments{
\item{epochs}{: number of step}

\item{ATCtree}{: ATC tree with upper bound of the DFS (without the root)}

\item{temperature}{: starting temperature, default = 1}

\item{observation}{: real observation of the ADR based on the medications of each real patients
(a DataFrame containing the medication on the first column and the ADR (boolean) on the second)}

\item{p_type1:}{probability to operate type1 mutation. Note :
the probability to operate the type 2 mutation is then 1 - P_type1. P_type1 must be in [0;1].}

\item{alpha}{: a hyperparameter allowing us to manage to probability of adding a drug to the cocktail. The probability
to add a drug to the cocktail is the following : \eqn{\frac12}{\alpha/n} Where n is the original size of the cocktail. 1 is the default value.}
}
\value{
if no problem return an array of the approximation of the RR distribution : the distribution of RR we've met; Otherwise the list is empty
}
\description{
The Evolutionary MCMC method that runs the random walk on a single cocktail
}
