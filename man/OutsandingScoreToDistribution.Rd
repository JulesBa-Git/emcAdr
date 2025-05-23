% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{OutsandingScoreToDistribution}
\alias{OutsandingScoreToDistribution}
\title{Output the outstanding score (Outstanding_score) outputed by the MCMC algorithm
in a special format}
\usage{
OutsandingScoreToDistribution(outstanding_score, max_score)
}
\arguments{
\item{outstanding_score}{: Outstanding_score outputed by MCMC algorithm to be converted
to the ScoreDistribution format}

\item{max_score}{: max_score parameter used during the MCMC algorithm}
}
\value{
outstanding_score in a format compatible with MCMC algorithm output
}
\description{
Output the outstanding score (Outstanding_score) outputed by the MCMC algorithm
in a special format
}
\examples{
\donttest{
 data("ATC_Tree_UpperBound_2024")
 data("FAERS_myopathy")

  DistributionApproximationResults = DistributionApproximation(epochs = 10,
            ATCtree = ATC_Tree_UpperBound_2024, observations = FAERS_myopathy)
  OutsandingScoreToDistribution(DistributionApproximationResults$Outstanding_score, max_score = 100)
}
}
