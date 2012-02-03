\name{resight.matrix}
\alias{resight.matrix}
\title{Various utility functions}
\usage{
  resight.matrix(x) naive.survival(x,...)
}
\arguments{
  \item{x}{processed data list - result from process.data
  in marked or real estimates from fitted model}

  \item{age}{at which Phi or p should be shown across time}

  \item{time}{at which Phi or p should be shown across
  ages}
}
\value{
  matrix of values with cohort and year labels
}
\description{
  Several functions have been added to help visualize data
  including resight.matrix which provides for each cohort
  the number of releases on the diagonal and the number
  resighted on each occasion in the upper-triangular
  matrix. naive.survival provides a naive survival estimate
  at each time for each cohort.  The estimate for time i is
  the number resighted at time i+1 or later divided by the
  number seen at time i or later. These values should be
  interpreted cautiously because they are influenced by
  capture probability in year i but it is useful to
  identify particularly high or low survival values.
  Functions Phi.mean and p.mean compute average real
  parameter values Phi or p across time for a single age or
  across ages for a single time.
}
\author{
  Jeff Laake
}
\keyword{utility}

